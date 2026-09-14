//! High-performance parallel pipeline for BAM processing.

use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::mpsc::sync_channel;
use std::thread;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{self, Read};

use crate::cigar::build_cigar;
use crate::consensus::{call_consensus, BaseAccumulator, ConsensusTables};
use crate::interval::{Interval, IntervalSet};
use crate::isoform::{
    build_isoform_indices_from_gtf, dump_isoform_json, find_compatible_transcripts,
    load_isoform_json_with_kind, GeneIsoformMap, IndexKind,
};
use crate::matrix::{IsoformQuantifier, MoleculeSummary};

fn temporary_sibling(path: &Path, suffix: &str) -> PathBuf {
    let name = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("output");
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|value| value.as_nanos())
        .unwrap_or_default();
    path.with_file_name(format!(
        ".{}.{}-{}-{}",
        name,
        suffix,
        std::process::id(),
        nonce
    ))
}

fn prepare_bam_output(input: &Path, output: &Path) -> Result<()> {
    if let Some(parent) = output.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("Creating BAM output directory: {:?}", parent))?;
        }
    }

    let input_path = fs::canonicalize(input)
        .with_context(|| format!("Resolving input BAM path: {:?}", input))?;
    let output_path = if output.exists() {
        fs::canonicalize(output)
            .with_context(|| format!("Resolving output BAM path: {:?}", output))?
    } else {
        let parent = output
            .parent()
            .filter(|path| !path.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let file_name = output
            .file_name()
            .context("BAM output path has no file name")?;
        fs::canonicalize(parent)
            .with_context(|| format!("Resolving BAM output directory: {:?}", parent))?
            .join(file_name)
    };

    if input_path == output_path {
        bail!("Input and output BAM paths must be different: {:?}", input);
    }
    Ok(())
}

fn remove_path_if_present(path: &Path) {
    let _ = fs::remove_file(path);
}

#[derive(Debug, Clone)]
pub struct GeneInfo {
    pub gene_id: String,
    pub seqid: String,
    pub start: i64,
    pub end: i64,
    pub tid: u32,
}

#[derive(Debug, Clone)]
pub struct PipelineConfig {
    pub input_bam: PathBuf,
    pub output_bam: PathBuf,
    pub gtf_file: PathBuf,
    pub isoform_file: Option<PathBuf>,
    pub junction_file: Option<PathBuf>,
    pub dump_index: Option<PathBuf>,
    pub matrix_out: Option<PathBuf>,
    pub counts_tsv: Option<PathBuf>,
    pub molecules_tsv: Option<PathBuf>,
    pub threads: usize,
    pub single_end: bool,
    pub skip_iso: bool,
    pub umi_tag: String,
    pub cell_tag: String,
    pub cells_file: Option<PathBuf>,
    pub genes_file: Option<PathBuf>,
    pub contig: Option<String>,
    pub gene_identifier: String,
}

#[derive(Debug, Clone)]
struct AlignedBase {
    pos: i64,
    base: u8,
    qual: u8,
}

struct ReadData {
    qname_id: u32,
    is_read1: bool,
    is_reverse: bool,
    is_exonic: bool,
    is_intronic: bool,
    aligned_bases: Vec<AlignedBase>,
    skipped_intervals: Vec<Interval>,
}

type GeneBatch = (Vec<bam::Record>, Vec<MoleculeSummary>);
type ClusterResult = (usize, Vec<bam::Record>, Vec<MoleculeSummary>);

struct GeneProcessingContext<'a> {
    config: &'a PipelineConfig,
    cell_set: Option<&'a HashSet<String>>,
    tables: &'a ConsensusTables,
    isoform_map: Option<&'a GeneIsoformMap>,
    refskip_map: Option<&'a GeneIsoformMap>,
    collect_quantification: bool,
}

/// Load cell whitelist
fn load_cell_set<P: AsRef<Path>>(path: P) -> Result<HashSet<String>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut set = HashSet::new();
    for line in reader.lines() {
        let l = line?.trim().to_string();
        if !l.is_empty() {
            set.insert(l);
        }
    }
    Ok(set)
}

/// Load gene whitelist
fn load_gene_set<P: AsRef<Path>>(path: P) -> Result<HashSet<String>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut set = HashSet::new();
    for line in reader.lines() {
        let l = line?.trim().to_string();
        if !l.is_empty() {
            set.insert(l);
        }
    }
    Ok(set)
}

/// Extract attribute value by key from GTF line
fn extract_gtf_attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    for part in attrs.split(';') {
        let part = part.trim();
        if let Some(rem) = part.strip_prefix(key) {
            if rem.starts_with(' ') || rem.starts_with('\t') || rem.starts_with('=') {
                let val = rem.trim().trim_matches('"');
                return Some(val);
            }
        }
    }
    None
}

/// Parse genes from GTF file
pub fn load_genes_from_gtf(
    gtf_path: &Path,
    header: &bam::HeaderView,
    contig_filter: Option<&str>,
    gene_set: Option<&HashSet<String>>,
    gene_identifier: &str,
) -> Result<Vec<GeneInfo>> {
    let file = File::open(gtf_path).with_context(|| format!("Opening GTF: {:?}", gtf_path))?;
    let reader = BufReader::new(file);

    let mut genes_map: HashMap<String, GeneInfo> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        let feature = fields[2];
        if feature != "gene" {
            continue;
        }
        let seqid = fields[0];
        if let Some(c) = contig_filter {
            if seqid != c {
                continue;
            }
        }

        // Check if contig exists in BAM header
        let tid = match header.tid(seqid.as_bytes()) {
            Some(t) => t,
            None => continue,
        };

        let start: i64 = match fields[3].parse() {
            Ok(s) => s,
            Err(_) => continue,
        };
        let end: i64 = match fields[4].parse() {
            Ok(e) => e,
            Err(_) => continue,
        };

        let gene_id = match extract_gtf_attr(fields[8], gene_identifier) {
            Some(id) => id.to_string(),
            None => {
                // Fallback to gene_id
                match extract_gtf_attr(fields[8], "gene_id") {
                    Some(id) => id.to_string(),
                    None => continue,
                }
            }
        };

        if let Some(set) = gene_set {
            if !set.contains(&gene_id) {
                continue;
            }
        }

        genes_map.entry(gene_id.clone()).or_insert(GeneInfo {
            gene_id,
            seqid: seqid.to_string(),
            start,
            end,
            tid,
        });
    }

    let mut genes: Vec<GeneInfo> = genes_map.into_values().collect();
    genes.sort_by(|a, b| {
        (a.tid, a.start, a.end, a.gene_id.as_str()).cmp(&(
            b.tid,
            b.start,
            b.end,
            b.gene_id.as_str(),
        ))
    });
    Ok(genes)
}

/// Dual-cursor alignment extractor: synchronizes CIGAR operations with query sequence and quality scores.
/// Correctly handles SoftClip (S) and Insertion (I) without coordinate-sequence desynchronization,
/// and handles Deletion (D) and RefSkip (N) without misaligning query bases.
fn extract_alignment(record: &bam::Record) -> (Vec<AlignedBase>, Vec<Interval>) {
    let mut aligned = Vec::with_capacity(record.seq_len());
    let mut skipped = Vec::new();

    let seq = record.seq().as_bytes();
    let qual = record.qual();
    let mut curr_ref = record.pos();
    let mut read_idx = 0usize;

    for op in record.cigar().iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let l = *len as usize;
                for _ in 0..l {
                    if read_idx < seq.len() {
                        aligned.push(AlignedBase {
                            pos: curr_ref,
                            base: seq[read_idx],
                            qual: qual[read_idx],
                        });
                        read_idx += 1;
                        curr_ref += 1;
                    }
                }
            }
            Cigar::Ins(len) => {
                // Consumes query sequence only
                read_idx += *len as usize;
            }
            Cigar::SoftClip(len) => {
                // Consumes query sequence only
                read_idx += *len as usize;
            }
            Cigar::Del(len) => {
                // Consumes reference coordinate only
                curr_ref += *len as i64;
            }
            Cigar::RefSkip(len) => {
                // Consumes reference coordinate only (intron junction)
                let l = *len as i64;
                if l > 0 {
                    skipped.push(Interval::new(curr_ref, curr_ref + l - 1));
                    curr_ref += l;
                }
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    (aligned, skipped)
}

/// Process a single gene's reads and return stitched BAM records and molecule summaries
fn process_single_gene(
    reader: &mut bam::IndexedReader,
    gene: &GeneInfo,
    context: &GeneProcessingContext<'_>,
) -> Result<GeneBatch> {
    let config = context.config;
    let start_pos = (gene.start - 1).max(0);
    let end_pos = gene.end;

    reader
        .fetch((gene.tid, start_pos, end_pos))
        .with_context(|| {
            format!(
                "Failed to fetch BAM region for gene {} (tid={}, {}-{})",
                gene.gene_id, gene.tid, start_pos, end_pos
            )
        })?;

    let cell_tag_bytes = config.cell_tag.as_bytes();
    let umi_tag_bytes = config.umi_tag.as_bytes();

    let mut read_dict: HashMap<(String, String), Vec<ReadData>> = HashMap::new();
    let mut qname_ids: HashMap<Vec<u8>, u32> = HashMap::new();

    let mut record = bam::Record::new();
    while let Some(res) = reader.read(&mut record) {
        res.with_context(|| format!("Error reading BAM record in gene {}", gene.gene_id))?;

        if record.is_unmapped() {
            continue;
        }

        // Match the Python/zUMIs counting contract: secondary and supplementary
        // alignments are not independent molecules. Duplicate/QC flags are left
        // untouched because upstream filtering owns that policy.
        if record.is_secondary() || record.is_supplementary() {
            continue;
        }

        if !config.single_end
            && (!record.is_paired() || record.is_mate_unmapped() || !record.is_proper_pair())
        {
            continue;
        }

        // Cell tag check
        let cell = match record.aux(cell_tag_bytes) {
            Ok(Aux::String(s)) => s.to_string(),
            _ => continue,
        };

        if let Some(set) = context.cell_set {
            if !set.contains(&cell) {
                continue;
            }
        }

        // UMI tag check
        let umi = match record.aux(umi_tag_bytes) {
            Ok(Aux::String(s)) => {
                let s = s.trim();
                if s.is_empty() {
                    continue;
                }
                s.to_string()
            }
            _ => continue,
        };

        // Gene assignment check
        let (is_exonic, gene_exon) = match record.aux(b"GE") {
            Ok(Aux::String(s)) => (true, Some(s)),
            _ => match (record.aux(b"GX"), record.aux(b"RE")) {
                (Ok(Aux::String(gx)), Ok(Aux::Char(b'E') | Aux::String("E"))) => (true, Some(gx)),
                _ => (false, None),
            },
        };

        let (is_intronic, gene_intron) = match record.aux(b"GI") {
            Ok(Aux::String(s)) => (true, Some(s)),
            _ => match (record.aux(b"GX"), record.aux(b"RE")) {
                (Ok(Aux::String(gx)), Ok(Aux::Char(b'N') | Aux::String("N"))) => (true, Some(gx)),
                _ => (false, None),
            },
        };

        let assigned_gene = match (gene_exon, gene_intron) {
            (Some(e), Some(i)) => {
                if e == i {
                    Some(e)
                } else {
                    None
                }
            }
            (Some(e), None) => Some(e),
            (None, Some(i)) => Some(i),
            (None, None) => None,
        };

        let assigned = match assigned_gene {
            Some(g) => g,
            None => continue,
        };

        if assigned != gene.gene_id {
            continue;
        }

        let (aligned_bases, skipped) = extract_alignment(&record);
        if aligned_bases.is_empty() {
            continue;
        }

        let qname_id = match qname_ids.get(record.qname()) {
            Some(&id) => id,
            None => {
                let id = qname_ids.len() as u32;
                qname_ids.insert(record.qname().to_vec(), id);
                id
            }
        };

        let r_data = ReadData {
            qname_id,
            is_read1: record.is_first_in_template(),
            is_reverse: record.is_reverse(),
            is_exonic,
            is_intronic,
            aligned_bases,
            skipped_intervals: skipped,
        };

        read_dict.entry((cell, umi)).or_default().push(r_data);
    }

    // ReadData keeps only compact IDs after ingestion; release the original
    // qname strings before the per-molecule consensus pass.
    drop(qname_ids);

    let mut out_records = Vec::new();
    let mut out_summaries = Vec::new();

    for ((cell, umi), mol) in read_dict {
        if mol.is_empty() {
            continue;
        }

        let n_read1 = mol.iter().filter(|r| r.is_read1).count();
        if !config.single_end && n_read1 == 0 {
            continue;
        }

        // Keep one qname map with two classification bits instead of three
        // separate HashSets.  ReadData stores compact per-gene qname IDs, so
        // this pass hashes integers rather than cloning or hashing qname text.
        let mut qname_flags: HashMap<u32, u8> = HashMap::new();
        let mut rev_count = 0;
        let mut fwd_count = 0;

        let mut all_pos_set: HashSet<i64> = HashSet::new();
        let mut all_skipped: Vec<Interval> = Vec::new();

        for r in &mol {
            let flags = qname_flags.entry(r.qname_id).or_insert(0);
            if r.is_exonic {
                *flags |= 0b01;
            }
            if r.is_intronic {
                *flags |= 0b10;
            }
            if config.single_end || r.is_read1 {
                if r.is_reverse {
                    rev_count += 1;
                } else {
                    fwd_count += 1;
                }
            }
            for b in &r.aligned_bases {
                all_pos_set.insert(b.pos);
            }
            all_skipped.extend_from_slice(&r.skipped_intervals);
        }

        if all_pos_set.is_empty() {
            continue;
        }

        let n_reads = qname_flags.len() as i32;
        let n_exonic = qname_flags
            .values()
            .filter(|flags| **flags & 0b01 != 0)
            .count() as i32;
        let n_intronic = qname_flags
            .values()
            .filter(|flags| **flags & 0b10 != 0)
            .count() as i32;
        let is_reverse = rev_count >= fwd_count;

        let mut sorted_ref: Vec<i64> = all_pos_set.into_iter().collect();
        sorted_ref.sort_unstable();

        // Map position to index for accumulator
        let mut accumulators = vec![BaseAccumulator::new(); sorted_ref.len()];

        for r in &mol {
            for b in &r.aligned_bases {
                if let Ok(idx) = sorted_ref.binary_search(&b.pos) {
                    accumulators[idx].add_read_base(b.base, b.qual, context.tables);
                }
            }
        }

        let (seq, qual) = call_consensus(&accumulators);

        let ref_intervals = IntervalSet::from_positions(&sorted_ref);
        let skipped_intervals = IntervalSet::from_intervals(all_skipped);

        let cigar_res = match build_cigar(&ref_intervals, &skipped_intervals) {
            Some(c) => c,
            None => continue,
        };

        // Construct htslib CIGAR
        let mut hts_cigar = Vec::with_capacity(cigar_res.ops.len());
        for op in &cigar_res.ops {
            match op.op {
                'M' => hts_cigar.push(Cigar::Match(op.len)),
                'N' => hts_cigar.push(Cigar::RefSkip(op.len)),
                'D' => hts_cigar.push(Cigar::Del(op.len)),
                _ => {}
            }
        }
        let cigar_string = CigarString(hts_cigar);

        let mut out_rec = bam::Record::new();
        let qname = format!("{}:{}:{}", cell, gene.gene_id, umi);
        let pos_0based = cigar_res.pos_0based;

        out_rec.set(qname.as_bytes(), Some(&cigar_string), &seq, &qual);

        out_rec.set_tid(gene.tid as i32);
        out_rec.set_pos(pos_0based);
        out_rec.set_mapq(255);
        if is_reverse {
            out_rec.set_flags(16);
        } else {
            out_rec.set_flags(0);
        }

        // Tags
        let _ = out_rec.push_aux(b"NR", Aux::I32(n_reads));
        let _ = out_rec.push_aux(b"ER", Aux::I32(n_exonic));
        let _ = out_rec.push_aux(b"IR", Aux::I32(n_intronic));
        let _ = out_rec.push_aux(b"CB", Aux::String(&cell));
        let _ = out_rec.push_aux(b"GX", Aux::String(&gene.gene_id));
        let _ = out_rec.push_aux(umi_tag_bytes, Aux::String(&umi));

        if cigar_res.conflict {
            let _ = out_rec.push_aux(b"NC", Aux::I32(cigar_res.n_conflict as i32));
        }

        let mut ct_opt = None;
        if !config.skip_iso {
            if let Some(ct) = find_compatible_transcripts(
                &gene.gene_id,
                &ref_intervals,
                &skipped_intervals,
                context.isoform_map,
                context.refskip_map,
            ) {
                let _ = out_rec.push_aux(b"CT", Aux::String(&ct));
                ct_opt = Some(ct);
            }
        }

        out_records.push(out_rec);
        if context.collect_quantification {
            out_summaries.push(MoleculeSummary {
                cell: cell.clone(),
                gene: gene.gene_id.clone(),
                umi: umi.clone(),
                transcripts: ct_opt.map(|s| {
                    s.split(',')
                        .map(|x| x.trim().to_string())
                        .filter(|x| !x.is_empty())
                        .collect()
                }),
                reads: n_reads,
            });
        }
    }

    Ok((out_records, out_summaries))
}

fn gene_cluster_ranges(genes: &[GeneInfo]) -> Vec<(usize, usize)> {
    if genes.is_empty() {
        return Vec::new();
    }

    let mut ranges = Vec::new();
    let mut cluster_start = 0;
    let mut cluster_end = genes[0].end;

    for (index, gene) in genes.iter().enumerate().skip(1) {
        if gene.start > cluster_end.saturating_add(1) {
            ranges.push((cluster_start, index));
            cluster_start = index;
            cluster_end = gene.end;
        } else {
            cluster_end = cluster_end.max(gene.end);
        }
    }
    ranges.push((cluster_start, genes.len()));
    ranges
}

/// Run the full pipeline
pub fn run_pipeline(config: PipelineConfig) -> Result<()> {
    let start_time = Instant::now();

    prepare_bam_output(&config.input_bam, &config.output_bam)?;

    if config.matrix_out.is_some() && config.skip_iso {
        bail!("--matrix-out requires isoform calling; remove --skip-iso");
    }

    if !config.skip_iso {
        match (&config.isoform_file, &config.junction_file) {
            (Some(_), None) | (None, Some(_)) => {
                bail!("--isoform and --junction must be provided together");
            }
            _ => {}
        }
    }

    let cell_set = if let Some(ref path) = config.cells_file {
        println!("Loading cell barcodes from {:?}", path);
        Some(load_cell_set(path)?)
    } else {
        None
    };

    let gene_set = if let Some(ref path) = config.genes_file {
        println!("Loading gene filter list from {:?}", path);
        Some(load_gene_set(path)?)
    } else {
        None
    };

    let header = {
        let bam = bam::Reader::from_path(&config.input_bam)
            .with_context(|| format!("Opening input BAM: {:?}", config.input_bam))?;
        bam::Header::from_template(bam.header())
    };

    // Fail before creating any output when the required BAM index is missing
    // or unreadable. Each worker relies on IndexedReader::fetch below.
    bam::IndexedReader::from_path(&config.input_bam)
        .with_context(|| format!("Opening BAM index for input: {:?}", config.input_bam))?;

    println!("Reading gene coordinates from {:?}", config.gtf_file);
    let genes = {
        let header_view = bam::HeaderView::from_header(&header);
        load_genes_from_gtf(
            &config.gtf_file,
            &header_view,
            config.contig.as_deref(),
            gene_set.as_ref(),
            &config.gene_identifier,
        )?
    };

    if genes.is_empty() {
        bail!("No matching genes found from GTF and BAM contigs.");
    }
    println!("Total genes scheduled for stitching: {}", genes.len());

    let (isoform_map, refskip_map) = if !config.skip_iso {
        if config.isoform_file.is_some() && config.junction_file.is_some() {
            println!("Loading isoform dictionary from JSON...");
            let iso = load_isoform_json_with_kind(
                config.isoform_file.as_ref().unwrap(),
                IndexKind::Exon,
            )?;
            println!("Loading junction dictionary from JSON...");
            let jun = load_isoform_json_with_kind(
                config.junction_file.as_ref().unwrap(),
                IndexKind::RefSkip,
            )?;
            (Some(iso), Some(jun))
        } else {
            println!("Building in-memory isoform and junction index directly from GTF (Zero-friction on-the-fly mode)...");
            let index_start = Instant::now();
            let (iso, jun) = build_isoform_indices_from_gtf(
                &config.gtf_file,
                config.contig.as_deref(),
                gene_set.as_ref(),
                &config.gene_identifier,
            )?;
            println!(
                "Indexed {} genes for isoform calling in {:.2}s.",
                iso.len(),
                index_start.elapsed().as_secs_f64()
            );

            if let Some(ref prefix) = config.dump_index {
                let iso_path = format!("{}.intervals.json.gz", prefix.display());
                let jun_path = format!("{}.refskip.json.gz", prefix.display());
                println!(
                    "Saving generated indices to {:?} and {:?}",
                    iso_path, jun_path
                );
                dump_isoform_json(&iso, &iso_path)?;
                dump_isoform_json(&jun, &jun_path)?;
            }

            (Some(iso), Some(jun))
        }
    } else {
        println!("Skipping isoform calling (--skip-iso).");
        (None, None)
    };

    // Group genes by chromosome tid to guarantee coordinate-sorted output
    let mut genes_by_tid: HashMap<u32, Vec<GeneInfo>> = HashMap::new();
    for gene in genes {
        genes_by_tid.entry(gene.tid).or_default().push(gene);
    }
    for genes_of_chr in genes_by_tid.values_mut() {
        genes_of_chr.sort_unstable_by(|a, b| {
            (a.start, a.end, a.gene_id.as_str()).cmp(&(b.start, b.end, b.gene_id.as_str()))
        });
    }

    let mut tids: Vec<u32> = genes_by_tid.keys().copied().collect();
    tids.sort_unstable();

    let num_threads = config.threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .context("Failed to build Rayon thread pool")?;

    // Write to a sibling temporary BAM. The final BAM is published only after
    // all workers finish and the index has been built successfully.
    let staging_bam = temporary_sibling(&config.output_bam, "bam");
    let staging_bai = PathBuf::from(format!("{}.bai", staging_bam.display()));
    // Keep only a small number of completed clusters in flight.  A single
    // high-depth cluster can be large, so a fixed queue of 64 items could
    // otherwise turn writer backpressure into a large memory spike.
    let channel_capacity = num_threads.saturating_mul(2).clamp(2, 16);
    let (tx, rx) = sync_channel::<GeneBatch>(channel_capacity);

    let staging_bam_for_writer = staging_bam.clone();
    let header_clone = header.clone();
    let collect_quantification = config.matrix_out.is_some()
        || config.counts_tsv.is_some()
        || config.molecules_tsv.is_some();
    let collect_matrix = config.matrix_out.is_some();
    let collect_counts = config.counts_tsv.is_some();
    let record_molecules_detail = config.molecules_tsv.is_some();

    let writer_thread = thread::spawn(move || -> Result<(usize, usize, IsoformQuantifier)> {
        let mut writer =
            bam::Writer::from_path(&staging_bam_for_writer, &header_clone, bam::Format::Bam)
                .with_context(|| {
                    format!(
                        "Creating temporary output BAM: {:?}",
                        staging_bam_for_writer
                    )
                })?;
        if num_threads > 1 {
            writer
                .set_threads(num_threads.min(4))
                .context("Enabling multi-threaded BAM output compression")?;
        }

        let mut written = 0;
        let mut batches = 0;
        let mut quantifier = IsoformQuantifier::with_options(
            collect_matrix,
            collect_counts,
            record_molecules_detail,
        );

        while let Ok((records, summaries)) = rx.recv() {
            for rec in records {
                writer.write(&rec)?;
                written += 1;
            }
            for sum in summaries {
                quantifier.add_molecule(sum, record_molecules_detail);
            }
            batches += 1;
        }

        Ok((written, batches, quantifier))
    });

    let tables = ConsensusTables::new();
    let config_ref = &config;
    let cell_set_ref = cell_set.as_ref();
    let iso_ref = isoform_map.as_ref();
    let jun_ref = refskip_map.as_ref();
    let header_view = bam::HeaderView::from_header(&header);
    let gene_context = GeneProcessingContext {
        config: config_ref,
        cell_set: cell_set_ref,
        tables: &tables,
        isoform_map: iso_ref,
        refskip_map: jun_ref,
        collect_quantification,
    };

    println!(
        "Stitching reads across {} contig(s) using {} thread(s)...",
        tids.len(),
        num_threads
    );

    let processing_result: Result<()> = (|| {
        for tid in tids {
            let genes_of_chr = &genes_by_tid[&tid];
            let chr_name = String::from_utf8_lossy(header_view.tid2name(tid)).to_string();
            // A stitched molecule can start before the GTF interval used to
            // fetch its gene (for example, when an alignment spans a long
            // intron).  Therefore, sorting inside a gene cluster is not
            // sufficient to guarantee coordinate order.  Retain one contig
            // at a time, sort all emitted molecules by BAM position, and then
            // hand the contig to the writer.  This bounds the extra memory by
            // the largest contig rather than the whole BAM.
            let mut contig_records = Vec::new();
            let mut contig_summaries = Vec::new();

            let cluster_ranges = gene_cluster_ranges(genes_of_chr);
            // Most genes form singleton clusters.  Processing a bounded batch
            // of clusters in parallel exposes useful parallelism.  Results
            // are accumulated and sorted per contig below.
            let cluster_batch_size = num_threads.saturating_mul(2).clamp(1, 8);
            for batch_start in (0..cluster_ranges.len()).step_by(cluster_batch_size) {
                let batch_end =
                    std::cmp::min(batch_start + cluster_batch_size, cluster_ranges.len());
                // Collect the Result wrapper inside Vec rather than collecting
                // Result<Vec<_>>.  Rayon may use an unordered `while_some`
                // path for the latter, which would allow clusters to reach
                // the writer out of genomic order.
                let gene_results: Vec<Result<ClusterResult>> = pool.install(|| {
                    cluster_ranges[batch_start..batch_end]
                        .par_iter()
                        .enumerate()
                        .map(|(batch_offset, &(cluster_start, cluster_end))| {
                            let cluster_index = batch_start + batch_offset;
                            let gene_cluster = &genes_of_chr[cluster_start..cluster_end];
                            let batches: Result<Vec<GeneBatch>> = gene_cluster
                                .par_iter()
                                .map_init(
                                    || bam::IndexedReader::from_path(&config_ref.input_bam),
                                    |reader_res, gene| {
                                        let reader = reader_res.as_mut().map_err(|e| {
                                            anyhow::anyhow!(
                                                "Failed to open BAM in worker thread: {}",
                                                e
                                            )
                                        })?;
                                        process_single_gene(reader, gene, &gene_context)
                                    },
                                )
                                .collect();
                            let batches = batches.with_context(|| {
                                format!(
                                    "Error processing contig {} gene cluster {}",
                                    chr_name, cluster_index
                                )
                            })?;

                            let mut cluster_records = Vec::new();
                            let mut cluster_summaries = Vec::new();
                            for (recs, sums) in batches {
                                cluster_records.extend(recs);
                                cluster_summaries.extend(sums);
                            }

                            Ok((cluster_index, cluster_records, cluster_summaries))
                        })
                        .collect()
                });

                for result in gene_results {
                    let (_cluster_index, cluster_records, cluster_summaries) = result
                        .with_context(|| {
                            format!(
                                "Error processing contig {} cluster batch starting at {}",
                                chr_name, batch_start
                            )
                        })?;
                    contig_records.extend(cluster_records);
                    contig_summaries.extend(cluster_summaries);
                }
            }

            contig_records.sort_unstable_by(|a, b| {
                a.pos().cmp(&b.pos()).then_with(|| a.qname().cmp(b.qname()))
            });
            if !contig_records.is_empty() {
                tx.send((contig_records, contig_summaries))
                    .map_err(|_| anyhow::anyhow!("Writer thread disconnected"))?;
            }
        }
        Ok(())
    })();

    // Drop sender to signal writer thread to finish
    drop(tx);

    let writer_result = writer_thread
        .join()
        .map_err(|_| anyhow::anyhow!("Writer thread panicked"))?;

    let (written, quantifier) = match (processing_result, writer_result) {
        (Ok(()), Ok((written, _, quantifier))) => (written, quantifier),
        (Err(process_error), Ok(_)) => {
            remove_path_if_present(&staging_bam);
            remove_path_if_present(&staging_bai);
            return Err(process_error);
        }
        (Ok(()), Err(writer_error)) => {
            remove_path_if_present(&staging_bam);
            remove_path_if_present(&staging_bai);
            return Err(writer_error);
        }
        (Err(process_error), Err(_writer_error)) => {
            remove_path_if_present(&staging_bam);
            remove_path_if_present(&staging_bai);
            return Err(process_error);
        }
    };

    // Build the index before publishing the staged BAM.
    if written > 0 {
        println!("Building BAM index (.bai) for {:?}...", staging_bam);
        if let Err(error) = bam::index::build(
            &staging_bam,
            Some(&staging_bai),
            bam::index::Type::Bai,
            num_threads as u32,
        ) {
            remove_path_if_present(&staging_bam);
            remove_path_if_present(&staging_bai);
            return Err(error).with_context(|| {
                format!("Building BAI index for output BAM {:?}", config.output_bam)
            });
        }
    }

    if let Err(error) = fs::rename(&staging_bam, &config.output_bam) {
        remove_path_if_present(&staging_bam);
        remove_path_if_present(&staging_bai);
        return Err(error)
            .with_context(|| format!("Publishing output BAM {:?}", config.output_bam));
    }

    let final_bai = PathBuf::from(format!("{}.bai", config.output_bam.display()));
    let final_csi = PathBuf::from(format!("{}.csi", config.output_bam.display()));
    if written > 0 {
        if let Err(error) = fs::rename(&staging_bai, &final_bai) {
            remove_path_if_present(&staging_bai);
            return Err(error).with_context(|| format!("Publishing BAM index {:?}", final_bai));
        }
        remove_path_if_present(&final_csi);
        println!("Successfully created index {:?}", final_bai);
    } else {
        // Prevent a stale index from being paired with a newly written empty BAM.
        remove_path_if_present(&staging_bai);
        remove_path_if_present(&final_bai);
        remove_path_if_present(&final_csi);
    }

    // Export quantification matrices if requested
    if let Some(ref mex_dir) = config.matrix_out {
        quantifier.export_mex(mex_dir)?;
    }
    if let Some(ref counts_path) = config.counts_tsv {
        quantifier.export_counts_tsv(counts_path)?;
    }
    if let Some(ref mol_path) = config.molecules_tsv {
        quantifier.export_molecules_tsv(mol_path)?;
    }

    let elapsed = start_time.elapsed();
    println!(
        "Finished writing {} stitched molecules to {:?}, took {:.2}s.",
        written,
        config.output_bam,
        elapsed.as_secs_f64()
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gene_clusters_are_non_overlapping_and_ordered() {
        let genes = vec![
            GeneInfo {
                gene_id: "G1".to_string(),
                seqid: "chr1".to_string(),
                start: 10,
                end: 20,
                tid: 0,
            },
            GeneInfo {
                gene_id: "G2".to_string(),
                seqid: "chr1".to_string(),
                start: 21,
                end: 30,
                tid: 0,
            },
            GeneInfo {
                gene_id: "G3".to_string(),
                seqid: "chr1".to_string(),
                start: 100,
                end: 110,
                tid: 0,
            },
        ];

        assert_eq!(gene_cluster_ranges(&genes), vec![(0, 2), (2, 3)]);
        assert!(gene_cluster_ranges(&[]).is_empty());
    }

    #[test]
    fn test_prepare_bam_output_rejects_same_path() {
        let path =
            std::env::temp_dir().join(format!("mfs_stitcher_same_bam_{}.bam", std::process::id()));
        let _ = fs::remove_file(&path);
        File::create(&path).unwrap();

        let error = prepare_bam_output(&path, &path).unwrap_err();
        assert!(error.to_string().contains("must be different"));

        let _ = fs::remove_file(path);
    }
}
