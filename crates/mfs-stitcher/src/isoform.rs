//! Isoform annotation and compatible transcript assignment.

use anyhow::{bail, Context, Result};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read};
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::interval::{Interval, IntervalSet};

const INDEX_FORMAT_VERSION: u32 = 1;
const BAM_COORDINATE_SYSTEM: &str = "bam-0based-inclusive";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexKind {
    Exon,
    RefSkip,
}

#[derive(Serialize)]
struct IndexDocument<'a> {
    format_version: u32,
    coordinate_system: &'static str,
    entries: &'a BTreeMap<String, BTreeMap<String, String>>,
}

fn temporary_json_path(path: &Path) -> std::path::PathBuf {
    let name = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("index.json.gz");
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|value| value.as_nanos())
        .unwrap_or_default();
    path.with_file_name(format!(".{}.tmp-{}-{}", name, std::process::id(), nonce))
}

/// Parsed interval set with its associated transcripts.
#[derive(Debug, Clone)]
pub struct AnnotatedIntervalEntry {
    pub intervals: IntervalSet,
    pub transcripts: HashSet<String>,
}

pub type GeneIsoformMap = HashMap<String, Vec<AnnotatedIntervalEntry>>;

/// Parse portion string like "[11869,12009] | [12058,12178]" into IntervalSet.
pub fn parse_portion_interval_string(s: &str) -> IntervalSet {
    let mut ivs = Vec::new();
    for part in s.split('|') {
        let part = part.trim();
        let stripped = if part.starts_with('[') && part.ends_with(']') {
            &part[1..part.len() - 1]
        } else {
            continue;
        };
        let mut nums = stripped.split(',');
        if let (Some(s_str), Some(e_str)) = (nums.next(), nums.next()) {
            if let (Ok(start), Ok(end)) = (s_str.trim().parse::<i64>(), e_str.trim().parse::<i64>())
            {
                if start <= end {
                    ivs.push(Interval::new(start, end));
                }
            }
        }
    }
    IntervalSet::from_intervals(ivs)
}

/// Format IntervalSet to portion string representation: "[11869,12009] | [12058,12178]"
pub fn format_interval_set_to_string(iv_set: &IntervalSet) -> String {
    let parts: Vec<String> = iv_set
        .intervals
        .iter()
        .map(|iv| format!("[{},{}]", iv.start, iv.end))
        .collect();
    parts.join(" | ")
}

/// Load an isoform JSON using exon-style legacy coordinate conversion.
///
/// New files written by `dump_isoform_json` carry an explicit BAM coordinate
/// metadata envelope. Legacy files from the Python stitcher are raw GTF-style
/// maps, so this compatibility wrapper treats them as exon intervals.
pub fn load_isoform_json<P: AsRef<Path>>(path: P) -> Result<GeneIsoformMap> {
    load_isoform_json_with_kind(path, IndexKind::Exon)
}

/// Load an isoform or refskip JSON (plain or .gz).
///
/// Legacy Python exon intervals are converted from 1-based inclusive GTF
/// coordinates. Legacy refskip intervals were written from exon-end to the
/// next-exon-start; they are converted to the corresponding BAM N interval.
pub fn load_isoform_json_with_kind<P: AsRef<Path>>(
    path: P,
    kind: IndexKind,
) -> Result<GeneIsoformMap> {
    let path_ref = path.as_ref();
    let file = File::open(path_ref).with_context(|| format!("Failed to open {:?}", path_ref))?;
    let mut reader: Box<dyn Read> = if path_ref.to_string_lossy().ends_with(".gz") {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(BufReader::new(file))
    };

    let document: serde_json::Value = serde_json::from_reader(&mut reader)
        .with_context(|| format!("Failed to parse JSON from {:?}", path_ref))?;

    let (raw, is_legacy) = if document.get("format_version").is_some() {
        let version = document
            .get("format_version")
            .and_then(serde_json::Value::as_u64)
            .unwrap_or_default();
        let coordinate_system = document
            .get("coordinate_system")
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        if version != INDEX_FORMAT_VERSION as u64 || coordinate_system != BAM_COORDINATE_SYSTEM {
            bail!(
                "Unsupported isoform index format in {:?}: version={}, coordinates={:?}",
                path_ref,
                version,
                coordinate_system
            );
        }
        let entries = document
            .get("entries")
            .cloned()
            .context("Isoform index metadata is missing entries")?;
        let raw: HashMap<String, HashMap<String, String>> = serde_json::from_value(entries)
            .with_context(|| format!("Invalid entries in isoform index {:?}", path_ref))?;
        (raw, false)
    } else {
        let raw: HashMap<String, HashMap<String, String>> = serde_json::from_value(document)
            .with_context(|| format!("Invalid legacy isoform index {:?}", path_ref))?;
        (raw, true)
    };

    let mut result: GeneIsoformMap = HashMap::with_capacity(raw.len());

    for (gene_id, entries) in raw {
        let mut entry_list = Vec::with_capacity(entries.len());
        for (interval_str, tx_str) in entries {
            let iv_set = parse_portion_interval_string(&interval_str);
            let iv_set = if is_legacy {
                convert_legacy_intervals(&iv_set, kind)
            } else {
                iv_set
            };
            if iv_set.is_empty() {
                continue;
            }
            let transcripts: HashSet<String> = tx_str
                .split(',')
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect();
            entry_list.push(AnnotatedIntervalEntry {
                intervals: iv_set,
                transcripts,
            });
        }
        result.insert(gene_id, entry_list);
    }

    Ok(result)
}

fn convert_legacy_intervals(intervals: &IntervalSet, kind: IndexKind) -> IntervalSet {
    let converted = intervals
        .intervals
        .iter()
        .filter_map(|interval| {
            let (start, end) = match kind {
                IndexKind::Exon => (interval.start - 1, interval.end - 1),
                IndexKind::RefSkip => (interval.start, interval.end - 2),
            };
            (start <= end).then(|| Interval::new(start, end))
        })
        .collect();
    IntervalSet::from_intervals(converted)
}

/// Save an in-memory GeneIsoformMap to compressed JSON (.json.gz).
pub fn dump_isoform_json<P: AsRef<Path>>(map: &GeneIsoformMap, path: P) -> Result<()> {
    let path_ref = path.as_ref();
    if let Some(parent) = path_ref.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("Creating index directory: {:?}", parent))?;
        }
    }

    let mut raw: BTreeMap<String, BTreeMap<String, String>> = BTreeMap::new();
    for (gene_id, entries) in map {
        let mut sub_map = BTreeMap::new();
        for entry in entries {
            let iv_str = format_interval_set_to_string(&entry.intervals);
            let mut tx_list: Vec<String> = entry.transcripts.iter().cloned().collect();
            tx_list.sort();
            sub_map.insert(iv_str, tx_list.join(","));
        }
        raw.insert(gene_id.clone(), sub_map);
    }

    let document = IndexDocument {
        format_version: INDEX_FORMAT_VERSION,
        coordinate_system: BAM_COORDINATE_SYSTEM,
        entries: &raw,
    };

    let temporary = temporary_json_path(path_ref);
    let result = (|| {
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
            .with_context(|| format!("Creating temporary index: {:?}", temporary))?;
        let mut gz = GzEncoder::new(BufWriter::new(file), Compression::default());
        serde_json::to_writer(&mut gz, &document)
            .with_context(|| format!("Writing isoform index: {:?}", path_ref))?;
        gz.finish()
            .with_context(|| format!("Finishing isoform index: {:?}", path_ref))?;
        fs::rename(&temporary, path_ref)
            .with_context(|| format!("Replacing isoform index: {:?}", path_ref))?;
        Ok(())
    })();

    if result.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    result
}

/// Sweep-line interval discretization algorithm (matching gtf_to_json.py linear time logic)
pub fn create_interval_entries_sweep_line(
    isoform_dict: &HashMap<String, IntervalSet>,
) -> Vec<AnnotatedIntervalEntry> {
    #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
    struct Event<'a> {
        pos: i64,
        delta: i32, // +1 for start, -1 for end + 1
        transcript: &'a str,
    }

    let mut events = Vec::new();
    for (transcript, iv_set) in isoform_dict {
        for iv in &iv_set.intervals {
            events.push(Event {
                pos: iv.start,
                delta: 1,
                transcript: transcript.as_str(),
            });
            events.push(Event {
                pos: iv.end + 1,
                delta: -1,
                transcript: transcript.as_str(),
            });
        }
    }

    if events.is_empty() {
        return Vec::new();
    }

    // Sort by pos
    events.sort_by_key(|e| e.pos);

    let mut active: BTreeSet<&str> = BTreeSet::new();
    let mut last_pos = events[0].pos;
    let mut segs_by_set: HashMap<Vec<String>, Vec<Interval>> = HashMap::new();

    let n_events = events.len();
    let mut i = 0;
    while i < n_events {
        let pos = events[i].pos;
        if pos > last_pos && !active.is_empty() {
            let key: Vec<String> = active.iter().map(|s| s.to_string()).collect();
            segs_by_set
                .entry(key)
                .or_default()
                .push(Interval::new(last_pos, pos - 1));
        }

        while i < n_events && events[i].pos == pos {
            if events[i].delta > 0 {
                active.insert(events[i].transcript);
            } else {
                active.remove(events[i].transcript);
            }
            i += 1;
        }
        last_pos = pos;
    }

    let mut entries = Vec::with_capacity(segs_by_set.len());
    for (tr_list, segs) in segs_by_set {
        let merged_set = IntervalSet::from_intervals(segs);
        let transcripts: HashSet<String> = tr_list.into_iter().collect();
        entries.push(AnnotatedIntervalEntry {
            intervals: merged_set,
            transcripts,
        });
    }

    entries
}

/// Helper to extract attribute value from GTF line (strictly matching key followed by delimiter)
pub fn extract_attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    for part in attrs.split(';') {
        let part = part.trim();
        if let Some(rest) = part.strip_prefix(key) {
            if rest.starts_with(|c: char| c.is_whitespace() || c == '=') {
                let val = rest.trim_start_matches(|c: char| c.is_whitespace() || c == '=');
                let val = val.trim().trim_matches('"');
                if !val.is_empty() {
                    return Some(val);
                }
            }
        }
    }
    None
}

/// Build both isoform exonic intervals and refskip junction indices directly from a GTF file in memory.
/// Converts GTF 1-based coordinates to 0-based coordinates to match BAM reference coordinates.
pub fn build_isoform_indices_from_gtf(
    gtf_path: &Path,
    contig_filter: Option<&str>,
    gene_set: Option<&HashSet<String>>,
    gene_identifier: &str,
) -> Result<(GeneIsoformMap, GeneIsoformMap)> {
    let file = File::open(gtf_path).with_context(|| format!("Opening GTF: {:?}", gtf_path))?;
    let reader = BufReader::new(file);

    // Group exons by gene_id -> transcript_id -> list of (start, end)
    let mut gene_tx_exons: HashMap<String, HashMap<String, Vec<Interval>>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        if fields[2] != "exon" {
            continue;
        }

        let seqid = fields[0];
        if let Some(c) = contig_filter {
            if seqid != c {
                continue;
            }
        }

        let start_1based: i64 = match fields[3].parse() {
            Ok(s) => s,
            Err(_) => continue,
        };
        let end_1based: i64 = match fields[4].parse() {
            Ok(e) => e,
            Err(_) => continue,
        };
        if start_1based > end_1based {
            continue;
        }

        // Convert GTF 1-based inclusive coordinates to BAM 0-based inclusive coordinates
        let start_0based = (start_1based - 1).max(0);
        let end_0based = (end_1based - 1).max(0);

        let attrs = fields[8];
        let gene_id = match extract_attr(attrs, gene_identifier) {
            Some(g) => g.to_string(),
            None => match extract_attr(attrs, "gene_id") {
                Some(g) => g.to_string(),
                None => continue,
            },
        };

        if let Some(set) = gene_set {
            if !set.contains(&gene_id) {
                continue;
            }
        }

        let transcript_id = match extract_attr(attrs, "transcript_id") {
            Some(t) => t.to_string(),
            None => continue,
        };

        gene_tx_exons
            .entry(gene_id)
            .or_default()
            .entry(transcript_id)
            .or_default()
            .push(Interval::new(start_0based, end_0based));
    }

    // Parallel process each gene using Rayon
    let genes: Vec<(String, HashMap<String, Vec<Interval>>)> = gene_tx_exons.into_iter().collect();

    let computed: Vec<(
        String,
        Vec<AnnotatedIntervalEntry>,
        Vec<AnnotatedIntervalEntry>,
    )> = genes
        .into_par_iter()
        .map(|(gene_id, tx_exons)| {
            let mut isoform_dict: HashMap<String, IntervalSet> =
                HashMap::with_capacity(tx_exons.len());
            let mut refskip_dict: HashMap<String, IntervalSet> =
                HashMap::with_capacity(tx_exons.len());

            for (tx_id, mut exons) in tx_exons {
                exons.sort_by_key(|e| e.start);
                let exon_set = IntervalSet::from_intervals(exons.clone());
                isoform_dict.insert(tx_id.clone(), exon_set);

                let mut introns = Vec::new();
                if exons.len() > 1 {
                    for i in 0..exons.len() - 1 {
                        let int_start = exons[i].end + 1;
                        let int_end = exons[i + 1].start - 1;
                        if int_start <= int_end {
                            introns.push(Interval::new(int_start, int_end));
                        }
                    }
                }
                let refskip_set = IntervalSet::from_intervals(introns);
                refskip_dict.insert(tx_id, refskip_set);
            }

            let iso_entries = create_interval_entries_sweep_line(&isoform_dict);
            let refskip_entries = create_interval_entries_sweep_line(&refskip_dict);

            (gene_id, iso_entries, refskip_entries)
        })
        .collect();

    let mut isoform_map: GeneIsoformMap = HashMap::with_capacity(computed.len());
    let mut refskip_map: GeneIsoformMap = HashMap::with_capacity(computed.len());

    for (gene_id, iso_entries, refskip_entries) in computed {
        if !iso_entries.is_empty() {
            isoform_map.insert(gene_id.clone(), iso_entries);
        }
        if !refskip_entries.is_empty() {
            refskip_map.insert(gene_id, refskip_entries);
        }
    }

    Ok((isoform_map, refskip_map))
}

/// Find compatible transcripts for a stitched molecule.
pub fn find_compatible_transcripts(
    gene_id: &str,
    ref_intervals: &IntervalSet,
    skipped_intervals: &IntervalSet,
    isoform_map: Option<&GeneIsoformMap>,
    refskip_map: Option<&GeneIsoformMap>,
) -> Option<String> {
    let mut tx_candidates: Option<HashSet<String>> = None;

    if let Some(iso_map) = isoform_map {
        if let Some(entries) = iso_map.get(gene_id) {
            let mut matching_sets: Vec<&HashSet<String>> = Vec::new();
            for entry in entries {
                let inter = entry.intervals.intersection(ref_intervals);
                if inter.total_len() > 4 {
                    matching_sets.push(&entry.transcripts);
                }
            }

            if !matching_sets.is_empty() {
                let mut current: HashSet<String> = matching_sets[0].clone();
                for next_set in &matching_sets[1..] {
                    current.retain(|x| next_set.contains(x));
                }
                tx_candidates = Some(current);
            }
        }
    }

    if let Some(rs_map) = refskip_map {
        if let Some(entries) = rs_map.get(gene_id) {
            let mut matching_sets: Vec<&HashSet<String>> = Vec::new();
            for entry in entries {
                let inter = entry.intervals.intersection(skipped_intervals);
                if inter.total_len() > 4 {
                    matching_sets.push(&entry.transcripts);
                }
            }

            if !matching_sets.is_empty() {
                let mut current: HashSet<String> = matching_sets[0].clone();
                for next_set in &matching_sets[1..] {
                    current.retain(|x| next_set.contains(x));
                }

                if let Some(ref mut existing) = tx_candidates {
                    existing.retain(|x| current.contains(x));
                } else {
                    tx_candidates = Some(current);
                }
            }
        }
    }

    let candidates = tx_candidates?;
    if candidates.is_empty() {
        return None;
    }

    let mut sorted: Vec<String> = candidates.into_iter().collect();
    sorted.sort();
    Some(sorted.join(","))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_portion_string() {
        let s = "[11869,12009] | [12058,12178]";
        let set = parse_portion_interval_string(s);
        assert_eq!(set.intervals.len(), 2);
        assert_eq!(set.intervals[0], Interval::new(11869, 12009));
        assert_eq!(set.intervals[1], Interval::new(12058, 12178));
        assert_eq!(format_interval_set_to_string(&set), s);
    }

    #[test]
    fn test_sweep_line_discretization() {
        // T1: 100..200
        // T2: 150..250
        let mut d = HashMap::new();
        d.insert(
            "T1".to_string(),
            IntervalSet::from_intervals(vec![Interval::new(100, 200)]),
        );
        d.insert(
            "T2".to_string(),
            IntervalSet::from_intervals(vec![Interval::new(150, 250)]),
        );

        let entries = create_interval_entries_sweep_line(&d);
        assert_eq!(entries.len(), 3);

        let t1_only = entries
            .iter()
            .find(|e| e.transcripts.len() == 1 && e.transcripts.contains("T1"))
            .unwrap();
        assert_eq!(t1_only.intervals.intervals, vec![Interval::new(100, 149)]);

        let both = entries.iter().find(|e| e.transcripts.len() == 2).unwrap();
        assert_eq!(both.intervals.intervals, vec![Interval::new(150, 200)]);

        let t2_only = entries
            .iter()
            .find(|e| e.transcripts.len() == 1 && e.transcripts.contains("T2"))
            .unwrap();
        assert_eq!(t2_only.intervals.intervals, vec![Interval::new(201, 250)]);
    }

    #[test]
    fn test_index_roundtrip_preserves_bam_coordinates() {
        let path = std::env::temp_dir().join(format!(
            "mfs_stitcher_index_roundtrip_{}.json.gz",
            std::process::id()
        ));
        let _ = fs::remove_file(&path);

        let mut map = GeneIsoformMap::new();
        let mut transcripts = HashSet::new();
        transcripts.insert("TX1".to_string());
        map.insert(
            "GENE1".to_string(),
            vec![AnnotatedIntervalEntry {
                intervals: IntervalSet::from_intervals(vec![Interval::new(99, 198)]),
                transcripts,
            }],
        );

        dump_isoform_json(&map, &path).unwrap();
        let loaded = load_isoform_json_with_kind(&path, IndexKind::Exon).unwrap();
        assert_eq!(
            loaded["GENE1"][0].intervals.intervals,
            vec![Interval::new(99, 198)]
        );
        assert!(loaded["GENE1"][0].transcripts.contains("TX1"));

        let _ = fs::remove_file(path);
    }

    #[test]
    fn test_legacy_index_coordinates_are_converted_by_kind() {
        let exon_path = std::env::temp_dir().join(format!(
            "mfs_stitcher_legacy_exon_{}.json",
            std::process::id()
        ));
        let refskip_path = std::env::temp_dir().join(format!(
            "mfs_stitcher_legacy_refskip_{}.json",
            std::process::id()
        ));
        let _ = fs::remove_file(&exon_path);
        let _ = fs::remove_file(&refskip_path);

        let raw_exon = serde_json::json!({"GENE1": {"[100,199]": "TX1"}});
        serde_json::to_writer(File::create(&exon_path).unwrap(), &raw_exon).unwrap();
        let exon = load_isoform_json_with_kind(&exon_path, IndexKind::Exon).unwrap();
        assert_eq!(
            exon["GENE1"][0].intervals.intervals,
            vec![Interval::new(99, 198)]
        );

        let raw_refskip = serde_json::json!({"GENE1": {"[199,300]": "TX1"}});
        serde_json::to_writer(File::create(&refskip_path).unwrap(), &raw_refskip).unwrap();
        let refskip = load_isoform_json_with_kind(&refskip_path, IndexKind::RefSkip).unwrap();
        assert_eq!(
            refskip["GENE1"][0].intervals.intervals,
            vec![Interval::new(199, 298)]
        );

        let _ = fs::remove_file(exon_path);
        let _ = fs::remove_file(refskip_path);
    }
}
