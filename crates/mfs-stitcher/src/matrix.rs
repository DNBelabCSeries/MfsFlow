//! Direct single-cell transcript isoform quantification and matrix export (MEX & TSV).

use std::collections::{BTreeSet, HashMap};
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;

fn temporary_output_path(path: &Path) -> PathBuf {
    let name = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("output");
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|value| value.as_nanos())
        .unwrap_or_default();
    path.with_file_name(format!(".{}.tmp-{}-{}", name, std::process::id(), nonce))
}

fn write_gzip_atomically<F>(path: &Path, write_contents: F) -> Result<()>
where
    F: FnOnce(&mut GzEncoder<BufWriter<File>>) -> std::io::Result<()>,
{
    let temporary = temporary_output_path(path);
    let result = (|| {
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
            .with_context(|| format!("Creating temporary output: {:?}", temporary))?;
        let mut encoder = GzEncoder::new(BufWriter::new(file), Compression::default());
        write_contents(&mut encoder).with_context(|| format!("Writing {:?}", path))?;
        encoder
            .finish()
            .with_context(|| format!("Finishing gzip output: {:?}", path))?;
        fs::rename(&temporary, path)
            .with_context(|| format!("Replacing output atomically: {:?}", path))?;
        Ok(())
    })();

    if result.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    result
}

/// Molecule summary information used for isoform quantification
#[derive(Debug, Clone)]
pub struct MoleculeSummary {
    pub cell: String,
    pub gene: String,
    pub umi: String,
    pub transcripts: Option<Vec<String>>,
    pub reads: i32,
}

/// Accumulator for single-cell isoform expression counts
#[derive(Default)]
pub struct IsoformQuantifier {
    // (transcript_id, cell_barcode) -> unique molecule count
    pub unique_counts: HashMap<(String, String), u32>,
    // (transcript_id, cell_barcode) -> fractional molecule count
    pub fractional_counts: HashMap<(String, String), f64>,
    // (transcript_id, cell_barcode) -> total reads
    pub read_counts: HashMap<(String, String), i32>,
    // transcript_id -> gene_id
    pub transcript_to_gene: HashMap<String, String>,
    // All observed cells
    pub cells: BTreeSet<String>,
    // All observed transcripts
    pub transcripts: BTreeSet<String>,
    // Stored molecules (optional for molecules.tsv)
    pub molecules: Vec<MoleculeSummary>,
    collect_matrix: bool,
    collect_counts: bool,
    record_molecules: bool,
}

impl IsoformQuantifier {
    pub fn new() -> Self {
        Self::with_options(true, true, true)
    }

    /// Create a quantifier with only the accumulators required by requested
    /// outputs.  Keeping the legacy `new()` behavior makes the public helper
    /// convenient in tests and downstream callers, while the pipeline can
    /// avoid large maps when it only needs a stitched BAM.
    pub fn with_options(
        collect_matrix: bool,
        collect_counts: bool,
        record_molecules: bool,
    ) -> Self {
        Self {
            unique_counts: HashMap::new(),
            fractional_counts: HashMap::new(),
            read_counts: HashMap::new(),
            transcript_to_gene: HashMap::new(),
            cells: BTreeSet::new(),
            transcripts: BTreeSet::new(),
            molecules: Vec::new(),
            collect_matrix,
            collect_counts,
            record_molecules,
        }
    }

    pub fn add_molecule(&mut self, mol: MoleculeSummary, record_molecule_list: bool) {
        if self.collect_matrix {
            self.cells.insert(mol.cell.clone());
        }

        if let Some(ref tx_list) = mol.transcripts {
            if !tx_list.is_empty() {
                let k = tx_list.len();
                let is_unique = k == 1;
                let frac = 1.0 / (k as f64);

                for tx in tx_list {
                    if self.collect_matrix {
                        self.transcripts.insert(tx.clone());
                    }
                    if self.collect_matrix || self.collect_counts {
                        self.transcript_to_gene
                            .entry(tx.clone())
                            .or_insert_with(|| mol.gene.clone());
                    }

                    if self.collect_matrix || self.collect_counts {
                        let key = (tx.clone(), mol.cell.clone());
                        if is_unique {
                            *self.unique_counts.entry(key.clone()).or_insert(0) += 1;
                        }
                        if self.collect_counts {
                            *self.fractional_counts.entry(key.clone()).or_insert(0.0) += frac;
                            *self.read_counts.entry(key).or_insert(0) += mol.reads;
                        }
                    }
                }
            }
        }

        if record_molecule_list && self.record_molecules {
            self.molecules.push(mol);
        }
    }

    /// Export a standard 10x-compatible MEX directory.
    ///
    /// Matrix Market uses integer values, so this export intentionally contains
    /// only molecules with a unique transcript assignment. Ambiguous molecules
    /// remain available through `fractional_counts` and `export_counts_tsv`.
    pub fn export_mex<P: AsRef<Path>>(&self, out_dir: P) -> Result<()> {
        let dir = out_dir.as_ref();
        fs::create_dir_all(dir).with_context(|| format!("Creating MEX directory: {:?}", dir))?;

        let barcodes_path = dir.join("barcodes.tsv.gz");
        let features_path = dir.join("features.tsv.gz");
        let matrix_path = dir.join("matrix.mtx.gz");

        // 1. Write barcodes.tsv.gz
        let cell_list: Vec<&String> = self.cells.iter().collect();
        let mut cell_to_idx: HashMap<&String, usize> = HashMap::with_capacity(cell_list.len());
        {
            for (idx, cell) in cell_list.iter().enumerate() {
                cell_to_idx.insert(cell, idx + 1); // 1-based index
            }
            write_gzip_atomically(&barcodes_path, |gz| {
                for cell in &cell_list {
                    writeln!(gz, "{}", cell)?;
                }
                Ok(())
            })?;
        }

        // 2. Write features.tsv.gz
        let tx_list: Vec<&String> = self.transcripts.iter().collect();
        let mut tx_to_idx: HashMap<&String, usize> = HashMap::with_capacity(tx_list.len());
        {
            for (idx, tx) in tx_list.iter().enumerate() {
                tx_to_idx.insert(tx, idx + 1); // 1-based index
            }
            write_gzip_atomically(&features_path, |gz| {
                for tx in &tx_list {
                    let gene = self
                        .transcript_to_gene
                        .get(*tx)
                        .map(|s| s.as_str())
                        .unwrap_or("Unknown");
                    // 10x format: <feature_id>\t<feature_name>\t<feature_type>
                    writeln!(gz, "{}\t{}\tGene Expression", tx, gene)?;
                }
                Ok(())
            })?;
        }

        // 3. Write matrix.mtx.gz (using unique molecule counts)
        {
            // Collect and sort coordinate entries by (feature_idx, barcode_idx)
            let mut entries: Vec<(usize, usize, u32)> =
                Vec::with_capacity(self.unique_counts.len());
            for ((tx, cell), count) in &self.unique_counts {
                if *count > 0 {
                    if let (Some(&f_idx), Some(&b_idx)) = (tx_to_idx.get(tx), cell_to_idx.get(cell))
                    {
                        entries.push((f_idx, b_idx, *count));
                    }
                }
            }
            entries.sort_unstable_by_key(|e| (e.0, e.1));

            write_gzip_atomically(&matrix_path, |gz| {
                writeln!(gz, "%%MatrixMarket matrix coordinate integer general")?;
                writeln!(
                    gz,
                    "%metadata_json: {{\"software_version\": \"mfs-stitcher-{}\"}}",
                    env!("CARGO_PKG_VERSION")
                )?;
                writeln!(
                    gz,
                    "{} {} {}",
                    tx_list.len(),
                    cell_list.len(),
                    entries.len()
                )?;

                for (f_idx, b_idx, count) in entries {
                    writeln!(gz, "{} {} {}", f_idx, b_idx, count)?;
                }
                Ok(())
            })?;
        }

        println!(
            "Exported 10x-compatible MEX matrix to {:?} ({} transcripts x {} cells)",
            dir,
            tx_list.len(),
            cell_list.len()
        );

        Ok(())
    }

    /// Export isoform counts summary TSV (gzipped)
    pub fn export_counts_tsv<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path_ref = path.as_ref();
        if let Some(parent) = path_ref.parent() {
            fs::create_dir_all(parent)?;
        }

        // Sort borrowed keys instead of copying every row into a BTreeMap.
        // The values remain in the existing hash maps, which keeps the peak
        // memory proportional to the number of output rows with less node
        // allocation overhead.
        let mut rows: Vec<(&(String, String), &f64)> = self.fractional_counts.iter().collect();
        rows.sort_unstable_by(|(left, _), (right, _)| {
            left.1.cmp(&right.1).then_with(|| left.0.cmp(&right.0))
        });

        write_gzip_atomically(path_ref, |gz| {
            writeln!(
                gz,
                "cell\tgene_id\ttranscript_id\tunique_counts\tfractional_counts\ttotal_reads"
            )?;

            for (key, frac) in rows {
                let tx = &key.0;
                let cell = &key.1;
                let uniq = self.unique_counts.get(key).copied().unwrap_or(0);
                let reads = self.read_counts.get(key).copied().unwrap_or(0);
                let gene = self
                    .transcript_to_gene
                    .get(tx)
                    .map(|s| s.as_str())
                    .unwrap_or("Unknown");
                writeln!(
                    gz,
                    "{}\t{}\t{}\t{}\t{:.6}\t{}",
                    cell, gene, tx, uniq, frac, reads
                )?;
            }
            Ok(())
        })?;
        println!("Exported isoform counts TSV to {:?}", path_ref);
        Ok(())
    }

    /// Export molecule-level audit detail TSV (gzipped)
    pub fn export_molecules_tsv<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path_ref = path.as_ref();
        if let Some(parent) = path_ref.parent() {
            fs::create_dir_all(parent)?;
        }

        write_gzip_atomically(path_ref, |gz| {
            writeln!(
                gz,
                "cell\tgene_id\tumi\ttranscripts\tn_transcripts\treads\tis_unique"
            )?;

            for mol in &self.molecules {
                let (tx_str, n_tx, is_uniq) = match &mol.transcripts {
                    Some(list) if !list.is_empty() => {
                        let s = list.join(",");
                        let n = list.len();
                        let u = if n == 1 { "1" } else { "0" };
                        (s, n, u)
                    }
                    _ => ("Unassigned".to_string(), 0, "0"),
                };

                writeln!(
                    gz,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    mol.cell, mol.gene, mol.umi, tx_str, n_tx, mol.reads, is_uniq
                )?;
            }
            Ok(())
        })?;
        println!("Exported molecule detail TSV to {:?}", path_ref);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::read::GzDecoder;
    use std::io::Read;

    fn read_gzip(path: &Path) -> String {
        let file = File::open(path).unwrap();
        let mut decoder = GzDecoder::new(file);
        let mut contents = String::new();
        decoder.read_to_string(&mut contents).unwrap();
        contents
    }

    #[test]
    fn test_isoform_quantification() {
        let mut quant = IsoformQuantifier::new();

        // Molecule 1: CellA, GENE1, unique to T1
        quant.add_molecule(
            MoleculeSummary {
                cell: "CellA".to_string(),
                gene: "GENE1".to_string(),
                umi: "UMI1".to_string(),
                transcripts: Some(vec!["T1".to_string()]),
                reads: 10,
            },
            true,
        );

        // Molecule 2: CellA, GENE1, ambiguous between T1 and T2
        quant.add_molecule(
            MoleculeSummary {
                cell: "CellA".to_string(),
                gene: "GENE1".to_string(),
                umi: "UMI2".to_string(),
                transcripts: Some(vec!["T1".to_string(), "T2".to_string()]),
                reads: 5,
            },
            true,
        );

        // CellA, T1 unique count should be 1
        assert_eq!(
            quant
                .unique_counts
                .get(&("T1".to_string(), "CellA".to_string())),
            Some(&1)
        );

        // CellA, T2 unique count should be None/0
        assert_eq!(
            quant
                .unique_counts
                .get(&("T2".to_string(), "CellA".to_string())),
            None
        );

        // Fractional counts: T1 gets 1.0 + 0.5 = 1.5
        let t1_frac = quant
            .fractional_counts
            .get(&("T1".to_string(), "CellA".to_string()))
            .unwrap();
        assert!((t1_frac - 1.5).abs() < 1e-6);

        // T2 gets 0.5
        let t2_frac = quant
            .fractional_counts
            .get(&("T2".to_string(), "CellA".to_string()))
            .unwrap();
        assert!((t2_frac - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_mex_export_contains_only_unique_assignments() {
        let out_dir =
            std::env::temp_dir().join(format!("mfs_stitcher_matrix_test_{}", std::process::id()));
        let _ = fs::remove_dir_all(&out_dir);

        let mut quant = IsoformQuantifier::new();
        quant.add_molecule(
            MoleculeSummary {
                cell: "CellA".to_string(),
                gene: "GENE1".to_string(),
                umi: "UMI1".to_string(),
                transcripts: Some(vec!["T1".to_string()]),
                reads: 3,
            },
            false,
        );
        quant.add_molecule(
            MoleculeSummary {
                cell: "CellA".to_string(),
                gene: "GENE1".to_string(),
                umi: "UMI2".to_string(),
                transcripts: Some(vec!["T1".to_string(), "T2".to_string()]),
                reads: 2,
            },
            false,
        );
        quant.add_molecule(
            MoleculeSummary {
                cell: "CellB".to_string(),
                gene: "GENE1".to_string(),
                umi: "UMI3".to_string(),
                transcripts: Some(vec!["T2".to_string()]),
                reads: 4,
            },
            false,
        );

        quant.export_mex(&out_dir).unwrap();

        let barcodes = read_gzip(&out_dir.join("barcodes.tsv.gz"));
        assert_eq!(barcodes, "CellA\nCellB\n");

        let features = read_gzip(&out_dir.join("features.tsv.gz"));
        assert_eq!(
            features,
            "T1\tGENE1\tGene Expression\nT2\tGENE1\tGene Expression\n"
        );

        let matrix = read_gzip(&out_dir.join("matrix.mtx.gz"));
        assert!(matrix.contains("2 2 2\n"));
        assert!(matrix.contains("1 1 1\n"));
        assert!(matrix.contains("2 2 1\n"));
        assert!(!matrix.contains("2 1 1\n"));

        let _ = fs::remove_dir_all(out_dir);
    }

    #[test]
    fn test_counts_export_preserves_sorted_rows_and_values() {
        let path = std::env::temp_dir().join(format!(
            "mfs_stitcher_counts_test_{}.tsv.gz",
            std::process::id()
        ));
        let _ = fs::remove_file(&path);

        let mut quant = IsoformQuantifier::new();
        quant.add_molecule(
            MoleculeSummary {
                cell: "CellB".to_string(),
                gene: "GENE1".to_string(),
                umi: "UMI1".to_string(),
                transcripts: Some(vec!["T2".to_string()]),
                reads: 4,
            },
            false,
        );
        quant.add_molecule(
            MoleculeSummary {
                cell: "CellA".to_string(),
                gene: "GENE1".to_string(),
                umi: "UMI2".to_string(),
                transcripts: Some(vec!["T2".to_string()]),
                reads: 2,
            },
            false,
        );
        quant.add_molecule(
            MoleculeSummary {
                cell: "CellA".to_string(),
                gene: "GENE1".to_string(),
                umi: "UMI3".to_string(),
                transcripts: Some(vec!["T1".to_string()]),
                reads: 3,
            },
            false,
        );

        quant.export_counts_tsv(&path).unwrap();
        let contents = read_gzip(&path);
        assert_eq!(
            contents,
            "cell\tgene_id\ttranscript_id\tunique_counts\tfractional_counts\ttotal_reads\n"
                .to_string()
                + "CellA\tGENE1\tT1\t1\t1.000000\t3\n"
                + "CellA\tGENE1\tT2\t1\t1.000000\t2\n"
                + "CellB\tGENE1\tT2\t1\t1.000000\t4\n"
        );

        let _ = fs::remove_file(path);
    }
}
