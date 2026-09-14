use anyhow::{bail, Result};
use clap::{Args, Parser};
use std::path::PathBuf;
use std::time::Instant;

use mfs_stitcher::isoform::{build_isoform_indices_from_gtf, dump_isoform_json};
use mfs_stitcher::pipeline::{run_pipeline, PipelineConfig};

#[derive(Parser, Debug)]
#[command(
    name = "mfs-stitcher",
    version = "2.0.0",
    about = "High-performance Smart-seq3 molecule stitching and transcript isoform caller (Rust edition)",
    after_help = "Common run:\n  mfs-stitcher -i input.bam -o stitched.bam -g genes.gtf -t 20\n\nIsoform indexes are built in memory automatically. Use --isoform and --junction only when reusing precomputed indexes."
)]
struct Cli {
    #[command(flatten)]
    core: CoreArgs,

    #[command(flatten)]
    outputs: OutputArgs,

    #[command(flatten)]
    annotation: AnnotationArgs,

    #[command(flatten)]
    runtime: RuntimeArgs,
}

#[derive(Args, Debug)]
#[command(next_help_heading = "Core inputs")]
struct CoreArgs {
    /// Coordinate-sorted and indexed input BAM file
    #[arg(
        short = 'i',
        long = "input",
        visible_alias = "input-bam",
        value_name = "BAM"
    )]
    input: Option<PathBuf>,

    /// GTF file containing gene and transcript annotation
    #[arg(short = 'g', long = "gtf", value_name = "GTF")]
    gtf: PathBuf,
}

#[derive(Args, Debug)]
#[command(next_help_heading = "Outputs")]
struct OutputArgs {
    /// Output BAM file containing stitched molecules
    #[arg(
        short = 'o',
        long = "output",
        visible_alias = "output-bam",
        value_name = "BAM"
    )]
    output: Option<PathBuf>,

    /// Output 10x-compatible MEX matrix directory
    #[arg(long = "matrix-out", value_name = "DIR")]
    matrix_out: Option<PathBuf>,

    /// Output isoform counts TSV (usually .tsv.gz)
    #[arg(long = "counts-tsv", value_name = "FILE")]
    counts_tsv: Option<PathBuf>,

    /// Output molecule-level audit TSV (usually .tsv.gz)
    #[arg(long = "molecules-tsv", value_name = "FILE")]
    molecules_tsv: Option<PathBuf>,
}

#[derive(Args, Debug)]
#[command(next_help_heading = "Annotation and indexing")]
struct AnnotationArgs {
    /// Reuse precomputed exon interval index
    #[arg(long = "isoform", value_name = "JSON.GZ")]
    isoform: Option<PathBuf>,

    /// Reuse precomputed exon-junction index
    #[arg(long = "junction", value_name = "JSON.GZ")]
    junction: Option<PathBuf>,

    /// Save generated indexes as <prefix>.intervals.json.gz and <prefix>.refskip.json.gz
    #[arg(long = "dump-index", value_name = "PREFIX")]
    dump_index: Option<PathBuf>,

    /// Build indexes only; PREFIX is the output filename prefix
    #[arg(long = "index-only", value_name = "PREFIX")]
    index_only: Option<PathBuf>,

    /// Skip transcript isoform compatibility calling
    #[arg(long = "skip-iso")]
    skip_iso: bool,

    /// Restrict processing/indexing to one chromosome or contig
    #[arg(long = "contig", value_name = "NAME")]
    contig: Option<String>,

    /// GTF attribute used as the gene identifier
    #[arg(
        long = "gene-identifier",
        default_value = "gene_id",
        value_name = "ATTR"
    )]
    gene_identifier: String,
}

#[derive(Args, Debug)]
#[command(next_help_heading = "Read and runtime options")]
struct RuntimeArgs {
    /// Number of worker threads
    #[arg(short = 't', long = "threads", default_value_t = 4, value_name = "N")]
    threads: usize,

    /// Treat input as single-end instead of paired-end
    #[arg(long = "single-end")]
    single_end: bool,

    /// UMI tag in BAM file
    #[arg(
        long = "umi-tag",
        visible_alias = "UMI-tag",
        default_value = "UB",
        value_name = "TAG"
    )]
    umi_tag: String,

    /// Cell barcode tag in BAM file
    #[arg(long = "cell-tag", default_value = "CB", value_name = "TAG")]
    cell_tag: String,

    /// Keep only cell barcodes listed in this file
    #[arg(long = "cells", value_name = "FILE")]
    cells: Option<PathBuf>,

    /// Process only gene IDs/names listed in this file
    #[arg(long = "genes", value_name = "FILE")]
    genes: Option<PathBuf>,
}

fn validate_cli(cli: &Cli) -> Result<()> {
    let core = &cli.core;
    let outputs = &cli.outputs;
    let annotation = &cli.annotation;
    let runtime = &cli.runtime;

    if runtime.threads == 0 {
        bail!("--threads must be at least 1");
    }

    if annotation.index_only.is_some() {
        if core.input.is_some() || outputs.output.is_some() {
            bail!("--index-only cannot be combined with --input or --output");
        }
        if outputs.matrix_out.is_some()
            || outputs.counts_tsv.is_some()
            || outputs.molecules_tsv.is_some()
        {
            bail!("--index-only cannot be combined with output matrix/count options");
        }
        if annotation.isoform.is_some()
            || annotation.junction.is_some()
            || annotation.dump_index.is_some()
            || annotation.skip_iso
        {
            bail!("--index-only only needs --gtf, --contig, and --gene-identifier");
        }
        if runtime.single_end || runtime.cells.is_some() || runtime.genes.is_some() {
            bail!("--index-only cannot be combined with read or barcode filtering options");
        }
        return Ok(());
    }

    if core.input.is_none() {
        bail!("Missing required argument: --input (-i)");
    }
    if outputs.output.is_none() {
        bail!("Missing required argument: --output (-o)");
    }

    match (&annotation.isoform, &annotation.junction) {
        (Some(_), None) | (None, Some(_)) => {
            bail!("--isoform and --junction must be provided together")
        }
        _ => {}
    }

    if annotation.skip_iso
        && (annotation.isoform.is_some()
            || annotation.junction.is_some()
            || annotation.dump_index.is_some())
    {
        bail!("--skip-iso cannot be combined with --isoform, --junction, or --dump-index");
    }
    if annotation.skip_iso && outputs.matrix_out.is_some() {
        bail!("--matrix-out requires isoform calling; remove --skip-iso");
    }

    Ok(())
}

fn main() -> Result<()> {
    let args = Cli::parse();
    validate_cli(&args)?;

    // Mode 1: Index-only mode (replaces gtf_to_json.py)
    if let Some(prefix) = args.annotation.index_only.as_ref() {
        println!("Running in index-only mode for {:?}", args.core.gtf);
        let start = Instant::now();
        let (iso, jun) = build_isoform_indices_from_gtf(
            &args.core.gtf,
            args.annotation.contig.as_deref(),
            None,
            &args.annotation.gene_identifier,
        )?;

        let iso_path = format!("{}.intervals.json.gz", prefix.display());
        let jun_path = format!("{}.refskip.json.gz", prefix.display());
        println!("Saving indices to {:?} and {:?}...", iso_path, jun_path);
        dump_isoform_json(&iso, &iso_path)?;
        dump_isoform_json(&jun, &jun_path)?;

        println!(
            "Indexed {} genes successfully in {:.2}s.",
            iso.len(),
            start.elapsed().as_secs_f64()
        );
        return Ok(());
    }

    // Mode 2: Stitching pipeline mode
    let config = PipelineConfig {
        input_bam: args.core.input.expect("validated --input"),
        output_bam: args.outputs.output.expect("validated --output"),
        gtf_file: args.core.gtf,
        isoform_file: args.annotation.isoform,
        junction_file: args.annotation.junction,
        dump_index: args.annotation.dump_index,
        matrix_out: args.outputs.matrix_out,
        counts_tsv: args.outputs.counts_tsv,
        molecules_tsv: args.outputs.molecules_tsv,
        threads: args.runtime.threads,
        single_end: args.runtime.single_end,
        skip_iso: args.annotation.skip_iso,
        umi_tag: args.runtime.umi_tag,
        cell_tag: args.runtime.cell_tag,
        cells_file: args.runtime.cells,
        genes_file: args.runtime.genes,
        contig: args.annotation.contig,
        gene_identifier: args.annotation.gene_identifier,
    };

    run_pipeline(config)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn common_arguments_are_valid() {
        let cli = Cli::try_parse_from([
            "mfs-stitcher",
            "--input",
            "input.bam",
            "--output",
            "output.bam",
            "--gtf",
            "genes.gtf",
            "--threads",
            "8",
        ])
        .unwrap();

        validate_cli(&cli).unwrap();
    }

    #[test]
    fn index_only_rejects_stitching_outputs() {
        let cli = Cli::try_parse_from([
            "mfs-stitcher",
            "--gtf",
            "genes.gtf",
            "--index-only",
            "reference",
            "--matrix-out",
            "matrix",
        ])
        .unwrap();

        assert!(validate_cli(&cli).is_err());
    }

    #[test]
    fn precomputed_isoform_indexes_must_be_a_pair() {
        let cli = Cli::try_parse_from([
            "mfs-stitcher",
            "--input",
            "input.bam",
            "--output",
            "output.bam",
            "--gtf",
            "genes.gtf",
            "--isoform",
            "isoform.json.gz",
        ])
        .unwrap();

        assert!(validate_cli(&cli).is_err());
    }
}
