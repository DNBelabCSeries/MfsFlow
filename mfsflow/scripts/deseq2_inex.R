#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))

usage <- function() {
    cat(
        "Usage:\n",
        "  Rscript deseq2_inex.R --matrix-dir SAMPLE.inex.umi \\\n",
        "    --metadata groups.tsv --group-a control --group-b treatment \\\n",
        "    --output deseq2_results\n\n",
        "Recommended input: the UMI (molecule count) matrix, e.g. SAMPLE.inex.umi.\n",
        "DESeq2 models molecule counts; read-count matrices carry PCR/sequencing-depth\n",
        "noise and are not recommended for differential expression testing.\n\n",
        "Metadata file (--metadata, e.g. groups.tsv) format:\n",
        "  - Tab-separated, with a header row.\n",
        "  - Exactly two columns named 'sample' and 'condition' (case-sensitive).\n",
        "  - One row per matrix column (one well/sample). 'sample' values must exactly\n",
        "    match entries in barcodes.tsv.gz (trimmed, no header).\n",
        "  - 'condition' is the group label for each sample; only the two values passed\n",
        "    to --group-a and --group-b are used, all other rows are ignored.\n",
        "  - 'sample' values must be unique and non-empty; each group must contain at\n",
        "    least one sample (>=2 per group to run DESeq2 with p-values).\n",
        "  Example:\n",
        "    sample\tcondition\n",
        "    WELL_01\tcontrol\n",
        "    WELL_02\tcontrol\n",
        "    WELL_03\ttreatment\n",
        "    WELL_04\ttreatment\n\n",
        "Sample values must exactly match entries in barcodes.tsv.gz. Positive\n",
        "log2FoldChange means higher expression in group-b relative to group-a.\n\n",
        "Options:\n",
        "  --min-total-count N  Keep genes with at least N total reads (default: 10)\n",
        "  --help               Show this message\n",
        sep = ""
    )
}

parse_args <- function(args) {
    if (length(args) == 0 || "--help" %in% args) {
        usage()
        quit(status = if (length(args) == 0) 1 else 0)
    }

    value_options <- c(
        "--matrix-dir", "--metadata", "--group-a", "--group-b",
        "--output", "--min-total-count"
    )
    parsed <- list(min_total_count = 10L)
    i <- 1L
    while (i <= length(args)) {
        key <- args[[i]]
        if (!(key %in% value_options)) {
            stop("Unknown argument: ", key, call. = FALSE)
        }
        if (i == length(args)) {
            stop("Missing value for ", key, call. = FALSE)
        }
        name <- gsub("-", "_", sub("^--", "", key))
        parsed[[name]] <- args[[i + 1L]]
        i <- i + 2L
    }

    required <- c("matrix_dir", "metadata", "group_a", "group_b", "output")
    missing <- required[!vapply(required, function(x) {
        !is.null(parsed[[x]]) && nzchar(parsed[[x]])
    }, logical(1))]
    if (length(missing) > 0) {
        stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    parsed$min_total_count <- suppressWarnings(as.integer(parsed$min_total_count))
    if (is.na(parsed$min_total_count) || parsed$min_total_count < 0) {
        stop("--min-total-count must be a non-negative integer", call. = FALSE)
    }
    if (identical(parsed$group_a, parsed$group_b)) {
        stop("--group-a and --group-b must be different", call. = FALSE)
    }
    parsed
}

require_file <- function(path) {
    if (!file.exists(path) || isTRUE(file.info(path)$size == 0)) {
        stop("Missing or empty input file: ", path, call. = FALSE)
    }
}

read_gzip_table <- function(path) {
    connection <- gzfile(path, "rt")
    on.exit(close(connection), add = TRUE)
    read.delim(
        connection, header = FALSE, stringsAsFactors = FALSE,
        colClasses = "character", quote = "", comment.char = ""
    )
}

read_mex <- function(matrix_dir) {
    paths <- file.path(matrix_dir, c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"))
    invisible(lapply(paths, require_file))

    matrix_connection <- gzfile(paths[[1]], "rt")
    on.exit(close(matrix_connection), add = TRUE)
    counts <- readMM(matrix_connection)
    features <- read_gzip_table(paths[[2]])
    barcodes <- read_gzip_table(paths[[3]])

    if (nrow(features) != nrow(counts)) {
        stop("features.tsv.gz rows do not match matrix rows", call. = FALSE)
    }
    if (nrow(barcodes) != ncol(counts)) {
        stop("barcodes.tsv.gz rows do not match matrix columns", call. = FALSE)
    }
    if (ncol(features) < 1 || ncol(barcodes) < 1) {
        stop("Invalid feature or barcode table", call. = FALSE)
    }

    gene_id <- trimws(features[[1]])
    gene_name <- if (ncol(features) >= 2) trimws(features[[2]]) else gene_id
    sample_names <- trimws(barcodes[[1]])
    if (any(!nzchar(gene_id)) || anyDuplicated(gene_id)) {
        stop("Gene IDs must be non-empty and unique", call. = FALSE)
    }
    if (any(!nzchar(sample_names)) || anyDuplicated(sample_names)) {
        stop("Barcode/sample names must be non-empty and unique", call. = FALSE)
    }
    if (length(counts@x) > 0 && (any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-8))) {
        stop("DESeq2 requires non-negative integer raw counts", call. = FALSE)
    }

    counts@x <- round(counts@x)
    rownames(counts) <- gene_id
    colnames(counts) <- sample_names
    list(counts = counts, gene_id = gene_id, gene_name = gene_name)
}

read_groups <- function(path, available_samples, group_a, group_b) {
    require_file(path)
    metadata <- read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("sample", "condition")
    if (!all(required %in% colnames(metadata))) {
        stop("Metadata must contain columns: sample and condition", call. = FALSE)
    }
    metadata$sample <- trimws(as.character(metadata$sample))
    metadata$condition <- trimws(as.character(metadata$condition))
    if (any(!nzchar(metadata$sample)) || anyDuplicated(metadata$sample)) {
        stop("Metadata sample values must be non-empty and unique", call. = FALSE)
    }
    metadata <- metadata[metadata$condition %in% c(group_a, group_b), , drop = FALSE]
    if (nrow(metadata) == 0) {
        stop("Metadata contains no samples from the requested groups", call. = FALSE)
    }
    absent <- setdiff(metadata$sample, available_samples)
    if (length(absent) > 0) {
        stop("Metadata samples not found in barcodes.tsv.gz: ", paste(absent, collapse = ", "), call. = FALSE)
    }
    group_counts <- table(factor(metadata$condition, levels = c(group_a, group_b)))
    if (any(group_counts == 0)) {
        stop("Both requested groups must contain at least one sample", call. = FALSE)
    }
    metadata$condition <- factor(metadata$condition, levels = c(group_a, group_b))
    rownames(metadata) <- metadata$sample
    metadata
}

median_ratio_normalize <- function(counts) {
    dense <- as.matrix(counts)
    positive <- rowSums(dense > 0) == ncol(dense)
    if (any(positive)) {
        geometric_means <- exp(rowMeans(log(dense[positive, , drop = FALSE])))
        ratios <- sweep(dense[positive, , drop = FALSE], 1, geometric_means, "/")
        size_factors <- apply(ratios, 2, median, na.rm = TRUE)
    } else {
        size_factors <- rep(NA_real_, ncol(dense))
    }
    if (any(!is.finite(size_factors)) || any(size_factors <= 0)) {
        library_sizes <- colSums(dense)
        if (any(library_sizes <= 0)) {
            stop("Cannot normalize a sample with zero total counts", call. = FALSE)
        }
        size_factors <- library_sizes / exp(mean(log(library_sizes)))
    }
    list(counts = sweep(dense, 2, size_factors, "/"), size_factors = size_factors)
}

write_tsv_gz <- function(data, path, row_names = FALSE) {
    connection <- gzfile(path, "wt")
    on.exit(close(connection), add = TRUE)
    write.table(data, connection, sep = "\t", quote = FALSE, row.names = row_names)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
mex <- read_mex(args$matrix_dir)
metadata <- read_groups(args$metadata, colnames(mex$counts), args$group_a, args$group_b)
counts <- mex$counts[, metadata$sample, drop = FALSE]

keep <- Matrix::rowSums(counts) >= args$min_total_count
if (!any(keep)) {
    stop("No genes remain after --min-total-count filtering", call. = FALSE)
}
counts <- counts[keep, , drop = FALSE]
gene_names <- mex$gene_name[keep]
dir.create(args$output, recursive = TRUE, showWarnings = FALSE)

replicates <- table(metadata$condition)
has_replicates <- all(replicates >= 2)
result_path <- file.path(args$output, "differential_expression.tsv.gz")
normalized_path <- file.path(args$output, "normalized_counts.tsv.gz")
summary_path <- file.path(args$output, "analysis_summary.txt")

if (has_replicates) {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop(
            "DESeq2 is required for replicated analysis. Install it with: ",
            "BiocManager::install('DESeq2')", call. = FALSE
        )
    }
    integer_counts <- round(as.matrix(counts))
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = integer_counts,
        colData = metadata,
        design = ~ condition
    )
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    deseq_result <- DESeq2::results(
        dds,
        contrast = c("condition", args$group_b, args$group_a),
        independentFiltering = TRUE
    )
    result <- as.data.frame(deseq_result)
    result$gene_id <- rownames(result)
    result$gene_name <- gene_names[match(result$gene_id, rownames(counts))]
    result$analysis_mode <- "DESeq2"
    result <- result[, c(
        "gene_id", "gene_name", "baseMean", "log2FoldChange", "lfcSE",
        "stat", "pvalue", "padj", "analysis_mode"
    )]
    result <- result[order(result$padj, na.last = TRUE), , drop = FALSE]
    normalized <- DESeq2::counts(dds, normalized = TRUE)

    pdf(file.path(args$output, "ma_plot.pdf"), width = 7, height = 6)
    DESeq2::plotMA(deseq_result, alpha = 0.05)
    dev.off()
    mode_message <- "DESeq2 differential test with biological replication"
} else {
    normalized_data <- median_ratio_normalize(counts)
    normalized <- normalized_data$counts
    group_a_mean <- rowMeans(normalized[, metadata$condition == args$group_a, drop = FALSE])
    group_b_mean <- rowMeans(normalized[, metadata$condition == args$group_b, drop = FALSE])
    result <- data.frame(
        gene_id = rownames(counts),
        gene_name = gene_names,
        baseMean = rowMeans(normalized),
        log2FoldChange = log2((group_b_mean + 0.5) / (group_a_mean + 0.5)),
        lfcSE = NA_real_,
        stat = NA_real_,
        pvalue = NA_real_,
        padj = NA_real_,
        analysis_mode = "descriptive_no_replication",
        stringsAsFactors = FALSE
    )
    result <- result[order(abs(result$log2FoldChange), decreasing = TRUE), , drop = FALSE]
    mode_message <- paste(
        "Descriptive comparison only: at least one group has fewer than two",
        "independent samples, so p-values and adjusted p-values were not calculated."
    )
}

normalized_table <- data.frame(
    gene_id = rownames(normalized),
    gene_name = gene_names[match(rownames(normalized), rownames(counts))],
    normalized,
    check.names = FALSE
)
write_tsv_gz(result, result_path)
write_tsv_gz(normalized_table, normalized_path)
writeLines(c(
    paste("Mode:", mode_message),
    paste("Matrix:", normalizePath(args$matrix_dir)),
    paste("Comparison:", args$group_b, "vs", args$group_a),
    paste("Positive log2FoldChange:", args$group_b, "higher than", args$group_a),
    paste("Samples:", paste(names(replicates), as.integer(replicates), sep = "=", collapse = ", ")),
    paste("Genes tested/reported:", nrow(result)),
    paste("Minimum total count:", args$min_total_count)
), summary_path)

message(mode_message)
message("Results: ", result_path)
message("Normalized counts: ", normalized_path)
