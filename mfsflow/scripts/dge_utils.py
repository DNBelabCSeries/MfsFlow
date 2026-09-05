"""Dependency-light helpers for DGE analysis."""


def resolve_worker_count(requested, task_count, performance_opts=None):
    """Resolve a bounded worker count for DGE passes.

    ``max_dge_workers`` is intentionally opt-in so existing runs keep their
    historical parallelism. It provides a simple guardrail for hosts where
    the command-line thread count exceeds the available memory.
    """
    requested = max(1, int(requested))
    task_count = max(1, int(task_count))
    options = performance_opts or {}
    cap = options.get("max_dge_workers")
    if cap not in (None, ""):
        try:
            cap = int(cap)
        except (TypeError, ValueError) as exc:
            raise ValueError("performance_opts.max_dge_workers must be an integer") from exc
        if cap < 1:
            raise ValueError("performance_opts.max_dge_workers must be positive")
        requested = min(requested, cap)
    return min(requested, task_count)


def workload_order(workloads):
    """Return barcode ids in deterministic, heaviest-first order."""
    return [
        barcode
        for barcode, _weight in sorted(
            workloads,
            key=lambda item: (-int(item[1]), str(item[0])),
        )
    ]


def dynamic_chunksize(task_count, worker_count, max_chunks_per_worker=4, max_chunksize=20):
    """Choose a small pool chunksize while retaining enough work stealing."""
    task_count = max(1, int(task_count))
    worker_count = max(1, int(worker_count))
    target_chunks = max(1, worker_count * int(max_chunks_per_worker))
    return max(1, min(int(max_chunksize), (task_count + target_chunks - 1) // target_chunks))


def balance_reference_chunks(references, mapped_counts, worker_count):
    """Balance references by indexed mapped-read counts using greedy bin packing."""
    references = list(references)
    if not references:
        return []
    worker_count = max(1, min(int(worker_count), len(references)))
    order = {reference: index for index, reference in enumerate(references)}
    weighted = sorted(
        references,
        key=lambda reference: (-max(0, int(mapped_counts.get(reference, 0))), order[reference]),
    )
    chunks = [[] for _ in range(worker_count)]
    loads = [0] * worker_count
    for reference in weighted:
        target = min(range(worker_count), key=lambda index: (loads[index], len(chunks[index]), index))
        chunks[target].append(reference)
        loads[target] += max(0, int(mapped_counts.get(reference, 0)))
    return [chunk for chunk in chunks if chunk]


def summarize_exon_intron_counts(exon_counts, intron_counts):
    """Build per-cell totals matching the historical stats table semantics."""
    exon_counts = exon_counts or {}
    intron_counts = intron_counts or {}

    def summarize_one(counts):
        result = {}
        for barcode, gene_counts in counts.items():
            positive = {gene: int(value) for gene, value in gene_counts.items() if int(value) > 0}
            result[barcode] = {"umis": sum(positive.values()), "genes": len(positive)}
        return result

    exon_summary = summarize_one(exon_counts)
    intron_summary = summarize_one(intron_counts)
    combined = {}
    for barcode in set(exon_counts) | set(intron_counts):
        exon_genes = {gene for gene, value in exon_counts.get(barcode, {}).items() if int(value) > 0}
        intron_genes = {gene for gene, value in intron_counts.get(barcode, {}).items() if int(value) > 0}
        combined[barcode] = {
            # Preserve the existing stats.tsv definition: exon total + intron total.
            "umis": exon_summary.get(barcode, {}).get("umis", 0) + intron_summary.get(barcode, {}).get("umis", 0),
            "genes": len(exon_genes | intron_genes),
        }
    return {"exon": exon_summary, "intron": intron_summary, "inex": combined}
