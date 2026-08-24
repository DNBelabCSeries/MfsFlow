"""Dependency-light helpers for DGE analysis."""


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
