"""Dependency-light helpers for DGE analysis."""

import os
import pickle
import sqlite3
import tempfile
from collections import Counter, defaultdict


def open_pass1_store(out_dir, project, tmp_root=None):
    """Create a disk-backed store for one DGE pass-1 result.

    Keeping one serialized record per barcode/chunk avoids retaining the full
    raw UMI graph in the parent process between the counting and clustering
    passes. The caller owns the returned connection and path.
    """
    store_dir = os.fspath(tmp_root) if tmp_root else os.path.join(out_dir, "intermediate")
    os.makedirs(store_dir, exist_ok=True)
    fd, path = tempfile.mkstemp(
        prefix=f".{project}.dge-pass1-",
        suffix=".sqlite3",
        dir=store_dir,
    )
    os.close(fd)
    try:
        connection = sqlite3.connect(path)
        connection.execute(
            "CREATE TABLE pass1_chunks ("
            "kind TEXT NOT NULL, "
            "ftype TEXT NOT NULL, "
            "barcode TEXT NOT NULL, "
            "weight INTEGER NOT NULL, "
            "payload BLOB NOT NULL)"
        )
        # WAL is unnecessary for one writer and would leave sidecar files in
        # the temporary directory. The pass is rebuilt on every run.
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        return connection, path
    except Exception:
        try:
            os.unlink(path)
        except OSError:
            pass
        raise


def store_pass1_result(connection, partial_read, partial_umi, partial_global):
    """Persist one count-worker result and release it from the caller."""
    rows = []

    for ftype in ("exon", "intron"):
        for barcode, genes in partial_read.get(ftype, {}).items():
            if genes:
                rows.append(
                    (
                        "read",
                        ftype,
                        str(barcode),
                        sum(int(count) for count in genes.values()),
                        sqlite3.Binary(pickle.dumps(dict(genes), protocol=pickle.HIGHEST_PROTOCOL)),
                    )
                )
        for barcode, genes in partial_umi.get(ftype, {}).items():
            if genes:
                payload = {gene: dict(umis) for gene, umis in genes.items()}
                rows.append(
                    (
                        "umi",
                        ftype,
                        str(barcode),
                        sum(len(umis) for umis in payload.values()),
                        sqlite3.Binary(pickle.dumps(payload, protocol=pickle.HIGHEST_PROTOCOL)),
                    )
                )

    for barcode, umis in (partial_global or {}).items():
        if umis:
            rows.append(
                (
                    "global",
                    "global",
                    str(barcode),
                    len(umis),
                    sqlite3.Binary(pickle.dumps(dict(umis), protocol=pickle.HIGHEST_PROTOCOL)),
                )
            )

    if rows:
        connection.executemany(
            "INSERT INTO pass1_chunks(kind, ftype, barcode, weight, payload) "
            "VALUES (?, ?, ?, ?, ?)",
            rows,
        )
    connection.commit()


def finalize_pass1_store(connection):
    """Add the lookup index after pass-1 insertion is complete."""
    connection.execute(
        "CREATE INDEX pass1_chunks_lookup "
        "ON pass1_chunks(kind, ftype, barcode)"
    )
    connection.commit()


def close_pass1_store(connection, path):
    """Close and remove a temporary pass-1 store."""
    close_error = None
    if connection is not None:
        try:
            connection.close()
        except Exception as exc:  # pragma: no cover - sqlite close is normally infallible
            close_error = exc
    if path:
        try:
            os.unlink(path)
        except FileNotFoundError:
            pass
    if close_error is not None:
        raise close_error


def pass1_barcode_workloads(connection, include_global):
    """Return barcode workloads using the same unique-UMI weights as before."""
    kinds = ("umi", "global") if include_global else ("umi",)
    placeholders = ",".join("?" for _ in kinds)
    rows = connection.execute(
        "SELECT barcode, SUM(weight) FROM pass1_chunks "
        f"WHERE kind IN ({placeholders}) GROUP BY barcode",
        kinds,
    )
    return [(str(barcode), int(weight or 0)) for barcode, weight in rows]


def pass1_barcodes(connection, kind):
    """Return distinct barcodes stored for a pass-1 result kind."""
    rows = connection.execute(
        "SELECT DISTINCT barcode FROM pass1_chunks WHERE kind = ? ORDER BY barcode",
        (kind,),
    )
    return [str(row[0]) for row in rows]


def _load_payloads(connection, kind, ftype, barcode):
    return connection.execute(
        "SELECT payload FROM pass1_chunks "
        "WHERE kind = ? AND ftype = ? AND barcode = ?",
        (kind, ftype, barcode),
    )


def load_pass1_read_counts(connection, barcode, ftype):
    """Merge one barcode/type's per-gene read counts from disk."""
    merged = defaultdict(int)
    for (payload,) in _load_payloads(connection, "read", ftype, barcode):
        for gene, count in pickle.loads(payload).items():
            merged[gene] += int(count)
    return dict(merged)


def load_pass1_umi_counts(connection, barcode, ftype):
    """Merge one barcode/type's per-gene UMI counts from disk."""
    merged = defaultdict(Counter)
    for (payload,) in _load_payloads(connection, "umi", ftype, barcode):
        for gene, umis in pickle.loads(payload).items():
            merged[gene].update(umis)
    return {gene: Counter(umis) for gene, umis in merged.items()}


def load_pass1_global_counts(connection, barcode):
    """Merge one barcode's global UMI counts from disk."""
    merged = Counter()
    for (payload,) in _load_payloads(connection, "global", "global", barcode):
        merged.update(pickle.loads(payload))
    return merged


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
