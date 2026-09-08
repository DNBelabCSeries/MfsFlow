"""
UMI utility functions for standalone script execution.

This module provides UMI (Unique Molecular Identifier) processing utilities,
including UMI clustering, deduplication, and analysis for single-cell
RNA sequencing data processing.
"""

import itertools
from collections import Counter


def cluster_umis(umis, threshold=1):
    """
    Deterministically map low-count UMI neighbors to higher-count parents.
    Uses mask-based indexing for faster Hamming distance computation.
    """
    if not umis:
        return {}

    counts = Counter()
    if isinstance(umis, (dict, Counter)):
        counts.update(umis)
    else:
        counts.update(umis)

    unique_umis = sorted(counts.keys(), key=lambda x: (-counts[x], x))
    parent_map = {u: u for u in unique_umis}

    if threshold == 0:
        return parent_map

    if len(unique_umis) < 100:
        return _cluster_umis_small(unique_umis, counts, threshold)

    return _cluster_umis_indexed(unique_umis, counts, threshold)


def _cluster_umis_small(unique_umis, counts, threshold):
    """Optimized clustering for small UMI sets."""
    parent_map = {u: u for u in unique_umis}
    umi_set = set(unique_umis)
    visited = set()
    bases = ("A", "C", "G", "T", "N")

    for parent in unique_umis:
        if parent in visited:
            continue

        parent_len = len(parent)
        p_chars = list(parent)

        for dist in range(1, threshold + 1):
            for positions in itertools.combinations(range(parent_len), dist):
                originals = [p_chars[pos] for pos in positions]
                replacements = []
                for orig in originals:
                    replacements.append([b for b in bases if b != orig])

                for repl_tuple in itertools.product(*replacements):
                    for pos, repl in zip(positions, repl_tuple):
                        p_chars[pos] = repl

                    child = "".join(p_chars)
                    if child in umi_set and child not in visited:
                        parent_map[child] = parent
                        visited.add(child)

                    for pos, orig in zip(positions, originals):
                        p_chars[pos] = orig

    return parent_map


def _cluster_umis_indexed(unique_umis, counts, threshold):
    """Cluster larger UMI sets without changing the distance semantics.

    The previous implementation built a position/base index but never queried
    it, and its neighbor search stopped at distance two.  Use the same exact
    substitution enumeration as the small-set path so ``threshold`` has one
    meaning on both sides of the 100-UMI dispatch boundary.
    """
    parent_map = {u: u for u in unique_umis}
    umi_set = set(unique_umis)
    visited = set()

    if not unique_umis:
        return parent_map

    umi_len = len(unique_umis[0])
    for parent in unique_umis:
        if parent in visited:
            continue

        neighbors = _find_hamming_neighbors(parent, umi_len, threshold, umi_set, visited)
        for neighbor in neighbors:
            if neighbor != parent:
                parent_map[neighbor] = parent
                visited.add(neighbor)

    return parent_map


def _find_hamming_neighbors(umi, umi_len, threshold, umi_set, visited):
    """Find all present UMIs within the requested Hamming distance."""
    neighbors = set()

    chars = list(umi)
    bases = ("A", "C", "G", "T", "N")
    for distance in range(1, threshold + 1):
        for positions in itertools.combinations(range(umi_len), distance):
            originals = [chars[pos] for pos in positions]
            replacements = [
                [base for base in bases if base != original]
                for original in originals
            ]
            for replacement_tuple in itertools.product(*replacements):
                for pos, replacement in zip(positions, replacement_tuple):
                    chars[pos] = replacement
                candidate = "".join(chars)
                if candidate in umi_set and candidate not in visited:
                    neighbors.add(candidate)
                for pos, original in zip(positions, originals):
                    chars[pos] = original

    return neighbors
