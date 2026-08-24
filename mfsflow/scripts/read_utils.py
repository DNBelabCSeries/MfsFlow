"""Helpers for applying read-pair aware counting rules."""


def is_pair_representative(read_obj):
    """Return whether a SAM/BAM record represents its read pair.

    Secondary and supplementary alignments are never representatives. For PE
    data, the primary R1 record represents the fragment, matching zUMIs'
    first-mate counting convention. R2 is never promoted when R1 is unmapped.

    Unpaired records are always representatives, which preserves SE behavior
    and makes the helper tolerant of orphan records.
    """
    if isinstance(read_obj, str):
        fields = read_obj.split("\t", 9)
        if len(fields) < 2:
            return False
        flag = int(fields[1])
        if flag & 0x900:
            return False
        if not (flag & 0x1):
            return True

        return bool(flag & 0x40)

    if (
        getattr(read_obj, "is_secondary", False)
        or getattr(read_obj, "is_supplementary", False)
    ):
        return False
    if not getattr(read_obj, "is_paired", False):
        return True

    return bool(getattr(read_obj, "is_read1", False))
