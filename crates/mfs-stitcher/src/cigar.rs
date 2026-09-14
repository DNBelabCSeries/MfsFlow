//! CIGAR reconstruction from exonic blocks, spliced junctions, and deletions.

use crate::interval::{Interval, IntervalSet};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CigarOp {
    pub op: char, // 'M', 'N', 'D'
    pub len: u32,
}

#[derive(Debug, Clone)]
pub struct StitchedCigar {
    pub pos: i64,        // 1-based start position
    pub pos_0based: i64, // 0-based start position
    pub cigar_string: String,
    pub ops: Vec<CigarOp>,
    pub conflict: bool,
    pub n_conflict: i64,
    pub conflict_intervals: Vec<i64>,
}

/// Construct 1-based POS, CIGAR, and conflict metadata from ref and skipped intervals.
pub fn build_cigar(
    ref_intervals: &IntervalSet,
    skipped_intervals: &IntervalSet,
) -> Option<StitchedCigar> {
    if ref_intervals.is_empty() {
        return None;
    }

    let mut current_skipped = skipped_intervals.clone();
    let intersect = ref_intervals.intersection(&current_skipped);

    let conflict = !intersect.is_empty();
    let n_conflict = intersect.total_len();
    let mut conflict_intervals = Vec::new();

    if conflict {
        current_skipped = current_skipped.difference(&intersect);
        for iv in &intersect.intervals {
            conflict_intervals.push(iv.start);
            conflict_intervals.push(iv.end);
        }
    }

    let min_start = ref_intervals.min_start()?;
    let max_end = ref_intervals.max_end()?;

    // Ensure skipped intervals do not extend outside the bounds of aligned reference intervals
    let valid_range = IntervalSet::from_intervals(vec![Interval::new(min_start, max_end)]);
    current_skipped = current_skipped.intersection(&valid_range);

    // del_intervals are regions between min_start and max_end not covered by ref or skipped
    let combined = ref_intervals.union(&current_skipped);
    let del_intervals = combined.internal_gaps(min_start, max_end);

    // Collect all intervals tagged with operator
    let mut tagged: Vec<(Interval, char)> = Vec::new();
    for iv in &ref_intervals.intervals {
        tagged.push((*iv, 'M'));
    }
    for iv in &current_skipped.intervals {
        tagged.push((*iv, 'N'));
    }
    for iv in &del_intervals.intervals {
        tagged.push((*iv, 'D'));
    }

    // Sort by start position
    tagged.sort_by_key(|(iv, _)| iv.start);

    let mut ops: Vec<CigarOp> = Vec::new();
    let mut cigar_string = String::new();

    for (iv, op) in tagged {
        let len = (iv.end - iv.start + 1) as u32;
        if len == 0 {
            continue;
        }
        if let Some(last) = ops.last_mut() {
            if last.op == op {
                last.len += len;
                continue;
            }
        }
        ops.push(CigarOp { op, len });
    }

    for op in &ops {
        cigar_string.push_str(&format!("{}{}", op.len, op.op));
    }

    let pos_0based = min_start;
    let pos = min_start + 1;

    Some(StitchedCigar {
        pos,
        pos_0based,
        cigar_string,
        ops,
        conflict,
        n_conflict,
        conflict_intervals,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_exon_no_splice() {
        let r = IntervalSet::from_intervals(vec![Interval::new(100, 249)]);
        let s = IntervalSet::new();
        let res = build_cigar(&r, &s).unwrap();
        assert_eq!(res.pos, 101);
        assert_eq!(res.cigar_string, "150M");
        assert!(!res.conflict);
    }

    #[test]
    fn test_two_exons_with_splice() {
        let r = IntervalSet::from_intervals(vec![Interval::new(100, 249), Interval::new(500, 649)]);
        let s = IntervalSet::from_intervals(vec![Interval::new(250, 499)]);
        let res = build_cigar(&r, &s).unwrap();
        assert_eq!(res.pos, 101);
        assert_eq!(res.cigar_string, "150M250N150M");
        assert!(!res.conflict);
    }

    #[test]
    fn test_internal_deletion() {
        let r = IntervalSet::from_intervals(vec![Interval::new(100, 199), Interval::new(220, 299)]);
        let s = IntervalSet::new();
        let res = build_cigar(&r, &s).unwrap();
        assert_eq!(res.pos, 101);
        assert_eq!(res.cigar_string, "100M20D80M");
    }

    #[test]
    fn test_conflict_resolution() {
        // Exon covers 100..300, but splice says 250..400
        let r = IntervalSet::from_intervals(vec![Interval::new(100, 300), Interval::new(401, 500)]);
        let s = IntervalSet::from_intervals(vec![Interval::new(250, 400)]);
        let res = build_cigar(&r, &s).unwrap();
        assert!(res.conflict);
        assert_eq!(res.n_conflict, 51); // 250..300
                                        // Spliced interval 250..400 minus 250..300 -> 301..400 (100 bp)
        assert_eq!(res.cigar_string, "201M100N100M");
    }
}
