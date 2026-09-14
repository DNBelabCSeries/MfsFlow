//! Fast, zero-overhead closed integer interval set operations for genomic coordinates.

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Interval {
    pub start: i64,
    pub end: i64, // inclusive [start, end]
}

impl Interval {
    #[inline]
    pub fn new(start: i64, end: i64) -> Self {
        assert!(
            start <= end,
            "Invalid interval: start > end ({} > {})",
            start,
            end
        );
        Interval { start, end }
    }

    #[inline]
    pub fn len(&self) -> i64 {
        self.end - self.start + 1
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.start > self.end
    }

    #[inline]
    pub fn overlaps(&self, other: &Interval) -> bool {
        self.start <= other.end && other.start <= self.end
    }

    #[inline]
    pub fn contains_point(&self, point: i64) -> bool {
        point >= self.start && point <= self.end
    }
}

/// A set of disjoint, sorted closed intervals.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct IntervalSet {
    pub intervals: Vec<Interval>,
}

impl IntervalSet {
    pub fn new() -> Self {
        IntervalSet {
            intervals: Vec::new(),
        }
    }

    pub fn from_intervals(mut intervals: Vec<Interval>) -> Self {
        if intervals.is_empty() {
            return IntervalSet::new();
        }
        intervals.sort_by_key(|i| (i.start, i.end));
        let mut merged: Vec<Interval> = Vec::with_capacity(intervals.len());
        for iv in intervals {
            if iv.start > iv.end {
                continue;
            }
            if let Some(last) = merged.last_mut() {
                if iv.start <= last.end + 1 {
                    last.end = last.end.max(iv.end);
                } else {
                    merged.push(iv);
                }
            } else {
                merged.push(iv);
            }
        }
        IntervalSet { intervals: merged }
    }

    /// Construct IntervalSet from a list of integer positions (e.g. sorted 0-based coordinates).
    pub fn from_positions(positions: &[i64]) -> Self {
        if positions.is_empty() {
            return IntervalSet::new();
        }
        let mut intervals = Vec::new();
        let mut start = positions[0];
        let mut prev = positions[0];

        for &pos in &positions[1..] {
            if pos == prev + 1 {
                prev = pos;
            } else if pos > prev + 1 {
                intervals.push(Interval::new(start, prev));
                start = pos;
                prev = pos;
            }
        }
        intervals.push(Interval::new(start, prev));
        IntervalSet { intervals }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }

    #[inline]
    pub fn total_len(&self) -> i64 {
        self.intervals.iter().map(|iv| iv.len()).sum()
    }

    #[inline]
    pub fn min_start(&self) -> Option<i64> {
        self.intervals.first().map(|iv| iv.start)
    }

    #[inline]
    pub fn max_end(&self) -> Option<i64> {
        self.intervals.last().map(|iv| iv.end)
    }

    /// Union of two interval sets
    pub fn union(&self, other: &IntervalSet) -> IntervalSet {
        let mut all = Vec::with_capacity(self.intervals.len() + other.intervals.len());
        all.extend_from_slice(&self.intervals);
        all.extend_from_slice(&other.intervals);
        IntervalSet::from_intervals(all)
    }

    /// Intersection of two interval sets
    pub fn intersection(&self, other: &IntervalSet) -> IntervalSet {
        let mut res = Vec::new();
        let mut i = 0;
        let mut j = 0;

        while i < self.intervals.len() && j < other.intervals.len() {
            let a = &self.intervals[i];
            let b = &other.intervals[j];

            let start = a.start.max(b.start);
            let end = a.end.min(b.end);

            if start <= end {
                res.push(Interval::new(start, end));
            }

            if a.end < b.end {
                i += 1;
            } else {
                j += 1;
            }
        }

        IntervalSet { intervals: res }
    }

    /// Difference (self - other)
    pub fn difference(&self, other: &IntervalSet) -> IntervalSet {
        let mut res = Vec::new();

        for iv in &self.intervals {
            let mut curr_start = iv.start;
            let curr_end = iv.end;

            for o in &other.intervals {
                if o.end < curr_start {
                    continue;
                }
                if o.start > curr_end {
                    break;
                }
                if o.start > curr_start {
                    res.push(Interval::new(curr_start, o.start - 1));
                }
                curr_start = curr_start.max(o.end + 1);
                if curr_start > curr_end {
                    break;
                }
            }

            if curr_start <= curr_end {
                res.push(Interval::new(curr_start, curr_end));
            }
        }

        IntervalSet { intervals: res }
    }

    /// Internal gaps (complement between min_start and max_end)
    pub fn internal_gaps(&self, min_start: i64, max_end: i64) -> IntervalSet {
        let span = IntervalSet {
            intervals: vec![Interval::new(min_start, max_end)],
        };
        span.difference(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_positions() {
        let pos = vec![1, 2, 3, 5, 6, 10];
        let set = IntervalSet::from_positions(&pos);
        assert_eq!(
            set.intervals,
            vec![
                Interval::new(1, 3),
                Interval::new(5, 6),
                Interval::new(10, 10),
            ]
        );
    }

    #[test]
    fn test_union() {
        let a = IntervalSet::from_intervals(vec![Interval::new(1, 5), Interval::new(10, 15)]);
        let b = IntervalSet::from_intervals(vec![Interval::new(4, 8), Interval::new(20, 25)]);
        let u = a.union(&b);
        assert_eq!(
            u.intervals,
            vec![
                Interval::new(1, 8),
                Interval::new(10, 15),
                Interval::new(20, 25),
            ]
        );
    }

    #[test]
    fn test_intersection() {
        let a = IntervalSet::from_intervals(vec![Interval::new(1, 10), Interval::new(20, 30)]);
        let b = IntervalSet::from_intervals(vec![Interval::new(5, 15), Interval::new(25, 35)]);
        let inter = a.intersection(&b);
        assert_eq!(
            inter.intervals,
            vec![Interval::new(5, 10), Interval::new(25, 30)]
        );
    }

    #[test]
    fn test_difference() {
        let a = IntervalSet::from_intervals(vec![Interval::new(1, 20)]);
        let b = IntervalSet::from_intervals(vec![Interval::new(5, 10), Interval::new(15, 18)]);
        let diff = a.difference(&b);
        assert_eq!(
            diff.intervals,
            vec![
                Interval::new(1, 4),
                Interval::new(11, 14),
                Interval::new(19, 20),
            ]
        );
    }

    #[test]
    fn test_internal_gaps() {
        let a = IntervalSet::from_intervals(vec![Interval::new(10, 20), Interval::new(30, 40)]);
        let gaps = a.internal_gaps(10, 40);
        assert_eq!(gaps.intervals, vec![Interval::new(21, 29)]);
    }
}
