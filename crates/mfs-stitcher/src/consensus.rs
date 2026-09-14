//! Bayesian consensus nucleotide sequence and quality score calling.

use std::f64;

pub const MAX_Q: usize = 94;

pub struct ConsensusTables {
    pub ll_this: [f64; MAX_Q],
    pub ll_other: [f64; MAX_Q],
    pub ll_n: f64,
}

impl ConsensusTables {
    pub fn new() -> Self {
        let mut ll_this = [f64::NEG_INFINITY; MAX_Q];
        let mut ll_other = [f64::NEG_INFINITY; MAX_Q];
        let ln_3 = 3.0_f64.ln();
        let ln_10 = 10.0_f64.ln();

        // Q = 0
        ll_this[0] = f64::NEG_INFINITY;
        ll_other[0] = -ln_3;

        for q in 1..MAX_Q {
            let q_f = q as f64;
            let p_err = 10.0_f64.powf(-q_f / 10.0);
            ll_this[q] = (1.0 - p_err).ln();
            ll_other[q] = -(q_f * ln_10) / 10.0 - ln_3;
        }

        let ll_n = -(4.0_f64.ln());

        ConsensusTables {
            ll_this,
            ll_other,
            ll_n,
        }
    }
}

impl Default for ConsensusTables {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct BaseAccumulator {
    pub ll_a: f64,
    pub ll_t: f64,
    pub ll_c: f64,
    pub ll_g: f64,
}

impl BaseAccumulator {
    #[inline]
    pub fn new() -> Self {
        BaseAccumulator {
            ll_a: 0.0,
            ll_t: 0.0,
            ll_c: 0.0,
            ll_g: 0.0,
        }
    }

    #[inline]
    pub fn add_read_base(&mut self, base: u8, qual: u8, tables: &ConsensusTables) {
        let q = (qual as usize).min(MAX_Q - 1);
        match base {
            b'A' | b'a' => {
                let other = tables.ll_other[q];
                let delta = tables.ll_this[q] - other;
                self.ll_a += other + delta;
                self.ll_t += other;
                self.ll_c += other;
                self.ll_g += other;
            }
            b'T' | b't' => {
                let other = tables.ll_other[q];
                let delta = tables.ll_this[q] - other;
                self.ll_a += other;
                self.ll_t += other + delta;
                self.ll_c += other;
                self.ll_g += other;
            }
            b'C' | b'c' => {
                let other = tables.ll_other[q];
                let delta = tables.ll_this[q] - other;
                self.ll_a += other;
                self.ll_t += other;
                self.ll_c += other + delta;
                self.ll_g += other;
            }
            b'G' | b'g' => {
                let other = tables.ll_other[q];
                let delta = tables.ll_this[q] - other;
                self.ll_a += other;
                self.ll_t += other;
                self.ll_c += other;
                self.ll_g += other + delta;
            }
            _ => {
                // 'N' or unknown
                self.ll_a += tables.ll_n;
                self.ll_t += tables.ll_n;
                self.ll_c += tables.ll_n;
                self.ll_g += tables.ll_n;
            }
        }
    }

    #[inline]
    pub fn finalize(&self) -> (u8, u8) {
        let lls = [self.ll_a, self.ll_t, self.ll_c, self.ll_g];
        let bases = [b'A', b'T', b'C', b'G'];

        let mut max_ll = lls[0];
        let mut max_idx = 0;
        for (i, &ll) in lls.iter().enumerate().skip(1) {
            if ll > max_ll {
                max_ll = ll;
                max_idx = i;
            }
        }

        // logsumexp
        let sum_exp = (lls[0] - max_ll).exp()
            + (lls[1] - max_ll).exp()
            + (lls[2] - max_ll).exp()
            + (lls[3] - max_ll).exp();
        let full_ll = max_ll + sum_exp.ln();

        let prob_max = (max_ll - full_ll).exp();

        let call_base = if prob_max > 0.3 { bases[max_idx] } else { b'N' };

        let diff = (1.0 - prob_max + 1e-13).max(1e-13);
        let phred = (-10.0 * diff.log10()).round();
        let phred_clamped = phred.clamp(0.0, 93.0) as u8;

        (call_base, phred_clamped)
    }
}

/// Call consensus sequence and qualities from a vector of accumulators
pub fn call_consensus(accumulators: &[BaseAccumulator]) -> (Vec<u8>, Vec<u8>) {
    let mut seq = Vec::with_capacity(accumulators.len());
    let mut qual = Vec::with_capacity(accumulators.len());

    for acc in accumulators {
        let (b, q) = acc.finalize();
        seq.push(b);
        qual.push(q);
    }

    (seq, qual)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_high_quality_base() {
        let tables = ConsensusTables::new();
        let mut acc = BaseAccumulator::new();
        acc.add_read_base(b'A', 30, &tables);
        let (base, q) = acc.finalize();
        assert_eq!(base, b'A');
        assert!(q >= 29);
    }

    #[test]
    fn test_two_reads_agreement() {
        let tables = ConsensusTables::new();
        let mut acc = BaseAccumulator::new();
        acc.add_read_base(b'C', 20, &tables);
        acc.add_read_base(b'C', 30, &tables);
        let (base, q) = acc.finalize();
        assert_eq!(base, b'C');
        assert!(q >= 35);
    }

    #[test]
    fn test_conflict_resolution() {
        let tables = ConsensusTables::new();
        let mut acc = BaseAccumulator::new();
        acc.add_read_base(b'G', 35, &tables);
        acc.add_read_base(b'T', 10, &tables);
        let (base, _) = acc.finalize();
        assert_eq!(base, b'G');
    }
}
