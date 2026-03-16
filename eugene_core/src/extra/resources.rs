use std::ops::{Deref, Index};

use crate::abstract_plants::*;

pub struct RecRate {
    rec_between_loci: Vec<f64>,
}

impl RecRate {
    pub fn new(rec_between_loci: Vec<f64>) -> Self {
        Self { rec_between_loci }
    }

    pub fn probability_gamete<A, B>(&self, x: &A, gx: &B) -> f64
    where
        A: Diploid<B>,
        B: Haploid,
    {
        let v_xu = x.upper().alleles();
        let v_xl = x.lower().alleles();
        let v_gx = gx.alleles();
        assert!(!v_gx.is_empty());
        assert!(v_xu.len() == v_xl.len() && v_xu.len() == v_gx.len());
        //assert!((0..v_gx.len()).all(|i| v_gx[i] == v_xu[i] || v_gx[i] == v_xl[i]));
        // Dynamic programming to calculate the probability of gamete
        let mut dp = [
            if v_gx[0] == v_xu[0] { 0.5 } else { 0.0 }, // upper parent
            if v_gx[0] == v_xl[0] { 0.5 } else { 0.0 }, // lower parent
        ];
        for i in 0..v_gx.len() - 1 {
            // dp[c] is the probability of creating the gamete up to locus i having sourced locus i
            // from the upper (c=0) or lower (c=1) parent
            let r = self.rec_between_loci[i];
            let b_u = if v_gx[i + 1] == v_xu[i + 1] { 1.0 } else { 0.0 };
            let b_l = if v_gx[i + 1] == v_xl[i + 1] { 1.0 } else { 0.0 };

            dp = [
                b_u * (dp[0] * (1.0 - r) + dp[1] * r), // upper parent
                b_l * (dp[1] * (1.0 - r) + dp[0] * r), // lower parent
            ]
        }
        // The probability of creating the gamete is the sum of the probabilities of sourcing it
        // from either parent
        dp[0] + dp[1]
    }

    pub fn crossing_resources<A, B>(&self, gamma: f64, x: &A, y: &A, z: &A) -> f64
    where
        A: Diploid<B>,
        B: Haploid,
    {
        assert!((0.0..=1.0).contains(&gamma));
        let px = self.probability_gamete(x, &z.upper());
        let py = self.probability_gamete(y, &z.lower());
        cost_of_crossing(gamma, px, py)
    }
}

impl From<Vec<f64>> for RecRate {
    fn from(rec_between_loci: Vec<f64>) -> Self {
        assert!(rec_between_loci.iter().all(|&r| (0.0..=0.5).contains(&r)));
        RecRate { rec_between_loci }
    }
}

impl From<&[f64]> for RecRate {
    fn from(rec_between_loci: &[f64]) -> Self {
        RecRate::from(rec_between_loci.to_vec())
    }
}

impl RecRate {
    pub fn n_loci(&self) -> usize {
        self.rec_between_loci.len() + 1
    }
}

impl Deref for RecRate {
    type Target = [f64];

    fn deref(&self) -> &Self::Target {
        &self.rec_between_loci
    }
}

pub fn cost_of_crossing(gamma: f64, px: f64, py: f64) -> f64 {
    assert!((0.0..=1.0).contains(&gamma));
    assert!((0.0..=1.0).contains(&px));
    assert!((0.0..=1.0).contains(&px));
    let rho = px * py;
    match rho {
        0.0 => f64::INFINITY,
        1.0 => 1.0,
        _ => f64::log2(1.0 - gamma) / f64::log2(1.0 - rho),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rec_rate_probability_gamete_2_loci() {
        use crate::plants::bit_array::{SingleChromGamete, SingleChromGenotype};
        let rec_rate = RecRate::from(vec![0.1]);

        let x = SingleChromGenotype::from_str("01", "10");
        let gx = SingleChromGamete::from_str("01");
        assert_eq!(rec_rate.probability_gamete(&x, &gx), 0.45);
        let gx = SingleChromGamete::from_str("10");
        assert_eq!(rec_rate.probability_gamete(&x, &gx), 0.45);
        let gx = SingleChromGamete::from_str("00");
        assert_eq!(rec_rate.probability_gamete(&x, &gx), 0.05);
        let gx = SingleChromGamete::from_str("11");
        assert_eq!(rec_rate.probability_gamete(&x, &gx), 0.05);

        let x = SingleChromGenotype::from_str("10", "11");
        let gx = SingleChromGamete::from_str("01");
        assert_eq!(rec_rate.probability_gamete(&x, &gx), 0.0);
        let gx = SingleChromGamete::from_str("10");
        assert_eq!(rec_rate.probability_gamete(&x, &gx), 0.5);
        let gx = SingleChromGamete::from_str("00");
        assert_eq!(rec_rate.probability_gamete(&x, &gx), 0.0);
        let gx = SingleChromGamete::from_str("11");
        assert_eq!(rec_rate.probability_gamete(&x, &gx), 0.5);
    }
}
