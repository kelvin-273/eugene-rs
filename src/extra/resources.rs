use crate::abstract_plants::*;

pub struct RecRate {
    rec_between_loci: Vec<f64>,
}

impl RecRate {
    pub fn probability_gamete<A, B>(&self, x: &A, gx: &B) -> f64
    where
        A: Genotype<B> + Diploid<B>,
        B: Gamete<A> + Haploid
    {
        let v_xu = x.upper().alleles();
        let v_xl = x.lower().alleles();
        let v_gx = gx.alleles();
        assert!(!v_gx.is_empty());
        assert!(v_xu.len() == v_xl.len() && v_xu.len() == v_gx.len());
        assert!((0..v_gx.len()).all(|i| v_gx[i] == v_xu[i] || v_gx[i] == v_xl[i]));
        // Dynamic programming to calculate the probability of gamete
        let mut dp = [1.0; 2];
        for i in 0..v_gx.len() {
            // dp[c] is the probability of creating the gamete up to locus i having sourced locus i
            // from the upper (c=0) or lower (c=1) parent
            let r = self.rec_between_loci[i];
            let b_u = if v_gx[i] == v_xu[i] { 1.0 } else { 0.0 };
            let b_l = if v_gx[i] == v_xl[i] { 1.0 } else { 0.0 };

            dp = [
                dp[0] * (1.0 - r) * b_u + dp[1] * r * b_l, // upper parent
                dp[1] * (1.0 - r) * b_l + dp[0] * r * b_u // lower parent
            ]
        }
        // The probability of creating the gamete is the sum of the probabilities of sourcing it
        // from either parent
        dp[0] + dp[1]
    }

    pub fn crossing_resources<A, B>(&self, gamma: f64, x: &A, y: &A, z: &A) -> f64
    where
        A: Genotype<B> + Diploid<B>,
        B: Gamete<A> + Haploid
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

pub fn cost_of_crossing(gamma: f64, px: f64, py: f64) -> f64 {
    let rho = px * py;
    match rho {
        0.0 => f64::INFINITY,
        _ => f64::log2(1.0 - gamma) / f64::log2(1.0 - rho)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rec_rate_probability_gamete() {
        use crate::plants::bit_array::{SingleChromGamete, SingleChromGenotype};
        let rec_rate = RecRate {
            rec_between_loci: vec![0.1, 0.2, 0.3]
        };
        let x = SingleChromGenotype::from_str("011", "110");
        let gx = SingleChromGamete::from_str("011");
        let _prob = rec_rate.probability_gamete(&x, &gx);
        todo!("Implement test for probability_gamete");
    }
}
