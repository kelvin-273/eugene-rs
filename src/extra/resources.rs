use crate::abstract_plants::*;

struct RecRate {
    rec_between_loci: Vec<f64>
}

impl RecRate {
    pub fn probability_gamete<A, B>(&self, x: &A, gx: &B) -> f64
    where
        A: Genotype<B>,
        B: Gamete<A>
    {
        unimplemented!();
    }

    pub fn crossing_resources<A, B>(&self, gamma: f64, x: &A, y: &A, z: &A) -> f64
    where
        A: Genotype<B> + Diploid<B>,
        B: Gamete<A> + Haploid
    {
        assert!(0.0 <= gamma && gamma >= 1.0);
        let px = self.probability_gamete(x, &z.upper());
        let py = self.probability_gamete(y, &z.lower());
        let rho = px * py;
        match rho {
            0.0 => f64::INFINITY,
            _ => f64::log2(1.0 - gamma) / f64::log2(1.0 - rho)
        }
    }
}
