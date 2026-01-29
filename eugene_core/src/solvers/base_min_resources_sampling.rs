use crate::abstract_plants::{Crosspoint, WGam, WGen};
use crate::plants::bit_array::{CrosspointBitVec, SingleChromGamete, SingleChromGenotype};
use crate::plants::dist_array::DistArray;
use crate::extra::resources::{RecRate, cost_of_crossing};
use crate::solution::BaseSolution;
use rand::seq::SliceRandom;

type WGe = WGen<SingleChromGenotype, SingleChromGamete>;
type WGa = WGam<SingleChromGenotype, SingleChromGamete>;

pub fn breeding_program_distribute(xs: &DistArray, rec_rate: &RecRate, gamma: f64) -> Option<BaseSolution> {
    let n_loci = xs.n_loci();
    assert!(gamma >= 0.0 && gamma <= 1.0);
    let pop_0 = SingleChromGenotype::init_pop_distribute(xs);
    let ideotype = SingleChromGenotype::ideotype(n_loci);

    let mut total_cost = 0.0;
    let mut pop: Vec<WGe> = pop_0.into_iter().map(|x| WGen::new(x)).collect();
    let mut rng = rand::thread_rng();

    while !pop.iter().any(|wx| wx.genotype() == &ideotype) {
        let wx = pop.choose(&mut rng)?;
        let wy = pop.choose(&mut rng)?;

        let gx = {
            let crosspoint = CrosspointBitVec::random_crosspoint(&mut rng, n_loci, rec_rate);
            WGa::new_from_genotype(
                SingleChromGamete::from_crosspoint(
                    wx.genotype(),
                    wy.genotype(),
                    &crosspoint
                ),
                WGen::from_gametes(wx, wy)
            )
        };
        let gy = {
            let crosspoint = CrosspointBitVec::random(n_loci - 1, &mut rng);
            WGa::new_from_genotype(
                SingleChromGamete::from_crosspoint(
                    wy.genotype(),
                    wx.genotype(),
                    &crosspoint
                ),
                WGen::from_gametes(wy, wx)
            )
        };

        let px = rec_rate.probability_gamete(wx.genotype(), gx);
        let py = rec_rate.probability_gamete(wy.genotype(), gy);
        let cost = cost_of_crossing(gamma, px, py);
        total_cost += cost;

        let z = todo!();
        pop.push(z);
    }
    None
}
