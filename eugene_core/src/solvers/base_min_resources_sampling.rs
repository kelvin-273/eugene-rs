use std::collections::{HashMap, HashSet};

use crate::abstract_plants::{Crosspoint, WGam, WGen};
use crate::extra::resources::{cost_of_crossing, RecRate};
use crate::plants::bit_array::{CrosspointBitVec, SingleChromGamete, SingleChromGenotype};
use crate::plants::dist_array::DistArray;
use crate::solution::BaseSolution;
use rand::seq::SliceRandom;

type WGe = WGen<SingleChromGenotype, SingleChromGamete>;
type WGa = WGam<SingleChromGenotype, SingleChromGamete>;

/// Breeding program of a distribute instance of MinRes.
///
/// # Arguments
/// * `xs` - The distribute array representing the initial population
/// * `rec_rate` - The recombination rate between loci
/// * `gamma` - The threshold probability required a crossing to succeed
///
/// # Returns
/// * `Option<BaseSolution>` - The breeding program solution if found within the limits
///
/// # Panics
/// Panics if `gamma` is not in the range [0.0, 1.0]
pub fn breeding_program_distribute(
    xs: &DistArray,
    rec_rate: &RecRate,
    gamma: f64,
) -> Option<BaseSolution> {
    let n_loci = xs.n_loci();
    assert!(gamma >= 0.0 && gamma <= 1.0);
    let pop_0 = SingleChromGenotype::init_pop_distribute(xs);
    breeding_program(n_loci, &pop_0, rec_rate, gamma)
}

/// Breeding program of an instance of MinRes.
///
/// # Arguments
/// * `n_loci` - The number of loci in the genotype
/// * `pop_0` - The initial population of genotypes
/// * `rec_rate` - The recombination rate between loci
/// * `gamma` - The threshold probability required a crossing to succeed
///
/// # Returns
/// * `Option<BaseSolution>` - The breeding program solution if found within the limits
///
/// # Panics
/// Panics if `gamma` is not in the range [0.0, 1.0]
pub fn breeding_program(
    n_loci: usize,
    pop_0: &Vec<SingleChromGenotype>,
    rec_rate: &RecRate,
    gamma: f64,
) -> Option<BaseSolution> {
    let ideotype = SingleChromGenotype::ideotype(n_loci);

    let mut pop: Vec<WGe> = pop_0.iter().cloned().map(WGen::new).collect();
    let mut rng = rand::thread_rng();
    let mut h: HashMap<_, _> = pop
        .iter()
        .map(|wx| (wx.genotype().clone(), wx.clone()))
        .collect();

    while !pop.iter().any(|wx| wx.genotype() == &ideotype) {
        let wx = pop.choose(&mut rng)?;
        let wy = pop.choose(&mut rng)?;

        let kx = CrosspointBitVec::random_crosspoint(&mut rng, n_loci, rec_rate);
        let ky = CrosspointBitVec::random_crosspoint(&mut rng, n_loci, rec_rate);

        let wgy = wy.cross(kx);
        let wgx = wx.cross(ky);

        let wz = WGe::from_gametes(&wgx, &wgy);
        // If we have seen z before and we have now created it at a lower cost in way that is not
        // dependent on itself, then replace the parents of z previously seen with wx and wy
        if let Some(wz_prev) = h.get_mut(wz.genotype()) {
            if ancestors_contain(&wz) {
                continue;
            }
            // TODO: This can be its own function //
            if let Some((wgx_prev, wgy_prev)) = wz_prev.history() {
                let (wx_prev, wy_prev) = wz_prev.parents().unwrap();
                let cost_prev = cost_of_crossing(
                    gamma,
                    rec_rate.probability_gamete(wx_prev.genotype(), wgx_prev.gamete()),
                    rec_rate.probability_gamete(wy_prev.genotype(), wgy_prev.gamete()),
                );
                let cost_new = cost_of_crossing(
                    gamma,
                    rec_rate.probability_gamete(wx.genotype(), wgx.gamete()),
                    rec_rate.probability_gamete(wy.genotype(), wgy.gamete()),
                );
                if cost_new < cost_prev {
                    *wz_prev = wz;
                }
            }
        } else {
            h.insert(wz.genotype().clone(), wz.clone());
            pop.push(wz);
        }
    }
    let wz = h.get(&ideotype)?;
    Some(wz.to_base_sol(n_loci, crate::solution::Objective::Crossings))
}

/// Check if the ancestors of wz contain wz itself
///
/// # Arguments:
/// * `wz` - The genotype to check
/// * `h` - The hash map of genotypes to WGen
///
/// # Returns:
/// * `bool` - `true` if the ancestors of wz contain wz itself, `false` otherwise
fn ancestors_contain(wz: &WGe) -> bool {
    let mut stack: Vec<WGe> = vec![wz.clone()];
    let mut pushed: HashSet<SingleChromGenotype> = HashSet::new();

    while let Some(current) = stack.pop() {
        if let Some((wx, wy)) = current.parents() {
            for parent in &[wx, wy] {
                let genotype = parent.genotype().clone();
                if genotype == *wz.genotype() {
                    return true;
                }
                if !pushed.contains(&genotype) {
                    pushed.insert(genotype.clone());
                    stack.push(parent.clone());
                }
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ancestors_contain_test() {
        use crate::plants::bit_array::*;
        let wx = WGen::new(SingleChromGenotype::from_str("10", "10"));
        let wy = WGen::new(SingleChromGenotype::from_str("01", "01"));
        let gl = SingleChromGamete::from_str("10");
        let gr = SingleChromGamete::from_str("01");
        let wgl = WGam::new_from_genotype(gl, wx);
        let wgr = WGam::new_from_genotype(gr, wy);
        let wz = WGen::from_gametes(&wgl, &wgr);
        assert!(!ancestors_contain(&wz));
    }

    #[test]
    fn changeing_pointer_changes_subgraph() {
        use crate::plants::bit_array::*;
        let wx = WGen::new(SingleChromGenotype::from_str("10", "01"));
        let wy = WGen::new(SingleChromGenotype::from_str("01", "01"));
        let gl = SingleChromGamete::from_str("11");
        let gr = SingleChromGamete::from_str("01");
        let wgl = WGam::new_from_genotype(gl, wx);
        let wgr = WGam::new_from_genotype(gr, wy);
        let mut wz = WGen::from_gametes(&wgl, &wgr);

        let new_wx = WGen::new(SingleChromGenotype::from_str("11", "11"));
        let new_wgl = WGam::new_from_genotype(SingleChromGamete::from_str("11"), new_wx);
        let new_wz = WGen::from_gametes(&new_wgl, &wgr);

        assert_eq!(wz.genotype(), &SingleChromGenotype::from_str("11", "11"));
        assert_eq!(
            new_wz.genotype(),
            &SingleChromGenotype::from_str("11", "11")
        );

        wz = new_wz;
    }
}
