use crate::abstract_plants::{Crosspoint, Dominance, WGamS2, WGenS2};
use crate::plants::bit_array::*;
use crate::solution::{BaseSolution, Objective, PyBaseSolution};
use crate::solvers::base_min_generations_enumerator_dominance as bed;
use std::collections::HashSet;

/// Runs a breeding program given `n_loci` and `pop_0` where `pop_0` is a population of single
/// chromosome diploid genotypes with `n_loci` loci.
#[pyo3::pyfunction]
pub fn breeding_program_python(n_loci: usize, pop_0: Vec<Vec<Vec<bool>>>) -> PyBaseSolution {
    let pop_0 = pop_0
        .iter()
        .map(|x| {
            SingleChromGenotype::new(
                x[0].iter()
                    .zip(x[1].iter())
                    .map(|(a, b)| (*a, *b))
                    .collect(),
            )
        })
        .collect::<Vec<_>>();
    let res = breeding_program(n_loci, &pop_0);
    match res {
        None => Ok(None),
        Some(sol) => Ok(Some((
            sol.tree_data,
            sol.tree_type,
            sol.tree_left,
            sol.tree_right,
            sol.objective,
        ))),
    }
}

pub fn breeding_program(n_loci: usize, pop_0: &[SingleChromGenotype]) -> Option<BaseSolution> {
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    if pop_0.contains(&ideotype) {
        return Some(WGenS2::new(ideotype).to_base_sol(n_loci, Objective::Generations));
    }
    let pop_0: Vec<WGenS2<SingleChromGenotype, _>> =
        pop_0.iter().map(|x| WGenS2::new(x.clone())).collect();
    let mut h: HashSet<SingleChromGamete> = HashSet::new();
    let mut contained_gametes: Vec<WGamS2<SingleChromGenotype, SingleChromGamete>> = vec![];
    for wx in pop_0 {
        for gx in CrosspointBitVec::crosspoints(&n_loci).map(|k| k.cross(wx.genotype())) {
            if !h.contains(&gx) {
                let wgx = WGamS2::new_from_genotype(gx.clone(), wx.clone());
                h.insert(gx);
                contained_gametes.push(wgx);
            }
        }
    }
    contained_gametes = bed::filter_non_dominating_fn(contained_gametes, |wgx, wgy| {
        DomGamete::dom(wgx.gamete(), wgy.gamete())
    });
    let ideotype_gamete = SingleChromGamete::ideotype(n_loci);
    while !h.contains(&ideotype_gamete) {
        let mut v = vec![];
        for i in 0..contained_gametes.len() {
            for j in i + 1..contained_gametes.len() {
                let wgx = &contained_gametes[i];
                let wgy = &contained_gametes[j];
                let wz = WGenS2::from_gametes(wgx, wgy);
                for k in CrosspointBitVec::crosspoints(&n_loci) {
                    let gz = k.cross(wz.genotype());
                    if !h.contains(&gz) {
                        let wgz = WGamS2::new_from_genotype(gz.clone(), wz.clone());
                        h.insert(gz);
                        v.push(wgz);
                    }
                }
            }
        }
        contained_gametes =
            bed::filter_non_dominating_fn(v, |wgx, wgy| DomGamete::dom(wgx.gamete(), wgy.gamete()));
    }
    let wg_star = contained_gametes.first().expect("hash map doesn't have g*");
    let wx_star = WGenS2::from_gametes(wg_star, wg_star);
    Some(wx_star.to_base_sol(n_loci, Objective::Generations))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn enumerator_dominance_test() {
        use rand::prelude::*;
        let n_loci = 10;
        let n_pop = 5;
        let mut rng = thread_rng();
        for _ in 0..100 {
            let pop_0 = SingleChromGenotype::init_pop_random(&mut rng, n_loci, n_pop);
            let sol1 = breeding_program(n_loci, &pop_0);
            let sol2 = breeding_program(n_loci, &pop_0);
            assert!(sol1.is_some());
            assert!(sol2.is_some());
            assert_eq!(sol1.unwrap().objective, sol2.unwrap().objective);
        }
    }
}
