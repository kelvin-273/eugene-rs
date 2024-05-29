use crate::abstract_plants::{Crosspoint, WGamS2, WGenS2};
use crate::plants::bit_array::*;
use crate::solution::{BaseSolution, Objective, PyBaseSolution};
use std::collections::HashMap;

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
    let mut h: HashMap<SingleChromGamete, WGamS2<SingleChromGenotype, SingleChromGamete>> =
        HashMap::new();
    for wx in pop_0 {
        for gx in CrosspointBitVec::crosspoints(&n_loci).map(|k| k.cross(wx.genotype())) {
            if !h.contains_key(&gx) {
                let wgx = WGamS2::new_from_genotype(gx.clone(), wx.clone());
                h.insert(gx, wgx);
            }
        }
    }
    let ideotype_gamete = SingleChromGamete::ideotype(n_loci);
    while !h.contains_key(&ideotype_gamete) {
        let contained_gametes: Vec<SingleChromGamete> = h.keys().cloned().collect();
        for i in 0..contained_gametes.len() {
            for j in i + 1..contained_gametes.len() {
                let wgx = h.get(&contained_gametes[i])?;
                let wgy = h.get(&contained_gametes[j])?;
                let wz = WGenS2::from_gametes(wgx, wgy);
                for k in CrosspointBitVec::crosspoints(&n_loci) {
                    let gz = k.cross(wz.genotype());
                    if !h.contains_key(&gz) {
                        let wgz = WGamS2::new_from_genotype(gz.clone(), wz.clone());
                        h.insert(gz, wgz.clone());
                    }
                }
            }
        }
    }
    let wg_star = h.get(&ideotype_gamete).expect("hash map doesn't have g*");
    let wx_star = WGenS2::from_gametes(wg_star, wg_star);
    Some(wx_star.to_base_sol(n_loci, Objective::Generations))
}
