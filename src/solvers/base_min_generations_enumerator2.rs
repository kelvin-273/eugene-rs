use crate::abstract_plants::{Crosspoint, Genotype, WGamS, WGenS};
use crate::plants::bit_array::*;
use crate::solution::BaseSolution;
use std::collections::HashMap;
use std::rc::Rc;
use pyo3::prelude::*;

/// Runs a breeding program given `n_loci` and `pop_0` where `pop_0` is a population of single
/// chromosome diploid genotypes with `n_loci` loci.
#[pyo3::pyfunction]
pub fn breeding_program_python(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<bool>>>,
) -> PyResult<
    Option<(
        Vec<Vec<Vec<i32>>>,
        Vec<&'static str>,
        Vec<usize>,
        Vec<usize>,
        usize,
    )>,
> {
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
        .collect();
    let res = breeding_program(n_loci, &pop_0);
    match res {
        None => Ok(None),
        Some(sol) => {
            Ok(Some((
                sol.tree_data,
                sol.tree_type,
                sol.tree_left,
                sol.tree_right,
                sol.objective,
            )))
        }
    }
}

pub fn breeding_program(n_loci: usize, pop_0: &Vec<SingleChromGenotype>) -> Option<BaseSolution> {
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    if pop_0.contains(&ideotype) {
        return Some(BaseSolution::min_gen_from_wgens(n_loci, &WGenS::new(ideotype)));
    }
    let pop_0: Vec<Rc<WGenS<SingleChromGenotype, _>>> = pop_0
        .iter()
        .map(|x| Rc::new(WGenS::new(x.clone())))
        .collect();
    let mut h: HashMap<
        SingleChromGamete,
        Rc<WGamS<SingleChromGenotype, SingleChromGamete>>,
    > = HashMap::new();
    for wx in pop_0 {
        for gx in CrosspointBitVec::crosspoints(&n_loci).map(|k| k.cross(&wx.genotype)) {
            if !h.contains_key(&gx) {
                let wgx = Rc::new(WGamS {
                    gamete: gx.clone(),
                    history: Some(wx.clone()),
                });
                h.insert(gx, wgx);
            }
        }
    }
    let ideotype_gamete = SingleChromGamete::ideotype(n_loci);
    while !h.contains_key(&ideotype_gamete) {
        let contained_gametes: Vec<SingleChromGamete> = h.keys().map(|k| k.clone()).collect();
        for i in 0..contained_gametes.len() {
            for j in i+1..contained_gametes.len() {
                let wgx = h.get(&contained_gametes[i])?.clone();
                let wgy = h.get(&contained_gametes[j])?.clone();
                let z = SingleChromGenotype::from_gametes(&wgx.gamete, &wgy.gamete);
                let wz = Rc::new(WGenS { genotype: z, history: Some((wgx, wgy)) });
                for k in CrosspointBitVec::crosspoints(&n_loci) {
                    let gz = k.cross(&wz.genotype);
                    if !h.contains_key(&gz) {
                        let wgz = Rc::new(WGamS {
                            gamete: gz.clone(),
                            history: Some(wz.clone()),
                        });
                        h.insert(gz, wgz);
                    }
                }
            }
        }
    }
    let mut wx_star = WGenS::new(ideotype);
    let wgz = h.get(&ideotype_gamete)?;
    wx_star.history = Some((wgz.clone(), wgz.clone()));
    Some(BaseSolution::min_gen_from_wgens(n_loci, &wx_star))
}
