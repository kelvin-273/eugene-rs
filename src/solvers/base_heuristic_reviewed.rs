use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solution::{BaseSolution, Objective, PyBaseSolution};
use crate::solvers::base_min_generations_enumerator_dominance::filter_non_dominating_key;
use rand::prelude::*;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

///// Runs a breeding program given `n_loci` and `pop_0` where `pop_0` is a population of single
///// chromosome diploid genotypes with `n_loci` loci.
//#[pyo3::pyfunction]
//pub fn breeding_program_python(
//    n_loci: usize,
//    pop_0: Vec<Vec<Vec<bool>>>,
//    timeout: Option<u64>,
//) -> PyBaseSolution {
//    let pop_0 = pop_0
//        .iter()
//        .map(|x| {
//            SingleChromGenotype::new(
//                x[0].iter()
//                    .zip(x[1].iter())
//                    .map(|(a, b)| (*a, *b))
//                    .collect(),
//            )
//        })
//        .collect();
//    let (tx, rx) = mpsc::channel();
//    thread::spawn(move || {
//        let res = breeding_program(n_loci, &pop_0);
//        tx.send(res)
//    });
//    let res = rx
//        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
//        .ok()
//        .flatten();
//    match res {
//        None => Ok(None),
//        Some(sol) => Ok(Some((
//            sol.tree_data,
//            sol.tree_type,
//            sol.tree_left,
//            sol.tree_right,
//            sol.objective,
//        ))),
//    }
//}

pub fn breeding_program(
    n_loci: usize,
    pop_0: &Vec<SingleChromGenotype>,
    marker_mapping: &Vec<f64>,
    recomb_rates: &Vec<f64>,
    max_generations: usize,
    children_per_crossing: usize,
) -> Option<BaseSolution> {
    assert_eq!(marker_mapping.len(), n_loci);
    assert_eq!(recomb_rates.len(), n_loci - 1);
    for _ in 0..max_generations {
        todo!("choose parent");
        todo!("create children");
        todo!("profit");
    }
    None
}

type WGe = WGen<SingleChromGenotype, SingleChromGamete>;

pub fn trait_introgression_with_backcrossing(
    n_loci: usize,
    elite: &SingleChromGenotype,
    donor: &SingleChromGenotype,
    recomb_rates: &Vec<f64>,
    marker_mapping: &Vec<f64>,
    max_generations: usize,
    children_per_crossing: usize,
    rng: &mut impl Rng,
) -> BaseSolution {
    assert_eq!(marker_mapping.len(), n_loci);
    assert_eq!(recomb_rates.len(), n_loci - 1);
    assert!(recomb_rates.iter().all(|r| 0.0 <= *r && *r <= 1.0));
    let w_donor: WGe = WGen::new(donor.clone());
    let w_elite: WGe = WGen::new(elite.clone());
    let f = |wx: &WGe| {
        let upper_value: f64 = wx
            .genotype()
            .upper()
            .alleles()
            .iter()
            .zip(marker_mapping.iter())
            .map(|(a, c)| *c * Into::<bool>::into(*a) as i64 as f64)
            .sum();
        let lower_value: f64 = wx
            .genotype()
            .upper()
            .alleles()
            .iter()
            .zip(marker_mapping.iter())
            .map(|(a, c)| *c * Into::<bool>::into(*a) as i64 as f64)
            .sum();
        upper_value + lower_value
    };
    let mut cross_recw = |wx: &WGe, wy: &WGe| {
        WGen::from_gametes(
            &wx.cross(CrosspointBitVec::random_crosspoint(
                rng,
                &n_loci,
                recomb_rates,
            )),
            &wy.cross(CrosspointBitVec::random_crosspoint(
                rng,
                &n_loci,
                recomb_rates,
            )),
        )
    };
    let mut chosen_child = cross_recw(&w_elite, &w_donor);
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    let mut generation = 1;
    while generation < max_generations && chosen_child.genotype() != &ideotype {
        if let Some(max_child) = (0..children_per_crossing)
            .map(|_| cross_recw(&w_elite, &chosen_child))
            .max_by(|wx, wy| {
                f(&wx)
                    .partial_cmp(&f(&wy))
                    .expect("ranking function yielded a NaN")
            })
        {
            if f(&chosen_child) < f(&max_child) {
                chosen_child = max_child;
            }
        }
        generation += 1;
    }
    chosen_child.to_base_sol(n_loci, Objective::Crossings)
}
