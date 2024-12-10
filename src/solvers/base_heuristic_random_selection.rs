use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solution::{BaseSolution, Objective, PyBaseSolution};
use crate::solvers::base_min_generations_enumerator_dominance::filter_non_dominating_key;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

/// Runs a breeding program given `n_loci` and `pop_0` where `pop_0` is a population of single
/// chromosome diploid genotypes with `n_loci` loci.
#[pyo3::pyfunction]
pub fn breeding_program_python(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<bool>>>,
    timeout: Option<u64>,
) -> PyBaseSolution {
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
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = breeding_program(n_loci, &pop_0);
        tx.send(res)
    });
    let res = rx
        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
        .ok()
        .flatten();
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

pub fn breeding_program(n_loci: usize, pop_0: &Vec<SingleChromGenotype>) -> Option<BaseSolution> {
    let mut rng = rand::thread_rng();
    breeding_program_random(n_loci, pop_0, Some(20), &mut rng)
}

pub fn breeding_program_random(
    n_loci: usize,
    pop_0: &Vec<SingleChromGenotype>,
    k_cross: Option<usize>,
    rng: impl rand::Rng,
) -> Option<BaseSolution> {
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    // TODO: check for infeasibility <18-07-24> //
    repeated_breeding_random(n_loci, pop_0, k_cross, rng)
        .iter()
        .find(|wx| wx.genotype() == &ideotype)
        .map(|wx| wx.to_base_sol(n_loci, Objective::Crossings))
}

pub fn repeated_breeding_random(
    n_loci: usize,
    pop_0: &[SingleChromGenotype],
    k_cross: Option<usize>,
    mut rng: impl rand::Rng,
) -> Vec<WGen<SingleChromGenotype, SingleChromGamete>> {
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    let mut pop: Vec<WGen<SingleChromGenotype, SingleChromGamete>> =
        pop_0.iter().cloned().map(WGen::new).collect();
    let mut n_crossings = 0;
    while !pop.iter().any(|wx| wx.genotype() == &ideotype)
        && k_cross.map_or(true, |k_cross| n_crossings < k_cross)
    {
        let n_pop = pop.len();
        let i = rng.gen_range(0..n_pop);
        let j = rng.gen_range(i..n_pop);
        let wx = &pop[i];
        let wy = &pop[j];
        let wz = WGen::from_gametes(
            &wx.cross(CrosspointBitVec::new(
                Chrom::from(rng.gen::<bool>()),
                rng.gen_range(0..n_loci),
            )),
            &wy.cross(CrosspointBitVec::new(
                Chrom::from(rng.gen::<bool>()),
                rng.gen_range(0..n_loci),
            )),
        );
        pop.push(wz);
        n_crossings += 1;
    }
    pop
}

pub fn breeding_program_random_dominance(
    n_loci: usize,
    pop_0: &Vec<SingleChromGenotype>,
    k_cross: Option<usize>,
    rng: impl rand::Rng,
) -> Option<BaseSolution> {
    // TODO: check for infeasibility <18-07-24> //
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    repeated_breeding_random_dominance(n_loci, pop_0, k_cross, rng)
        .iter()
        .find(|wx| wx.genotype() == &ideotype)
        .map(|wx| wx.to_base_sol(n_loci, Objective::Crossings))
}

pub fn repeated_breeding_random_dominance(
    n_loci: usize,
    pop_0: &Vec<SingleChromGenotype>,
    k_cross: Option<usize>,
    mut rng: impl rand::Rng,
) -> Vec<WGen<SingleChromGenotype, SingleChromGamete>> {
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    let mut pop: Vec<WGen<SingleChromGenotype, SingleChromGamete>> =
        pop_0.iter().cloned().map(WGen::new).collect();
    let mut n_crossings = 0;
    while !pop.iter().any(|wx| wx.genotype() == &ideotype)
        && k_cross.map_or(true, |k_cross| n_crossings < k_cross)
    {
        let non_dominating_wgametes = filter_non_dominating_key::<_, _, _, DomGamete>(
            pop.iter()
                .flat_map(|wx| CrosspointBitVec::crosspoints(&n_loci).map(|k| wx.cross(k))),
            |wg| wg.gamete(),
        );
        let i = rng.gen_range(0..non_dominating_wgametes.len());
        let j = rng.gen_range(i..non_dominating_wgametes.len());
        pop.push(WGen::from_gametes(
            &non_dominating_wgametes[i],
            &non_dominating_wgametes[j],
        ));
        n_crossings += 1;
    }
    pop
}
