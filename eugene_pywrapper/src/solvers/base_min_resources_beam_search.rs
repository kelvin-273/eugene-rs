use eugene_core::abstract_plants::{WGam, WGen};
use eugene_core::plants::bit_array::{SingleChromGamete, SingleChromGenotype};
use eugene_core::solvers::base_min_resources_beam_search;
use crate::solution::PyBaseSolution;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

#[derive(pyo3::IntoPyObject, Debug, Clone, Default)]
pub struct Output {
    expansions: usize,
    pushed_nodes: usize,
    children_created_by_branching: usize,
    objective: Option<usize>,
}

type WGe = WGen<SingleChromGenotype, SingleChromGamete>;
type WGa = WGam<SingleChromGenotype, SingleChromGamete>;

/// Runs a breeding program given `n_loci` and `pop_0` where `pop_0` is a population of single
/// chromosome diploid genotypes with `n_loci` loci.
#[pyo3::pyfunction]
pub fn breeding_program_python(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<bool>>>,
    recombination_rates: Vec<f64>,
    gamma: f64,
    timeout: Option<u64>,
) -> PyBaseSolution {
    assert_eq!(recombination_rates.len(), n_loci - 1);
    assert!(recombination_rates.iter().all(|&r| (0.0..=0.5).contains(&r)));
    assert!((0.0..=1.0).contains(&gamma));

    let pop_0: Vec<_> = pop_0
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
        let res = base_min_resources_beam_search::breeding_program(n_loci, &pop_0, &recombination_rates, gamma);
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
