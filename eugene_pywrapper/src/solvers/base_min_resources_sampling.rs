use std::sync::mpsc;
use std::thread;
use std::time::Duration;

use eugene_core::extra::resources::RecRate;
use eugene_core::plants::bit_array::SingleChromGenotype;
use eugene_core::solution::CrossingSchedule;
use eugene_core::solvers::base_min_resources_sampling::{
    breeding_program, breeding_program_distribute,
};
use pyo3::prelude::*;

use crate::solution::PyCrossingSchedule;

/// Heuristically computes a crossing schedule that minimises the required number of resources by
/// repeatedly sampling crossings according to the distribution of offspring from a uniformly
/// chosen pair of parents given a number `n_loci` of loci, an initial population `pop_0`, a vector
/// `rec_rat` of recombination rates, and a threshold probability `gamma`.
#[pyfunction]
pub fn breeding_program_python(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<bool>>>,
    rec_rate: Vec<f64>,
    gamma: f64,
    timeout: Option<u64>,
) -> PyResult<Option<PyCrossingSchedule>> {
    let rec_rate = RecRate::new(rec_rate);
    let pop_0: Vec<SingleChromGenotype> = pop_0
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
        let res: Option<CrossingSchedule> = breeding_program(n_loci, &pop_0, &rec_rate, gamma);
        tx.send(res)
    });
    let res = rx
        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
        .ok()
        .flatten();
    Ok(res.map(PyCrossingSchedule::new))
}

/// Heuristically computes a crossing schedule that minimises the required number of resources by
/// repeatedly sampling crossings according to the distribution of offspring from a uniformly
/// chosen pair of parents given a distribute array `xs`, a vector `rec_rate` of recombination
/// rates, and a threshold probability `gamma`.
#[pyfunction]
pub fn breeding_program_distribute_python(
    xs: Vec<usize>,
    rec_rate: Vec<f64>,
    gamma: f64,
    timeout: Option<u64>,
) -> PyResult<Option<PyCrossingSchedule>> {
    let rec_rate = RecRate::new(rec_rate);
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res: Option<CrossingSchedule> =
            breeding_program_distribute(&xs.into(), &rec_rate, gamma);
        tx.send(res)
    });
    let res = rx
        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
        .ok()
        .flatten();
    Ok(res.map(PyCrossingSchedule::new))
}
