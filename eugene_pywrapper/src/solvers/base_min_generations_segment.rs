use eugene_core::plants::bit_array::SingleChromGenotype;
use eugene_core::solution::Objective;
use eugene_core::solvers::base_min_generations_segment;
use crate::solution::PyBaseSolution;
use pyo3::PyResult;
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
        let res = base_min_generations_segment::breeding_program(n_loci, &pop_0)
            .map(|x_star| x_star.to_base_sol(n_loci, Objective::Generations));
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

#[pyo3::pyfunction]
#[pyo3(name = "mingen_answer")]
pub fn mingen_answer_segment(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<bool>>>,
    timeout: Option<u64>,
) -> PyResult<Option<usize>> {
    //let instant = std::time::Instant::now();
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
    //println!("Rust's half of the serialisation: {:#?}", instant.elapsed());
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = base_min_generations_segment::min_generations(n_loci, &pop_0);
        tx.send(res)
    });
    let n_gen = rx.recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0));
    Ok(n_gen.ok())
}
