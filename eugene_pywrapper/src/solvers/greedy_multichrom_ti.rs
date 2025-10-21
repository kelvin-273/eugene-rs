use eugene_core::plants::bit_array::SingleChromGenotype;
use eugene_core::solvers::greedy_multichrom_ti;
use eugene_core::solution::Objective;
use crate::solution::PyBaseSolution;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

/// Runs a breeding program given `n_loci` and `pop_0` where `pop_0` is a population of single
/// chromosome diploid genotypes with `n_loci` loci.
#[pyo3::pyfunction]
pub fn breeding_program_python(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<(bool, bool)>>>,
    timeout: Option<u64>,
) -> PyBaseSolution {
    let pop_0: Vec<_> = pop_0.iter().map(|v| {
        let single_chrom_genotypes: Vec<SingleChromGenotype> = v.iter().map(|chrom| {
            SingleChromGenotype::from_str(
                &chrom.iter().map(|(a, _)| if *a { '1' } else { '0' }).collect::<String>(),
                &chrom.iter().map(|(_, b)| if *b { '1' } else { '0' }).collect::<String>(),
            )
        }).collect();
        greedy_multichrom_ti::MultiChromGenotype::new(single_chrom_genotypes)
    }).collect();
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = greedy_multichrom_ti::breeding_program(&pop_0)
            .map(|x_star| x_star.to_base_sol(n_loci, Objective::Crossings));
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
