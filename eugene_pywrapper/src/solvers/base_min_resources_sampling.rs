use eugene_core::extra::resources::RecRate;
use eugene_core::plants::bit_array::SingleChromGenotype;
use eugene_core::plants::dist_array::DistArray;
use eugene_core::solvers::base_min_resources_sampling::{
    breeding_program, breeding_program_distribute,
};
use pyo3::prelude::*;

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
) -> PyResult<
    Option<(
        Vec<Vec<Vec<i32>>>,
        Vec<&'static str>,
        Vec<usize>,
        Vec<usize>,
        usize,
    )>,
> {
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
    Ok(
        breeding_program(n_loci, &pop_0, &rec_rate, gamma).map(|res| {
            let (tree_data, tree_type, tree_left, tree_right) = (&res).into();
            (
                tree_data,
                tree_type,
                tree_left,
                tree_right,
                res.resources(&rec_rate, gamma),
            )
        }),
    )
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
) -> PyResult<
    Option<(
        Vec<Vec<Vec<i32>>>,
        Vec<&'static str>,
        Vec<usize>,
        Vec<usize>,
        usize,
    )>,
> {
    let rec_rate = RecRate::new(rec_rate);
    Ok(
        breeding_program_distribute(&DistArray::from(xs), &rec_rate, gamma).map(|res| {
            let (tree_data, tree_type, tree_left, tree_right) = (&res).into();
            (
                tree_data,
                tree_type,
                tree_left,
                tree_right,
                res.resources(&rec_rate, gamma),
            )
        }),
    )
}
