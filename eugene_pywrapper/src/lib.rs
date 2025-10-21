//! Plant breeding program modelling and solving library.
//!
//! Components of the plant breeding program are separated by module.
//! Typical use is as follows:
//! ```
//! use eugene::plants::bit_array::*;
//! use eugene::solvers::base_min_generations_segment::*;
//!
//! use rand::prelude::*;
//! let mut rng = thread_rng();
//!
//! let n_loci = 10;
//! let n_pop = 6;
//! let pop_0 = SingleChromGenotype::init_pop_random(&mut rng, n_loci, n_pop);
//!
//! let result = breeding_program(n_loci, pop_0);
//! ```
//#![deny(missing_docs)]

use pyo3::prelude::*;

/// Exports the type for representing breeding programs
pub mod solution;
/// Exports the algorithms for solving the plant breeding program problem
pub mod solvers;

#[pymodule]
#[pyo3(name = "eugene_pywrapper")]
fn eugene_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let min_gen = PyModule::new(m.py(), "min_gen")?;
    let min_cross = PyModule::new(m.py(), "min_cross")?;

    let mg_naive1 = PyModule::new(m.py(), "naive1")?;
    let mg_naive2 = PyModule::new(m.py(), "naive2")?;
    let mg_segment = PyModule::new(m.py(), "segment")?;

    mg_naive1.add_function(wrap_pyfunction!(
        solvers::base_min_generations_enumerator::breeding_program_python,
        &mg_naive1
    )?)?;

    mg_naive1.add_function(wrap_pyfunction!(
        solvers::base_min_generations_enumerator::mingen_answer_enumerator,
        &mg_naive1
    )?)?;

    mg_naive2.add_function(wrap_pyfunction!(
        solvers::base_min_generations_enumerator_dominance::breeding_program_python,
        &mg_naive2
    )?)?;

    mg_naive2.add_function(wrap_pyfunction!(
        solvers::base_min_generations_enumerator_dominance::mingen_answer_dominance,
        &mg_naive2
    )?)?;

    mg_segment.add_function(wrap_pyfunction!(
        solvers::base_min_generations_segment::breeding_program_python,
        &mg_segment
    )?)?;

    mg_segment.add_function(wrap_pyfunction!(
        solvers::base_min_generations_segment::mingen_answer_segment,
        &mg_segment
    )?)?;

    min_gen.add_submodule(&mg_naive1)?;
    min_gen.add_submodule(&mg_naive2)?;
    min_gen.add_submodule(&mg_segment)?;

    let mc_astar = PyModule::new(m.py(), "astar")?;
    //let mc_mip = PyModule::new_bound(m.py(), "mip")?;
    let mc_distribute_astar = PyModule::new(m.py(), "distribute_astar")?;

    mc_astar.add_function(wrap_pyfunction!(
        solvers::base_min_crossings_astar::breeding_program_python,
        &mc_astar
    )?)?;

    mc_distribute_astar.add_function(wrap_pyfunction!(
        solvers::base_min_crossings_distribute_astar::breeding_program_distribute_general_python,
        &mc_distribute_astar
    )?)?;

    mc_distribute_astar.add_function(wrap_pyfunction!(
        solvers::base_min_crossings_distribute_astar::breeding_program_distribute_python,
        &mc_distribute_astar
    )?)?;

    mc_distribute_astar.add_function(wrap_pyfunction!(
        solvers::base_min_crossings_distribute_astar::breeding_program_distribute_diving_python,
        &mc_distribute_astar
    )?)?;

    mc_distribute_astar.add_function(wrap_pyfunction!(
        solvers::base_min_crossings_distribute_astar::breeding_program_distribute_no_full_join_python,
        &mc_distribute_astar
    )?)?;

    //mc_mip.add_function(wrap_pyfunction!(
    //    solvers::base_min_crossings_mip::breeding_program_python_mip,
    //    &mc_mip
    //)?)?;

    min_cross.add_submodule(&mc_astar)?;
    //min_cross.add_submodule(&mc_mip)?;
    min_cross.add_submodule(&mc_distribute_astar)?;

    m.add_submodule(&min_gen)?;
    m.add_submodule(&min_cross)?;
    Ok(())
}
