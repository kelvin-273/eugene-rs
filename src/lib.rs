//! Plant breeding program modelling and solving library.
//!
//! Components are of the plant breeding programs are separated modules.
//! Typical use is as follows:
//! ```
//! use eugene::plants::bit_array::*;
//! use eugene::solvers::greedy_base::*;
//!
//! use rand::prelude::*;
//! let mut rng = thread_rng();
//!
//! let n_loci = 10;
//! let n_pop = 6;
//! let pop_0 = SingleChromGenotype::init_pop_random(&mut rng, n_loci, n_pop);
//! let ideotype = SingleChromGenotype::ideotype(n_loci);
//!
//! let result = breeding_program::<
//!     SingleChromGenotype,
//!     SingleChromGamete,
//!     CrosspointBitVec,
//!     SegmentBitVec
//!     >(n_loci, pop_0, ideotype);
//! ```

use pyo3::prelude::*;

/// Exports the traits used by the algorithms in `solvers`
pub mod abstract_plants;
pub mod extra;
pub mod plants;
pub mod solution;
pub mod solvers;

mod play {
    #[cfg(test)]
    mod tests {}
}

#[pymodule]
//#[pyo3(name = "eugene_rs")]
fn eugene_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let min_gen = PyModule::new_bound(m.py(), "min_gen")?;
    let min_cross = PyModule::new_bound(m.py(), "min_cross")?;

    let mg_naive1 = PyModule::new_bound(m.py(), "naive1")?;
    let mg_naive12 = PyModule::new_bound(m.py(), "naive12")?;
    let mg_naive2 = PyModule::new_bound(m.py(), "naive2")?;
    let mg_segment = PyModule::new_bound(m.py(), "segment")?;

    mg_naive1.add_function(wrap_pyfunction!(
        solvers::base_min_generations_enumerator::breeding_program_python,
        &mg_naive1
    )?)?;

    mg_naive12.add_function(wrap_pyfunction!(
        solvers::base_min_generations_enumerator2::breeding_program_python,
        &mg_naive12
    )?)?;

    mg_naive2.add_function(wrap_pyfunction!(
        solvers::base_min_generations_enumerator_dominance::breeding_program_python,
        &mg_naive2
    )?)?;

    mg_segment.add_function(wrap_pyfunction!(
        solvers::base_min_generations_segment::breeding_program_python,
        &mg_segment
    )?)?;

    min_gen.add_submodule(&mg_naive1)?;
    min_gen.add_submodule(&mg_naive12)?;
    min_gen.add_submodule(&mg_naive2)?;
    min_gen.add_submodule(&mg_segment)?;

    let mc_astar = PyModule::new_bound(m.py(), "astar")?;
    //let mc_mip = PyModule::new_bound(m.py(), "mip")?;

    mc_astar.add_function(wrap_pyfunction!(
        solvers::base_min_crossings_astar::breeding_program_python,
        &mc_astar
    )?)?;

    //mc_mip.add_function(wrap_pyfunction!(
    //    solvers::base_min_crossings_mip::breeding_program_python_mip,
    //    &mc_mip
    //)?)?;

    min_cross.add_submodule(&mc_astar)?;
    //min_cross.add_submodule(&mc_mip)?;

    m.add_submodule(&min_gen)?;
    m.add_submodule(&min_cross)?;
    Ok(())
}
