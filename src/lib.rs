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
    let mg_naive2 = PyModule::new_bound(m.py(), "naive2")?;
    let mg_segment = PyModule::new_bound(m.py(), "segment")?;

    mg_naive1.add_function(wrap_pyfunction!(
        solvers::base_min_generations_enumerator::breeding_program_python,
        &mg_naive1
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

#[cfg(test)]
mod tests {
    use crate::extra::*;
    use crate::solvers::*;
    use rand::prelude::*;
    use std::env;
    use std::io;

    fn _main() -> io::Result<()> {
        env::set_var("RUST_BACKTRACE", "1");
        _main7()
    }

    fn _main7() -> io::Result<()> {
        let mut lines = io::stdin().lines();
        while let Some(Ok(s)) = lines.next() {
            println!("string: {}", &s);
            let (n_loci, pop0) = instance_generators::parse_homozygous(&s);
            dbg!(&n_loci);
            dbg!(&pop0);
            let t = base_min_crossings_astar::breeding_program(n_loci, &pop0).unwrap();
            let n_seg = base_min_generations_segment::min_covering_segments(n_loci, &pop0).len();
            let n_crs = t.crossings();
            dbg!(&t);
            println!("Segments: {}\tCrossings: {}", n_seg, n_crs);
        }
        Ok(())
    }

    fn _main6() -> io::Result<()> {
        use crate::plants::bit_array::*;
        let mut rng = thread_rng();
        let n_loci = 18;
        let n_pop = 6;
        for _ in 0..1 {
            let pop0 = SingleChromGenotype::init_pop_random(&mut rng, n_loci, n_pop);
            dbg!(&pop0);
            let wt = base_min_generations_segment::breeding_program(n_loci, &pop0).unwrap();
            let n_seg = base_min_generations_segment::min_covering_segments(n_loci, &pop0).len();
            let n_crs = wt.crossings();
            dbg!(&wt);
            println!("Segments: {}\tCrossings: {}", n_seg, n_crs);
        }
        Ok(())
    }

    fn _main5() -> io::Result<()> {
        use crate::plants::bit_array::*;
        let n_loci = 16;
        let n_pop = 2;
        let _res = base_min_crossings_astar::breeding_program(
            n_loci,
            &SingleChromGenotype::init_pop_random(&mut thread_rng(), n_loci, n_pop),
        );
        Ok(())
    }

    fn _main4() -> io::Result<()> {
        use crate::plants::bit_array::*;
        let n_loci = 20;
        let n_pop = 1024;
        let _res = base_min_generations_segment::breeding_program(
            n_loci,
            &SingleChromGenotype::init_pop_random(&mut thread_rng(), n_loci, n_pop),
        );
        Ok(())
    }

    fn _main1() -> io::Result<()> {
        use crate::plants::bit_array::*;
        let n_loci = 8;
        let pop0: Vec<SingleChromGenotype> =
            vec![SingleChromGenotype::from_str("01010101", "10101010")];
        if let Some(res) =
            base_min_generations_enumerator_dominance::breeding_program(n_loci, &pop0)
        {
            print!("breeding_program successful");

            visualisation::draw_tree_genotypes::<SingleChromGenotype, SingleChromGamete>(&res);
        }
        Ok(())
    }

    fn _main0() -> io::Result<()> {
        use crate::plants::bit_array::*;
        let n_loci = 10;
        let pop0: Vec<SingleChromGenotype> = vec![
            SingleChromGenotype::from_str("0000011000", "0000000000"),
            SingleChromGenotype::from_str("0000000100", "0000000000"),
            SingleChromGenotype::from_str("0000000101", "0000000010"),
            SingleChromGenotype::from_str("1100000000", "0000000000"),
            SingleChromGenotype::from_str("0010000000", "0000000000"),
            SingleChromGenotype::from_str("0010100000", "0001000000"),
        ];
        if let Some(res) =
            base_min_generations_enumerator_dominance::breeding_program(n_loci, &pop0)
        {
            print!("breeding_program successful");

            visualisation::draw_tree_genotypes::<SingleChromGenotype, SingleChromGamete>(&res);
        }
        Ok(())
    }
}
