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

/// Exports the traits used by the algorithms in `solvers`
pub mod abstract_plants;
/// Exports the utility functions for parsing and generating instances
pub mod extra;
/// Exports the data structures for plant genotypes and gametes
pub mod plants;
/// Exports the type for representing breeding programs
pub mod solution;
/// Exports the algorithms for solving the plant breeding program problem
pub mod solvers;

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

            let _ = visualisation::draw_tree_genotypes::<SingleChromGenotype, SingleChromGamete>(&res);
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

            let _ = visualisation::draw_tree_genotypes::<SingleChromGenotype, SingleChromGamete>(&res);
        }
        Ok(())
    }
}
