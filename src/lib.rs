#![allow(unused)]

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

/// Exports the traits used by the algorithms in `solvers`
pub mod abstract_plants;
pub mod plants;
pub mod solvers;
pub mod visualisation;

mod play {
    #[cfg(test)]
    mod tests {}
}
