use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solution::*;
use crate::solvers::base_min_generations_segment;
use std::rc::Rc;
use std::sync::mpsc;
use std::{panic, thread};
use std::time::Duration;
use std::collections::HashMap;

struct MultiChromGenotype {
    arr: Vec<SingleChromGenotype>,
}

impl MultiChromGenotype {
    /// Returns the number of chromosomes in this genotype.
    ///
    /// ```
    /// let genotype = MultiChromGenotype {
    ///     arr: vec![
    ///         SingleChromGenotype::from_str("00111000110", "00111000110"),
    ///         SingleChromGenotype::from_str("11000111000", "11000111000"),
    ///     ],
    /// };
    /// assert_eq!(2, genotype.n_chrom());
    /// ```
    #[inline]
    fn n_chrom(&self) -> usize {
        self.arr.len()
    }
}

impl WGen<MultiChromGenotype, MultiChromGamete> {
    fn to_base_sol(&self, n_loci: usize, objective: Objective) -> BaseSolution {
        unimplemented!()
    }
}

struct MultiChromGamete {
    arr: Vec<SingleChromGamete>,
}

type WGe = WGen<MultiChromGenotype, MultiChromGamete>;
type WGa = WGam<MultiChromGenotype, MultiChromGamete>;

/// A triple containing a start-point `s`, end-point `e`, and gamete `g` used to represent a
/// contiguous sequence of favourable alleles on `g` from index `s` to `e` inclusive.
#[derive(Debug, Clone)]
pub struct Segment<B> {
    /// Start point
    s: usize,
    /// End point
    e: usize,
    /// Gamete
    g: B,
}

impl<B> Segment<B> {
    /// Returns the start point of the segment.
    ///
    /// ```
    /// let seg = Segment { s: 2, e, 4, g: SingleChromGamete::from_str("00111000110") };
    /// assert_eq!(2, seg.start())
    /// ```
    pub fn start(&self) -> usize {
        self.s
    }

    /// Returns the end point of the segment.
    ///
    /// ```
    /// let seg = Segment { s: 2, e, 4, g: SingleChromGamete::from_str("00111000110") };
    /// assert_eq!(4, seg.start())
    /// ```
    pub fn end(&self) -> usize {
        self.e
    }

    pub fn gamete(&self) -> &B {
        &self.g
    }
}

type Seg = Segment<SingleChromGamete>;
type SegW = Segment<WGam<SingleChromGenotype, SingleChromGamete>>;

impl Seg {
    pub fn join(&self, other: &Self) -> Self {
        assert!(self.s < other.s);
        // TODO: complete adjacency condition <27-05-24> //
        Self {
            s: self.s,
            e: other.e,
            g: {
                CrosspointBitVec::new(Chrom::Upper, other.s)
                    .cross(&SingleChromGenotype::from_gametes(&self.g, &other.g))
            },
        }
    }
}

impl SegW {
    pub fn join(&self, other: &Self) -> Self {
        assert!(self.s < other.s);
        // TODO: complete adjacency condition <27-05-24> //
        Self {
            s: self.s,
            e: other.e,
            g: {
                let crosspoint_bit_vec = &CrosspointBitVec::new(Chrom::Upper, other.s);
                let wz = WGen::from_gametes(&self.g, &other.g);
                wz.cross_fn(|z| crosspoint_bit_vec.cross(z))
            },
        }
    }

    pub fn join_with_store(&self, other: &Self, h: &mut HashMap<SingleChromGamete, WGa>) -> Self {
        assert!(self.s < other.s);
        // TODO: complete adjacency condition <27-05-24> //
        Self {
            s: self.s,
            e: other.e,
            g: {
                let z = SingleChromGenotype::from_gametes(self.g.gamete(), other.g.gamete());
                let wz: WGen<SingleChromGenotype, SingleChromGamete> = WGen::new(z);
                let crosspoint_bit_vec = &CrosspointBitVec::new(Chrom::Upper, other.s);
                let gz = crosspoint_bit_vec.cross(wz.genotype());
                todo!()
            },
        }
    }
}

enum Segment2<B> {
    SegNone(B, usize),
    SegSome(B, usize, usize, usize)
}

impl<B> Segment2<B> {
    pub fn chromosome(&self) -> usize {
        match self {
            Segment2::SegNone(_, chrom) => *chrom,
            Segment2::SegSome(_, chrom, _, _) => *chrom,
        }
    }
    pub fn gamete(&self) -> &B {
        match self {
            Segment2::SegNone(g, _) => g,
            Segment2::SegSome(g, _, _, _) => g,
        }
    }
}

/// Runs a breeding program given `n_loci` and `pop_0` where `pop_0` is a population of single
/// chromosome diploid genotypes with `n_loci` loci.
#[pyo3::pyfunction]
pub fn breeding_program_python(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<bool>>>,
    timeout: Option<u64>,
) -> PyBaseSolution {
    let pop_0 = todo!("Convert pop_0 to Vec<MultiChromGenotype> from Vec<Vec<Vec<bool>>>");
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = breeding_program(&pop_0)
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

fn breeding_program(pop_0: &Vec<MultiChromGenotype>) -> Option<WGe> {
    unimplemented!()
}

fn is_ti_instance(pop_0: &Vec<MultiChromGenotype>) -> bool {
    if pop_0.len() != 2 {
        return false;
    }
    let first = &pop_0[0];
    if pop_0.iter().any(|x| x.n_chrom() != first.n_chrom()) {
        return false;
    }
    if pop_0.iter().any(|x| x.arr.len() != x.n_chrom()) {
        panic!("All genotypes must have the same number of chromosomes");
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multi_chrom_min_generations_infeasible() {
        let pop_0 = vec![
            MultiChromGenotype {
                arr: vec![
                    SingleChromGenotype::from_str("00111000110", "00111000110"),
                    SingleChromGenotype::from_str("00111000110", "00111000110"),
                ],
            },
            MultiChromGenotype {
                arr: vec![
                    SingleChromGenotype::from_str("11000111000", "11000111000"),
                    SingleChromGenotype::from_str("11000111000", "11000111000"),
                ],
            },
        ];
        let res = breeding_program(&pop_0);
        assert!(res.is_none());
    }

    #[test]
    fn test_multi_chrom_min_generations_feasible() {
        let pop_0 = vec![
            MultiChromGenotype {
                arr: vec![
                    SingleChromGenotype::from_str("00111000110", "00111000111"),
                    SingleChromGenotype::from_str("00111000110", "00111000111"),
                ],
            },
            MultiChromGenotype {
                arr: vec![
                    SingleChromGenotype::from_str("11000111000", "11000111000"),
                    SingleChromGenotype::from_str("11000111000", "11000111000"),
                ],
            },
        ];
        let res = breeding_program(&pop_0);
        assert!(res.is_none());
    }

}
