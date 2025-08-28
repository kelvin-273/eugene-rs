use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solution::*;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;
use std::collections::HashMap;

pub struct MultiChromGenotype {
    arr: Vec<SingleChromGenotype>,
}

impl MultiChromGenotype {
    /// Creates a new `MultiChromGenotype` from a vector of `SingleChromGenotype`.
    /// ```
    /// let genotype = MultiChromGenotype {
    ///     arr: vec![
    ///         SingleChromGenotype::from_str("00111000110", "00111000110"),
    ///         SingleChromGenotype::from_str("10111000", "11011000"),
    ///     ],
    /// };
    /// assert_eq!(2, genotype.arr.len());
    /// ```
    pub fn new(arr: Vec<SingleChromGenotype>) -> Self {
        Self { arr }
    }

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
    fn to_base_sol(&self, _n_loci: usize, _objective: Objective) -> BaseSolution {
        unimplemented!()
    }
}

pub struct MultiChromGamete {
    _arr: Vec<SingleChromGamete>,
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
        let result = Self {
            s: self.s,
            e: other.e,
            g: {
                let crosspoint_bit_vec = &CrosspointBitVec::new(Chrom::Upper, other.s);
                let wz = WGen::from_gametes(&self.g, &other.g);
                wz.cross_fn(|z| crosspoint_bit_vec.cross(z))
            },
        };
        // Store the result in the hash map
        if !h.contains_key(result.g.gamete()) {
            todo!("Insert new segment into hash map");
        }
        result
    }
}

enum _Segment2<B> {
    SegNone(B, usize),
    SegSome(B, usize, usize, usize)
}

impl<B> _Segment2<B> {
    fn _chromosome(&self) -> usize {
        match self {
            _Segment2::SegNone(_, chrom) => *chrom,
            _Segment2::SegSome(_, chrom, _, _) => *chrom,
        }
    }
    fn _gamete(&self) -> &B {
        match self {
            _Segment2::SegNone(g, _) => g,
            _Segment2::SegSome(g, _, _, _) => g,
        }
    }
}

/// Runs a breeding program given `n_loci` and `pop_0` where `pop_0` is a population of single
/// chromosome diploid genotypes with `n_loci` loci.
#[pyo3::pyfunction]
pub fn breeding_program_python(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<(bool, bool)>>>,
    timeout: Option<u64>,
) -> PyBaseSolution {
    let pop_0 = pop_0.iter().map(|v| {
        let single_chrom_genotypes: Vec<SingleChromGenotype> = v.iter().map(|chrom| {
            SingleChromGenotype::from_str(
                &chrom.iter().map(|(a, _)| if *a { '1' } else { '0' }).collect::<String>(),
                &chrom.iter().map(|(_, b)| if *b { '1' } else { '0' }).collect::<String>(),
            )
        }).collect();
        MultiChromGenotype::new(single_chrom_genotypes)
    }).collect();
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
    // Check if the instance is a valid TI instance
    if !is_ti_instance(pop_0) {
        return None;
    }

    // Construct a set for each chromosome that contains the segments between the two genotypes
    let mut segments: Vec<Vec<SegW>> = Vec::with_capacity(pop_0[0].n_chrom());
    for i in 0..pop_0[0].n_chrom() {
        let elite_chrom = &pop_0[0].arr[i];
        let donor_chrom = &pop_0[1].arr[i];
        let mut segs: Vec<SegW> = Vec::new();
        // Find segments in the elite chromosome that are not in the donor chromosome
        // We can assume that elite_allele is the opposite of donor_allele at each locus
        let mut _start = 0;
        for j in 0..elite_chrom.n_loci() {
            let _elite_allele = elite_chrom.get(true, j).unwrap();
            let _donor_allele = donor_chrom.get(true, j).unwrap();
            // Start a new segment if we are at the start or the previous segment ended
            if segs.is_empty() || segs.last().unwrap().end() < j {
                todo!("Push new segment");
            }
            // Extend the last segment
            let last_seg = segs.last_mut().unwrap();
            last_seg.e = j;
        }
        segs.retain(|seg| seg.e >= seg.s); // Remove empty segments
        segs.sort_by_key(|seg| seg.s); // Sort segments by start point
        segments.push(segs);
    }
    unimplemented!()
}

/// Returns true if the given population is a valid TI instance.
/// A valid TI instance must have exactly two genotypes, each with the same number of chromosomes,
/// all chromosomes must be of the same length, and all chromosomes must be homozygous.
/// Additionally, the genotypes must be complementary, meaning that for each chromosome,
/// the alleles in the first genotype must be the opposite of the alleles in the second genotype.
/// Note, that any valid TI is also feasible.
fn is_ti_instance(pop_0: &Vec<MultiChromGenotype>) -> bool {
    if pop_0.len() != 2 {
        return false;
    }
    let (elite, donor) = (&pop_0[0], &pop_0[1]);
    if elite.n_chrom() != donor.n_chrom() {
        return false;
    }
    // Check if all chromosomes are of the same length
    let chrom_length = elite.arr[0].n_loci();
    if pop_0.iter().any(|x| x.arr.iter().any(|chrom| chrom.n_loci() != chrom_length)) {
        return false;
    }
    // Check if all chromosomes are homozygous
    if pop_0.iter().any(|x| {
        let first_chrom = &x.arr[0];
        x.arr.iter().any(|chrom| {
            chrom != first_chrom
        })
    }) {
        return false;
    }
    // Check if genotypes are complimentary
    for i in 0..elite.n_chrom() {
        let elite_chrom = &elite.arr[i];
        let donor_chrom = &donor.arr[i];
        if !elite_chrom.is_complementary(donor_chrom) {
            return false;
        }
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
        assert!(res.is_some());
    }

    #[test]
    fn test_multi_chrom_min_generations_different_chromosomes() {
        let pop_0 = vec![
            MultiChromGenotype {
                arr: vec![
                    SingleChromGenotype::from_str("00000111111", "00000111111"),
                    SingleChromGenotype::from_str("00011110000", "00011110000"),
                ],
            },
            MultiChromGenotype {
                arr: vec![
                    SingleChromGenotype::from_str("11111000000", "11111000000"),
                    SingleChromGenotype::from_str("11100001111", "11100001111"),
                ],
            },
        ];
        let res = breeding_program(&pop_0);
        assert!(res.is_some());
        let wge = res.unwrap();
        assert_eq!(3, wge.generations());
    }
}
