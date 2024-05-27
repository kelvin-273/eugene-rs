use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solution::Objective;
use pyo3::prelude::*;
use std::collections::HashMap;

/// Runs a breeding program given `n_loci` and `pop_0` where `pop_0` is a population of single
/// chromosome diploid genotypes with `n_loci` loci.
#[pyo3::pyfunction]
pub fn breeding_program_python(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<bool>>>,
) -> PyResult<
    Option<(
        Vec<Vec<Vec<i32>>>,
        Vec<&'static str>,
        Vec<usize>,
        Vec<usize>,
        usize,
    )>,
> {
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
    let res = breeding_program(n_loci, pop_0);
    match res {
        None => Ok(None),
        Some(x_star) => {
            let sol = x_star.to_base_sol(n_loci, Objective::Generations);
            Ok(Some((
                sol.tree_data,
                sol.tree_type,
                sol.tree_left,
                sol.tree_right,
                sol.objective,
            )))
        }
    }
}

#[derive(Debug, Clone)]
pub struct Segment<B> {
    pub s: usize,
    pub e: usize,
    pub g: B,
}

impl<B> Segment<B> {
    pub fn start(&self) -> usize {
        self.s
    }

    pub fn end(&self) -> usize {
        self.e
    }

    pub fn gamete(&self) -> &B {
        &self.g
    }
}

type Seg = Segment<SingleChromGamete>;
type WGa = WGamS2<SingleChromGenotype, SingleChromGamete>;
type SegW = Segment<WGamS2<SingleChromGenotype, SingleChromGamete>>;

impl Seg {
    pub fn join(&self, other: &Self) -> Self {
        assert!(self.s < other.s);
        // TODO: complete adjacency condition <27-05-24> //
        Self {
            s: self.s,
            e: other.e,
            g: {
                CrosspointBitVec::new(false, other.s)
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
                let crosspoint_bit_vec = &CrosspointBitVec::new(false, other.s);
                let wz = WGenS2::from_gametes(&self.g, &other.g);
                wz.cross(|z| crosspoint_bit_vec.cross(&z))
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
                let z = SingleChromGenotype::from_gametes(&self.g.gamete(), &other.g.gamete());
                let wz = WGenS2::new(z);
                let crosspoint_bit_vec = &CrosspointBitVec::new(false, other.s);
                let gz = crosspoint_bit_vec.cross(wz.genotype());
                h.entry(gz.clone())
                    .or_insert_with(|| WGamS2::new_from_genotype(gz, wz))
                    .clone()
            },
        }
    }
}

/// Constraints:
/// - has to be able to produce segments
/// - has to have two genotypes
pub fn breeding_program(
    n_loci: usize,
    pop_0: Vec<SingleChromGenotype>,
) -> Option<WGenS2<SingleChromGenotype, SingleChromGamete>> {
    // Generate minimal ordered set of segments
    let min_segments = min_covering_segments(n_loci, &pop_0);

    // Construct crossing tree from segments
    None
}

pub fn min_covering_segments(n_loci: usize, pop_0: &Vec<SingleChromGenotype>) -> Vec<SegW> {
    let mut segment_pigeonholes: Vec<Option<SegW>> = vec![None; n_loci];
    for x in pop_0 {
        for c in generate_segments_genotype_diploid(n_loci, x) {
            let s = c.s;
            let e = c.e;
            if segment_pigeonholes[s].is_none() {
                segment_pigeonholes[s] = Some(c);
            } else {
                segment_pigeonholes[s].as_mut().map(|c_old| {
                    if c_old.e < e {
                        *c_old = c
                    } else {
                    }
                });
            }
        }
    }
    segment_pigeonholes
        .iter()
        .scan(None, |state, c| {
            c.as_ref().map(|c| match state {
                Some((ref mut s, ref mut e)) => {
                    assert!(*s < c.s);
                    if *e < c.e {
                        *s = c.s;
                        *e = c.e;
                        Some(c)
                    } else {
                        None
                    }
                }
                None => {
                    *state = Some((c.s, c.e));
                    Some(c)
                }
            })
        })
        .filter_map(|c| c.map(|c| c.clone()))
        .collect()
}

fn generate_segments_genotype_diploid(n_loci: usize, x: &SingleChromGenotype) -> Vec<SegW> {
    let wx: WGenS2<SingleChromGenotype, SingleChromGamete> = WGenS2::new(x.clone());
    let q1_true = generate_segments_gamete_haploid::<SingleChromGamete>(x.upper());
    let q2_true = generate_segments_gamete_haploid::<SingleChromGamete>(x.lower());
    let mut q = (&q1_true, &q2_true);

    let mut used1_true = vec![false; q1_true.len()];
    let mut used2_true = vec![false; q2_true.len()];
    let mut used1 = &mut used1_true;
    let mut used2 = &mut used2_true;

    let mut i = 0;
    let mut j = 0;

    let mut out = vec![];

    // oh dear god, does this even work?
    while i < q.0.len() && j < q.1.len() {
        let c1 = q.0[i].clone();
        let c2 = q.1[j].clone();
        let (s1, e1) = (c1.s, c1.e);
        let (s2, e2) = (c2.s, c2.e);

        if s1 > s2 || s1 == s2 && e1 < e2 {
            q = (q.1, q.0);
            (used1, used2) = (used2, used1);
            (i, j) = (j, i);
        } else if s1 <= s2 && e1 >= e2 {
            j += 1;
        } else if e1 + 1 < s2 {
            assert!(s1 < e1 + 1 && e1 + 1 < s2 && s2 <= e2);
            if !used1[i] {
                out.push(c1);
                used1[i] = true;
            }
            i += 1;
        } else {
            out.push(Seg::join(&c1, &c2));
            used1[i] = true;
            used2[j] = true;
            // break early if possible
            if e2 == n_loci - 1 {
                i = q.0.len();
                j = q.1.len();
            } else {
                i += 1;
            }
        }
    }

    out.into_iter()
        .chain(q.0[i..].iter().map(|c| c.clone()))
        .chain(q.1[j..].iter().map(|c| c.clone()))
        .map(|c| SegW {
            s: c.s,
            e: c.e,
            g: WGamS2::new_from_genotype(c.g.clone(), wx.clone()),
        })
        .collect()
    //q.0[i..].iter().for_each(|c| out.push(c.clone()));
    //q.1[j..].iter().for_each(|c| out.push(c.clone()));
    //out
}

fn generate_segments_gamete_haploid<B>(gx: B) -> Vec<Segment<B>>
where
    B: Haploid + Clone,
{
    let alleles = gx.alleles();
    let n_loci = alleles.len();
    let mut out = vec![];
    let mut i = 0;
    let mut j = 0;
    while i < n_loci && j < n_loci {
        if alleles[i] == Allele::Z {
            assert_eq!(i, j);
            i += 1;
        } else if alleles[j] == Allele::Z {
            out.push(Segment {
                s: i,
                e: j - 1,
                g: gx.clone(),
            });
            i = j + 1;
        }
        j += 1;
    }
    if i < j {
        out.push(Segment {
            s: i,
            e: j - 1,
            g: gx,
        });
    }
    out
}

/// Finds the number of generations required to construct the ideotype from pop_0.
/// Assumes that the ideotype is not in pop_0.
pub fn min_generations(n_loci: usize, pop_0: &Vec<SingleChromGenotype>) -> usize {
    let n_segments = min_covering_segments(n_loci, pop_0).len();
    (n_segments as f64).log2().ceil() as usize + 1
}

#[cfg(test)]
mod tests {
    use super::*;

    impl SingleChrom for u32 {}

    impl SingleChrom for u64 {}

    #[test]
    fn generate_segments_gamete_haploid_test() {
        let f = |gx: u32| {
            generate_segments_gamete_haploid::<u32>(gx)
                .iter()
                .map(|c| (c.start(), c.end()))
                .collect::<Vec<(usize, usize)>>()
        };

        assert_eq!(f(0x0000_0f0f), vec![(20, 23), (28, 31)]);
        assert_eq!(f(0x0000_0f01), vec![(20, 23), (31, 31)]);
        assert_eq!(f(0x1000_0f01), vec![(3, 3), (20, 23), (31, 31)]);
        assert_eq!(f(0x8000_0f01), vec![(0, 0), (20, 23), (31, 31)]);
    }

    #[test]
    fn bitvec_breeding_program_test() {
        use crate::plants::bit_array::*;
        let n_loci = 10;
        let n_pop = 6;
        let mut rng = rand::thread_rng();

        breeding_program(
            n_loci,
            SingleChromGenotype::init_pop_random(&mut rng, n_loci, n_pop),
        );
    }
}
