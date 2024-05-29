use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solution::{Objective, PyBaseSolution};
use std::collections::HashMap;

/// Runs a breeding program given `n_loci` and `pop_0` where `pop_0` is a population of single
/// chromosome diploid genotypes with `n_loci` loci.
#[pyo3::pyfunction]
pub fn breeding_program_python(n_loci: usize, pop_0: Vec<Vec<Vec<bool>>>) -> PyBaseSolution {
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
    let res = breeding_program(n_loci, &pop_0);
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

/// A triple containing a start-point `s`, end-point `e`, and gamete `g` used to represent a
/// contiguous sequence of favourable alleles on `g` from index `s` to `e` inclusive.
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
    pop_0: &Vec<SingleChromGenotype>,
) -> Option<WGenS2<SingleChromGenotype, SingleChromGamete>> {
    // TODO: Test for feasibility <27-05-24> //

    // Generate minimal ordered set of segments
    let mut min_segments = min_covering_segments(n_loci, pop_0);

    // Construct crossing tree from segments
    while min_segments.len() > 1 {
        let n_segments = min_segments.len();
        let mut new_segments = Vec::with_capacity(n_segments.div_ceil(2));
        for i in 0..n_segments / 2 {
            let cx = &min_segments[2 * i];
            let cy = &min_segments[2 * i + 1];
            let cz = cx.join(cy);
            new_segments.push(cz);
        }
        if n_segments % 2 == 1 {
            new_segments.push(min_segments[n_segments - 1].clone());
        }
        min_segments = new_segments;
    }
    Some(WGenS2::from_gametes(&min_segments[0].g, &min_segments[0].g))
}

/// Given a vector of segments, returns the minimum cardinality subset of segments that covers
/// {0..n_loci-1}.
pub fn min_covering_segments(n_loci: usize, pop_0: &Vec<SingleChromGenotype>) -> Vec<SegW> {
    let mut segment_pigeonholes: Vec<Option<SegW>> = vec![None; n_loci];

    // Fill each hole with the largest segment that starts at s
    for x in pop_0 {
        for Segment { s, e, g } in generate_segments_genotype_diploid(n_loci, x) {
            match segment_pigeonholes[s].as_mut() {
                None => {
                    segment_pigeonholes[s] = Some(Segment { s, e, g });
                }
                Some(c_old) => {
                    if c_old.e < e {
                        *c_old = Segment { s, e, g };
                    }
                }
            }
        }
    }

    // Filter non-dominated segments
    let mut i = 1;
    assert!(segment_pigeonholes[0].is_some());
    let mut j = 1;
    let mut e_covered = segment_pigeonholes[0].as_ref().unwrap().e;
    while i < n_loci && e_covered < n_loci - 1 {
        // INV: segment_pigeonholes[j-1].e == e_covered.
        // INV: segment_pigeonholes[..j] is the minimum cardinality subset of segments in
        //      segment_pigeonholes[..i] that covers {0..e_covered}.
        if let Some(mut c_next) = segment_pigeonholes[i].clone() {
            while i < n_loci && i <= e_covered + 1 {
                // INV: Let c_prev = segment_pigeonholes[j-1].
                //      c_next is the segment in segment_pigeonholes[..i] with the largest endpoint
                //      such that c_prev.s < c_next.s <= e_covered + 1
                // EXC: i == n_loci \/ i > e_covered + 1
                if let Some(c) = segment_pigeonholes[i].as_ref() {
                    if c.e > c_next.e {
                        c_next = c.clone();
                    }
                }
                i += 1;
            }
            e_covered = c_next.e;
            segment_pigeonholes[j] = Some(c_next);
            j += 1;
        } else {
            i += 1;
        }
    }

    segment_pigeonholes.into_iter().take(j).flatten().collect()
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
        .chain(q.0[i..].iter().cloned())
        .chain(q.1[j..].iter().cloned())
        .map(|c| SegW {
            s: c.s,
            e: c.e,
            g: WGamS2::new_from_genotype(c.g.clone(), wx.clone()),
        })
        .collect()
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
    fn generate_segments_genotype_diploid3_test() {
        let segments =
            generate_segments_genotype_diploid(3, &SingleChromGenotype::from_str("010", "010"));
        assert_eq!(
            vec![(1, 1)],
            segments
                .iter()
                .map(|c| (c.s, c.e))
                .collect::<Vec<(usize, usize)>>()
        );
        let segments =
            generate_segments_genotype_diploid(3, &SingleChromGenotype::from_str("101", "101"));
        assert_eq!(
            vec![(0, 0), (2, 2)],
            segments
                .iter()
                .map(|c| (c.s, c.e))
                .collect::<Vec<(usize, usize)>>()
        );
    }

    #[test]
    fn generate_segments_genotype_diploid4_test() {
        let segments =
            generate_segments_genotype_diploid(4, &SingleChromGenotype::from_str("0101", "1010"));
        assert_eq!(
            vec![(0, 1), (1, 2), (2, 3)],
            segments
                .iter()
                .map(|c| (c.s, c.e))
                .collect::<Vec<(usize, usize)>>()
        );
    }

    #[test]
    fn min_segments_test() {
        let segments = min_covering_segments(
            3,
            &vec![
                SingleChromGenotype::from_str("010", "010"),
                SingleChromGenotype::from_str("101", "101"),
            ],
        );
        assert_eq!(
            vec![(0, 0), (1, 1), (2, 2)],
            segments
                .iter()
                .map(|c| (c.s, c.e))
                .collect::<Vec<(usize, usize)>>()
        );

        let segments =
            min_covering_segments(4, &vec![SingleChromGenotype::from_str("0101", "1010")]);
        assert_eq!(
            vec![(0, 1), (2, 3)],
            segments
                .iter()
                .map(|c| (c.s, c.e))
                .collect::<Vec<(usize, usize)>>()
        );
    }

    #[test]
    fn bitvec_breeding_program_4_zigzag_test() {
        use crate::plants::bit_array::*;
        let n_loci = 4;
        let pop_0 = vec![SingleChromGenotype::from_str("0101", "1010")];
        assert_eq!(
            2,
            breeding_program(n_loci, &pop_0)
                .expect("infeasible")
                .to_base_sol(n_loci, Objective::Generations)
                .objective
        );
    }

    #[test]
    fn bitvec_breeding_program_4_dist_zigzag_test() {
        use crate::plants::bit_array::*;
        let n_loci = 4;
        let pop_0 = vec![
            SingleChromGenotype::from_str("0101", "0101"),
            SingleChromGenotype::from_str("1010", "1010"),
        ];
        let op_sol = breeding_program(n_loci, &pop_0)
            .map(|sol| sol.to_base_sol(n_loci, Objective::Generations));
        dbg!(&op_sol);
        assert_eq!(3, op_sol.expect("infeasible").objective);
    }

    #[test]
    fn bitvec_breeding_program_dist_test() {
        use crate::extra::instance_generators;
        fn tester(actual: usize, v: &[usize]) {
            let pop_0 = instance_generators::distarray_to_homo(v);
            let n_loci = v.len();
            let op_sol = breeding_program(n_loci, &pop_0)
                .map(|sol| sol.to_base_sol(n_loci, Objective::Generations));
            dbg!(&op_sol);
            assert_eq!(actual, op_sol.expect("infeasible").objective);
        }
        tester(2, &[0, 1]);
        tester(3, &[0, 1, 0]);
        tester(3, &[0, 1, 2, 0]);
        tester(4, &[0, 1, 0, 2, 1]);
    }

    #[test]
    fn bitvec_breeding_program_random_test() {
        use crate::plants::bit_array::*;
        use crate::solvers::base_min_generations_enumerator_dominance2;
        let n_loci = 17;
        let n_pop = 6;
        let mut rng = rand::thread_rng();
        for _ in 0..20 {
            let pop_0 = SingleChromGenotype::init_pop_random(&mut rng, n_loci, n_pop);

            let sol_seg = breeding_program(n_loci, &pop_0);
            let sol_dom =
                base_min_generations_enumerator_dominance2::breeding_program(n_loci, &pop_0);
            assert_eq!(
                sol_dom.map(|sol| sol.objective),
                sol_seg.map(|sol| sol.to_base_sol(n_loci, Objective::Generations).objective)
            );
        }
    }
}
