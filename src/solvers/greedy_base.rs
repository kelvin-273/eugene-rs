use crate::abstract_plants::*;
use std::rc::Rc;

/// Constraints:
/// - has to be able to produce segments
/// - has to have two genotypes
pub fn breeding_program<A, B, K, S>(n_loci: usize, pop_0: Vec<A>, ideotype: A) -> Option<WGen<A, B>>
where
    A: Genotype<B> + Diploid<B> + SingleChrom,
    B: Gamete<A> + Haploid + SingleChrom,
    S: HaploidSegment<A, B> + Clone,
{
    // Generate minimal ordered set of segments
    let min_segments = min_covering_segments::<A, B, S>(n_loci, &pop_0);

    // Construct crossing tree from segments
    let c_star = join_segments(&min_segments, 0, min_segments.len());
    Some(WGen::new(A::from_gametes(
        &c_star.gamete().gamete,
        &c_star.gamete().gamete,
    )))
}

fn join_segments<A, B, S>(segments_s0: &Vec<S>, i: usize, j: usize) -> S
where
    A: Genotype<B> + Diploid<B>,
    B: Gamete<A> + Haploid,
    S: HaploidSegment<A, B> + Clone,
{
    // The base case can't return the raw segments as they are owned by the vector.
    // Although, we could clone the base segments.
    assert!(segments_s0.len() > 0);
    assert!(i < j && j <= segments_s0.len());
    if j == i + 1 {
        return segments_s0[i].clone();
    } else {
        let mid = (i + j) >> 1;
        let res1 = join_segments::<A, B, S>(segments_s0, i, mid);
        let res2 = join_segments::<A, B, S>(segments_s0, mid, j);
        S::join(&res1, &res2)
    }
}

pub fn min_covering_segments<A, B, S>(n_loci: usize, pop_0: &Vec<A>) -> Vec<S>
where
    A: Genotype<B> + Diploid<B> + SingleChrom,
    B: Gamete<A> + Haploid + SingleChrom,
    S: HaploidSegment<A, B> + Clone,
{
    let mut segment_pigeonholes: Vec<Option<S>> = vec![None; n_loci];
    for x in pop_0 {
        for c in generate_segments_genotype_diploid::<A, B, S>(n_loci, x) {
            let s = c.start();
            let e = c.end();
            if segment_pigeonholes[s].is_none() {
                segment_pigeonholes[s] = Some(c);
            } else {
                segment_pigeonholes[s].as_mut().map(|c_old| {
                    if c_old.end() < e {
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
                Some((s, e)) => {
                    assert!(*s < c.start());
                    if *e < c.end() {
                        *s = c.start();
                        *e = c.end();
                        Some(c)
                    } else {
                        None
                    }
                }
                None => {
                    *state = Some((c.start(), c.end()));
                    Some(c)
                }
            })
        })
        .filter_map(|c| c.map(|c| c.clone()))
        .collect()
}

fn generate_segments_genotype_diploid<A, B, S>(n_loci: usize, x: &A) -> Vec<S>
where
    A: Genotype<B> + Diploid<B> + SingleChrom,
    B: Gamete<A> + Haploid + SingleChrom,
    S: HaploidSegment<A, B> + Clone,
{
    let q1_true = generate_segments_gamete_haploid::<A, B, S>(&x.upper());
    let q2_true = generate_segments_gamete_haploid::<A, B, S>(&x.lower());
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
        let (s1, e1) = (c1.start(), c1.end());
        let (s2, e2) = (c2.start(), c2.end());

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
            out.push(S::join(&c1, &c2));
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
    q.0[i..].iter().for_each(|c| out.push(c.clone()));
    q.1[j..].iter().for_each(|c| out.push(c.clone()));
    out
}

fn generate_segments_gamete_haploid<A, B, S>(gx: &B) -> Vec<S>
where
    A: Genotype<B> + Diploid<B> + SingleChrom,
    B: Gamete<A> + Haploid + SingleChrom,
    S: HaploidSegment<A, B>,
{
    let alleles = gx.alleles();
    let n_loci = alleles.len();
    let wgx = Rc::new(gx.lift_b());
    let mut out = vec![];
    let mut i = 0;
    let mut j = 0;
    while i < n_loci && j < n_loci {
        if alleles[i] == Allele::Z {
            assert_eq!(i, j);
            i += 1;
        } else if alleles[j] == Allele::Z {
            out.push(S::from_start_end_gamete(i, j - 1, wgx.clone()));
            i = j + 1;
        }
        j += 1;
    }
    if i < j {
        out.push(S::from_start_end_gamete(i, j - 1, wgx.clone()));
    }
    out
}

/// Finds the number of generations required to construct the ideotype from pop_0.
/// Assumes that the ideotype is not in pop_0.
pub fn min_generations<A, B, S>(n_loci: usize, pop_0: &Vec<A>) -> usize
where
    A: Genotype<B> + Diploid<B> + SingleChrom,
    B: Gamete<A> + Haploid + SingleChrom,
    S: HaploidSegment<A, B> + Clone,
{
    let n_segments = min_covering_segments::<A, B, S>(n_loci, pop_0).len();
    (n_segments as f64).log2().ceil() as usize + 1
}

#[cfg(test)]
mod tests {
    use super::*;

    impl SingleChrom for u32 {}

    impl SingleChrom for u64 {}

    #[derive(Debug, PartialEq, Eq, Clone)]
    struct SegmentU32 {
        start: usize,
        end: usize,
        wgam: Rc<WGam<u64, u32>>,
    }

    impl Segment<u64, u32> for SegmentU32 {}

    impl HaploidSegment<u64, u32> for SegmentU32 {
        fn from_start_end_gamete(s: usize, e: usize, g: Rc<WGam<u64, u32>>) -> Self {
            Self {
                start: s,
                end: e,
                wgam: g,
            }
        }

        fn start(&self) -> usize {
            self.start
        }

        fn end(&self) -> usize {
            self.end
        }

        fn gamete(&self) -> Rc<WGam<u64, u32>> {
            self.wgam.clone()
        }

        fn join(&self, other: &Self) -> Self {
            self.clone()
            //SegmentU32 {
            //    start: self.start(),
            //    end: other.end(),
            //    wgam: ()
            //}
        }
    }

    #[test]
    fn generate_segments_gamete_haploid_test() {
        let f = |gx: u32| {
            generate_segments_gamete_haploid::<u64, u32, SegmentU32>(&gx)
                .iter()
                .map(|c| (c.start, c.end))
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

        breeding_program::<SingleChromGenotype, SingleChromGamete, CrosspointBitVec, SegmentBitVec>(
            n_loci,
            SingleChromGenotype::init_pop_random(&mut rng, n_loci, n_pop),
            SingleChromGenotype::ideotype(n_loci),
        );
    }
}
