use crate::abstract_plants::*;
use std::rc::Rc;

/// Constraints:
/// - has to be able to produce segments
/// - has to have two genotypes
fn breeding_program<A, B, K, S>(pop_0: Vec<A>, ideotype: A) -> Option<WGen<A, B>>
where
    A: Genotype<B> + Diploid<B> + SingleChrom,
    B: Gamete<A> + Haploid + SingleChrom,
    S: HaploidSegment<A, B> + Copy,
{
    // Generate ordered set of segments
    let segments = generate_segments::<A, B, K, S>(&pop_0);
    // Choose min set cover
    let min_segments = min_1d_geometric_set_cover::<A, B, S>(segments);
    // Construct crossing tree from segments
    None
}

fn generate_segments<A, B, K, S>(pop_0: &Vec<A>) -> Vec<S>
where
    A: Genotype<B> + Diploid<B> + SingleChrom,
    B: Gamete<A> + Haploid + SingleChrom,
    S: HaploidSegment<A, B> + Copy,
{
    pop_0
        .iter()
        .flat_map(|x| generate_segments_genotype_diploid::<A, B, S>(x).into_iter())
        .collect()
}

fn generate_segments_genotype_diploid<A, B, S>(x: &A) -> Vec<S>
where
    A: Genotype<B> + Diploid<B> + SingleChrom,
    B: Gamete<A> + Haploid + SingleChrom,
    S: HaploidSegment<A, B> + Copy,
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
        let c1 = q.0[i];
        let c2 = q.1[i];
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
        } else if e2 == 2 - 1 {
            panic!("That 2 should be an instance of n_loci");
            out.push(S::join(&c1, &c2));
            used1[i] = true;
            used2[j] = true;
            i += 1;
        }
    }
    out.extend_from_slice(&q.0[i..]);
    out.extend_from_slice(&q.1[j..]);
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

fn min_1d_geometric_set_cover<A, B, S>(mut segments: Vec<S>) -> Vec<S>
where
    A: Genotype<B> + Diploid<B>,
    B: Gamete<A> + Haploid,
    S: HaploidSegment<A, B>,
{
    segments.sort_by(|x, y| match x.start().partial_cmp(&y.start()).unwrap() {
        std::cmp::Ordering::Equal => y.end().partial_cmp(&x.end()).unwrap(),
        res => res,
    });
    vec![]
}

trait SingleChrom {}

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
}
