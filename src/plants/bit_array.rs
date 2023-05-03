use crate::abstract_plants::*;
use bit_vec::BitVec;
use num;
use rand::prelude::*;

#[derive(Debug, Clone, PartialEq, PartialOrd, Ord, Eq, Hash)]
pub struct SingleChromGamete {
    n_loci: usize,
    gamete: BitVec,
}

impl IndexAllele<usize> for SingleChromGamete {
    fn index(&self, idx: usize) -> Allele {
        Allele::from(self.gamete[idx])
    }
}

impl Haploid for SingleChromGamete {
    fn alleles(&self) -> Vec<Allele> {
        self.gamete.iter().map(|x| Allele::from(x)).collect()
    }
}

impl BioSize for SingleChromGamete {
    fn get_sizes(&self) -> Vec<(usize, usize)> {
        vec![(1, self.n_loci)]
    }
}

impl Gamete<SingleChromGenotype> for SingleChromGamete {
    fn lift_b(&self) -> WGam<SingleChromGenotype, Self> {
        WGam::new(self.clone())
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd, Ord, Eq)]
pub struct SingleChromGenotype {
    n_loci: usize,
    chrom1: BitVec,
    chrom2: BitVec,
}

impl SingleChromGenotype {
    pub fn new(v: Vec<(bool, bool)>) -> Self {
        Self {
            n_loci: v.len(),
            chrom1: BitVec::from_fn(v.len(), |i| v[i].0),
            chrom2: BitVec::from_fn(v.len(), |i| v[i].1),
        }
    }

    pub fn from_str(s1: &str, s2: &str) -> Self {
        assert_eq!(s1.len(), s2.len());
        let f = |c: char| match c {
            '0' => false,
            '1' => true,
            _ => panic!("invalid character"),
        };
        Self {
            n_loci: s1.len(),
            chrom1: BitVec::from_iter(s1.chars().map(f)),
            chrom2: BitVec::from_iter(s2.chars().map(f)),
        }
    }

    pub fn ideotype(n_loci: usize) -> Self {
        Self {
            n_loci,
            chrom1: BitVec::from_elem(n_loci, true),
            chrom2: BitVec::from_elem(n_loci, true),
        }
    }

    fn random_genotype<R>(rng: &mut R, n_loci: usize) -> Self
    where
        R: Rng + ?Sized,
    {
        Self {
            n_loci,
            chrom1: BitVec::from_fn(n_loci, |_| rng.gen()),
            chrom2: BitVec::from_fn(n_loci, |_| rng.gen()),
        }
    }

    pub fn init_pop_random<R>(rng: &mut R, n_loci: usize, n_pop: usize) -> Vec<SingleChromGenotype>
    where
        R: rand::Rng + ?Sized,
    {
        // TODO: This is naive, do masking instead <26-04-23> //
        let mut pop_0 = (0..n_pop)
            .map(|_| SingleChromGenotype::random_genotype::<R>(rng, n_loci))
            .collect();
        while !SingleChromGenotype::is_feasible(&n_loci, &pop_0) {
            pop_0 = (0..n_pop)
                .map(|_| SingleChromGenotype::random_genotype::<R>(rng, n_loci))
                .collect();
        }
        pop_0
    }
}

impl BioSize for SingleChromGenotype {
    fn get_sizes(&self) -> Vec<(usize, usize)> {
        vec![(2, self.n_loci)]
    }
}

impl Diploid<SingleChromGamete> for SingleChromGenotype {
    fn upper(&self) -> SingleChromGamete {
        unimplemented!();
    }

    fn lower(&self) -> SingleChromGamete {
        unimplemented!();
    }
}

impl Genotype<SingleChromGamete> for SingleChromGenotype {
    fn lift_a(&self) -> WGen<Self, SingleChromGamete> {
        WGen::new(self.clone())
    }

    fn from_gametes(gx: &SingleChromGamete, gy: &SingleChromGamete) -> Self {
        assert_eq!(gx.n_loci, gy.n_loci);
        Self {
            n_loci: gx.n_loci,
            chrom1: gx.gamete.clone().into(),
            chrom2: gy.gamete.clone().into(),
        }
    }
}

impl Feasible<usize> for SingleChromGenotype {
    fn is_feasible(n_loci: &usize, pop: &Vec<Self>) -> bool {
        let mut total = BitVec::from_elem(*n_loci, false);
        for x in pop {
            total.or(&x.chrom1);
            total.or(&x.chrom2);
        }
        total == BitVec::from_elem(*n_loci, true)
    }
}

impl SingleChromGamete {
    pub fn from_str(s1: &str) -> Self {
        let f = |c: char| match c {
            '0' => false,
            '1' => true,
            _ => panic!("invalid character"),
        };
        Self {
            n_loci: s1.len(),
            gamete: BitVec::from_iter(s1.chars().map(f)),
        }
    }

    pub fn ideotype(n_loci: usize) -> Self {
        Self {
            n_loci,
            gamete: BitVec::from_elem(n_loci, true),
        }
    }
}

pub struct DomGamete {}

impl Dominance<SingleChromGamete> for DomGamete {
    fn dom(x: &SingleChromGamete, y: &SingleChromGamete) -> bool {
        x.gamete
            .iter()
            .zip(y.gamete.iter())
            .all(|(xi, yi)| xi >= yi)
    }
}

#[derive(Debug)]
pub struct CrosspointBitVec {
    start: bool,
    head: usize,
}

impl CrosspointBitVec {
    pub fn new(start: bool, head: usize) -> Self {
        Self { start, head }
    }
}

impl Crosspoint<SingleChromGenotype, SingleChromGamete, usize> for CrosspointBitVec {
    fn cross(self, x: &SingleChromGenotype) -> SingleChromGamete {
        let mut c1 = &x.chrom1;
        let mut c2 = &x.chrom2;
        if self.start {
            (c1, c2) = (c2, c1);
        }
        SingleChromGamete {
            n_loci: x.n_loci,
            gamete: {
                let mut gx = BitVec::from_elem(x.n_loci, false);
                for i in (0..self.head) {
                    let b = c1[i];
                    gx.set(i, b);
                }
                for i in (self.head..x.n_loci) {
                    let b = c2[i];
                    gx.set(i, b);
                }
                gx
            },
        }
    }

    fn crosspoints(n_loci: &usize) -> Box<dyn std::iter::Iterator<Item = Self>> {
        let n_loci = *n_loci;
        Box::new([false, true].iter().flat_map(move |start| {
            (0..n_loci).map(|head| CrosspointBitVec {
                start: *start,
                head,
            })
        }))
    }
}

fn test_bitvec() {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cross_test() {
        macro_rules! f {
            ($start:expr, $head:expr, $s1:expr, $s2:expr, $s3:expr) => {
                assert_eq!(
                    CrosspointBitVec::new($start, $head)
                        .cross(&SingleChromGenotype::from_str($s1, $s2)),
                    SingleChromGamete::from_str($s3)
                )
            };
        }
        f!(false, 0, "010010", "101001", "101001");
        f!(false, 1, "010010", "101001", "001001");
        f!(false, 2, "010010", "101001", "011001");
        f!(false, 3, "010010", "101001", "010001");
        f!(false, 4, "010010", "101001", "010001");
        f!(false, 5, "010010", "101001", "010011");
        f!(true, 0, "010010", "101001", "010010");
        f!(true, 1, "010010", "101001", "110010");
        f!(true, 2, "010010", "101001", "100010");
        f!(true, 3, "010010", "101001", "101010");
        f!(true, 4, "010010", "101001", "101010");
        f!(true, 5, "010010", "101001", "101000");
    }

    #[test]
    fn feasibility_test() {
        assert!(SingleChromGenotype::is_feasible(
            &5,
            &vec![SingleChromGenotype::from_str("01010", "10101")]
        ));
        assert!(SingleChromGenotype::is_feasible(
            &5,
            &vec![
                SingleChromGenotype::from_str("01010", "01010"),
                SingleChromGenotype::from_str("10101", "10101")
            ]
        ));
        assert!(SingleChromGenotype::is_feasible(
            &5,
            &vec![
                SingleChromGenotype::from_str("01010", "01010"),
                SingleChromGenotype::from_str("10101", "00101")
            ]
        ));
        assert!(!SingleChromGenotype::is_feasible(
            &5,
            &vec![
                SingleChromGenotype::from_str("01010", "01010"),
                SingleChromGenotype::from_str("00101", "00101")
            ]
        ));
    }

    #[test]
    fn feasibility_random_population_test() {
        let n_loci = 13;
        let mut rng = thread_rng();
        for _ in (0..100) {
            assert!(SingleChromGenotype::is_feasible(
                &n_loci,
                &SingleChromGenotype::init_pop_random(&mut rng, n_loci, 6)
            ));
        }
    }
}
