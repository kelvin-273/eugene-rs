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
        self.gamete
            .iter()
            .map(|x| Allele::from(x))
            .collect()
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

impl SingleChromGenotype {
    pub fn from_str(s1: &str, s2: &str) -> Self {
        assert_eq!(s1.len(), s2.len());
        let f = |c: char| { match c {
            '0' => false,
            '1' => true,
            _ => panic!("invalid character")
        } };
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

    fn random_genotype<R>(rng: &mut R, n_loci: usize) -> Self where
        R: Rng + ?Sized
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
        let mut pop_0 = (0..n_pop)
            .map(|_| SingleChromGenotype::random_genotype::<R>(rng, n_loci))
            .collect();
        pop_0
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

pub struct CrosspointBitVec {
    start: bool,
    head: usize,
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
                for (i, b) in (0..self.head).zip(c1.iter()) {
                    gx.set(i, b);
                }
                for (i, b) in (self.head..x.n_loci).zip(c1.iter()) {
                    gx.set(i, b);
                }
                gx
            },
        }
    }

    fn crosspoints(n_loci: &usize) -> Box<dyn std::iter::Iterator<Item = Self>> {
        Box::new([false, true]
            .iter()
            .zip(0..*n_loci)
            .map(|(&start, head)| CrosspointBitVec { start, head }))
    }
}

fn test_bitvec() {}
