use bitvec::prelude::*;
use rand::prelude::*;
use crate::abstract_plants::*;

type MyBits1 = BitArray<[u8; 5], Lsb0>;
type MyBits2 = BitVec;

#[derive(Debug, Clone, PartialEq, PartialOrd, Ord, Eq)]
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
        self.gamete.iter().by_vals().map(|x| Allele::from(x)).collect()
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

impl SingleChromGenotype {
    pub fn init_pop_random(n_loci: usize, n_pop: usize) -> Vec<SingleChromGenotype> {
        unimplemented!()
    }
}
