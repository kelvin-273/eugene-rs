use crate::abstract_plants::*;

impl Genotype<u32> for u64 {
    fn lift_a(&self) -> WGen<Self, u32> {
        WGen {
            genotype: *self,
            history: None,
        }
    }

    fn from_gametes(gx: &u32, gy: &u32) -> Self {
        ((gx << 32) as u64) + (*gy as u64)
    }
}

impl BioSize for u64 {
    fn get_sizes(&self) -> Vec<(usize, usize)> {
        vec![(2, 32)]
    }
}

impl BioSize for u32 {
    fn get_sizes(&self) -> Vec<(usize, usize)> {
        vec![(1, 32)]
    }
}

impl IndexAllele<usize> for u32 {
    fn index(&self, index: usize) -> Allele {
        match self >> (index as u32) & 1 {
            0 => Allele::Z,
            1 => Allele::O,
            _ => panic!("Bit-shift returning non-{{0,1}} value"),
        }
    }
}

impl Haploid for u32 {
    fn alleles(&self) -> Vec<Allele> {
        (0..32)
            .rev()
            .map(|i| match (self >> i) & 1 {
                0 => Allele::Z,
                1 => Allele::O,
                _ => panic!("should be a single digit"),
            })
            .collect()
    }
}

impl Diploid<u32> for u64 {
    fn upper(&self) -> u32 {
        (*self >> 32) as u32
    }

    fn lower(&self) -> u32 {
        (*self & ((1 << 32) - 1)) as u32
    }
}

impl Gamete<u64> for u32 {
    fn lift_b(&self) -> WGam<u64, Self> {
        WGam {
            gamete: *self,
            history: vec![],
        }
    }
}

impl Crosspoint<u64, u32, ()> for u32 {
    fn cross(self, x: &u64) -> u32 {
        ((*x as u32) & self) | (((*x >> 32) as u32) & (!self))
    }

    fn crosspoints(_data: &()) -> Box<dyn Iterator<Item = u32> + 'static> {
        Box::new(u32::MIN..u32::MAX)
    }
}

#[derive(Debug)]
struct CrosspointSingleU64U32 {
    start: bool,
    len_prefix: usize,
}

impl Crosspoint<u64, u32, ()> for CrosspointSingleU64U32 {
    fn cross(self, x: &u64) -> u32 {
        if self.len_prefix > 31 {
            panic!("consuming too many bits for n_loci = 32")
        };
        let mask_u32 = (1 << 32) - 1;
        let mask_shift = 31 - self.len_prefix;
        let mask_pref: u32 = ((mask_u32 >> (mask_shift)) << (mask_shift)) as u32;
        let head: u32 = (x >> 32) as u32;
        let tail: u32 = (x & mask_u32) as u32;
        match self.start {
            true => (head & mask_pref) | (tail & !mask_pref),
            false => (head & !mask_pref) | (tail & mask_pref),
        }
    }

    fn crosspoints(_data: &()) -> Box<dyn std::iter::Iterator<Item = Self>> {
        Box::new((0..=1).flat_map(|start| {
            (0..32).map(move |len_prefix| CrosspointSingleU64U32 {
                start: start != 0,
                len_prefix,
            })
        }))
    }
}
