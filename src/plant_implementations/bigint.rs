pub struct SingleChromGamete {
    n_loci: usize,
    array: Vec<u8>
}

impl SingleChromGamete {
    pub fn new(n_loci: usize, array: Vec<u8>) -> Self {
        Self { n_loci, array }
    }

    pub fn ideotype(n_loci: usize) -> Self {
        let n_array = (n_loci + 7).div_euclid(8);
        let mut array = vec![u8::MAX; n_array];
        if n_loci % 8 != 0 {
            // TODO: set the last bits that are not in [n_loci] to 0 <06-02-23> //

        }
        Self { n_loci, array }
    }
}

pub struct SingleChromGenotype {
    n_loci: usize,
    
}
