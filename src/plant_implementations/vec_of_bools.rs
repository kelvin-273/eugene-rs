use crate::abstract_plants::*;
use crate::visualisation;
use crate::visualisation::Draw;
use rand::distributions::{Distribution, Uniform};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SingleChromGamete {
    pub array: Vec<bool>,
}

impl SingleChromGamete {
    pub fn new(array: Vec<bool>) -> Self {
        Self { array }
    }

    pub fn ideotype(n_loci: usize) -> Self {
        Self {
            array: vec![true; n_loci],
        }
    }

    pub fn from_str(s: &str) -> Self {
        Self {
            array: s
                .chars()
                .map(|c| match c {
                    '0' => false,
                    '1' => true,
                    _ => panic!("non-01 character entered!"),
                })
                .collect(),
        }
    }
}

impl BioSize for SingleChromGamete {
    fn get_sizes(&self) -> Vec<(usize, usize)> {
        vec![(1, self.array.len())]
    }
}

impl Gamete<SingleChromGenotype> for SingleChromGamete {
    fn lift_b(&self) -> WGam<SingleChromGenotype, Self> {
        WGam {
            gamete: self.clone(),
            history: vec![],
        }
    }
}

impl IndexAllele<usize> for SingleChromGamete {
    fn index(&self, idx: usize) -> Allele {
        Allele::from(self.array[idx])
    }
}

impl Haploid for SingleChromGamete {
    fn alleles(&self) -> Vec<Allele> {
        self.array.iter().map(|x| Allele::from(*x)).collect()
    }
}

impl Draw for SingleChromGamete {
    fn view_box_size(&self) -> Option<(usize, usize)> {
        Some((
            self.get_n_loci(0) * visualisation::BLOCKSIZE + 1,
            self.get_n_chrom(0) * visualisation::BLOCKSIZE + 1,
        ))
    }

    fn draw(&self) -> svg::node::element::Group {
        visualisation::draw_gam_base(self, (0, 0))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SingleChromGenotype {
    pub array: Vec<(bool, bool)>,
}

impl SingleChromGenotype {
    pub fn new(array: Vec<(bool, bool)>) -> Self {
        Self { array }
    }

    pub fn ideotype(n_loci: usize) -> Self {
        Self {
            array: vec![(true, true); n_loci],
        }
    }

    pub fn from_str(s1: &str, s2: &str) -> Self {
        assert_eq!(s1.len(), s2.len());
        Self {
            array: SingleChromGamete::from_str(s1)
                .array
                .into_iter()
                .zip(SingleChromGamete::from_str(s2).array.into_iter())
                .collect(),
        }
    }
}

impl BioSize for SingleChromGenotype {
    fn get_sizes(&self) -> Vec<(usize, usize)> {
        vec![(2, self.array.len())]
    }
}
impl Genotype<SingleChromGamete> for SingleChromGenotype {
    fn lift_a(&self) -> WGen<Self, SingleChromGamete> {
        WGen {
            genotype: self.clone(),
            history: None,
        }
    }

    fn from_gametes(gx: &SingleChromGamete, gy: &SingleChromGamete) -> Self {
        assert_eq!(gx.get_sizes(), gy.get_sizes());
        Self {
            array: gx
                .array
                .iter()
                .zip(gy.array.iter())
                .map(|(a, b)| (*a, *b))
                .collect(),
        }
    }
}

impl Diploid<SingleChromGamete> for SingleChromGenotype {
    fn upper(&self) -> SingleChromGamete {
        SingleChromGamete::new((0..self.get_n_loci(0)).map(|i| self.array[i].0).collect())
    }

    fn lower(&self) -> SingleChromGamete {
        SingleChromGamete::new((0..self.get_n_loci(0)).map(|i| self.array[i].1).collect())
    }
}

impl Draw for SingleChromGenotype {
    fn view_box_size(&self) -> Option<(usize, usize)> {
        Some((
            self.get_n_loci(0) * visualisation::BLOCKSIZE + 1,
            self.get_n_chrom(0) * visualisation::BLOCKSIZE + 1,
        ))
    }

    fn draw(&self) -> svg::node::element::Group {
        visualisation::draw_gen_base(self, (0, 0))
    }
}

pub struct CrosspointSingleVob {
    start: bool,
    len_prefix: usize,
}

impl CrosspointSingleVob {
    pub fn new(start: bool, len_prefix: usize) -> Self {
        Self { start, len_prefix }
    }
}

impl Crosspoint<SingleChromGenotype, SingleChromGamete, usize> for CrosspointSingleVob {
    fn cross(self, x: &SingleChromGenotype) -> SingleChromGamete {
        SingleChromGamete::new(
            (0..x.get_n_loci(0))
                .map(|i| match (self.start, i < self.len_prefix) {
                    (true, true) => x.array[i].0,
                    (true, false) => x.array[i].1,
                    (false, true) => x.array[i].1,
                    (false, false) => x.array[i].0,
                })
                .collect(),
        )
    }

    fn crosspoints(data: &usize) -> Box<dyn std::iter::Iterator<Item = Self>> {
        let _data = data.clone();
        Box::new([false, true].iter().flat_map(move |start| {
            return (0.._data).map(move |len_prefix| CrosspointSingleVob {
                start: *start,
                len_prefix,
            });
        }))
    }
}

pub struct CrosspointMultiVob {
    start: bool,
    len_prefix: Vec<usize>,
}

impl CrosspointMultiVob {
    pub fn new(start: bool, len_prefix: Vec<usize>) -> Self {
        Self { start, len_prefix }
    }
}

impl Crosspoint<SingleChromGenotype, SingleChromGamete, (usize, usize)> for CrosspointMultiVob {
    fn cross(self, x: &SingleChromGenotype) -> SingleChromGamete {
        let n_loci = x.array.len();
        let j_max = self.len_prefix.len();
        let mut j = 0;
        let mut v = Vec::with_capacity(n_loci);
        let mut current_chrom = self.start;
        for i in (0..n_loci) {
            while j < self.len_prefix.len() && i >= self.len_prefix[j] {
                current_chrom = !current_chrom;
                j += 1;
            }
            v[i] = match current_chrom {
                true => x.array[i].0,
                false => x.array[i].1,
            }
        }
        SingleChromGamete::new(v)
    }

    fn crosspoints(data: &(usize, usize)) -> Box<dyn std::iter::Iterator<Item = Self>> {
        let k = data.0;
        let n_loci = data.1;
        unimplemented!()
    }
}

pub struct WeakDomSingle {}

impl Dominance<SingleChromGenotype> for WeakDomSingle {
    fn dom(x: &SingleChromGenotype, y: &SingleChromGenotype) -> bool {
        assert_eq!(x.array.len(), y.array.len());
        x.array
            .iter()
            .zip(y.array.iter())
            .all(|((aa, ab), (ba, bb))| (aa >= ba && ab >= bb) || (aa >= bb && ab >= ba))
    }
}

impl Dominance<SingleChromGamete> for WeakDomSingle {
    fn dom(x: &SingleChromGamete, y: &SingleChromGamete) -> bool {
        assert_eq!(x.array.len(), y.array.len());
        x.array.iter().zip(y.array.iter()).all(|(a, b)| a >= b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vob_cross_test() {
        /// Want to test that all of the crossings work
        let x = SingleChromGenotype::from_str("00101", "10011");
        assert_eq!(CrosspointSingleVob::new(true, 0).cross(&x), SingleChromGamete::from_str("10011"));
        assert_eq!(CrosspointSingleVob::new(true, 1).cross(&x), SingleChromGamete::from_str("00011"));
        assert_eq!(CrosspointSingleVob::new(true, 2).cross(&x), SingleChromGamete::from_str("00011"));
        assert_eq!(CrosspointSingleVob::new(true, 3).cross(&x), SingleChromGamete::from_str("00111"));
        assert_eq!(CrosspointSingleVob::new(true, 4).cross(&x), SingleChromGamete::from_str("00101"));
        assert_eq!(CrosspointSingleVob::new(false, 0).cross(&x), SingleChromGamete::from_str("00101"));
        assert_eq!(CrosspointSingleVob::new(false, 1).cross(&x), SingleChromGamete::from_str("10101"));
        assert_eq!(CrosspointSingleVob::new(false, 2).cross(&x), SingleChromGamete::from_str("10101"));
        assert_eq!(CrosspointSingleVob::new(false, 3).cross(&x), SingleChromGamete::from_str("10001"));
        assert_eq!(CrosspointSingleVob::new(false, 4).cross(&x), SingleChromGamete::from_str("10011"));
    }
}
