use crate::abstract_plants::*;
use crate::extra::visualisation;
use crate::solution::{BaseSolution, Objective};
use bit_vec::BitVec;
use rand::prelude::*;
use std::collections::{HashMap, VecDeque};

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
        self.gamete.iter().map(Allele::from).collect()
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

#[derive(Debug, Clone, PartialEq, PartialOrd, Ord, Eq, Hash)]
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

    pub fn blank(n_loci: usize) -> Self {
        Self {
            n_loci,
            chrom1: BitVec::from_elem(n_loci, false),
            chrom2: BitVec::from_elem(n_loci, false),
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

    pub fn get(&self, chrom: bool, locus: usize) -> Option<bool> {
        match chrom {
            true => self.chrom1.get(locus),
            false => self.chrom2.get(locus),
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
        SingleChromGamete {
            n_loci: self.n_loci,
            gamete: self.chrom1.clone(),
        }
    }

    fn lower(&self) -> SingleChromGamete {
        SingleChromGamete {
            n_loci: self.n_loci,
            gamete: self.chrom2.clone(),
        }
    }
}

impl IndexAllele<(Chrom, usize)> for SingleChromGenotype {
    fn index(&self, idx: (Chrom, usize)) -> Allele {
        match idx.0 {
            Chrom::Upper => self.get(true, idx.1).expect("Index out of bounds").into(),
            Chrom::Lower => self.get(false, idx.1).expect("Index out of bounds").into(),
        }
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
            chrom1: gx.gamete.clone(),
            chrom2: gy.gamete.clone(),
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
    pub fn bool_array(v: &Vec<bool>) -> Self {
        Self {
            n_loci: v.len(),
            gamete: BitVec::from_iter(v.iter().copied()),
        }
    }

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

#[derive(Debug, Clone)]
pub struct CrosspointBitVec {
    start: Chrom,
    head: usize,
}

impl CrosspointBitVec {
    pub fn new(start: Chrom, head: usize) -> Self {
        Self { start, head }
    }

    pub fn random_crosspoint_uniform(rng: &mut impl Rng, n_loci: &usize) -> Self {
        Self { start: rng.gen::<bool>().into(), head: rng.gen_range(0..*n_loci) }
    }

    /// Generates a crosspoint with a uniform randomly starting chromosome and a head that is
    /// chosen with respect to a vector of recombination rates. The recombination rates are
    /// marginalised so that the crosspoints are single-point recombinations.
    pub fn random_crosspoint(rng: &mut impl Rng, n_loci: &usize, recomb_rates: &[f64]) -> Self {
        let start = rng.gen::<bool>().into();
        if *n_loci <= 1 {
            return Self { start, head: 0 }
        }
        let total: f64 = recomb_rates.iter().sum();
        let pinpoint: f64 = rng.gen();
        let mut j = 0;
        let mut running_sum = 0.0;
        while j < n_loci - 1 &&  running_sum / total < pinpoint {
            running_sum += recomb_rates[j];
            j += 1;
        }
        Self { start, head: j }
    }
}

impl Crosspoint<SingleChromGenotype, SingleChromGamete, usize> for CrosspointBitVec {
    fn cross(&self, x: &SingleChromGenotype) -> SingleChromGamete {
        let mut c1 = &x.chrom1;
        let mut c2 = &x.chrom2;
        if self.start == Chrom::Lower {
            (c1, c2) = (c2, c1);
        }
        SingleChromGamete {
            n_loci: x.n_loci,
            gamete: {
                let mut gx = BitVec::from_elem(x.n_loci, false);
                for i in 0..self.head {
                    let b = c1[i];
                    gx.set(i, b);
                }
                for i in self.head..x.n_loci {
                    let b = c2[i];
                    gx.set(i, b);
                }
                gx
            },
        }
    }

    fn crosspoints(n_loci: &usize) -> Box<dyn std::iter::Iterator<Item = Self>> {
        let n_loci = *n_loci;
        Box::new(
            [Chrom::Upper, Chrom::Lower]
                .into_iter()
                .flat_map(move |start| {
                    (0..n_loci).map(move |head| CrosspointBitVec { start, head })
                }),
        )
    }
}

impl SingleChrom for SingleChromGenotype {}

impl SingleChrom for SingleChromGamete {}

impl WGen<SingleChromGenotype, SingleChromGamete> {
    pub fn to_base_sol(&self, n_loci: usize, obj: Objective) -> BaseSolution {
        // get the number of crossings and cells required
        let mut q_node = VecDeque::new();
        let mut n_cells = 1;
        let mut n_cross = 0;

        /*
         * map every genotype to a unique index in the order they were seen starting from 0
         */

        // create index_map
        let mut index_map: HashMap<SingleChromGenotype, usize> = HashMap::new();
        index_map.insert(self.genotype().clone(), 0);

        q_node.push_back(self.clone());
        if self.history().is_some() {
            n_cross += 1;
        }

        // using q_node to process all genotypes for the counting
        while let Some(z) = q_node.pop_front() {
            if let Some((wgx, wgy)) = z.history() {
                // QUESTION: how is it that a single & here
                // makes all the difference?
                let x = wgx.history().expect("gamete wgx should have a history");
                if !index_map.contains_key(x.genotype()) {
                    index_map.insert(x.genotype().clone(), n_cells);
                    n_cells += 1;
                    q_node.push_back(x.clone());
                    if x.history().is_some() {
                        n_cross += 1;
                    }
                }

                let y = wgy.history().expect("gamete wgy should have a history");
                if !index_map.contains_key(y.genotype()) {
                    index_map.insert(y.genotype().clone(), n_cells);
                    n_cells += 1;
                    q_node.push_back(y.clone());
                    if y.history().is_some() {
                        n_cross += 1;
                    }
                }
            }
        }
        index_map.drain();

        // initialise arrays
        let mut tree_data = vec![vec![vec![0; n_loci]; 2]; n_cells];
        let mut tree_type = vec!["Null"; n_cells];
        let mut tree_left = vec![0; n_cells];
        let mut tree_right = vec![0; n_cells];

        // flush index_map because nodes and leaves have mixed order
        index_map.drain();
        index_map.insert(self.genotype().clone(), 0);

        q_node.push_back(self.clone());
        let mut i_node = 0;
        let mut i_leaf = n_cross;

        // using q_node to process all genotypes for the counting
        while let Some(z) = q_node.pop_front() {
            let &i_z = index_map
                .get(z.genotype())
                .expect("z.genotype is not in the index_map");

            // fill in data
            for (j, a) in z.genotype().upper().alleles().iter().enumerate() {
                tree_data[i_z][0][j] = match a {
                    Allele::Z => 0,
                    Allele::O => 1,
                };
            }
            for (j, a) in z.genotype().lower().alleles().iter().enumerate() {
                tree_data[i_z][1][j] = match a {
                    Allele::Z => 0,
                    Allele::O => 1,
                };
            }

            if let Some((wgx, wgy)) = z.history() {
                tree_type[i_z] = "Node";

                let x = wgx
                    .history()
                    .expect("wgx should have a history (in the queue)");
                if !index_map.contains_key(x.genotype()) {
                    if x.history().is_some() {
                        i_node += 1;
                        index_map.insert(x.genotype().clone(), i_node);
                    } else {
                        index_map.insert(x.genotype().clone(), i_leaf);
                        i_leaf += 1;
                    }
                    q_node.push_back(x.clone());
                }
                let &i_x = index_map
                    .get(x.genotype())
                    .expect("x.genotype should be in the index_map");
                tree_left[i_z] = i_x + 1;

                let y = wgy
                    .history()
                    .expect("wgy should have a history (in the queue)");
                if !index_map.contains_key(y.genotype()) {
                    if y.history().is_some() {
                        i_node += 1;
                        index_map.insert(y.genotype().clone(), i_node);
                    } else {
                        index_map.insert(y.genotype().clone(), i_leaf);
                        i_leaf += 1;
                    }
                    q_node.push_back(y.clone());
                }
                let &i_y = index_map
                    .get(y.genotype())
                    .expect("y.genotype should be in the index_map");
                tree_right[i_z] = i_y + 1;
            } else {
                tree_type[i_z] = "Leaf";
            }
        }

        // compute generations if required
        let objective = match obj {
            Objective::Crossings => n_cross,
            Objective::Generations => {
                let mut generations_dp = vec![-1; n_cells];
                fn aux(
                    tree_type: &Vec<&str>,
                    tree_left: &Vec<usize>,
                    tree_right: &Vec<usize>,
                    generations_dp: &mut Vec<i32>,
                    i: usize,
                ) -> i32 {
                    if generations_dp[i] < 0 {
                        generations_dp[i] = match tree_type[i] {
                            "Leaf" => 0,
                            "Node" => {
                                let gen_l = aux(
                                    tree_type,
                                    tree_left,
                                    tree_right,
                                    generations_dp,
                                    tree_left[i] - 1,
                                );
                                let gen_r = aux(
                                    tree_type,
                                    tree_left,
                                    tree_right,
                                    generations_dp,
                                    tree_right[i] - 1,
                                );
                                gen_l.max(gen_r) + 1
                            }
                            _ => panic!("tree type is neither Leaf nor Node"),
                        }
                    }
                    generations_dp[i]
                }
                aux(&tree_type, &tree_left, &tree_right, &mut generations_dp, 0) as usize
            }
        };

        BaseSolution {
            tree_data,
            tree_type,
            tree_left,
            tree_right,
            objective,
        }
    }
}

impl visualisation::Draw for SingleChromGenotype {
    fn view_box_size(&self) -> Option<(usize, usize)> {
        Some((
            self.get_n_loci(0) * visualisation::BLOCKSIZE + 1,
            self.get_n_chrom() * visualisation::BLOCKSIZE + 1,
        ))
    }

    fn draw(&self) -> svg::node::element::Group {
        visualisation::draw_gen_base(self, (0, 0))
    }
}

impl visualisation::Draw for SingleChromGamete {
    fn view_box_size(&self) -> Option<(usize, usize)> {
        Some((
            self.get_n_loci(0) * visualisation::BLOCKSIZE + 1,
            self.get_n_chrom() * visualisation::BLOCKSIZE + 1,
        ))
    }

    fn draw(&self) -> svg::node::element::Group {
        visualisation::draw_gam_base(self, (0, 0))
    }
}

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
        f!(Chrom::Upper, 0, "010010", "101001", "101001");
        f!(Chrom::Upper, 1, "010010", "101001", "001001");
        f!(Chrom::Upper, 2, "010010", "101001", "011001");
        f!(Chrom::Upper, 3, "010010", "101001", "010001");
        f!(Chrom::Upper, 4, "010010", "101001", "010001");
        f!(Chrom::Upper, 5, "010010", "101001", "010011");
        f!(Chrom::Lower, 0, "010010", "101001", "010010");
        f!(Chrom::Lower, 1, "010010", "101001", "110010");
        f!(Chrom::Lower, 2, "010010", "101001", "100010");
        f!(Chrom::Lower, 3, "010010", "101001", "101010");
        f!(Chrom::Lower, 4, "010010", "101001", "101010");
        f!(Chrom::Lower, 5, "010010", "101001", "101000");
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
        for _ in 0..100 {
            assert!(SingleChromGenotype::is_feasible(
                &n_loci,
                &SingleChromGenotype::init_pop_random(&mut rng, n_loci, 6)
            ));
        }
    }
}
