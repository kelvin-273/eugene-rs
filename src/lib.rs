mod base {
    use bitvec::prelude::*;
    use std::rc::Rc;

    pub struct Gamete {
        _g: BitArray,
    }

    pub struct Genotype {
        _c1: Rc<Gamete>,
        _c2: Rc<Gamete>,
    }

    pub struct TreeGamete {
        gamete: Gamete,
        history: Vec<TreeGenotype>,
    }

    impl TreeGamete {
        pub fn new(gamete: Gamete) -> Self {
            Self { 
                gamete,
                history: Vec::new()
            }
        }
    }

    pub struct TreeGenotype {
        genotype: Genotype,
        history: Option<(TreeGamete, TreeGamete)>,
    }

    impl TreeGenotype {
        pub fn new(genotype: Genotype) -> Self {
            Self { 
                genotype,
                history: None
            }
        }
    }

    pub struct Results {}

    fn breeding_program(
        n_loci: usize,
        ideotype: Genotype,
        pop0: &Vec<Genotype>,
    ) -> Option<Results> {
        // test feasibility
        if !test_feasibility(n_loci, pop0) {
            return None;
        }
        // enumerate
        None
    }

    fn test_feasibility(n_loci: usize, pop0: &Vec<Genotype>) -> bool {
        true
    }
}

mod full {
    use bitvec;
    use bitvec::prelude::*;

    /// Gametes and Genotypes are:
    /// - sets of
    /// - homologous groups of
    /// - chromosomes of
    /// - alleles
    ///
    /// Genotypes have 2x the number of homologous groups as the Gametes
    pub struct Gamete {
        g: Vec<Vec<BitVec>>,
    }

    #[derive(Clone)]
    pub struct Genotype {
        g: Vec<Vec<BitVec>>,
    }

    pub struct TreeGamete {
        gamete: Gamete,
        history: Vec<TreeGenotype>,
    }

    impl TreeGamete {
        pub fn new(gamete: Gamete) -> Self {
            Self { 
                gamete,
                history: Vec::new()
            }
        }
    }

    pub struct TreeGenotype {
        genotype: Genotype,
        history: Option<(Box<TreeGamete>, Box<TreeGamete>)>,
    }

    struct NodeGenotype {}

    struct TreeGenotype2 {
        head: NodeGenotype
    }

    pub struct Config {
        n_chrm: usize,
        n_homs: usize,
        n_loci: Vec<usize>
    }

    impl TreeGenotype {
        pub fn new(genotype: Genotype) -> Self {
            Self { 
                genotype,
                history: None
            }
        }

        pub fn lift2(f: &dyn Fn(Genotype, Genotype) -> Genotype, m: &Self, n: &Self) -> Self {
            let h1 = m.history;
            let h2 = n.history;
            TreeGenotype {
                genotype: f(m.genotype, n.genotype),
                history: Some((Box::new(*m), Box::new(n)))
            }
        }
    }

    pub struct Results {}

    fn breeding_program(
        n_homs: usize,
        n_chrm: usize,
        n_loci: &Vec<usize>,
        ideotype: &Genotype,
        pop0: &Vec<Genotype>,
    ) -> Option<Results> {
        // test feasibility
        if !test_feasibility(n_homs, n_chrm, n_loci, ideotype, pop0) {
            return None;
        }
        // clone initial population and lift to TreeGenotype
        let mut pop: Vec<TreeGenotype> = pop0.iter().map(|x| TreeGenotype::new(x.clone())).collect();
        // enumerate
        None
    }

    fn test_feasibility(
        n_homs: usize,
        n_chrm: usize,
        n_loci: &Vec<usize>,
        ideotype: &Genotype,
        pop0: &Vec<Genotype>,
    ) -> bool {
        // Verifying all lengths are the same
        if !(pop0.iter().all(|x| x.g.len() == n_homs)) {
            return false;
        }
        if !(pop0.iter().all(|x| x.g.iter().all(|h| h.len() == n_chrm))) {
            return false;
        }
        if !(pop0.iter().all(|x| {
            x.g.iter()
                .zip(n_loci)
                .all(|(h, &nl)| h.iter().all(|c| c.len() == nl))
        })) {
            return false;
        }
        // taking the union of the chromosomes
        let mut union: Vec<BitVec> = n_loci.iter().map(|&i| BitVec::with_capacity(i)).collect();
        for x in pop0 {
            for i in 0..n_homs {
                for j in 0..n_chrm {
                    union[i] = bitvec_union(&union[i], &x.g[i][j])
                }
            }
        }
        if union.iter().any(|x| x.first_zero().is_some()) { return false }
        // finally
        true
    }

    fn bitvec_union(x: &BitVec, y: &BitVec) -> BitVec {
        assert_eq!(x.len(), y.len());
        BitVec::from_iter(x.iter().zip(y.iter()).map(|(a, b)| *a | *b))
    }
} /* Full */
