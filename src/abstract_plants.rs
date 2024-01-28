use std::rc::Rc;

pub trait BioSize {
    /// Creates an array of sizes where each size is a tuple (n_chrom_i, n_loci_i) is a tuple of
    /// the number of chromosomes in the homologous group and the number of QTL on those
    /// chromosomes
    fn get_sizes(&self) -> Vec<(usize, usize)>;

    fn get_n_loci(&self, chromosome_i: usize) -> usize {
        self.get_sizes()[chromosome_i].1
    }

    fn get_n_chrom(&self) -> usize {
        self.get_sizes().len()
    }

    fn get_n_ploid(&self, chromosome_i: usize) -> usize {
        self.get_sizes()[chromosome_i].0
    }
}

pub trait IndexAllele<Idx> {
    fn index(&self, idx: Idx) -> Allele;
}

pub trait Genotype<B: Gamete<Self>>: Sized + BioSize + std::fmt::Debug {
    fn lift_a(&self) -> WGen<Self, B>;

    fn from_gametes(gx: &B, gy: &B) -> Self;
}

pub trait Gamete<A: Genotype<Self>>: Sized + BioSize + std::fmt::Debug {
    fn lift_b(&self) -> WGam<A, Self>;
}

pub trait Crosspoint<A: Genotype<B>, B: Gamete<A>, Data>: std::fmt::Debug {
    fn cross(self, x: &A) -> B;

    fn crosspoints(data: &Data) -> Box<dyn std::iter::Iterator<Item = Self>>;
}

#[derive(Debug, PartialEq, Eq)]
pub struct WGen<A, B> {
    pub genotype: A,
    pub history: Option<(Rc<WGam<A, B>>, Rc<WGam<A, B>>)>,
}

impl<A: Genotype<B>, B: Gamete<A>> WGen<A, B> {
    fn from_gametes(gx: Rc<WGam<A, B>>, gy: Rc<WGam<A, B>>) -> WGen<A, B> {
        WGen {
            genotype: Genotype::from_gametes(&gx.gamete, &gy.gamete),
            history: Some((gx, gy)),
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct WGam<A, B> {
    pub gamete: B,
    pub history: Vec<Rc<WGen<A, B>>>,
}

impl<A, B> WGam<A, B> {
    pub fn new(gamete: B) -> Self {
        Self {
            gamete,
            history: vec![],
        }
    }

    pub fn extract_first_from_ref(self: &Self) -> WGamS<A, B>
    where
        A: Clone,
        B: Clone,
    {
        WGamS {
            gamete: self.gamete.clone(),
            history: self
                .history
                .first()
                .map(|x| x.clone().extract_first())
        }
    }

    /// Given an Rc pointer, returns a new Rc pointer to the first extractable wrapped gamete.
    ///
    /// # Panics
    ///
    /// Panics if the list of source genotypes is empty.
    pub fn extract_first(self: Rc<Self>) -> Rc<WGamS<A, B>>
    where
        A: Clone,
        B: Clone,
    {
        Rc::new(WGamS {
            gamete: self.gamete.clone(),
            history: self
                .history
                .first()
                .map(|x| x.clone().extract_first())
        })
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct WGenS<A, B> {
    pub genotype: A,
    pub history: Option<(Rc<WGamS<A, B>>, Rc<WGamS<A, B>>)>,
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct WGamS<A, B> {
    pub gamete: B,
    pub history: Option<Rc<WGenS<A, B>>>,
}

impl<A, B> WGen<A, B> {
    pub fn new(genotype: A) -> Self {
        Self {
            genotype,
            history: None,
        }
    }

    pub fn extract_first(self: Rc<Self>) -> Rc<WGenS<A, B>>
    where
        A: Clone,
        B: Clone,
    {
        Rc::new(WGenS {
            genotype: self.genotype.clone(),
            history: self
                .history
                .as_ref()
                .map(|(lg, rg)| (lg.clone().extract_first(), rg.clone().extract_first())),
        })
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Allele {
    Z,
    O,
}

impl From<bool> for Allele {
    fn from(value: bool) -> Self {
        match value {
            true => Allele::O,
            false => Allele::Z,
        }
    }
}

impl Into<bool> for Allele {
    fn into(self) -> bool {
        match self {
            Allele::Z => false,
            Allele::O => true,
        }
    }
}

pub trait Segment<A: Genotype<B>, B: Gamete<A>> {}

pub trait HaploidSegment<A, B>: Segment<A, B>
where
    A: Genotype<B> + Diploid<B>,
    B: Gamete<A> + Haploid,
{
    fn from_start_end_gamete(s: usize, e: usize, g: Rc<WGam<A, B>>) -> Self
    where
        B: Haploid;

    fn start(&self) -> usize;

    fn end(&self) -> usize;

    fn gamete(&self) -> Rc<WGam<A, B>>;

    fn join(&self, other: &Self) -> Self;
}

pub trait Haploid: IndexAllele<usize> {
    fn alleles(&self) -> Vec<Allele>;
}

pub trait Diploid<B: Haploid> {
    fn upper(&self) -> B;

    fn lower(&self) -> B;
}

pub trait SingleChrom {}

pub trait Dominance<T> {
    fn dom(x: &T, y: &T) -> bool;
}

pub trait Traverse<A, B, S> {
    fn f_gen(x: &Rc<WGen<A, B>>, state: &mut S);

    fn f_gam(gx: &Rc<WGam<A, B>>, state: &mut S);

    fn traverse_gen_preorder(x: &Rc<WGen<A, B>>, state: &mut S) {
        Self::f_gen(x, state);
        match &x.history {
            None => {}
            Some((gx, gy)) => {
                Self::traverse_gam_preorder(gx, state);
                Self::traverse_gam_preorder(gy, state);
            }
        }
    }

    fn traverse_gam_preorder(gx: &Rc<WGam<A, B>>, state: &mut S) {
        Self::f_gam(gx, state);
        gx.history
            .iter()
            .for_each(|x| Self::traverse_gen_preorder(x, state));
    }

    fn traverse_gen_postorder(x: &Rc<WGen<A, B>>, state: &mut S) {
        match &x.history {
            None => {}
            Some((gx, gy)) => {
                Self::traverse_gam_postorder(gx, state);
                Self::traverse_gam_postorder(gy, state);
            }
        }
        Self::f_gen(x, state);
    }

    fn traverse_gam_postorder(gx: &Rc<WGam<A, B>>, state: &mut S) {
        gx.history
            .iter()
            .for_each(|x| Self::traverse_gen_postorder(x, state));
        Self::f_gam(gx, state);
    }
}

pub trait Feasible<Data>: Sized {
    fn is_feasible(data: &Data, pop: &Vec<Self>) -> bool;
}
