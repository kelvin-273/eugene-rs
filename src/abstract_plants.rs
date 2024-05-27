use std::collections::HashMap;
use std::hash::Hash;
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
    fn cross(&self, x: &A) -> B;

    fn crosspoints(data: &Data) -> Box<dyn std::iter::Iterator<Item = Self>>;
}

#[derive(Debug, PartialEq, Eq)]
pub struct WGen<A, B> {
    pub genotype: A,
    pub history: Option<(Rc<WGam<A, B>>, Rc<WGam<A, B>>)>,
}

impl<A: Genotype<B>, B: Gamete<A>> WGen<A, B> {
    pub fn from_gametes(gx: Rc<WGam<A, B>>, gy: Rc<WGam<A, B>>) -> WGen<A, B> {
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

    pub fn gamete(&self) -> &B {
        &self.gamete
    }

    pub fn extract_first_from_ref(self: &Self) -> WGamS<A, B>
    where
        A: Clone,
        B: Clone,
    {
        WGamS {
            gamete: self.gamete.clone(),
            history: self.history.first().map(|x| x.clone().extract_first()),
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
            history: self.history.first().map(|x| x.clone().extract_first()),
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

    pub fn genotype(&self) -> &A {
        &self.genotype
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

impl<A, B> WGenS<A, B> {
    pub fn new(genotype: A) -> Self {
        Self {
            genotype,
            history: None,
        }
    }

    pub fn genotype(&self) -> &A {
        &self.genotype
    }
}

impl<A, B> WGamS<A, B> {
    pub fn new(gamete: B) -> Self {
        Self {
            gamete,
            history: None,
        }
    }

    pub fn gamete(&self) -> &B {
        &self.gamete
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

///////////////////////////////
//  New version of WGen API  //
///////////////////////////////

#[derive(Debug, PartialEq, Eq)]
pub struct WGenS2<A, B> {
    head: Rc<WGenSCell<A, B>>,
}

#[derive(Debug, PartialEq, Eq)]
pub struct WGamS2<A, B> {
    head: Rc<WGamSCell<A, B>>,
}

#[derive(Debug, PartialEq, Eq)]
struct WGenSCell<A, B> {
    genotype: A,
    history: Option<(Rc<WGamSCell<A, B>>, Rc<WGamSCell<A, B>>)>,
}

#[derive(Debug, PartialEq, Eq)]
struct WGamSCell<A, B> {
    gamete: B,
    history: Option<Rc<WGenSCell<A, B>>>,
}

impl<A, B> WGenS2<A, B> {
    pub fn new(x: A) -> Self {
        Self {
            head: Rc::new(WGenSCell {
                genotype: x,
                history: None,
            }),
        }
    }

    pub fn genotype(&self) -> &A {
        &self.head.genotype
    }

    pub fn history(&self) -> Option<(WGamS2<A, B>, WGamS2<A, B>)> {
        self.head
            .history
            .as_ref()
            .map(|(wgl, wgr)| (WGamS2 { head: wgl.clone() }, WGamS2 { head: wgr.clone() }))
    }

    pub fn from_gametes(wgl: &WGamS2<A, B>, wgr: &WGamS2<A, B>) -> Self
    where
        A: Genotype<B>,
        B: Gamete<A>,
    {
        Self {
            head: Rc::new(WGenSCell {
                genotype: A::from_gametes(wgl.gamete(), wgr.gamete()),
                history: Some((wgl.head.clone(), wgr.head.clone())),
            }),
        }
    }

    pub fn from_gametes_with_store(
        h: &mut HashMap<A, WGenS2<A, B>>,
        wgl: &WGamS2<A, B>,
        wgr: &WGamS2<A, B>,
    ) -> Self
    where
        A: Genotype<B> + Sized + PartialEq + Eq + Hash + Clone,
        B: Gamete<A> + Clone,
    {
        let x = A::from_gametes(wgl.gamete(), wgr.gamete());
        h.entry(x.clone())
            .or_insert_with(|| Self {
                head: Rc::new(WGenSCell {
                    genotype: x,
                    history: Some((wgl.head.clone(), wgr.head.clone())),
                }),
            })
            .clone()
    }

    pub fn cross<K: Fn(&A) -> B>(&self, k: K) -> WGamS2<A, B> {
        WGamS2 {
            head: Rc::new(WGamSCell {
                gamete: k(self.genotype()),
                history: Some(self.head.clone()),
            }),
        }
    }

    pub fn cross_with_store(
        &self,
        h: &mut HashMap<B, WGamS2<A, B>>,
        k: &dyn Fn(&A) -> B,
    ) -> WGamS2<A, B>
    where
        B: Sized + PartialEq + Eq + Hash + Clone,
    {
        let g = k(self.genotype());
        h.entry(g.clone())
            .or_insert_with(|| WGamS2 {
                head: Rc::new(WGamSCell {
                    gamete: g,
                    history: Some(self.head.clone()),
                }),
            })
            .clone()
    }

    pub fn generations(&self) -> usize {
        unimplemented!()
    }

    pub fn crossings(&self) -> usize {
        unimplemented!()
    }
}

impl<A, B> WGamS2<A, B> {
    pub fn new(gx: B) -> Self {
        Self {
            head: Rc::new(WGamSCell {
                gamete: gx,
                history: None,
            }),
        }
    }

    pub fn gamete(&self) -> &B {
        &self.head.gamete
    }

    pub fn history(&self) -> Option<WGenS2<A, B>> {
        self.head
            .history
            .as_ref()
            .map(|wx| WGenS2 { head: wx.clone() })
    }

    pub fn new_from_genotype(gx: B, wx: WGenS2<A, B>) -> Self {
        Self {
            head: Rc::new(WGamSCell {
                gamete: gx,
                history: Some(wx.head.clone()),
            }),
        }
    }

    pub fn generations(&self) -> usize {
        unimplemented!()
    }

    pub fn crossings(&self) -> usize {
        unimplemented!()
    }
}

impl<A, B> Clone for WGenS2<A, B> {
    fn clone(&self) -> Self {
        Self {
            head: self.head.clone(),
        }
    }
}

impl<A, B> Clone for WGamS2<A, B> {
    fn clone(&self) -> Self {
        Self {
            head: self.head.clone(),
        }
    }
}
