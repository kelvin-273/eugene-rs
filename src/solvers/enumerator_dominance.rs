use crate::abstract_plants::*;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::rc::Rc;

/// Runs the enumeration with dominance pruning from the pop_0 to the ideotype.
/// Dominance pruning is performed on gametes
/// Returns the breeding_program rooted at the ideotype
pub fn breeding_program<A, B, K, D, Data>(
    pop_0: Vec<A>,
    ideotype: A,
    data: Data,
) -> Option<WGen<A, B>>
where
    A: Genotype<B> + PartialEq + Clone + Feasible<Data>,
    B: Gamete<A> + Hash + Eq + Clone,
    K: Crosspoint<A, B, Data>,
    D: Dominance<B>,
{
    breeding_program_timeout_gametes::<A, B, K, D, Data>(pop_0, ideotype, None, data)
}

pub fn breeding_program_timeout_genotypes<A, B, K, D, Data>(
    pop_0: Vec<A>,
    ideotype: A,
    gen_max_opt: Option<usize>,
    data: Data,
) -> Option<WGen<A, B>>
where
    A: Genotype<B> + PartialEq + Clone + Feasible<Data>,
    B: Gamete<A> + Hash + Eq + Clone,
    K: Crosspoint<A, B, Data>,
    D: Dominance<A>,
{
    if !A::is_feasible(&data, &pop_0) {
        return None;
    }

    let mut pop: Vec<Rc<WGen<A, B>>> = filter_non_dominating::<A, D>(pop_0)
        .iter()
        .map(|x| Rc::new(x.lift_a()))
        .collect();
    let mut gen_i = 0;
    let timeout_check = |t: usize| gen_max_opt.map_or(true, |t_max| t < t_max);

    while !pop.iter().any(|x| x.genotype == ideotype) && timeout_check(gen_i) {
        // filter filter_non_dominating genotypes
        pop = filter_non_dominating_key::<Rc<WGen<A, B>>, A, _, D>(pop, |wx| wx.genotype.clone());

        // create all gametes
        let pop_g: Vec<Rc<WGam<A, B>>> = pop
            .iter()
            .flat_map(|wx| {
                K::crosspoints(&data).map(|k| {
                    Rc::new(WGam {
                        gamete: k.cross(&wx.genotype),
                        history: vec![wx.clone()],
                    })
                })
            })
            .collect();

        // create all genotypes from genotypes
        pop = (0..pop_g.len())
            .flat_map(|i| (i..pop_g.len()).map(move |j| (i, j)))
            //.flat_map(|i| (0..pop_g.len()).map(move |j| (i, j)))
            .map(|(i, j)| {
                let wgx = &pop_g[i];
                let wgy = &pop_g[j];
                Rc::new(WGen {
                    genotype: A::from_gametes(&wgx.gamete.clone(), &wgy.gamete.clone()),
                    history: Some((wgx.clone(), wgy.clone())),
                })
            })
            .collect();
        gen_i += 1;
    }

    pop.into_iter()
        .find(|x| x.genotype == ideotype)
        .map(|x| WGen {
            genotype: x.genotype.clone(),
            history: x.history.to_owned(),
        })
}

pub fn breeding_program_timeout_gametes<A, B, K, D, Data>(
    pop_0: Vec<A>,
    ideotype: A,
    gen_max_opt: Option<usize>,
    data: Data,
) -> Option<WGen<A, B>>
where
    A: Genotype<B> + PartialEq + Clone + Feasible<Data>,
    B: Gamete<A> + Hash + Eq + Clone,
    K: Crosspoint<A, B, Data>,
    D: Dominance<B>,
{
    if !A::is_feasible(&data, &pop_0) {
        return None;
    }

    let mut pop: Vec<Rc<WGen<A, B>>> = pop_0.iter().map(|x| Rc::new(x.lift_a())).collect();
    let mut gen_i = 0;
    let timeout_check = |t: usize| gen_max_opt.map_or(true, |t_max| t < t_max);

    while !pop.iter().any(|x| x.genotype == ideotype) && timeout_check(gen_i) {
        let mut h: HashMap<B, Vec<Rc<WGen<A, B>>>> = HashMap::new();

        //K::crosspoints(&data).for_each(|kx| println!("k: {:?}", kx));
        // collect all possible wrapped gametes
        for x in &pop {
            let s: HashSet<B> = K::crosspoints(&data)
                .map(|k| k.cross(&x.genotype))
                .collect();
            for g in filter_non_dominating::<B, D>(s) {
                h.entry(g.clone())
                    .and_modify(|v| (*v).push(x.clone()))
                    .or_insert(vec![x.clone()]);
            }
        }
        // remove dominated gametes
        let next_gametes: Vec<Rc<WGam<A, B>>> =
            filter_non_dominating_key::<_, _, _, D>(h, |(g, v)| g.clone())
                .into_iter()
                .map(|(g, v)| {
                    Rc::new(WGam {
                        gamete: g,
                        history: v,
                    })
                })
                .collect();
        //println!("{:?}", next_gametes.iter().map(|wgx| &wgx.gamete).collect::<Vec<&B>>());
        // construct next population from collected gametes
        pop = {
            next_gametes
                .iter()
                .enumerate()
                .flat_map(|(i, wg1)| {
                    next_gametes[i..].iter().enumerate().map(|(j, wg2)| {
                        Rc::new(WGen {
                            genotype: A::from_gametes(&wg1.gamete, &wg2.gamete),
                            history: Some((wg1.clone(), wg2.clone())),
                        })
                    })
                })
                .collect()
        };
        println!("gen_i: {},\tpop.len: {}", gen_i, pop.len());
        gen_i += 1;
    }
    pop.into_iter()
        .find(|x| x.genotype == ideotype)
        .map(|x| WGen {
            genotype: x.genotype.clone(),
            history: x.history.to_owned(),
        })
}

pub fn filter_non_dominating<T, D>(s: impl IntoIterator<Item = T>) -> Vec<T>
where
    D: Dominance<T>,
{
    let v: Vec<T> = s.into_iter().collect();
    let mut keeps: Vec<bool> = vec![true; v.len()];
    for i in 0..v.len() {
        let x = &v[i];
        for j in i + 1..v.len() {
            let y = &v[j];
            if D::dom(&x, &y) {
                keeps[j] = false;
            } else if D::dom(y, x) {
                keeps[i] = false;
                break;
            }
        }
    }
    v.into_iter()
        .enumerate()
        .filter_map(|(i, x)| match keeps[i] {
            true => Some(x),
            false => None,
        })
        .collect()
}

pub fn filter_non_dominating_key<T, U, F, D>(s: impl IntoIterator<Item = T>, key: F) -> Vec<T>
where
    F: Fn(&T) -> U,
    D: Dominance<U>,
{
    let v: Vec<T> = s.into_iter().collect();
    let mut keeps: Vec<bool> = vec![true; v.len()];
    for i in 0..v.len() {
        let x = &v[i];
        for j in i + 1..v.len() {
            let y = &v[j];
            let _x = key(x);
            let _y = key(y);
            if D::dom(&_x, &_y) {
                keeps[j] = false;
            } else if D::dom(&_y, &_x) {
                keeps[i] = false;
                break;
            }
        }
    }
    v.into_iter()
        .enumerate()
        .filter_map(|(i, x)| match keeps[i] {
            true => Some(x),
            false => None,
        })
        .collect()
}

pub fn filter_non_dominating_fn<T>(
    s: impl IntoIterator<Item = T>,
    func: impl Fn(&T, &T) -> bool,
) -> Vec<T> {
    let v: Vec<T> = s.into_iter().collect();
    let mut keeps: Vec<bool> = vec![true; v.len()];
    for i in 0..v.len() {
        let x = &v[i];
        for j in i + 1..v.len() {
            let y = &v[j];
            if func(x, y) {
                keeps[j] = false;
            } else if func(y, x) {
                keeps[i] = false;
                break;
            }
        }
    }
    v.into_iter()
        .enumerate()
        .filter_map(|(i, x)| match keeps[i] {
            true => Some(x),
            false => None,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plants::bit_array::*;
    use rand::prelude::*;

    #[test]
    fn filter_non_dominating_test() {
        use rand::prelude::*;

        impl Dominance<i32> for i32 {
            fn dom(x: &i32, y: &i32) -> bool {
                x >= y
            }
        }

        for _ in (0..100) {
            let n: usize = random::<usize>() % 100 + 1;
            let v: Vec<i32> = (0..n).map(|_| random()).collect();
            let res = filter_non_dominating::<i32, i32>(v.to_owned());
            assert_eq!(res.len(), 1);
            assert_eq!(v.iter().max(), res.first());
        }
    }

    struct Pair {
        a: i32,
        b: i32,
    }
    impl Dominance<Pair> for Pair {
        fn dom(x: &Pair, y: &Pair) -> bool {
            x.a >= y.a && x.b >= y.b
        }
    }
    impl Pair {
        pub fn new(a: i32, b: i32) -> Self {
            Self { a, b }
        }
        fn random() -> Self {
            Self {
                a: random(),
                b: random(),
            }
        }
    }

    #[test]
    fn filter_non_dominating_pair_test() {
        for _ in (0..100) {
            let n: usize = random::<usize>() % 100 + 1;
            let v: Vec<Pair> = (0..n).map(|_| Pair::random()).collect();
            let res = filter_non_dominating::<Pair, Pair>(v);
            assert!(res.len() > 0);
            for i in (0..res.len()) {
                for j in (0..i) {
                    let x = &res[i];
                    let y = &res[j];
                    assert!((x.a < y.a || x.b < y.b) && (x.a > y.a || x.b > y.b));
                }
            }
        }
    }

    #[test]
    fn simplest_test() {
        let x = SingleChromGenotype::from_str("010", "101");
        assert_ne!(
            breeding_program_timeout_gametes::<
                SingleChromGenotype,
                SingleChromGamete,
                CrosspointBitVec,
                DomGamete,
                usize,
            >(vec![x], SingleChromGenotype::ideotype(3), Some(10), 3),
            None
        );
    }

    #[test]
    fn simple_test() {
        let x = SingleChromGenotype::from_str("01010101", "10101010");
        assert_ne!(
            breeding_program_timeout_gametes::<
                SingleChromGenotype,
                SingleChromGamete,
                CrosspointBitVec,
                DomGamete,
                usize,
            >(vec![x], SingleChromGenotype::ideotype(8), Some(10), 8),
            None
        );
    }

    #[test]
    fn time_test() {
        let n_loci = 8;
        let pop_0 = (0..n_loci)
            .map(|i| {
                SingleChromGenotype::new(
                    (0..n_loci)
                        .map(|j| if i == j { (true, true) } else { (false, false) })
                        .collect(),
                )
            })
            .collect();
        assert_ne!(
            breeding_program_timeout_gametes::<
                SingleChromGenotype,
                SingleChromGamete,
                CrosspointBitVec,
                DomGamete,
                usize,
            >(
                pop_0,
                SingleChromGenotype::ideotype(n_loci),
                Some(4),
                n_loci
            ),
            None
        );
    }
}
