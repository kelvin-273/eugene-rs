use crate::abstract_plants::*;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::rc::Rc;

fn breeding_program<A, B, K, D>(pop_0: Vec<A>, ideotype: A) -> Option<WGen<A, B>>
where
    A: Genotype<B> + PartialEq + Clone,
    B: Gamete<A> + Hash + Eq + Clone,
    K: Crosspoint<A, B>,
    D: Dominance<B>,
{
    let mut pop: Vec<Rc<WGen<A, B>>> = pop_0.iter().map(|x| Rc::new(x.lift_a())).collect();
    if !is_feasible(pop_0) {
        return None;
    }
    while !pop.iter().any(|x| x.genotype == ideotype) {
        let mut h: HashMap<B, Vec<Rc<WGen<A, B>>>> = HashMap::new();

        // collect all possible wrapped gametes
        for x in &pop {
            let s: HashSet<B> = K::crosspoints().map(|k| k.cross(&x.genotype)).collect();
            for g in filter_non_dominating_gametes::<B, D>(s) {
                h.entry(g.clone())
                    .and_modify(|v| (*v).push(x.clone()))
                    .or_insert(vec![x.clone()]);
            }
        }
        let next_gametes: Vec<Rc<WGam<A, B>>> = h
            .into_iter()
            .map(|(g, v)| {
                Rc::new(WGam {
                    gamete: g,
                    history: v,
                })
            })
            .collect();
        // construct next population
        pop = next_gametes
            .iter()
            .flat_map(|wg1| {
                next_gametes.iter().map(|wg2| {
                    Rc::new(WGen {
                        genotype: A::from_gametes(&wg1.gamete, &wg2.gamete),
                        history: Some((wg1.clone(), wg2.clone())),
                    })
                })
            })
            .collect()
    }
    pop.into_iter()
        .find(|x| x.genotype == ideotype)
        .map(|x| WGen {
            genotype: x.genotype.clone(),
            history: x.history.to_owned(),
        })
}

fn is_feasible<T>(x: T) -> bool {
    unimplemented!()
}

fn filter_non_dominating_gametes<T, D>(s: impl IntoIterator<Item = T>) -> Vec<T>
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
