use crate::abstract_plants::*;
use std::hash::Hash;
use std::rc::Rc;

pub fn breeding_program<A, B, K, Data>(pop_0: Vec<A>, ideotype: A, data: Data) -> Option<WGen<A, B>>
where
    A: Genotype<B> + PartialEq + Clone + Feasible<Data>,
    B: Gamete<A> + Hash + Eq + Clone,
    K: Crosspoint<A, B, Data>,
{
    breeding_program_gamete::<A, B, K, Data>(pop_0, ideotype, data)
}

pub fn breeding_program_gamete<A, B, K, Data>(
    pop_0: Vec<A>,
    ideotype: A,
    data: Data,
) -> Option<WGen<A, B>>
where
    A: Genotype<B> + PartialEq + Clone + Feasible<Data>,
    B: Gamete<A> + Hash + Eq + Clone,
    K: Crosspoint<A, B, Data>,
{
    use std::collections::{HashMap, HashSet};

    let mut pop: Vec<Rc<WGen<A, B>>> = pop_0.iter().map(|x| Rc::new(x.lift_a())).collect();
    if !A::is_feasible(&data, &pop_0) {
        return None;
    }
    while !pop.iter().any(|x| x.genotype == ideotype) {
        let mut h: HashMap<B, Vec<Rc<WGen<A, B>>>> = HashMap::new();

        // collect all possible wrapped gametes
        for x in &pop {
            let s: HashSet<B> = K::crosspoints(&data)
                .map(|k| k.cross(&x.genotype))
                .collect();
            for g in &s {
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

pub fn breeding_program_naive<A, B, K, Data>(
    pop_0: Vec<A>,
    ideotype: A,
    data: Data,
) -> Option<WGen<A, B>>
where
    A: Genotype<B> + PartialEq + Clone + Feasible<Data>,
    B: Gamete<A> + Hash + Eq + Clone,
    K: Crosspoint<A, B, Data>,
{
    let mut pop: Vec<Rc<WGen<A, B>>> = pop_0.iter().map(|x| Rc::new(x.lift_a())).collect();
    let mut gen_i = 0;
    while !pop.iter().any(|x| x.genotype == ideotype) {
        pop = pop
            .iter()
            .zip(pop.iter())
            .zip(K::crosspoints(&data).zip(K::crosspoints(&data)))
            .map(|((x, y), (kx, ky))| {
                Rc::new(A::from_gametes(&kx.cross(&x.genotype), &ky.cross(&y.genotype)).lift_a())
            })
            .collect();
        gen_i += 1;
        println!("gen_i: {},\tpop.len: {}", gen_i, pop.len());
    }
    Some(ideotype.lift_a())
}
