use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solution::BaseSolution;
use pyo3::prelude::*;
use std::hash::Hash;
use std::rc::Rc;

/// Runs a breeding program given `n_loci` and `pop_0` where `pop_0` is a population of single
/// chromosome diploid genotypes with `n_loci` loci.
#[pyo3::pyfunction]
pub fn breeding_program_python(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<bool>>>,
) -> PyResult<
    Option<(
        Vec<Vec<Vec<i32>>>,
        Vec<&'static str>,
        Vec<usize>,
        Vec<usize>,
        usize,
    )>,
> {
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    let pop_0 = pop_0
        .iter()
        .map(|x| {
            SingleChromGenotype::new(
                x[0].iter()
                    .zip(x[1].iter())
                    .map(|(a, b)| (*a, *b))
                    .collect(),
            )
        })
        .collect();
    let res = breeding_program::<SingleChromGenotype, SingleChromGamete, CrosspointBitVec, usize>(
        pop_0, ideotype, n_loci,
    );
    match res {
        None => Ok(None),
        Some(x_star) => {
            let sol = BaseSolution::min_gen_from_wgen(n_loci, &x_star);
            Ok(Some((
                sol.tree_data,
                sol.tree_type,
                sol.tree_left,
                sol.tree_right,
                sol.objective,
            )))
        }
    }
}

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plants::u64_and_u32;

    #[test]
    fn function_name_test() {
        let _res = breeding_program::<u64, u32, u64_and_u32::CrosspointSingleU64U32, usize>(vec![
            (0b001 << 32) | 0b001,
            (0b010 << 32) | 0b010,
            (0b100 << 32) | 0b100,
        ],
        (0b111 << 32) | 0b111,
        3);
    }
}
