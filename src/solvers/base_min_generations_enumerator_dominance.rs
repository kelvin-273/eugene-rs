use crate::abstract_plants::{Crosspoint, Dominance, WGam, WGen};
use crate::plants::bit_array::*;
use crate::solution::{BaseSolution, Objective, PyBaseSolution};
use std::collections::HashSet;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

/// Runs a breeding program given `n_loci` and `pop_0` where `pop_0` is a population of single
/// chromosome diploid genotypes with `n_loci` loci.
#[pyo3::pyfunction]
pub fn breeding_program_python(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<bool>>>,
    timeout: Option<u64>,
) -> PyBaseSolution {
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
        .collect::<Vec<_>>();
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = breeding_program(n_loci, &pop_0);
        tx.send(res)
    });
    let res = rx
        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
        .ok()
        .flatten();
    match res {
        None => Ok(None),
        Some(sol) => Ok(Some((
            sol.tree_data,
            sol.tree_type,
            sol.tree_left,
            sol.tree_right,
            sol.objective,
        ))),
    }
}

pub fn breeding_program(n_loci: usize, pop_0: &[SingleChromGenotype]) -> Option<BaseSolution> {
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    if pop_0.contains(&ideotype) {
        return Some(WGen::new(ideotype).to_base_sol(n_loci, Objective::Generations));
    }
    let pop_0: Vec<WGen<SingleChromGenotype, _>> =
        pop_0.iter().map(|x| WGen::new(x.clone())).collect();
    let mut h: HashSet<SingleChromGamete> = HashSet::new();
    let mut contained_gametes: Vec<WGam<SingleChromGenotype, SingleChromGamete>> = vec![];
    for wx in pop_0 {
        for gx in CrosspointBitVec::crosspoints(&n_loci).map(|k| k.cross(wx.genotype())) {
            if !h.contains(&gx) {
                let wgx = WGam::new_from_genotype(gx.clone(), wx.clone());
                h.insert(gx);
                contained_gametes.push(wgx);
            }
        }
    }
    contained_gametes = filter_non_dominating_fn(contained_gametes, |wgx, wgy| {
        DomGamete::dom(wgx.gamete(), wgy.gamete())
    });
    let ideotype_gamete = SingleChromGamete::ideotype(n_loci);
    while !h.contains(&ideotype_gamete) {
        let mut v = vec![];
        for i in 0..contained_gametes.len() {
            for j in i + 1..contained_gametes.len() {
                let wgx = &contained_gametes[i];
                let wgy = &contained_gametes[j];
                let wz = WGen::from_gametes(wgx, wgy);
                for k in CrosspointBitVec::crosspoints(&n_loci) {
                    let gz = k.cross(wz.genotype());
                    if !h.contains(&gz) {
                        let wgz = WGam::new_from_genotype(gz.clone(), wz.clone());
                        h.insert(gz);
                        v.push(wgz);
                    }
                }
            }
        }
        contained_gametes =
            filter_non_dominating_fn(v, |wgx, wgy| DomGamete::dom(wgx.gamete(), wgy.gamete()));
    }
    let wg_star = contained_gametes.first().expect("hash map doesn't have g*");
    let wx_star = WGen::from_gametes(wg_star, wg_star);
    Some(wx_star.to_base_sol(n_loci, Objective::Generations))
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
            if D::dom(x, y) {
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
    F: Fn(&T) -> &U,
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
            if D::dom(_x, _y) {
                keeps[j] = false;
            } else if D::dom(_y, _x) {
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
    use rand::prelude::*;

    #[test]
    fn filter_non_dominating_test() {
        impl Dominance<i32> for i32 {
            fn dom(x: &i32, y: &i32) -> bool {
                x >= y
            }
        }

        for _ in 0..100 {
            let n: usize = random::<usize>() % 100 + 1;
            let v: Vec<i32> = (0..n).map(|_| random()).collect();
            let res = filter_non_dominating::<i32, i32>(v.to_owned());
            assert_eq!(res.len(), 1);
            assert_eq!(v.iter().max(), res.first());
        }
    }

    #[derive(Debug, PartialEq, Eq)]
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
        let pairs = vec![
            Pair::new(2, 3),
            Pair::new(4, 5),
            Pair::new(0, 9),
            Pair::new(5, 1),
            Pair::new(5, 0),
            Pair::new(2, 5),
            Pair::new(9, 7),
            Pair::new(5, 4),
            Pair::new(0, 8),
            Pair::new(5, 5),
            Pair::new(6, 2),
            Pair::new(4, 4),
            Pair::new(4, 3),
            Pair::new(2, 2),
            Pair::new(9, 1),
            Pair::new(8, 8),
            Pair::new(4, 3),
            Pair::new(2, 1),
            Pair::new(0, 4),
        ];
        assert_eq!(
            vec![Pair::new(0, 9), Pair::new(9, 7), Pair::new(8, 8),],
            filter_non_dominating::<Pair, Pair>(pairs)
        );

        let pairs = vec![
            Pair::new(2, 3),
            Pair::new(4, 5),
            Pair::new(0, 9),
            Pair::new(5, 1),
            Pair::new(5, 0),
            Pair::new(2, 5),
            Pair::new(5, 4),
            Pair::new(0, 8),
            Pair::new(5, 5),
            Pair::new(6, 2),
            Pair::new(4, 4),
            Pair::new(4, 3),
            Pair::new(2, 2),
            Pair::new(9, 1),
            Pair::new(4, 3),
            Pair::new(2, 1),
            Pair::new(0, 4),
        ];
        assert_eq!(
            vec![
                Pair::new(0, 9),
                Pair::new(5, 5),
                Pair::new(6, 2),
                Pair::new(9, 1),
            ],
            filter_non_dominating::<Pair, Pair>(pairs)
        );

        for _ in 0..100 {
            let n: usize = random::<usize>() % 100 + 1;
            let v: Vec<Pair> = (0..n).map(|_| Pair::random()).collect();
            let res = filter_non_dominating::<Pair, Pair>(v);
            assert!(res.len() > 0);
            for i in 0..res.len() {
                for j in 0..i {
                    let x = &res[i];
                    let y = &res[j];
                    assert!((x.a < y.a || x.b < y.b) && (x.a > y.a || x.b > y.b));
                }
            }
        }
    }

    #[test]
    fn bitvec_breeding_program_4_zigzag_test() {
        use crate::plants::bit_array::*;
        let n_loci = 4;
        let pop_0 = vec![SingleChromGenotype::from_str("0101", "1010")];
        assert_eq!(
            2,
            breeding_program(n_loci, &pop_0)
                .expect("infeasible")
                .objective
        );
    }

    #[test]
    fn bitvec_breeding_program_4_dist_zigzag_test() {
        use crate::plants::bit_array::*;
        let n_loci = 4;
        let pop_0 = vec![
            SingleChromGenotype::from_str("0101", "0101"),
            SingleChromGenotype::from_str("1010", "1010"),
        ];
        let op_sol = breeding_program(n_loci, &pop_0);
        dbg!(&op_sol);
        assert_eq!(3, op_sol.expect("infeasible").objective);
    }

    #[test]
    fn bitvec_breeding_program_dist_test() {
        use crate::extra::instance_generators;
        fn tester(actual: usize, v: &[usize]) {
            let pop_0 = instance_generators::distarray_to_homo(v);
            let n_loci = v.len();
            let op_sol = breeding_program(n_loci, &pop_0);
            dbg!(&op_sol);
            assert_eq!(actual, op_sol.expect("infeasible").objective);
        }
        tester(2, &[0, 1]);
        tester(3, &[0, 1, 0]);
        tester(3, &[0, 1, 2, 0]);
        tester(4, &[0, 1, 0, 2, 1]);
    }

    #[test]
    fn enumerator_dominance_test() {
        use rand::prelude::*;
        let n_loci = 10;
        let n_pop = 5;
        let mut rng = thread_rng();
        for _ in 0..10 {
            let pop_0 = SingleChromGenotype::init_pop_random(&mut rng, n_loci, n_pop);
            let sol1 = breeding_program(n_loci, &pop_0);
            let sol2 =
                crate::solvers::base_min_generations_enumerator::breeding_program(n_loci, &pop_0);
            assert!(sol1.is_some());
            assert!(sol2.is_some());
            assert_eq!(sol1.unwrap().objective, sol2.unwrap().objective);
        }
    }
}
