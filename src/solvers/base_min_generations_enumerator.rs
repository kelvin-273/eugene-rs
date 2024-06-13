use crate::abstract_plants::{Crosspoint, WGam, WGen};
use crate::plants::bit_array::*;
use crate::solution::{BaseSolution, Objective, PyBaseSolution};
use std::collections::HashMap;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;
use pyo3::PyResult;

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

#[pyo3::pyfunction]
#[pyo3(name = "mingen_answer")]
pub fn mingen_answer_enumerator(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<bool>>>,
    timeout: Option<u64>,
) -> PyResult<Option<usize>> {
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
        let res = min_generations_enumerator(n_loci, &pop_0);
        tx.send(res)
    });
    let res = rx
        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
        .ok();
    Ok(res)
}

pub fn breeding_program(n_loci: usize, pop_0: &[SingleChromGenotype]) -> Option<BaseSolution> {
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    if pop_0.contains(&ideotype) {
        return Some(WGen::new(ideotype).to_base_sol(n_loci, Objective::Generations));
    }
    let pop_0: Vec<WGen<SingleChromGenotype, _>> =
        pop_0.iter().map(|x| WGen::new(x.clone())).collect();
    let mut h: HashMap<SingleChromGamete, WGam<SingleChromGenotype, SingleChromGamete>> =
        HashMap::new();
    for wx in pop_0 {
        for gx in CrosspointBitVec::crosspoints(&n_loci).map(|k| k.cross(wx.genotype())) {
            if !h.contains_key(&gx) {
                let wgx = WGam::new_from_genotype(gx.clone(), wx.clone());
                h.insert(gx, wgx);
            }
        }
    }
    let ideotype_gamete = SingleChromGamete::ideotype(n_loci);
    while !h.contains_key(&ideotype_gamete) {
        let contained_gametes: Vec<SingleChromGamete> = h.keys().cloned().collect();
        for i in 0..contained_gametes.len() {
            for j in i + 1..contained_gametes.len() {
                let wgx = h.get(&contained_gametes[i])?;
                let wgy = h.get(&contained_gametes[j])?;
                let wz = WGen::from_gametes(wgx, wgy);
                for k in CrosspointBitVec::crosspoints(&n_loci) {
                    let gz = k.cross(wz.genotype());
                    if !h.contains_key(&gz) {
                        let wgz = WGam::new_from_genotype(gz.clone(), wz.clone());
                        h.insert(gz, wgz.clone());
                    }
                }
            }
        }
    }
    let wg_star = h.get(&ideotype_gamete).expect("hash map doesn't have g*");
    let wx_star = WGen::from_gametes(wg_star, wg_star);
    Some(wx_star.to_base_sol(n_loci, Objective::Generations))
}

pub fn min_generations_enumerator(n_loci: usize, pop_0: &[SingleChromGenotype]) -> usize {
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    if pop_0.contains(&ideotype) {
        return 0;
    }
    let pop_0: Vec<WGen<SingleChromGenotype, _>> =
        pop_0.iter().map(|x| WGen::new(x.clone())).collect();
    let mut h: HashMap<SingleChromGamete, WGam<SingleChromGenotype, SingleChromGamete>> =
        HashMap::new();
    for wx in pop_0 {
        for gx in CrosspointBitVec::crosspoints(&n_loci).map(|k| k.cross(wx.genotype())) {
            if !h.contains_key(&gx) {
                let wgx = WGam::new_from_genotype(gx.clone(), wx.clone());
                h.insert(gx, wgx);
            }
        }
    }
    let mut gen = 0;
    let ideotype_gamete = SingleChromGamete::ideotype(n_loci);
    while !h.contains_key(&ideotype_gamete) {
        gen += 1;
        let contained_gametes: Vec<SingleChromGamete> = h.keys().cloned().collect();
        for i in 0..contained_gametes.len() {
            for j in i + 1..contained_gametes.len() {
                let wgx = h.get(&contained_gametes[i]).expect("gamete not found in hash map");
                let wgy = h.get(&contained_gametes[j]).expect("gamete not found in hash map");
                let wz = WGen::from_gametes(wgx, wgy);
                for k in CrosspointBitVec::crosspoints(&n_loci) {
                    let gz = k.cross(wz.genotype());
                    if !h.contains_key(&gz) {
                        let wgz = WGam::new_from_genotype(gz.clone(), wz.clone());
                        h.insert(gz, wgz.clone());
                    }
                }
            }
        }
    }
    gen + 1
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
