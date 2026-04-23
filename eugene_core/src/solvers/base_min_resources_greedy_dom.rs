#![deny(unused)]
use std::collections::{BTreeMap, BTreeSet};

use crate::abstract_plants::{Crosspoint, Dominance, WGam, WGen};
use crate::extra::resources::{cost_of_crossing, RecRate};
use crate::plants::bit_array::{
    CrosspointBitVec, DomGamete, SingleChromGamete, SingleChromGenotype,
};
use crate::plants::dist_array::DistArray;
use crate::solution::CrossingSchedule;

type WGe = WGen<SingleChromGenotype, SingleChromGamete>;
type WGa = WGam<SingleChromGenotype, SingleChromGamete>;

pub fn breeding_program_distribute(
    xs: &DistArray,
    rec_rate: &RecRate,
    gamma: f64,
) -> Option<CrossingSchedule> {
    let n_loci = xs.n_loci();
    let pop_0 = SingleChromGenotype::init_pop_distribute(xs);
    breeding_program(n_loci, &pop_0, rec_rate, gamma)
}

pub fn breeding_program(
    n_loci: usize,
    pop_0: &[SingleChromGenotype],
    rec_rate: &RecRate,
    gamma: f64,
) -> Option<CrossingSchedule> {
    assert!((0.0..=1.0).contains(&gamma));
    assert!(!pop_0.is_empty());
    assert!(pop_0.iter().all(|x| x.n_loci() == n_loci));

    let target = SingleChromGenotype::ideotype(n_loci);
    if pop_0.contains(&target) {
        return Some(WGen::new(target).into());
    }
    let mut pop: Vec<WGe> = pop_0.iter().map(|x| WGen::new(x.clone())).collect();
    let mut d_gen: BTreeMap<SingleChromGenotype, WGe> =
        pop_0.iter().cloned().zip(pop.iter().cloned()).collect();
    let mut d_gam: BTreeMap<SingleChromGamete, (WGa, f64)> = BTreeMap::new();

    // Generate initial population of non-dominated gametes
    for wx in &pop {
        for k_x in CrosspointBitVec::crosspoints(&n_loci) {
            let gx = k_x.cross(wx.genotype());
            let rgx = rec_rate.probability_gamete(wx.genotype(), &gx);
            d_gam
                .entry(gx.clone())
                .and_modify(|(wgy, rgy)| {
                    if rgx > *rgy {
                        *wgy = WGa::new_from_genotype(gx.clone(), wx.clone());
                        *rgy = rgx;
                    }
                })
                .or_insert_with(|| (WGa::new_from_genotype(gx, wx.clone()), rgx));
        }
    }

    loop {
        //dbg!("Current population: {:?}", pop.iter().map(|w| w.genotype()).collect::<Vec<_>>());
        //dbg!("Current gametes: {:?}", d_gam.iter().map(|(k, (_, p))| (k, p)).collect::<Vec<_>>());
        let wz = min_cost_crossing(rec_rate, gamma, &d_gen, &d_gam)
            .expect("no new supporting crossing available");
        if wz.genotype() == &target {
            return Some(wz.into());
        }
        d_gen.insert(wz.genotype().clone(), wz.clone());
        pop.push(wz.clone());
        augment_d_gam(wz, rec_rate, &mut d_gam, n_loci);
    }
}

fn min_cost_crossing(
    rec_rate: &RecRate,
    gamma: f64,
    d_gen: &BTreeMap<SingleChromGenotype, WGe>,
    d_gam: &BTreeMap<SingleChromGamete, (WGa, f64)>,
) -> Option<WGe> {
    let mut candidates = vec![];

    for (wgx, rgx) in d_gam.values() {
        for (wgy, rgy) in d_gam.values() {
            if wgx.gamete() < wgy.gamete() {
                continue;
            }
            if wgx.history().is_none() || wgy.history().is_none() {
                continue;
            }

            let wz = WGe::from_gametes(wgx, wgy);
            if d_gen.contains_key(wz.genotype()) {
                continue;
            }

            let z_cost = cost_of_crossing(gamma, *rgx, *rgy);

            for (gz, rgz) in unique_gametes_with_probs(rec_rate, wz.genotype()) {
                let mut flag = true;
                for (gx, (_, rgx)) in d_gam {
                    if DomGamete::dom(gx, &gz) && *rgx >= rgz {
                        flag = false;
                        break;
                    }
                }
                if flag {
                    candidates.push((wz.clone(), z_cost));
                    break;
                }
            }
        }
    }

    candidates
        .into_iter()
        .min_by(|(_, x_cost), (_, y_cost)| x_cost.total_cmp(y_cost))
        .map(|(wz, _)| wz)
}

pub fn unique_gametes_with_probs(
    rec_rate: &RecRate,
    z: &SingleChromGenotype,
) -> Vec<(SingleChromGamete, f64)> {
    let n_loci = z.n_loci();
    CrosspointBitVec::crosspoints(&n_loci)
        .map(|k| k.cross(z))
        .collect::<BTreeSet<SingleChromGamete>>()
        .into_iter()
        .map(|gz| {
            let rgz = rec_rate.probability_gamete(z, &gz);
            (gz, rgz)
        })
        .collect()
}

fn augment_d_gam(
    wz: WGe,
    rec_rate: &RecRate,
    d_gam: &mut BTreeMap<SingleChromGamete, (WGa, f64)>,
    n_loci: usize,
) {
    for gz in CrosspointBitVec::crosspoints(&n_loci).map(|k| k.cross(wz.genotype())) {
        let rgz = rec_rate.probability_gamete(wz.genotype(), &gz);
        d_gam
            .entry(gz.clone())
            .and_modify(|(wgy, rgy)| {
                if rgz > *rgy {
                    *wgy = WGa::new_from_genotype(gz.clone(), wz.clone());
                    *rgy = rgz
                }
            })
            .or_insert_with(|| (WGa::new_from_genotype(gz, wz.clone()), rgz));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plants::dist_array::dist_array;

    #[test]
    fn breeding_program_trivial_test() {
        let n_loci = 1;
        let pop_0 = SingleChromGenotype::init_pop_distribute(&dist_array![0]);
        let rec_rate = RecRate::new(vec![]);
        let gamma = 0.99;
        let result = breeding_program(n_loci, &pop_0, &rec_rate, gamma)
            .expect("breeding_program returned None");
        assert_eq!(result.resources(&rec_rate, gamma), 0);
        assert_eq!(result.crossings(), 0)
    }

    #[test]
    fn breeding_program_2_loci_test() {
        let n_loci = 2;
        let pop_0 = SingleChromGenotype::init_pop_distribute(&dist_array![0, 1]);
        let rec_rate = RecRate::new(vec![0.1]);
        let gamma = 0.99;
        let result = breeding_program(n_loci, &pop_0, &rec_rate, gamma)
            .expect("breeding_program returned None");
        assert_eq!(result.resources(&rec_rate, gamma), 108);
        assert_eq!(result.crossings(), 3);
        assert_eq!(result.generations(), 3);
    }

    #[test]
    fn unique_gametes_with_probs_2_loci_test() {
        let z = SingleChromGenotype::from_str("10", "01");
        let rec_rate = RecRate::new(vec![0.1]);
        let mut results = unique_gametes_with_probs(&rec_rate, &z);
        results.sort_by(|x, y| x.0.cmp(&y.0));
        assert_eq!(
            results,
            Vec::from([
                (SingleChromGamete::from_str("00"), 0.05),
                (SingleChromGamete::from_str("01"), 0.45),
                (SingleChromGamete::from_str("10"), 0.45),
                (SingleChromGamete::from_str("11"), 0.05),
            ])
        )
    }

    #[test]
    fn min_cost_crossing_2_loci_test() {
        let rec_rate = RecRate::new(vec![0.1]);
        let gamma = 0.99;

        let x = SingleChromGenotype::from_str("10", "10");
        let y = SingleChromGenotype::from_str("01", "01");

        let gx = SingleChromGamete::from_str("10");
        let gy = SingleChromGamete::from_str("01");

        let rgx = rec_rate.probability_gamete(&x, &gx);
        let rgy = rec_rate.probability_gamete(&y, &gy);

        let wx = WGe::new(x.clone());
        let wy = WGe::new(y.clone());

        let wgx = WGa::new_from_genotype(gx.clone(), wx.clone());
        let wgy = WGa::new_from_genotype(gy.clone(), wy.clone());

        let d_gam = BTreeMap::from([(gx, (wgx.clone(), rgx)), (gy, (wgy.clone(), rgy))]);
        let d_gen = BTreeMap::from([(x, wx), (y, wy)]);
        let results = min_cost_crossing(&rec_rate, gamma, &d_gen, &d_gam);
        assert_eq!(results, Some(WGe::from_gametes(&wgx, &wgy)));
    }

    #[test]
    fn breeding_program_2_loci_2_test() {
        let n_loci = 2;
        let pop_0 = vec![
            SingleChromGenotype::from_str("01", "00"),
            SingleChromGenotype::from_str("01", "10"),
            SingleChromGenotype::from_str("00", "00"),
            SingleChromGenotype::from_str("00", "01"),
        ];
        let rec_rate = RecRate::new(vec![0.3]);
        let gamma = 0.7;
        let result = breeding_program(n_loci, &pop_0, &rec_rate, gamma)
            .expect("breeding_program_distribute returned None");
        assert_eq!(result.resources(&rec_rate, gamma), 13);
        assert_eq!(result.crossings(), 4);
        assert_eq!(result.generations(), 3);
    }
}
