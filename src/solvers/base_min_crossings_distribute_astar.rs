use crate::solution::BaseSolution;
use itertools::Itertools;
use pathfinding::prelude::*;

pub fn breeding_program(dist_array: Vec<u8>) -> BaseSolution {
    let _res = astar(&dist_array, successors, heuristic, success);
    unimplemented!()
}

fn heuristic(state: &Vec<u8>) -> u8 {
    get_n_pop(state) + get_n_loci(state)
}

fn successors(state: &Vec<u8>) -> Vec<(Vec<u8>, u8)> {
    let n_pop = get_n_pop(state);
    let mut out: Vec<(Vec<u8>, u8)> = vec![];
    for gx in 0..n_pop - 1 {
        for gy in gx + 1..n_pop {
            // Note where the segment joins are and which direction they're going
            let segjoins = state
                .iter()
                .tuple_windows()
                .enumerate()
                .filter(|(p, (&i, &j))| (i, j) == (gx, gy) || (i, j) == (gy, gx))
                .map(|(p, (&i, &j))| match (i, j) == (gx, gy) {
                    true => (p, CrossoverDirection::Forward),
                    false => (p, CrossoverDirection::Backward),
                })
                .collect_vec();

            if segjoins.is_empty() {
                todo!();
            } else {
                todo!();
            }
            todo!();
        }
    }
    unimplemented!()
}

fn success(state: &Vec<u8>) -> bool {
    state == &vec![0]
}

fn get_n_pop(state: &Vec<u8>) -> u8 {
    state.iter().max().map(|res| *res + 1).unwrap_or(0)
}

fn get_n_loci(state: &Vec<u8>) -> u8 {
    state.len() as u8
}

enum CrossoverDirection {
    Forward,
    Backward,
}
