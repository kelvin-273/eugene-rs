use pathfinding::prelude::*;

pub struct BaseSolution {}

pub fn breeding_program(dist_array: Vec<u8>) -> BaseSolution {
    let res = astar(&dist_array, successors, heuristic, success);
    unimplemented!()
}

fn heuristic(state: &Vec<u8>) -> u8 {
    get_n_pop(state) + get_n_loci(state)
}

fn successors(_state: &Vec<u8>) -> Vec<(Vec<u8>, u8)> {
    // 
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
