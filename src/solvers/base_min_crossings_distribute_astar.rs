use crate::solution::BaseSolution;
use itertools::Itertools;
use pathfinding::prelude::*;
use std::collections::HashSet;

pub fn breeding_program(dist_array: Vec<u8>) -> BaseSolution {
    let _res = astar(&dist_array, successors, heuristic, success);
    unimplemented!()
}

fn heuristic(state: &Vec<u8>) -> u8 {
    get_n_pop(state) + get_n_loci(state)
}

fn successors(state: &Vec<u8>) -> Vec<(Vec<u8>, u8)> {
    // first check for full joins and return only the first full join if any are available
    if let Some((gx, gy)) = detect_fulljoin(state) {
        let gz = gx.min(gy);
        return vec![
            (state
                .iter()
                .map(|&g| if g == gx || g == gy { gz } else { g })
                .collect_vec(), 2)
        ]
    }

    let n_pop = get_n_pop(state);
    let mut out: Vec<(Vec<u8>, u8)> = vec![];

    for gx in 0..n_pop - 1 {
        for gy in gx + 1..n_pop {
            // Record where the segment joins are and which direction they're going
            let segjoins = state
                .iter()
                .tuple_windows()
                .enumerate()
                .filter(|(_p, (&i, &j))| (i, j) == (gx, gy) || (i, j) == (gy, gx))
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

fn successors_given_chosen_segjoins(state: &[u8], (gx, gy): (usize, usize), chosen_segjoins: &Vec<(usize, CrossoverDirection)>) -> Vec<Vec<usize>> {
    unimplemented!()
}

fn success(state: &Vec<u8>) -> bool {
    *state == vec![0]
}

fn get_n_pop(state: &Vec<u8>) -> u8 {
    state.iter().max().map(|res| *res + 1).unwrap_or(0)
}

fn get_n_loci(state: &Vec<u8>) -> u8 {
    state.len() as u8
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
enum CrossoverDirection {
    Forward,
    Backward,
}

type RemainingRange = Option<(usize, usize)>;

fn successors_given_chosen_gametes(
    state: &[u8],
    gx: u8,
    gy: u8,
    chosen_segjoins: &Vec<(usize, CrossoverDirection)>,
) -> Vec<(Vec<u8>, u8)> {
    // Test if sorted
    assert_eq!(*chosen_segjoins, {
        let mut temp: Vec<(usize, CrossoverDirection)> = chosen_segjoins.clone();
        temp.sort();
        temp
    });

    let n_loci = state.len();

    let is_gxgy = (0..n_loci)
        .map(|i| i == gx as usize || i == gy as usize)
        .collect_vec();

    let is_on_chosen_segjoins = (0..n_loci)
        .map(|i| is_gxgy[i] && (chosen_segjoins.iter().any(|x| x.0 == i || x.0 + 1 == i)))
        .collect_vec();

    let assignees = (0..n_loci)
        .filter(|i| is_gxgy[*i] && !is_on_chosen_segjoins[*i])
        .collect_vec();

    let mut assignments = todo!();

    vec![]
}

fn remaining_ranges(
    n_loci: usize,
    chosen_segjoins: &Vec<(usize, CrossoverDirection)>,
) -> (RemainingRange, RemainingRange) {
    chosen_segjoins
        .iter()
        .map(|(p, d)| match d {
            CrossoverDirection::Forward => (Some((*p + 1, n_loci - 1)), Some((0, *p))),
            CrossoverDirection::Backward => (Some((0, *p)), Some((*p + 1, n_loci - 1))),
        })
        .reduce(|x, y| (join(x.0, y.0), join(x.1, y.1)))
        .unwrap_or((Some((0, n_loci - 1)), Some((0, n_loci - 1))))
}

#[inline]
fn join(x: RemainingRange, y: RemainingRange) -> RemainingRange {
    match (x, y) {
        (Some((xs, xe)), Some((ys, ye))) => {
            let zs = xs.max(ys);
            let ze = xe.min(ye);
            if zs > ze {
                None
            } else {
                Some((zs, ze))
            }
        }
        _ => None,
    }
}

fn detect_fulljoin(state: &[u8]) -> Option<(u8, u8)> {
    let n = state.len();
    let mut is_last = vec![false; n];

    let mut s = HashSet::new();
    let mut e = HashSet::new();

    for i in (0..n).rev() {
        let x = state[i];
        is_last[i] = !e.contains(&x);
        e.insert(x);
    }

    s.insert(state[0]);
    for i in 1..n {
        let x = state[i];
        if !s.contains(&x) && is_last[i-1] {
            return Some((state[i-1], state[i]));
        }
        s.insert(x);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn join_test() {
        assert_eq!(None, join(None, None));
        assert_eq!(None, join(None, Some((2, 6))));
        assert_eq!(None, join(Some((1, 7)), None));
        assert_eq!(None, join(Some((1, 3)), Some((4, 6))));
        assert_eq!(None, join(Some((4, 6)), Some((1, 3))));
        assert_eq!(Some((3, 3)), join(Some((1, 3)), Some((3, 6))));
        assert_eq!(Some((3, 5)), join(Some((1, 5)), Some((3, 6))));
        assert_eq!(Some((3, 7)), join(Some((3, 8)), Some((1, 7))));
    }

    #[test]
    fn remaining_ranges_test() {
        let n_loci = 10;
        let mut v = vec![];
        assert_eq!((Some((0, 9)), Some((0, 9))), remaining_ranges(n_loci, &v));
        v.push((4, CrossoverDirection::Forward));
        assert_eq!((Some((5, 9)), Some((0, 4))), remaining_ranges(n_loci, &v));
        v.push((3, CrossoverDirection::Forward));
        assert_eq!((Some((5, 9)), Some((0, 3))), remaining_ranges(n_loci, &v));
        v.push((7, CrossoverDirection::Backward));
        assert_eq!((Some((5, 7)), None), remaining_ranges(n_loci, &v));
        v.push((8, CrossoverDirection::Forward));
        assert_eq!((None, None), remaining_ranges(n_loci, &v));
    }
}
