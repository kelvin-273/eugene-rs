use crate::plants::bit_array::*;

pub fn min_crossings(n_loci: usize, pop_0: Vec<SingleChromGenotype>) -> Option<usize> {
    if pop_0.contains(&SingleChromGenotype::ideotype(n_loci)) {
        return Some(0);
    }
    unimplemented!();
}

pub fn min_crossings_distribtue(dist_array: &Vec<usize>) -> Option<usize> {
    // TODO: assert valid dist_array <19-09-23> //
    let n_loci = dist_array.len();
    let _n_pop = dist_array.iter().max().unwrap() + 1;
    let n_tree_cells = 2 * n_loci;
    let mut _tree_left = vec![0; n_tree_cells];
    let mut _tree_right = vec![0; n_tree_cells];
    let mut tree_xs = vec![vec![vec![0; n_loci]; 2]; n_tree_cells];
    // create each of the initial population
    for (i, x) in dist_array.iter().enumerate() {
        tree_xs[*x][0][i] = 1;
        tree_xs[*x][1][i] = 1;
    }
    None
}

pub fn enumerate_covering_subsets<S>(_subsets: Vec<S>) {
    unimplemented!();
}

/// Given a distribute array, returns the ranges, as tuples (s, e), where for any range (s, e) and
/// any i in [s, e] dist_array[i] cannot be found in dist_array outside of [s, e].
///
/// Example:
/// ```
/// use eugene::solvers::base_min_crossings_backtracking::subproblems_shallow;
///
/// assert_eq!(
///     subproblems_shallow(&vec![0, 1, 0, 1, 0, 1, 0, 2, 3, 2, 4, 2, 5]),
///     vec![(0, 6), (7, 11), (12, 12)]
/// );
/// ```
pub fn subproblems_shallow(dist_array: &Vec<usize>) -> Vec<(usize, usize)> {
    use std::collections::HashMap;
    let n_pop = dist_array.iter().max().unwrap_or(&0) + 1;
    // collect start and end points of each genotype
    let mut s = HashMap::new();
    let mut e = HashMap::new();
    for (i, x) in dist_array.iter().enumerate() {
        if !s.contains_key(x) {
            s.insert(*x, i);
        }
        e.insert(*x, i);
    }
    // union all intersecting ranges
    let mut ranges: Vec<(usize, usize)> = (0..n_pop)
        .map(|i| -> (usize, usize) {
            (
                *s.get(&i).expect("i not in (0..n_pop)"),
                *e.get(&i).expect("i not in (0..n_pop)"),
            )
        })
        .collect();
    ranges.sort();
    let mut chosen: Vec<bool> = (0..n_pop).map(|_| true).collect();
    for i in 0..n_pop - 1 {
        if ranges[i].1 > ranges[i + 1].0 {
            ranges[i + 1] = (ranges[i].0, ranges[i].1.max(ranges[i + 1].1));
            chosen[i] = false;
        }
    }
    ranges
        .iter()
        .enumerate()
        .filter_map(|(i, x)| match chosen[i] {
            true => Some((x.0, x.1)),
            false => None,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn subproblems_shallow_test() {
        assert_eq!(
            subproblems_shallow(&vec![1, 2, 1, 2, 1, 3, 4, 3, 5, 0]),
            vec![(0, 4), (5, 7), (8, 8), (9, 9)]
        );
        assert_eq!(
            subproblems_shallow(&vec![0, 1, 0, 1, 0, 1, 0, 2, 3, 2, 4, 2]),
            vec![(0, 6), (7, 11)]
        );
    }
}
