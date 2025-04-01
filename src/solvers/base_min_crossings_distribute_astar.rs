use crate::solution::{BaseSolution, PyBaseSolution};
use crate::solvers::base_min_generations_enumerator_dominance::filter_non_dominating_fn;
use core::cmp::Reverse;
use itertools::Itertools;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fs;
use std::io::Write;
use std::rc::Rc;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

/// Computes an optimal crossing schedule given a distribute array `xs`.
///
/// For now, only the objective is computed.
#[pyo3::pyfunction]
pub fn breeding_program_distribute_python(xs: DistArray, timeout: Option<u64>) -> PyBaseSolution {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = breeding_program_distribute(&xs);
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
pub fn breeding_program_distribute_diving_python(xs: DistArray, timeout: Option<u64>) -> PyBaseSolution {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = breeding_program_distribute_diving(&xs);
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

type DistArray = Vec<usize>;

#[inline]
fn distribute_n_pop(xs: &DistArray) -> usize {
    *xs.iter().max().expect("xs is empty") + 1
}

#[derive(Debug, PartialEq, Eq)]
struct AstarNodeBase {
    xs: DistArray,
    parent_node: Option<Rc<AstarNodeBase>>,
    parent_gametes: Option<(usize, usize)>,
    g: usize,
    n_pop: usize,
    n_segments: usize,
}

#[derive(Debug, PartialEq, Eq, Clone)]
struct AstarNode {
    head: Rc<AstarNodeBase>,
}

impl AstarNode {
    pub fn new(xs: &DistArray) -> Self {
        Self {
            head: Rc::new(AstarNodeBase::new(xs)),
        }
    }

    #[inline]
    pub fn dist_array(&self) -> &DistArray {
        &self.head.xs
    }

    pub fn success(&self) -> bool {
        let xs = self.dist_array();
        !xs.is_empty() && xs.iter().all(|x| *x == xs[0])
    }

    #[inline]
    fn parent_node(&self) -> Option<Self> {
        self.head.parent_node.as_ref().map(|node_base| AstarNode {
            head: node_base.clone(),
        })
    }

    #[inline]
    fn parent_gametes(&self) -> Option<(usize, usize)> {
        self.head.parent_gametes
    }

    #[inline]
    fn g(&self) -> usize {
        self.head.g
    }

    #[inline]
    fn h(&self) -> usize {
        self.n_segments() + self.n_pop()
    }

    #[inline]
    fn f(&self) -> f32 {
        self.g() as f32 + self.h() as f32 * 1.1
    }

    #[inline]
    fn n_segments(&self) -> usize {
        self.head.n_segments
    }

    #[inline]
    fn n_pop(&self) -> usize {
        self.head.n_pop
    }

    fn create_offspring(&self, zs: DistArray, gx: usize, gy: usize) -> Self {
        let n_segments = zs.len();
        let n_pop = distribute_n_pop(&zs);
        Self {
            head: Rc::new(AstarNodeBase {
                xs: zs.clone(),
                parent_node: Some(self.head.clone()),
                parent_gametes: Some((gx, gy)),
                g: self.g() + 2,
                n_segments,
                n_pop,
            }),
        }
    }
}

impl AstarNodeBase {
    pub fn new(xs: &DistArray) -> Self {
        Self {
            xs: xs.clone(),
            parent_node: None,
            parent_gametes: None,
            g: 0,
            n_pop: distribute_n_pop(xs),
            n_segments: xs.len(),
        }
    }
}

impl PartialOrd for AstarNode {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        (self.head.g + self.head.n_pop + self.head.n_segments)
            .partial_cmp(&(other.head.g + other.head.n_pop + other.head.n_segments))
    }
}

impl Ord for AstarNode {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.f()
            .partial_cmp(&(other.f()))
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

pub fn breeding_program_distribute(xs: &DistArray) -> Option<BaseSolution> {
    let path = astar(xs)?;
    let mut node_ref = path.clone();
    // TODO: convert path into BaseSolution <10-03-25> //
    let mut obj = 0;
    // Add the final selfing to the total
    if node_ref.parent_node().is_some() {
        obj += 1;
    }
    while node_ref.parent_node().is_some() {
        node_ref = node_ref.parent_node().unwrap();
        obj += 1;
    }
    Some(BaseSolution {
        tree_data: vec![],
        tree_type: vec![],
        tree_left: vec![],
        tree_right: vec![],
        objective: obj,
    })
}

pub fn breeding_program_distribute_diving(xs: &DistArray) -> Option<BaseSolution> {
    let path = diving(xs)?;
    let mut node_ref = path.clone();
    // TODO: convert path into BaseSolution <10-03-25> //
    let mut obj = 0;
    // Add the final selfing to the total
    if node_ref.parent_node().is_some() {
        obj += 1;
    }
    while node_ref.parent_node().is_some() {
        node_ref = node_ref.parent_node().unwrap();
        obj += 1;
    }
    Some(BaseSolution {
        tree_data: vec![],
        tree_type: vec![],
        tree_left: vec![],
        tree_right: vec![],
        objective: obj,
    })
}

struct ClosedList {
    //set: HashMap<String, AstarNode>,
    set: HashSet<String>,
}

impl ClosedList {
    pub fn new() -> Self {
        Self {
            //set: HashMap::new(),
            set: HashSet::new(),
        }
    }

    pub fn contains(&self, node: &AstarNode) -> bool {
        //self.set.contains_key(&format!("{:?}", node.dist_array()))
        self.set.contains(&format!("{:?}", node.dist_array()))
    }

    pub fn insert(&mut self, node: &AstarNode) -> bool {
        //self.set
        //    .insert(format!("{:?}", node.dist_array()), node.clone())
        self.set.insert(format!("{:?}", node.dist_array()))
    }
}

fn astar(xs: &DistArray) -> Option<AstarNode> {
    if xs.is_empty() {
        return None;
    }
    // open file to write posthoc-trace
    let mut file = fs::File::create(format!(
        "posthoc-trace/posthoc-{}.trace.yml",
        xs.iter().map(|x| x.to_string()).join("")
    ))
    .ok()
    .expect("couldn't create file");
    file.write(b"version: 1.4.0\n");
    //dbg!(xs);
    let mut open_list = BinaryHeap::from([Reverse(AstarNode::new(xs))]);
    let mut closed_list = ClosedList::new();
    //dbg!(&open_list);

    while let Some(Reverse(node)) = open_list.pop() {
        //dbg!(&node);
        if closed_list.contains(&node) {
            continue;
        }

        // write node
        write_node(&mut file, &node).expect("failed to write node");

        closed_list.insert(&node);
        if node.success() {
            return Some(node);
        }
        let children = branching(&node, None);
        for child in children {
            write_expansion(&mut file, &child).expect("failed to write successor");
            if !closed_list.contains(&child) {
                open_list.push(Reverse(child))
            }
        }
    }
    None
}

fn diving(xs: &DistArray) -> Option<AstarNode> {
    if xs.is_empty() {
        return None;
    }
    // TODO: Open file to write posthoc-trace <18-03-25> //
    let mut node = AstarNode::new(xs);
    while !node.success() {
        node = branching(&node, None).iter().min()?.clone();
    }
    Some(node)
}

fn write_node(file: &mut fs::File, node: &AstarNode) -> std::io::Result<usize> {
    file.write(
        format!(
            "- {{ type: node, id: {:?}, pId: {} }}\n",
            node.dist_array(),
            node.parent_node()
                .map_or("null".to_string(), |xs| format!("{:?}", xs.dist_array()))
        )
        .as_bytes(),
    )
}

fn write_expansion(file: &mut fs::File, child: &AstarNode) -> std::io::Result<usize> {
    let xs = child.dist_array();
    let zs = child
        .parent_node()
        .map_or(String::from("null"), |x| format!("{:?}", x.dist_array()));
    let (gx, gy) = child
        .parent_gametes()
        .expect("child not equipped with parent gametes");
    file.write(
        format!(
            "- {{ type: successor, id: {}, pId: {:?}, gx: {}, gy: {} }}\n",
            zs, xs, gx, gy
        )
        .as_bytes(),
    )
}

fn branching(node: &AstarNode, closed_list: Option<&ClosedList>) -> Vec<AstarNode> {
    if let Some((zs, gx, gy)) = first_full_join(node.dist_array()) {
        let zs = simplify_dist_array(&zs);
        return vec![node.create_offspring(zs, gx, gy)];
    }

    let solution_iter = generate_redistributions(node.dist_array())
        .into_iter()
        .filter(|(zs, _, _)| zs != node.dist_array());
    //let solution_iter = filter_non_dominating_fn(solution_iter, |(zs1, _, _), (zs2, _, _)| {
    //    dominates_gametewise(zs1, zs2)
    //});
    let solution_iter = solution_iter
        .into_iter()
        .map(|(zs, gx, gy)| (simplify_dist_array(&zs), gx, gy));
    //let solution_iter = filter_non_dominating_fn(solution_iter, |(zs1, _, _), (zs2, _, _)| {
    //    dominates_gametewise(zs1, zs2) || dominates_as_subsequence(zs1, zs2)
    //});
    let solution_iter = solution_iter
        .into_iter()
        .map(|(zs, gx, gy)| node.create_offspring(zs, gx, gy))
        .collect();
    solution_iter
}

enum Direction {
    Down,
    Up,
}

fn generate_redistributions_two_delta(xs: &DistArray) -> Vec<DistArray> {
    let n_loci = xs.len();
    let n_pop = distribute_n_pop(xs);
    assert!(n_loci > 1);

    #[inline]
    fn f(xs: &Vec<usize>, i: usize, gx: usize, gy: usize) -> bool {
        (xs[i], xs[i + 1]) == (gx, gy) || (xs[i], xs[i + 1]) == (gy, gx)
    }

    struct TwoDeltaState<'a> {
        xs: &'a DistArray,
        gx: usize,
        gy: usize,
        available_gz_values: &'a mut Vec<usize>,
        zs: &'a mut DistArray,
        out: &'a mut Vec<DistArray>,
        segjoin_choice: &'a mut Vec<bool>,
    }

    impl TwoDeltaState<'_> {
        #[inline]
        fn is_segjoin_index(&self, i: usize) -> bool {
            (self.xs[i], self.xs[i + 1]) == (self.gx, self.gy)
                || (self.xs[i], self.xs[i + 1]) == (self.gy, self.gx)
        }

        #[inline]
        fn is_gx_or_gy(&self, i: usize) -> bool {
            self.xs[i] == self.gx || self.xs[i] == self.gy
        }
    }

    fn bt_segjoin_choices(state: &mut TwoDeltaState, i: usize) {
        if i >= state.xs.len() - 1 {
            todo!("return redistribution")
        } else if !state.is_segjoin_index(i) {
            bt_segjoin_choices(state, i + 1);
        } else if i < state.xs.len() - 2 && state.is_segjoin_index(i + 1) {
            state.segjoin_choice[i] = true;
            bt_segjoin_choices(state, i + 2);
            state.segjoin_choice[i] = false;
            state.segjoin_choice[i + 1] = true;
            bt_segjoin_choices(state, i + 3);
            state.segjoin_choice[i + 1] = false;
        } else if i < state.xs.len() - 2 && !state.is_segjoin_index(i + 1) {
            state.segjoin_choice[i] = true;
            bt_segjoin_choices(state, i + 2);
            state.segjoin_choice[i] = false;
        } else {
            // i == xs.len()
            state.segjoin_choice[i] = true;
            bt_segjoin_choices(state, i + 1);
            state.segjoin_choice[i] = false;
        }
    }

    fn assign_given_segjoins(state: &mut TwoDeltaState) {
        let n_loci = state.xs.len();
        let ranges = remaining_ranges(state.xs, state.gx, state.gy, state.segjoin_choice);
        if (0..n_loci).any(|i| {
            state.is_gx_or_gy(i)
                && !state.segjoin_choice[i]
                && !state.segjoin_choice[i + 1]
                && ![(ranges[0][0]..ranges[0][1]), (ranges[1][0]..ranges[1][1])]
                    [(state.xs[i] > state.gx) as usize]
                    .contains(&i)
        }) {
            return;
        }
        fn bt_final_gamete_construction(
            state: &mut TwoDeltaState,
            i: usize,
        ) {
            let mut j = 0;
        }
    }

    let mut out = Vec::new();
    let mut segjoin_choice = vec![false; n_loci - 1];
    let mut zs = xs.clone();
    let mut available_gz_values: Vec<usize> = (n_pop - 2..n_loci).collect();
    for gx in 0..n_pop - 1 {
        for gy in gx + 1..n_pop {
            let mut state = TwoDeltaState {
                xs,
                gx,
                gy,
                available_gz_values: &mut available_gz_values,
                out: &mut out,
                segjoin_choice: &mut segjoin_choice,
                zs: &mut zs,
            };

            if (0..n_loci - 1).all(|i| !f(xs, i, gx, gy)) {
                continue;
            }
            state.available_gz_values[0] = gx;
            state.available_gz_values[1] = gy;
            bt_segjoin_choices(&mut state, 0);
        }
    }
    out
}

/// Returns ranges of loci on each input gamete that are not covered by any chosen segment join
/// Ranges are not inclusive of the endpoints.
///
/// ```
/// assert_eq!(remaining_ranges(
///     &vec![0, 1, 2, 0, 1, 3], 0, 1,
///     &vec![true, false, false, true, false]
/// ), [[4, 6], [0, 1]]);
///
/// assert_eq!(remaining_ranges(
///     &vec![0, 1, 2, 1, 0, 3], 0, 1,
///     &vec![true, false, false, true, false]
/// ), [[1, 4], [4, 1]]);
/// ```
fn remaining_ranges(
    xs: &Vec<usize>,
    gx: usize,
    gy: usize,
    segjoin_choice: &Vec<bool>,
) -> [[usize; 2]; 2] {
    let mut ranges_remaining = [[0, xs.len()], [0, xs.len()]];
    for i in 0..xs.len() - 1 {
        if segjoin_choice[i] && ((xs[i], xs[i + 1]) == (gx, gy)) {
            ranges_remaining[0][1] = ranges_remaining[0][1].min(i + 1);
            ranges_remaining[1][0] = ranges_remaining[1][0].max(i + 1);
        } else if segjoin_choice[i] && ((xs[i], xs[i + 1]) == (gy, gx)) {
            ranges_remaining[0][0] = ranges_remaining[0][0].min(i + 1);
            ranges_remaining[1][1] = ranges_remaining[1][1].max(i + 1);
        }
    }
    ranges_remaining
}

fn first_full_join(xs: &DistArray) -> Option<(DistArray, usize, usize)> {
    let n_loci = xs.len();
    let n_pop = distribute_n_pop(xs);

    let mut start_points = vec![n_loci; n_pop];
    let mut end_points = vec![0; n_pop];
    for (i, x) in xs.iter().enumerate() {
        end_points[*x] = i;
    }
    for (i, x) in xs.iter().enumerate().rev() {
        start_points[*x] = i;
    }

    // find the (gx,gy) with the smallest start_points[gx] such that
    // a full join occurs between gx and gy
    for (gy, sy) in start_points.iter().enumerate() {
        if *sy > 0 && end_points[xs[sy - 1]] == sy - 1 {
            let gx = xs[sy - 1];
            let mut out = xs.clone();
            for i in *sy..end_points[gy] + 1 {
                if xs[i] == gy {
                    out[i] = gx;
                }
            }
            return Some((out, gx, gy));
        }
    }
    None
}

fn simplify_dist_array(xs: &DistArray) -> DistArray {
    let mut zs = xs.clone();
    let n_loci = zs.len();
    let n_pop = distribute_n_pop(&zs);
    let mut mapping = vec![None; n_pop];
    mapping[zs[0]] = Some(0);
    let mut x_max = 0;
    let mut i = 1;
    let mut j = 1;
    while j < n_loci {
        if zs[j] != xs[j - 1] {
            if mapping[zs[j]].is_some() {
                zs[i] = mapping[zs[j]]
                    .expect("mapping[xs[j]] contains None despite .is_some() reporting otherwise");
            } else {
                x_max += 1;
                mapping[zs[j]] = Some(x_max);
                zs[i] = x_max;
            }
            i += 1;
        }
        j += 1;
    }
    zs.resize(i, 0);
    zs
}

fn requires_multipoint(zs: &DistArray, xs: &DistArray, g: usize, i: usize) -> bool {
    let mut swaps = -1;
    let mut prev_value = None;
    for (&gz, &gx) in zs.iter().take(i + 1).zip(xs.iter()) {
        if gz == g && Some(gx) != prev_value {
            prev_value = Some(gx);
            swaps += 1;
        }
        if swaps > 1 {
            return true;
        }
    }
    false
}

/// Creates redistributions of Distribute array `xs`.
/// Currently does not prune dominated redistributions.
fn generate_redistributions(xs: &DistArray) -> Vec<(DistArray, usize, usize)> {
    let mut zs = xs.clone();
    let n_loci = xs.len();
    let n_pop = distribute_n_pop(xs);
    if n_pop == 1 {
        return vec![];
    }
    let mut available_gz_values = (n_pop - 2..n_loci).collect::<Vec<_>>();
    let mut out = Vec::new();
    // Backtracking for raw redistributions
    // Assumes that gx and gy are already fixed
    fn bt(
        xs: &DistArray,
        out: &mut Vec<(DistArray, usize, usize)>,
        gx: usize,
        gy: usize,
        zs: &mut DistArray,
        i: usize,
        available_gz_values: &Vec<usize>,
        j_max: usize,
    ) {
        if i >= xs.len() {
            out.push((zs.clone(), gx, gy));
        } else if xs[i] != gx && xs[i] != gy {
            bt(xs, out, gx, gy, zs, i + 1, available_gz_values, j_max);
        } else {
            for (j, &gz) in available_gz_values.iter().take(j_max + 1).enumerate() {
                // set value
                zs[i] = gz;
                // check if set value would require more than 1 crossover
                if requires_multipoint(zs, xs, gz, i) {
                    continue;
                }
                // TODO: check if set value would create a dominated gamete
                // Backtrack
                if j == j_max {
                    bt(xs, out, gx, gy, zs, i + 1, available_gz_values, j_max + 1);
                } else {
                    bt(xs, out, gx, gy, zs, i + 1, available_gz_values, j_max);
                }
            }
        }
    }
    for gx in 0..n_pop - 1 {
        for gy in gx + 1..n_pop {
            for i in 0..n_loci {
                zs[i] = xs[i];
            }
            available_gz_values[0] = gx;
            available_gz_values[1] = gy;
            bt(xs, &mut out, gx, gy, &mut zs, 0, &available_gz_values, 0);
        }
    }
    out
}

/// Returns true if each gamete in `ys` are a subset of some gamete in `xs`
fn dominates_gametewise(xs: &DistArray, ys: &DistArray) -> bool {
    if xs.len() != ys.len() {
        return false;
    }
    let n_pop_ys = distribute_n_pop(ys);
    let mut d = vec![None; n_pop_ys];
    for (&gx, &gy) in xs.iter().zip(ys.iter()) {
        if d[gy].is_none() {
            d[gy] = Some(gx);
        } else if d[gy] != Some(gx) {
            return false;
        }
    }
    true
}

/// Returns true if and only if `ys` is a subsequence of `zs`
fn dominates_as_subsequence(zs: &DistArray, ys: &DistArray) -> bool {
    return false;
    let nz = zs.len();
    let ny = ys.len();
    if nz > ny {
        return false;
    }
    let mut dp = vec![vec![0; ny + 1]; nz + 1];
    for iz in 0..nz {
        for iy in 0..ny {
            dp[1 + iz][1 + iy] = (((zs[iz] == ys[iy]) as usize) + dp[iz][iy])
                .max(dp[iz + 1][iy])
                .max(dp[iz][iy + 1])
        }
    }
    dp[nz][ny] == nz
}

type PyPath = pyo3::PyResult<Option<(Vec<DistArray>, Vec<Option<(usize, usize)>>)>>;

/// Computes an optimal crossing schedule given a distribute array `xs`.
///
/// For now, only the objective is computed.
#[pyo3::pyfunction]
pub fn experiment_distribute_over_2delta_python(xs: DistArray, timeout: Option<u64>) -> PyPath {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        tx.send(astar(&xs).map(|path| {
            let mut node_ref = path.clone();
            let mut out_distarrays = vec![path.dist_array().clone()];
            let mut out_parents = vec![path.parent_gametes()];
            while node_ref.parent_node().is_some() {
                node_ref = node_ref.parent_node().unwrap();
                out_distarrays.push(node_ref.dist_array().clone());
                out_parents.push(node_ref.parent_gametes())
            }
            out_distarrays.reverse();
            out_parents.reverse();
            (out_distarrays, out_parents)
        }))
    });
    let res = rx
        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
        .ok()
        .flatten();
    Ok(res)
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! pretty_print {
        ($xs: expr) => {
            format!(
                "[\n{}\n]",
                $xs.iter().map(|zs| format!("\t{:?}", zs)).join("\n")
            )
        };
    }

    #[test]
    fn branching_failed_test() {
        macro_rules! f {
            ($xs: expr, $zs: expr) => {
                let node = AstarNode::new(&Vec::from($xs));
                let output: Vec<DistArray> = branching(&node, None)
                    .iter()
                    .map(|x| x.dist_array())
                    .cloned()
                    .collect();
                assert!(
                    output.contains(&Vec::from($zs)),
                    "{:?} not in {}",
                    $zs,
                    pretty_print!(output)
                );
            };
        }
        f!([0, 1, 0], [0, 1]);
        f!([0, 1, 0, 2, 1, 0, 2], [0, 1, 0, 1, 2]);
    }

    #[test]
    fn breeding_program_distribute_test() {
        macro_rules! f {
            ($xs: expr, $obj_check: expr) => {
                assert_eq!(
                    breeding_program_distribute(&Vec::from($xs)).map(|sol| sol.objective),
                    Some($obj_check)
                )
            };
        }
        f!([0], 0);
        f!([0, 1], 2);
        f!([0, 1, 0], 3);
        f!([0, 1, 2], 3);
        f!([0, 1, 0, 1], 3);
        f!([0, 1, 0, 2], 4);
        f!([0, 1, 2, 0], 4);
        f!([0, 1, 2, 1], 4);
        f!([0, 1, 0, 1, 0], 4);
        f!([0, 1, 2, 0, 1], 4);
        f!([0, 1, 2, 1, 0], 4);
        f!([0, 1, 0, 2, 0], 5);
        f!([0, 1, 0, 2, 0, 1, 0], 5);
        f!([0, 1, 0, 2, 1, 0, 2], 5);
    }

    #[test]
    fn simplify_dist_array_test() {
        assert_eq!(
            simplify_dist_array(&vec![0, 1, 0, 0, 1, 2, 2]),
            vec![0, 1, 0, 1, 2]
        )
    }

    #[test]
    fn generate_redistributions_test() {
        macro_rules! f {
            ($xs: expr, $zs: expr) => {
                let output = generate_redistributions(&Vec::from($xs));
                assert!(
                    output.contains(&(Vec::from($zs.0), $zs.1, $zs.2)),
                    "{:?} not in {}",
                    $zs,
                    pretty_print!(output)
                )
            };
        }
        f!([0, 1, 0, 2, 1, 0, 2], ([0, 1, 0, 0, 1, 2, 2], 0, 2));
    }

    #[test]
    fn dominates_gametewise_test() {
        macro_rules! f {
            ($xs: expr, $ys: expr) => {
                assert!(
                    dominates_gametewise(&Vec::from($xs), &Vec::from($ys)),
                    "dominates_gametewise({:?}, {:?}) == false",
                    $xs,
                    $ys
                )
            };
        }
        f!([0, 0, 1], [0, 1, 2]);
        f!([0, 1, 1], [0, 1, 2]);
    }

    #[test]
    fn branching_full_joins_test() {
        unimplemented!()
    }

    #[test]
    fn remaining_ranges_test() {
        assert_eq!(
            remaining_ranges(
                &vec![0, 1, 2, 0, 1, 3],
                0,
                1,
                &vec![true, false, false, true, false]
            ),
            [[4, 6], [0, 1]]
        );

        assert_eq!(
            remaining_ranges(
                &vec![0, 1, 2, 1, 0, 3],
                0,
                1,
                &vec![true, false, false, true, false]
            ),
            [[1, 4], [4, 1]]
        );
    }
}
