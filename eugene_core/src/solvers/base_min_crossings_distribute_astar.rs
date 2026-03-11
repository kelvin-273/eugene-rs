use crate::plants::dist_array::DistArray;
use crate::solution::BaseSolution;
use crate::solvers::base_min_generations_enumerator_dominance::IteratorNonDominating;
use core::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fs;
use std::io::Write;
use std::rc::Rc;

pub struct Config {
    full_join: bool,
    dominance: bool,
    diving: bool,
    debug_trace_file: Option<String>,
}

impl Config {
    pub fn new(
        full_join: bool,
        dominance: bool,
        diving: bool,
        debug_trace_file: Option<String>,
    ) -> Self {
        Self {
            full_join,
            dominance,
            diving,
            debug_trace_file,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct Output {
    expansions: usize,
    pushed_nodes: usize,
    children_created_by_branching: usize,
    objective: Option<usize>,
}

impl Output {
    pub fn expansions(&self) -> usize {
        self.expansions
    }

    pub fn pushed_nodes(&self) -> usize {
        self.pushed_nodes
    }

    pub fn children_created_by_branching(&self) -> usize {
        self.children_created_by_branching
    }

    pub fn objective(&self) -> Option<usize> {
        self.objective
    }
}

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
    fn crossings(&self) -> usize {
        debug_assert_eq!(self.g() & 1, 0);
        self.g() >> 1
    }

    #[inline]
    fn h(&self) -> usize {
        self.n_segments() + self.n_pop()
    }

    #[inline]
    fn f(&self) -> f32 {
        self.g() as f32 + self.h() as f32
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
    #[allow(clippy::non_canonical_partial_ord_impl)]
    #[inline]
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

#[inline]
fn compute_objective(path: &AstarNode) -> usize {
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
    obj
}

pub fn breeding_program_distribute_general(xs: &DistArray, config: &Config) -> Option<Output> {
    let mut output = Output {
        expansions: 0,
        pushed_nodes: 0,
        children_created_by_branching: 0,
        objective: None,
    };
    let path = astar_general(xs, config, &mut output)?;
    let _ = output.objective.insert(compute_objective(&path));
    Some(output)
}

fn astar_general(xs: &DistArray, config: &Config, output: &mut Output) -> Option<AstarNode> {
    if xs.is_empty() {
        return None;
    }
    let mut open_list = BinaryHeap::from([Reverse(AstarNode::new(xs))]);
    let mut closed_list = ClosedList::new();

    let mut file = config
        .debug_trace_file
        .as_ref()
        .map(|s| fs::File::create(s).expect("unable to create debug trace file"));
    if let Some(f) = file.as_mut() {
        let _ = f.write(b"version: 1.4.0\n");
        let _ = f.write(b"events:\n");
    }

    while let Some(Reverse(node)) = open_list.pop() {
        let _ = file.as_mut().map(|f| write_node(f, &node));
        output.expansions += 1;
        if closed_list.contains(&node) {
            continue;
        }

        closed_list.insert(&node);
        if node.success() {
            return Some(node);
        }
        let children = branching_general(&node, config);
        output.children_created_by_branching += children.len();
        for child in children {
            let _ = file.as_mut().map(|f| write_expansion(f, &node));
            if !closed_list.contains(&child) {
                open_list.push(Reverse(child));
                output.pushed_nodes += 1;
            }
        }
    }
    None
}

fn branching_general(node: &AstarNode, config: &Config) -> Vec<AstarNode> {
    //if config.full_join || config.diving {  // Adding config.diving passes the tests,
    // but why do the tests fail without it? Surely the full join should be in the set of solutions
    // returned by generate_redistributions
    if config.full_join {
        if let Some((zs, gx, gy)) = first_full_join(node.dist_array()) {
            let zs = simplify_dist_array(&zs);
            return vec![node.create_offspring(zs, gx, gy)];
        }
    }
    if config.diving {
        // TODO: return either the minimum one or the one that meets the lower bound
        let lower_bound = node.n_segments() + node.n_pop();
        assert!(lower_bound >= 2);

        let mut min_obj = usize::MAX;
        let mut min_dist_arr = (node.dist_array().clone(), 0, 0);

        for (zs, gx, gy) in generate_redistributions(node.dist_array()) {
            if &zs == node.dist_array() {
                continue;
            }
            let zs_simplified = simplify_dist_array(&zs);
            let n_pop = zs_simplified.iter().max().expect("empty dist_array") + 1;
            let n_segments = zs_simplified.len();
            if n_pop + n_segments == lower_bound - 2 {
                return vec![node.create_offspring(zs_simplified, gx, gy)];
            } else if n_pop + n_segments < min_obj {
                min_dist_arr = (zs_simplified, gx, gy);
                min_obj = n_pop + n_segments;
            }
        }
        let (zs, gx, gy) = min_dist_arr;
        vec![node.create_offspring(zs, gx, gy)]
    } else if config.dominance {
        generate_redistributions(node.dist_array())
            .into_iter()
            .filter(|(zs, _, _)| zs != node.dist_array())
            .filter_non_dominating_fn(|(xs, _, _), (ys, _, _)| dominates_gametewise(xs, ys))
            .map(|(zs, gx, gy)| (simplify_dist_array(&zs), gx, gy))
            .filter_non_dominating_fn(|(xs, _, _), (ys, _, _)| {
                dominates_gametewise(xs, ys) || dominates_as_subsequence(xs, ys)
            })
            .map(|(zs, gx, gy)| node.create_offspring(zs, gx, gy))
            .collect()
    } else {
        generate_redistributions(node.dist_array())
            .into_iter()
            .filter(|(zs, _, _)| zs != node.dist_array())
            .map(|(zs, gx, gy)| (simplify_dist_array(&zs), gx, gy))
            .map(|(zs, gx, gy)| node.create_offspring(zs, gx, gy))
            .collect()
    }
}

pub fn breeding_program_distribute(xs: &DistArray) -> Option<BaseSolution> {
    if xs.is_target() {
        return Some(BaseSolution {
            tree_data: vec![vec![vec![1; xs.n_loci()]; 2]],
            tree_type: vec!["Leaf"],
            tree_left: vec![0],
            tree_right: vec![0],
            objective: 0,
        });
    }

    let path = astar(xs)?;
    let obj = path.crossings() + (path.crossings() > 0) as usize;

    let mut sol = BaseSolution {
        tree_data: vec![],
        tree_type: vec![],
        tree_left: vec![],
        tree_right: vec![],
        objective: obj,
    };

    // Create population of gametes (as `Vec<i32>'s`)
    let popg: Vec<Vec<i32>> = (0..xs.n_pop())
        .map(|i| {
            (0..xs.n_loci())
                .map(|j| if xs[j] == i { 1 } else { 0 })
                .collect()
        })
        .collect();
    let _popg_idx = (0..popg.len()).collect::<Vec<_>>();

    for g in popg.iter() {
        sol.tree_data.push(vec![g.clone(), g.clone()]);
        sol.tree_type.push("Leaf");
        sol.tree_left.push(0);
        sol.tree_right.push(0);
    }

    // collect distribute arrays on the path
    let mut node_ref = path.clone();
    let mut path_dist_arrays = vec![node_ref.dist_array().clone()];
    while node_ref.parent_node().is_some() {
        node_ref = node_ref.parent_node().unwrap();
        path_dist_arrays.push(node_ref.dist_array().clone());
    }

    // Working backwards, convert the distribute arrays to unsimplifed form
    //todo!("Oh God!");

    Some(sol)
}

pub fn breeding_program_distribute_dominance(xs: &DistArray) -> Option<BaseSolution> {
    let path = astar_dominance(xs)?;
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

pub fn breeding_program_distribute_no_full_join(xs: &DistArray) -> Option<BaseSolution> {
    let path = astar_no_full_join(xs)?;
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

pub fn breeding_program_distribute_dominance_no_full_join(xs: &DistArray) -> Option<BaseSolution> {
    let path = astar_dominance_no_full_join(xs)?;
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
    let mut open_list = BinaryHeap::from([Reverse(AstarNode::new(xs))]);
    let mut closed_list = ClosedList::new();

    while let Some(Reverse(node)) = open_list.pop() {
        if closed_list.contains(&node) {
            continue;
        }

        closed_list.insert(&node);
        if node.success() {
            return Some(node);
        }
        let children = branching(&node);
        for child in children {
            if !closed_list.contains(&child) {
                open_list.push(Reverse(child))
            }
        }
    }
    None
}

fn astar_no_full_join(xs: &DistArray) -> Option<AstarNode> {
    if xs.is_empty() {
        return None;
    }
    let mut open_list = BinaryHeap::from([Reverse(AstarNode::new(xs))]);
    let mut closed_list = ClosedList::new();

    while let Some(Reverse(node)) = open_list.pop() {
        if closed_list.contains(&node) {
            continue;
        }

        closed_list.insert(&node);
        if node.success() {
            return Some(node);
        }
        let children = branching_no_full_join(&node);
        for child in children {
            if !closed_list.contains(&child) {
                open_list.push(Reverse(child))
            }
        }
    }
    None
}

fn astar_dominance(xs: &DistArray) -> Option<AstarNode> {
    if xs.is_empty() {
        return None;
    }
    let mut open_list = BinaryHeap::from([Reverse(AstarNode::new(xs))]);
    let mut closed_list = ClosedList::new();

    while let Some(Reverse(node)) = open_list.pop() {
        if closed_list.contains(&node) {
            continue;
        }

        closed_list.insert(&node);
        if node.success() {
            return Some(node);
        }
        let children = branching_dominance(&node);
        for child in children {
            if !closed_list.contains(&child) {
                open_list.push(Reverse(child))
            }
        }
    }
    None
}

fn astar_dominance_no_full_join(xs: &DistArray) -> Option<AstarNode> {
    if xs.is_empty() {
        return None;
    }
    let mut open_list = BinaryHeap::from([Reverse(AstarNode::new(xs))]);
    let mut closed_list = ClosedList::new();

    while let Some(Reverse(node)) = open_list.pop() {
        if closed_list.contains(&node) {
            continue;
        }

        closed_list.insert(&node);
        if node.success() {
            return Some(node);
        }
        let children = branching_dominance_no_full_join(&node);
        for child in children {
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
        node = branching(&node).iter().min()?.clone();
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

fn branching_no_full_join(node: &AstarNode) -> Vec<AstarNode> {
    generate_redistributions(node.dist_array())
        .into_iter()
        .filter(|(zs, _, _)| zs != node.dist_array())
        .map(|(zs, gx, gy)| (simplify_dist_array(&zs), gx, gy))
        .map(|(zs, gx, gy)| node.create_offspring(zs, gx, gy))
        .collect()
}

fn branching_dominance_no_full_join(node: &AstarNode) -> Vec<AstarNode> {
    generate_redistributions(node.dist_array())
        .into_iter()
        .filter(|(zs, _, _)| zs != node.dist_array())
        .filter_non_dominating_fn(|(zs1, _, _), (zs2, _, _)| dominates_gametewise(zs1, zs2))
        .map(|(zs, gx, gy)| (simplify_dist_array(&zs), gx, gy))
        .filter_non_dominating_fn(|(zs1, _, _), (zs2, _, _)| {
            dominates_gametewise(zs1, zs2) || dominates_as_subsequence(zs1, zs2)
        })
        .map(|(zs, gx, gy)| node.create_offspring(zs, gx, gy))
        .collect()
}

fn branching(node: &AstarNode) -> Vec<AstarNode> {
    if let Some((zs, gx, gy)) = first_full_join(node.dist_array()) {
        let zs = simplify_dist_array(&zs);
        return vec![node.create_offspring(zs, gx, gy)];
    }
    generate_redistributions(node.dist_array())
        .into_iter()
        .filter(|(zs, _, _)| zs != node.dist_array())
        .map(|(zs, gx, gy)| (simplify_dist_array(&zs), gx, gy))
        .map(|(zs, gx, gy)| node.create_offspring(zs, gx, gy))
        .collect()
}

fn branching_dominance(node: &AstarNode) -> Vec<AstarNode> {
    if let Some((zs, gx, gy)) = first_full_join(node.dist_array()) {
        let zs = simplify_dist_array(&zs);
        return vec![node.create_offspring(zs, gx, gy)];
    }
    generate_redistributions(node.dist_array())
        .into_iter()
        .filter(|(zs, _, _)| zs != node.dist_array())
        .filter_non_dominating_fn(|(zs1, _, _), (zs2, _, _)| dominates_gametewise(zs1, zs2))
        .map(|(zs, gx, gy)| (simplify_dist_array(&zs), gx, gy))
        .filter_non_dominating_fn(|(zs1, _, _), (zs2, _, _)| {
            dominates_gametewise(zs1, zs2) || dominates_as_subsequence(zs1, zs2)
        })
        .map(|(zs, gx, gy)| node.create_offspring(zs, gx, gy))
        .collect()
}

struct _TwoDeltaState<'a> {
    xs: &'a DistArray,
    gx: usize,
    gy: usize,
    available_gz_values: &'a mut Vec<usize>,
    zs: &'a mut DistArray,
    segjoin_choice: &'a mut Vec<bool>,
    out: &'a mut Vec<DistArray>,
}

impl _TwoDeltaState<'_> {
    #[inline]
    fn _is_segjoin_index(&self, i: usize) -> bool {
        (self.xs[i], self.xs[i + 1]) == (self.gx, self.gy)
            || (self.xs[i], self.xs[i + 1]) == (self.gy, self.gx)
    }

    #[inline]
    fn _is_gx_or_gy(&self, i: usize) -> bool {
        self.xs[i] == self.gx || self.xs[i] == self.gy
    }
}

fn _generate_redistributions_two_delta(xs: &DistArray) -> Vec<DistArray> {
    let n_loci = xs.len();
    let n_pop = distribute_n_pop(xs);
    assert!(n_loci > 1);

    #[inline]
    fn f(xs: &[usize], i: usize, gx: usize, gy: usize) -> bool {
        (xs[i], xs[i + 1]) == (gx, gy) || (xs[i], xs[i + 1]) == (gy, gx)
    }

    fn bt_segjoin_choices(state: &mut _TwoDeltaState, i: usize) {
        if i >= state.xs.len() - 1 {
            todo!("return redistribution")
        } else if !state._is_segjoin_index(i) {
            bt_segjoin_choices(state, i + 1);
        } else if i < state.xs.len() - 2 && state._is_segjoin_index(i + 1) {
            state.segjoin_choice[i] = true;
            bt_segjoin_choices(state, i + 2);
            state.segjoin_choice[i] = false;
            state.segjoin_choice[i + 1] = true;
            bt_segjoin_choices(state, i + 3);
            state.segjoin_choice[i + 1] = false;
        } else if i < state.xs.len() - 2 && !state._is_segjoin_index(i + 1) {
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

    fn assign_given_segjoins(state: &mut _TwoDeltaState) {
        let n_loci = state.xs.len();
        let ranges = _remaining_ranges(state.xs, state.gx, state.gy, state.segjoin_choice);
        if (0..n_loci).any(|i| {
            state._is_gx_or_gy(i)
                && !state.segjoin_choice[i]
                && !state.segjoin_choice[i + 1]  // i.e. we didn't choose the next segjoin given
                                                 // that we didn't choose the current segjoin
                && ![(ranges[0][0]..ranges[0][1]), (ranges[1][0]..ranges[1][1])]
                    [(state.xs[i] > state.gx) as usize]
                    .contains(&i)
        }) {
            return;
        }
        enum Direction {
            Down,
            Up,
        }

        /// We've chosen which segment joins we're going to use.
        /// We know that
        fn bt_final_gamete_construction(state: &mut _TwoDeltaState, _i: usize) {
            let [_x_rng, _y_rng] =
                _remaining_ranges(state.xs, state.gx, state.gy, state.segjoin_choice);
            todo!("Complete bt_final_gamete_construction");
        }
    }

    let mut out = Vec::new();
    let mut segjoin_choice = vec![false; n_loci - 1];
    let mut zs = xs.clone();
    let mut available_gz_values: Vec<usize> = (n_pop - 2..n_loci).collect();
    for gx in 0..n_pop - 1 {
        for gy in gx + 1..n_pop {
            let mut state = _TwoDeltaState {
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
fn _remaining_ranges(
    xs: &[usize],
    gx: usize,
    gy: usize,
    segjoin_choice: &[bool],
) -> [[usize; 2]; 2] {
    let mut ranges_remaining = [[0, xs.len()], [0, xs.len()]];
    for i in 0..xs.len() - 1 {
        if segjoin_choice[i] && ((xs[i], xs[i + 1]) == (gx, gy)) {
            ranges_remaining[0][0] = ranges_remaining[0][0].max(i + 1);
            ranges_remaining[1][1] = ranges_remaining[1][1].min(i + 1);
        } else if segjoin_choice[i] && ((xs[i], xs[i + 1]) == (gy, gx)) {
            ranges_remaining[1][0] = ranges_remaining[1][0].max(i + 1);
            ranges_remaining[0][1] = ranges_remaining[0][1].min(i + 1);
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

fn simplify_dist_array_new(xs: &DistArray) -> DistArray {
    let mut zs = xs.clone();
    let n_loci = xs.n_loci();
    let n_pop = xs.n_pop();
    let mut mapping = vec![None; n_pop];
    mapping[zs[0]] = Some(0);
    let mut x_max = 0;
    let mut i = 1;
    let mut j = 1;
    while j < n_loci {
        if xs[j] != xs[j - 1] {
            if mapping[zs[j]].is_none() {
                mapping[zs[j]] = Some(x_max);
                x_max += 1;
            }
            zs[i] = mapping[xs[j]]
                // TODO: is there a way to rewrite this so that we don't have to call expect?
                .expect("mapping[xs[j]] contains None despite .is_some() reporting otherwise");
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
        available_gz_values: &[usize],
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
            zs[..n_loci].copy_from_slice(&xs[..n_loci]);
            available_gz_values[0] = gx;
            available_gz_values[1] = gy;
            bt(xs, &mut out, gx, gy, &mut zs, 0, &available_gz_values, 0);
        }
    }
    if out.is_empty() {
        dbg!("No redistributions found for {:?}", xs);
    }
    out
}

/// Generates redistributions of Distribute array `xs` using brute force.
/// This is a very inefficient method and should only be used for testing purposes.
///
/// Returns a vector of tuples containing the redistribution, gx, and gy.
/// The gx and gy are the gametes that were used to create the redistribution.
///
/// # Examples
/// ```
/// let xs = vec![0, 1, 0, 2];
/// let redistributions = generate_redistributions_brute_force(&xs);
///
/// // TODO: Add assertions to check the correctness of the redistributions
/// ```
#[allow(dead_code)]
fn generate_redistributions_brute_force(xs: &DistArray) -> Vec<(DistArray, usize, usize)> {
    let mut zs = xs.clone();
    let n_pop = distribute_n_pop(xs);
    if n_pop == 1 {
        return vec![];
    }
    let mut out = Vec::new();
    // Backtracking for raw redistributions
    fn bt(xs: &DistArray, out: &mut Vec<(DistArray, usize, usize)>, zs: &mut DistArray, i: usize) {
        if i >= xs.len() {
            if let Some((gx, gy)) = find_parent_gametes(xs, zs) {
                out.push((zs.clone(), gx, gy));
            }
        } else {
            for gz in 0..i + 1 {
                zs[i] = gz;
                bt(xs, out, zs, i + 1);
            }
        }
    }
    bt(xs, &mut out, &mut zs, 0);
    out
}

/// Finds the parent gametes of a given redistribution `zs` of `xs`.
fn find_parent_gametes(xs: &DistArray, zs: &DistArray) -> Option<(usize, usize)> {
    let n_pop = distribute_n_pop(xs);

    #[inline]
    fn is_single_point(xs: &DistArray, zs: &DistArray, gx: usize, gy: usize) -> bool {
        xs.iter().zip(zs.iter()).all(|(&x, &z)| {
            if x == gx || x == gy {
                // Check that all of the indices j where zs[j] == z switch between gx and gy at
                // most once
                let mut current_value = None;
                let mut switched = false;
                for (&gz, &gx) in zs.iter().zip(xs.iter()) {
                    if gz != z || Some(gx) == current_value {
                        continue;
                    } else if current_value.is_none() {
                        current_value = Some(gx);
                    } else if current_value != Some(gx) {
                        if switched {
                            return false; // More than one switch
                        }
                        switched = true;
                        current_value = Some(gx);
                    }
                }
                true
            } else {
                z == x
            }
        })
    }

    #[inline]
    fn is_semi_canonical(xs: &DistArray, zs: &DistArray, gx: usize, gy: usize) -> bool {
        let n_pop = distribute_n_pop(xs);
        assert!(gx < gy);
        let mut gz_max = gx;
        for (&_, &z) in xs
            .iter()
            .zip(zs.iter())
            .filter(|(&x, &_)| x == gx || x == gy)
        {
            if z != gx && z != gy && z < n_pop {
                return false;
            }
            if z > gz_max {
                if gz_max == gx && z != gy {
                    return false;
                }
                if gz_max == gy && z != n_pop {
                    return false;
                }
                if gz_max >= n_pop && z != gz_max + 1 {
                    return false;
                }
                gz_max = z;
            }
        }
        true
    }

    for gx in 0..n_pop {
        for gy in gx + 1..n_pop {
            if is_single_point(xs, zs, gx, gy) && is_semi_canonical(xs, zs, gx, gy) {
                return Some((gx, gy));
            }
        }
    }
    None
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
fn dominates_as_subsequence(xs: &DistArray, ys: &DistArray) -> bool {
    let nx = xs.len();
    let ny = ys.len();
    if nx > ny {
        return false;
    }
    if nx == ny && xs != ys {
        return false;
    }
    // Constructing DP array where dp[ix][iy] stores the length of the longest common subsequence
    // between xs[..ix] and ys[..iy]
    let mut dp = vec![vec![0; ny + 1]; nx + 1];
    for (ix, x) in xs.iter().enumerate() {
        for (iy, y) in ys.iter().enumerate() {
            dp[ix + 1][iy + 1] = (dp[ix][iy] + ((x == y) as usize))
                .max(dp[ix + 1][iy])
                .max(dp[ix][iy + 1])
        }
    }
    dp[nx][ny] == nx
}

/// Constructs a `BaseSolution` from a list of `DistArray`s corresponding to a path from the
/// starting distribute array to the target distribute array. If the list of distribute arrays does
/// not correspond to a valid path, then `None` is returned.
fn distribute_arrays_to_solution_with_parents(
    xss: &[DistArray],
    gxys: &[(usize, usize)],
) -> Option<BaseSolution> {
    todo!("Implement distribute_arrays_to_solution_with_parents");
}

fn gametes_from_redistribution_forwards(
    xs: &DistArray,
    zs: &DistArray,
    gxs: &Vec<Vec<i32>>,
) -> Vec<Vec<i32>> {
    todo!("Implement gametes_from_redistribution_forwards");
}

fn gametes_from_redistribution_forwards_with_parents(
    xs: &DistArray,
    zs: &DistArray,
    gzs: &Vec<Vec<i32>>,
    gx: usize,
    gy: usize,
) -> Vec<Vec<i32>> {
    todo!("Implement gametes_from_redistribution_backwards");
}

fn gmaetes_used_in_redistribution(xs: &DistArray, zs: &DistArray) -> (usize, usize) {
    let mut d_forwards = HashMap::new();
    let mut d_backwards = HashMap::new();
    for (gx, gz) in xs.iter().zip(zs.iter()) {
        d_forwards
            .entry(*gx)
            .and_modify(|s: &mut HashSet<usize>| {
                s.insert(*gz);
            })
            .or_insert_with(|| {
                let mut s = HashSet::new();
                s.insert(*gz);
                s
            });
        d_backwards
            .entry(*gz)
            .and_modify(|s: &mut HashSet<usize>| {
                s.insert(*gx);
            })
            .or_insert_with(|| {
                let mut s = HashSet::new();
                s.insert(*gx);
                s
            });
    }
    debug_assert!(
        d_forwards.len() >= 2,
        "Less than two gametes are used in the redistribution"
    );
    debug_assert!(
        d_backwards.values().all(|s| s.len() >= 1),
        "Gamete with too few parents"
    );
    debug_assert!(
        d_backwards.values().any(|s| s.len() <= 2),
        "Gamete with too many parents"
    );
    debug_assert!(
        d_forwards.values().filter(|s| s.len() > 1).count() <= 2,
        "More than two gametes used in the redistribution"
    );
    for gxs in d_backwards.values() {
        if gxs.len() == 2 {
            let items: [usize; 2] = gxs.iter().cloned().collect::<Vec<_>>().try_into().unwrap();
            return (items[0], items[1]);
        }
    }
    let current_out: Vec<usize> = d_forwards
        .into_iter()
        .filter_map(|(gx, gzs)| if gzs.len() > 1 { Some(gx) } else { None })
        .collect();
    match current_out.len() {
        2 => (current_out[0], current_out[1]),
        1 => (current_out[0], if current_out[0] == 0 { 1 } else { 0 }),
        0 => (0, 1),
        _ => panic!("More than two parents got through this code."),
    }
}

#[cfg(test)]
mod tests {
    use crate::plants::dist_array::dist_array;
    use std::vec;

    use crate::{
        abstract_plants::{Allele, Haploid},
        plants::bit_array::CrosspointBitVec,
    };

    use super::*;

    macro_rules! pretty_print {
        ($xs: expr) => {
            format!(
                "[\n{}\n]",
                $xs.iter()
                    .map(|zs| format!("\t{:?}", zs))
                    .collect::<Vec<_>>()
                    .join("\n")
            )
        };
    }

    #[test]
    fn branching_failed_test() {
        macro_rules! f {
            ($xs: expr, $zs: expr) => {
                let node = AstarNode::new(&DistArray::from(Vec::from($xs)));
                let output: Vec<DistArray> = branching(&node)
                    .iter()
                    .map(|x| x.dist_array())
                    .cloned()
                    .collect();
                assert!(
                    output.contains(&DistArray::from(Vec::from($zs))),
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
                    breeding_program_distribute(&DistArray::from(Vec::from($xs)))
                        .map(|sol| sol.objective),
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
            simplify_dist_array(&dist_array![0, 1, 0, 0, 1, 2, 2]),
            dist_array![0, 1, 0, 1, 2]
        )
    }

    #[test]
    fn generate_redistributions_brute_force_test() {
        assert_eq!(
            vec![
                (dist_array![0, 0, 1], 0, 1),
                (dist_array![0, 1, 0], 0, 1),
                (dist_array![0, 1, 1], 0, 1),
                (dist_array![0, 1, 2], 0, 1),
            ],
            generate_redistributions_brute_force(&dist_array![0, 1, 0])
        );
        assert_eq!(
            vec![
                (dist_array![0, 0, 1, 2], 0, 1),
                (dist_array![0, 1, 0, 0], 0, 2),
                (dist_array![0, 1, 0, 1], 1, 2),
                (dist_array![0, 1, 0, 2], 0, 1),
                (dist_array![0, 1, 1, 2], 0, 1),
                (dist_array![0, 1, 2, 0], 0, 2),
                (dist_array![0, 1, 2, 2], 0, 2),
                (dist_array![0, 1, 2, 3], 0, 2),
            ],
            generate_redistributions_brute_force(&dist_array![0, 1, 0, 2])
        );
    }

    #[test]
    fn generate_redistributions_test() {
        macro_rules! f {
            ($xs: expr, $zs: expr) => {
                let output = generate_redistributions(&DistArray::from(Vec::from($xs)));
                assert!(
                    output.contains(&(DistArray::from(Vec::from($zs.0)), $zs.1, $zs.2)),
                    "{:?} not in {}",
                    $zs,
                    pretty_print!(output)
                )
            };
        }
        f!([0, 1, 0, 2, 1, 0, 2], ([0, 1, 0, 0, 1, 2, 2], 0, 2));

        use std::collections::HashMap;
        macro_rules! f2 {
            ($xs: expr) => {
                let guess = generate_redistributions(&DistArray::from(Vec::from($xs)));
                let deduped_guess: HashMap<_, _> = guess
                    .into_iter()
                    .map(|zss| (format!("{:?}", zss.0.clone()), zss.clone()))
                    .collect();
                let mut guess: Vec<_> = deduped_guess.into_values().collect();
                guess.sort();

                let mut check = generate_redistributions_brute_force(&DistArray::from(Vec::from($xs)));
                check.sort();
                assert!(guess.iter().zip(check.iter()).all(|(x, y)| x.0 == y.0), "assertion failed `guess == check` failed for xs = {:?}\n guess: {}\n check: {}\n",
                    $xs,
                    pretty_print!(guess),
                    pretty_print!(check)
                );
            };
        }
        f2!([0, 1]);
        f2!([0, 1, 0]);
        f2!([0, 1, 2]);
        f2!([0, 1, 0, 1]);
        f2!([0, 1, 0, 2]);
        f2!([0, 1, 2, 0]);
        f2!([0, 1, 2, 1]);
        f2!([0, 1, 0, 1, 0]);
        f2!([0, 1, 2, 0, 1]);
        f2!([0, 1, 2, 1, 0]);
        f2!([0, 1, 0, 2, 0]);
        f2!([0, 1, 0, 2, 0, 1, 0]);
        f2!([0, 1, 0, 2, 1, 0, 2]);
    }

    #[test]
    fn dominates_gametewise_test() {
        macro_rules! f {
            ($xs: expr, $ys: expr) => {
                assert!(
                    dominates_gametewise(
                        &DistArray::from(Vec::from($xs)),
                        &DistArray::from(Vec::from($ys))
                    ),
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
    fn remaining_ranges_test() {
        assert_eq!(
            _remaining_ranges(
                &vec![0, 1, 2, 0, 1, 3],
                0,
                1,
                &vec![true, false, false, true, false]
            ),
            [[4, 6], [0, 1]]
        );

        assert_eq!(
            _remaining_ranges(
                &vec![0, 1, 2, 1, 0, 3],
                0,
                1,
                &vec![true, false, false, true, false]
            ),
            [[1, 4], [4, 1]]
        );
    }

    #[test]
    fn branching_full_join_test() {
        unimplemented!()
    }

    #[test]
    fn path_to_crossing_schedule_test() {
        unimplemented!()
    }

    #[test]
    fn breeding_program_distribute_diving_test() {
        macro_rules! f {
            ($xs: expr, $obj_check: expr) => {
                assert_eq!(
                    //breeding_program_distribute_diving(&DistArray::from(Vec::from($xs))).map(|sol| sol.objective),
                    breeding_program_distribute_diving(&DistArray::from(Vec::from($xs))).map(|sol| sol.objective),
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
    fn breeding_program_distribute_diving_general_test() {
        macro_rules! f {
            ($xs: expr, 0, $ub_check: expr) => {
                if let Some(res) = breeding_program_distribute_general(&DistArray::from(Vec::from($xs)), &Config {
                    full_join: false,
                    dominance: false,
                    diving: true,
                    debug_trace_file: None,
                })
                    .map(|sol| sol.objective).flatten() {
                        assert!(res <= $ub_check, "breeding_program_distribute_general({:?}) returned a solution above the upper bound {}", $xs, $ub_check);
                } else {
                    panic!("breeding_program_distribute_general({:?}) returned None", $xs);
                }
            };
            ($xs: expr, $lb_check: expr, $ub_check: expr) => {
                if let Some(res) = breeding_program_distribute_general(&DistArray::from(Vec::from($xs)), &Config::new(false, false, true, None))
                    .map(|sol| sol.objective).flatten() {
                        assert!(res >= $lb_check, "breeding_program_distribute_general({:?}) returned a solution below the lower bound {}", $xs, $lb_check);
                        assert!(res <= $ub_check, "breeding_program_distribute_general({:?}) returned a solution above the upper bound {}", $xs, $ub_check);
                } else {
                    panic!("breeding_program_distribute_general({:?}) returned None", $xs);
                }
            };
        }
        f!([0], 0, 0);
        f!([0, 1], 2, 2);
        f!([0, 1, 0], 3, 3);
        f!([0, 1, 2], 3, 3);
        f!([0, 1, 0, 1], 3, 4);
        f!([0, 1, 0, 2], 4, 4);
        f!([0, 1, 2, 0], 4, 4);
        f!([0, 1, 2, 1], 4, 4);
        f!([0, 1, 0, 1, 0], 4, 5);
        f!([0, 1, 2, 0, 1], 4, 5);
        f!([0, 1, 2, 1, 0], 4, 5);
        f!([0, 1, 0, 2, 0], 5, 5);
        f!([0, 1, 0, 2, 0, 1, 0], 5, 7);
        f!([0, 1, 0, 2, 1, 0, 2], 5, 7);
    }

    #[test]
    fn generate_redistributions_2_test() {
        let xs = dist_array![0, 1, 0, 2];
        let output = generate_redistributions(&xs);
        assert!(
            vec![
                (dist_array![0, 0, 1, 2], 0, 1),
                (dist_array![0, 1, 0, 2], 0, 1),
                (dist_array![0, 1, 1, 2], 0, 1),
                (dist_array![0, 1, 0, 0], 0, 2),
                (dist_array![0, 1, 0, 1], 1, 2),
                (dist_array![0, 1, 2, 0], 0, 2),
            ]
            .into_iter()
            .all(|x| output.contains(&x)),
            "Output was: {}",
            pretty_print!(output)
        );
    }
}
