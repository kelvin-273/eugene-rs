use crate::solution::{BaseSolution, PyBaseSolution};
use crate::solvers::base_min_generations_enumerator_dominance::IteratorNonDominating;
use core::cmp::Reverse;
use std::collections::{BinaryHeap, HashSet};
use std::fs;
use std::io::Write;
use std::rc::Rc;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

struct Config {
    full_join: bool,
    dominance: bool,
    diving: bool,
}

#[derive(pyo3::IntoPyObject)]
pub struct Output {
    expansions: usize,
    pushed_nodes: usize,
    children_created_by_branching: usize,
    objective: Option<usize>,
}

impl Output {
    pub fn new() -> Self {
        Self {
            expansions: 0,
            pushed_nodes: 0,
            children_created_by_branching: 0,
            objective: None,
        }
    }
}

/// Returns the number of crossings required for distribute array `xs` along with several other
/// runtime statistics.
#[pyo3::pyfunction]
#[pyo3(signature = (xs, timeout=None, diving=false, full_join=false, dominance=false))]
pub fn breeding_program_distribute_general_python(
    xs: DistArray,
    timeout: Option<u64>,
    diving: bool,
    full_join: bool,
    dominance: bool,
) -> pyo3::PyResult<Option<Output>> {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = breeding_program_distribute_general(
            &xs,
            &Config {
                diving,
                full_join,
                dominance,
            },
        );
        tx.send(res)
    });
    let res = rx
        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
        .ok()
        .flatten();
    Ok(res)
}

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
pub fn breeding_program_distribute_no_full_join_python(
    xs: DistArray,
    timeout: Option<u64>,
) -> PyBaseSolution {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = breeding_program_distribute_no_full_join(&xs);
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
pub fn breeding_program_distribute_dominance_python(
    xs: DistArray,
    timeout: Option<u64>,
) -> PyBaseSolution {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = breeding_program_distribute_dominance(&xs);
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
pub fn breeding_program_distribute_dominance_no_full_join_python(
    xs: DistArray,
    timeout: Option<u64>,
) -> PyBaseSolution {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = breeding_program_distribute_dominance_no_full_join(&xs);
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
pub fn breeding_program_distribute_diving_python(
    xs: DistArray,
    timeout: Option<u64>,
) -> PyBaseSolution {
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

pub type DistArray = Vec<usize>;

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

fn breeding_program_distribute_general(xs: &DistArray, config: &Config) -> Option<Output> {
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

    while let Some(Reverse(node)) = open_list.pop() {
        output.expansions += 1;
        if closed_list.contains(&node) {
            continue;
        }

        closed_list.insert(&node);
        if node.success() {
            return Some(node);
        }
        let children = branching_general(&node, None, config);
        output.children_created_by_branching += children.len();
        for child in children {
            if !closed_list.contains(&child) {
                open_list.push(Reverse(child));
                output.pushed_nodes += 1;
            }
        }
    }
    None
}

fn branching_general(
    node: &AstarNode,
    closed_list: Option<&ClosedList>,
    config: &Config,
) -> Vec<AstarNode> {
    if config.full_join {
        if let Some((zs, gx, gy)) = first_full_join(node.dist_array()) {
            let zs = simplify_dist_array(&zs);
            return vec![node.create_offspring(zs, gx, gy)];
        }
    }
    if config.dominance {
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
    } else if config.diving {
        // TODO: return either the minimum one or the one that meets the lower bound
        let lower_bound = node.n_segments() + node.n_pop();

        let mut min_obj = usize::MAX;
        let mut min_dist_arr = (node.dist_array().clone(), 0, 0);

        for (zs, gx, gy) in generate_redistributions(node.dist_array()) {
            if &zs != node.dist_array() {
                let zs_simplified = simplify_dist_array(&zs);
                let n_pop = zs_simplified.iter().max().expect("empty dist_array") + 1;
                let n_segments = zs_simplified.len();
                if n_pop + n_segments == lower_bound {
                    return vec![node.create_offspring(zs_simplified, gx, gy)];
                } else if n_pop + n_segments < min_obj {
                    min_dist_arr = (zs_simplified, gx, gy);
                    min_obj = n_pop + n_segments;
                }
            }
        }
        let (zs, gx, gy) = min_dist_arr;
        vec![node.create_offspring(zs, gx, gy)]
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
        let children = branching(&node, None);
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
        let children = branching_no_full_join(&node, None);
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
        let children = branching_dominance(&node, None);
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
        let children = branching_dominance_no_full_join(&node, None);
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

fn branching_no_full_join(node: &AstarNode, closed_list: Option<&ClosedList>) -> Vec<AstarNode> {
    generate_redistributions(node.dist_array())
        .into_iter()
        .filter(|(zs, _, _)| zs != node.dist_array())
        .map(|(zs, gx, gy)| (simplify_dist_array(&zs), gx, gy))
        .map(|(zs, gx, gy)| node.create_offspring(zs, gx, gy))
        .collect()
}

fn branching_dominance_no_full_join(
    node: &AstarNode,
    closed_list: Option<&ClosedList>,
) -> Vec<AstarNode> {
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

fn branching(node: &AstarNode, closed_list: Option<&ClosedList>) -> Vec<AstarNode> {
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

fn branching_dominance(node: &AstarNode, closed_list: Option<&ClosedList>) -> Vec<AstarNode> {
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
            dominates_gametewise(&zs1, &zs2) || dominates_as_subsequence(&zs1, &zs2)
        })
        .map(|(zs, gx, gy)| node.create_offspring(zs, gx, gy))
        .collect()
}

fn nondom_test() {
    let u = vec![1, 2, 3];
    let v: Vec<i32> = u
        .into_iter()
        .filter_non_dominating_fn(|x, y| x > y)
        .map(|x| x)
        .collect();
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
        segjoin_choice: &'a mut Vec<bool>,
        out: &'a mut Vec<DistArray>,
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
        fn bt_final_gamete_construction(state: &mut TwoDeltaState, i: usize) {
            let [x_rng, y_rng] =
                remaining_ranges(state.xs, state.gx, state.gy, state.segjoin_choice);
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
    use std::{os::unix::raw::nlink_t, vec};

    use crate::{
        abstract_plants::{Allele, Crosspoint, Haploid, WGen}, extra::instance_generators::distarray_to_homo, plants::bit_array::CrosspointBitVec
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

    #[test]
    fn path_to_crossing_schedule_test() {
        use crate::abstract_plants::Chrom;
        use crate::plants::bit_array::SingleChromGenotype;
        use crate::abstract_plants::WGen;
        let dist_array = vec![0, 1];
        let n_loci = dist_array.len();
        let pop_0 = SingleChromGenotype::init_pop_distribute(&dist_array);

        // check that the initial population is valid
        assert_eq!(pop_0.len(), 2);
        assert!((0..pop_0.len())
            .all(|i| (0..n_loci).all(|j| { pop_0[i].get(true, j) == pop_0[i].get(false, j) })));
        assert!((0..pop_0.len()).all(|i| {
            (0..n_loci).all(|j| {
                pop_0[i]
                    .get(true, j)
                    .map(|b| b == (dist_array[j] == i))
                    .expect("get should not return None")
            })
        }));

        // construct the path using astar and check the output
        let path = astar(&dist_array).expect("astar failed");
        let mut node_ref = path.clone();
        let mut out_distarrays = vec![path.dist_array().clone()];
        let mut out_parents = vec![path.parent_gametes()];
        while node_ref.parent_node().is_some() {
            node_ref = node_ref.parent_node().unwrap();
            out_distarrays.push(node_ref.dist_array().clone());
            out_parents.push(node_ref.parent_gametes());
        }
        out_distarrays.reverse();
        out_parents.reverse();
        assert_eq!(out_distarrays, vec![vec![0, 1], vec![0],]);
        assert_eq!(out_parents, vec![None, Some((0, 1))]);

        // create a crossing schedule from the path

        let mut wgametes = pop_0
            .iter()
            .map(|x| Some(WGen::new(x.clone()).cross(CrosspointBitVec::new(Chrom::Upper, n_loci))))
            .collect::<Vec<_>>();
        wgametes.resize(n_loci, None);
        for (i, (dist_array, parent_gametes)) in
            out_distarrays.iter().zip(out_parents.iter()).enumerate()
        {
            if i == 0 {
                // first dist_array is the initial population
                continue;
            }
            let prev_dist_array = &out_distarrays[i - 1];
            let (gx, gy) = parent_gametes.expect("parent gametes should be set");
            let (wgx, wgy) = (
                wgametes[gx].as_ref().expect("gx should be set"),
                wgametes[gy].as_ref().expect("gy should be set"),
            );
            let wz = WGen::from_gametes(wgx, wgy);

            // infer the redistribution from the previous dist_array and the current dist_array
            let new_gametes = (0..dist_array.len())
                .filter(|j| prev_dist_array[*j] == gx || prev_dist_array[*j] == gy)
                .map(|j| dist_array[j])
                .collect::<HashSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            for new_gamete in new_gametes {
                // infer the crosspoint from the new gamete and the previous gametes
                let mut crosspoint = CrosspointBitVec::new(Chrom::Upper, n_loci);
                let mut first_source = None;
                for j in 0..dist_array.len() {
                    if dist_array[j] == gx || dist_array[j] == gy {
                        if first_source.is_none() {
                            first_source = Some(dist_array[j]);
                        } else if Some(dist_array[j]) == first_source {
                            continue;
                        } else {
                            let crosspoint_locus = j - 1;
                            let crosspoint_chrom = Chrom::from(Some(gx) == first_source);
                            crosspoint = CrosspointBitVec::new(crosspoint_chrom, crosspoint_locus);
                            break;
                        }
                    }
                }
                ///////////////
                //  BIG BUG  //
                ///////////////
                // The crosspoint is not always correct, because the size of the distribute array
                // decreases as we go through the path.
                // Solution: we need to figure out how to undo the simplification of the
                // between the current and previous distribute arrays and then apply that operation
                // to the remaining distribute arrays in the path.

                /////////////////////////
                //  THE OTHER BIG BUG  //
                /////////////////////////
                // The crosspoint is also not always correct because the simplification of the
                // distribute also relabels the gametes.
                // as well.

                // create a new WGam for the new gamete
                let wgz = wz.cross(crosspoint); 
                wgametes[new_gamete] = Some(wgz);
            }
        }
        // check that the first index of the wgametes is the target gamete
        assert!(
            wgametes[0]
                .as_ref()
                .expect("wgametes[0] should be set")
                .gamete()
                .alleles()
                .iter()
                .all(|&a| a == Allele::O),
        );
    }
}
