use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solution::*;
use crate::solvers::base_min_generations_enumerator_dominance::filter_non_dominating_key;
use crate::solvers::base_min_generations_segment;
use pathfinding::directed::astar::astar;
use std::rc::Rc;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

type WGe = WGen<SingleChromGenotype, SingleChromGamete>;
type WGa = WGam<SingleChromGenotype, SingleChromGamete>;

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
        .collect();
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = breeding_program(n_loci, &pop_0)
            .map(|x_star| x_star.to_base_sol(n_loci, Objective::Crossings));
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

pub fn breeding_program(n_loci: usize, pop_0: &Vec<SingleChromGenotype>) -> Option<WGe>
where
{
    // TODO: test for feasibility

    let wrapped_pop_0: Vec<WGe> = pop_0.iter().map(|x| WGen::new(x.clone())).collect();
    let state_0: State<WGe> = State::new(wrapped_pop_0);
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    let g_star = SingleChromGamete::ideotype(n_loci);

    let successors = |state: &State<WGe>| successors_single_node_extensions(n_loci, &g_star, state);
    let heuristic = |state: &State<WGe>| heuristic(n_loci, state);
    //let heuristic = |state: &State<WGe>| heuristic_baseline(n_loci, state);
    let success = |state: &State<WGe>| match state.head.as_ref() {
        StateData::Start(v) => v.iter().any(|wx| wx.genotype() == &ideotype),
        StateData::Next(x, _tail) => x.genotype() == &ideotype,
    };

    let (v, _obj) = astar(&state_0, successors, heuristic, success)?;

    // TODO: This is probably a confusing sequence of methods (last followed by first?)
    let state_star = v.last()?;
    state_star.first().cloned()
}

fn heuristic(n_loci: usize, state: &State<WGe>) -> usize {
    let n_segments = base_min_generations_segment::n_min_covering_segments(
        n_loci,
        &state.iter().map(|wgx| wgx.genotype()).cloned().collect(),
    );
    (n_segments + 1) >> 1
}

fn heuristic_baseline(_n_loci: usize, _state: &State<WGe>) -> usize {
    0
}

fn successors_single_node_extensions(
    n_loci: usize,
    g_star: &SingleChromGamete,
    state: &State<WGe>,
) -> impl IntoIterator<Item = (State<WGe>, usize)> {
    let current_cost = state.depth;

    // Generate all non dominated gametes from each genotype in state (wrapped)
    let gametes: Vec<WGa> = filter_non_dominating_key::<_, _, _, DomGamete>(
        state
            .iter()
            .flat_map(|wx| CrosspointBitVec::crosspoints(&n_loci).map(|k| wx.cross(k))),
        |wg: &WGa| wg.gamete(),
    );
    let n_gametes = gametes.len();

    // Return the ideotype if it can be created
    if let Some(wg) = gametes.iter().find(|wg| wg.gamete() == g_star) {
        let wx_star = WGen::from_gametes(wg, wg);
        let new_state = state.push(wx_star);
        return vec![(new_state, 2)];
    }

    // Return a vector of WGenS2 genotypes from all pairs of non-equivalent gametes up to symmetry.
    // Note: non-equivalent because a genotype (gx, gx) is dominated by all other (gx, gy) unless
    // gx is the g* gamete
    (0..n_gametes - 1)
        .flat_map(|i| (i + 1..n_gametes).map(move |j| (i, j)))
        .map(|(i, j)| {
            (
                state.push(WGen::from_gametes(&gametes[i], &gametes[j])),
                current_cost + 1,
            )
        })
        .collect::<Vec<_>>()
}

#[derive(Debug)]
enum StateData<A> {
    Start(Vec<A>),
    Next(A, Rc<StateData<A>>),
}

type Link<A> = Rc<StateData<A>>;

#[derive(Clone, Debug)]
struct State<A> {
    head: Link<A>,
    most_recent_parent: usize,
    depth: usize,
}

impl<A> State<A> {
    pub fn new(v: Vec<A>) -> Self {
        Self {
            head: Rc::new(StateData::Start(v)),
            most_recent_parent: 0,
            depth: 0,
        }
    }

    pub fn push(&self, x: A) -> Self {
        Self {
            most_recent_parent: 0,
            head: Rc::new(StateData::Next(x, self.head.clone())),
            depth: self.depth + 1,
        }
    }

    pub fn _count_nodes(&self) -> usize {
        let head = &self.head;
        fn aux<A>(head: &Link<A>) -> usize {
            match head.as_ref() {
                StateData::Start(_) => 0,
                StateData::Next(_, link) => 1 + aux(link),
            }
        }
        aux(head)
    }

    pub fn first(&self) -> Option<&A> {
        self.iter().next()
    }
}

impl<A> std::hash::Hash for State<A> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        state.write_usize(self.most_recent_parent)
        //state.write_u8(0)
    }
}

impl PartialEq for State<WGe> {
    fn eq(&self, other: &Self) -> bool {
        //match (self.head.as_ref(), other.head.as_ref()) {
        //    (StateData::Start(ref v1), StateData::Start(ref v2)) => true, // Assumes v1 == v2
        //    (StateData::Start(_), StateData::Next(_, _)) => false,
        //    (StateData::Next(_, _), StateData::Start(_)) => false,
        //    (StateData::Next(x1, _), StateData::Next(x2, _)) => x1.genotype == x2.genotype,
        //}
        let mut v1: Vec<_> = self.iter().map(|x| x.genotype()).collect();
        let mut v2: Vec<_> = other.iter().map(|x| x.genotype()).collect();
        v1.sort();
        v2.sort();
        v1 == v2
    }
}

//impl<A: PartialEq> Eq for State<A> {}

impl Eq for State<WGe> {}

struct StateIter<'a, A: 'a> {
    current: &'a StateData<A>,
    vector_idx: usize,
}

impl<'a, A> State<A> {
    fn iter(&'a self) -> StateIter<'a, A> {
        StateIter {
            current: self.head.as_ref(),
            vector_idx: 0,
        }
    }
}

impl<A: PartialEq> PartialEq for StateData<A> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (StateData::Start(v1), StateData::Start(v2)) => v1 == v2,
            (StateData::Next(x1, tail1), StateData::Next(x2, tail2)) => {
                x1 == x2 && tail1.as_ref() == tail2.as_ref()
            }
            (_, _) => false,
        }
    }
}

impl<A: Eq> Eq for StateData<A> {}

impl<'a, A> Iterator for StateIter<'a, A> {
    type Item = &'a A;

    fn next(&mut self) -> Option<Self::Item> {
        match self.current {
            StateData::Start(ref v) => {
                if self.vector_idx >= v.len() {
                    return None;
                }
                let x = &v[self.vector_idx];
                self.vector_idx += 1;
                Some(x)
            }
            StateData::Next(ref x, ref tail) => {
                self.current = tail.as_ref();
                Some(x)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn breeding_program_test() {
        let n_loci = 3;
        let res = breeding_program(n_loci, &vec![SingleChromGenotype::from_str("010", "101")]);
        assert_ne!(None, res);
        assert_eq!(
            Some(2),
            res.clone().map(|wx_star| wx_star
                .to_base_sol(n_loci, Objective::Generations)
                .objective)
        );
        assert_eq!(
            Some(2),
            res.map(|wx_star| wx_star.to_base_sol(n_loci, Objective::Crossings).objective)
        );
    }

    use crate::extra::instance_generators;
    macro_rules! distribute_tester {
        ($left:expr, $right:expr) => {
            let pop_0 = instance_generators::distarray_to_homo($right);
            let n_loci = $right.len();
            let op_sol = breeding_program(n_loci, &pop_0)
                .map(|sol| sol.to_base_sol(n_loci, Objective::Crossings));
            dbg!(&op_sol);
            assert_eq!($left, op_sol.expect("infeasible").objective);
        };
    }
    #[test]
    fn breeding_program_dist_test() {
        distribute_tester!(2, &[0, 1]);
        distribute_tester!(3, &[0, 1, 0]);
        distribute_tester!(3, &[0, 1, 2]);
        distribute_tester!(3, &[0, 1, 0, 1]);
        distribute_tester!(4, &[0, 1, 0, 2]);
        distribute_tester!(4, &[0, 1, 2, 0]);
        distribute_tester!(4, &[0, 1, 2, 1]);
        distribute_tester!(4, &[0, 1, 2, 3]);
        distribute_tester!(4, &[0, 1, 0, 1, 0]);
        distribute_tester!(4, &[0, 1, 0, 1, 2]);
        distribute_tester!(4, &[0, 1, 0, 2, 1]);
        distribute_tester!(5, &[0, 1, 0, 2, 0]);
    }
    #[test]
    #[ignore]
    fn breeding_program_dist_long_test() {
        distribute_tester!(5, &[0, 1, 0, 2, 0, 2]);
        distribute_tester!(6, &[0, 1, 2, 3, 4, 5]);
        distribute_tester!(5, &[0, 1, 2, 3, 2, 1]);
    }

    #[test]
    fn breeding_program_zigzag_test() {
        let res = breeding_program(
            4,
            &vec![
                SingleChromGenotype::from_str("1010", "1010"),
                SingleChromGenotype::from_str("0101", "0101"),
            ],
        );
        assert_ne!(None, res);
        let sol = res
            .expect("feasible instance returned None")
            .to_base_sol(4, Objective::Crossings);
        assert_eq!(3, sol.objective);
    }
}
