use crate::abstract_plants::*;
use crate::solvers::greedy_base::min_generations;
use pathfinding::directed::astar::astar;
use std::rc::Rc;

pub fn breeding_program<A, B, K, S>(n_loci: usize, pop_0: Vec<A>, ideotype: A) -> Option<WGen<A, B>>
where
    A: Genotype<B> + SingleChrom + Diploid<B> + Clone + PartialEq,
    B: Gamete<A> + SingleChrom + Haploid + PartialEq,
    K: Crosspoint<A, B, usize>,
    S: Segment<A, B> + HaploidSegment<A, B> + Clone,
{
    // test for feasibility
    // create new population / lift genotypes
    let start: Vec<WGen<A, B>> = pop_0.iter().map(|x| x.lift_a()).collect();
    let state_0: State<WGen<A, B>> = State::new(start);
    // Initialise priority queue

    let successors = |state: &State<WGen<A, B>>| unimplemented!();

    let heuristic = |state: &State<WGen<A, B>>| {
        min_generations::<A, B, K, S>(
            n_loci,
            &state.iter().map(|wx| wx.genotype.clone()).collect(),
        )
    };

    let success = |state: &State<WGen<A, B>>| match state.head.as_ref() {
        StateData::Start(ref v) => v.iter().any(|wx| wx.genotype == ideotype),
        StateData::Next(x, tail) => x.genotype == ideotype,
    };

    astar(&state_0, successors, heuristic, success);
    None
}

impl<A: Clone, B> Clone for WGen<A, B> {
    fn clone(&self) -> Self {
        WGen {
            genotype: self.genotype.clone(),
            history: self.history.clone(),
        }
    }
}

enum StateData<A> {
    Start(Vec<A>),
    Next(A, Rc<StateData<A>>),
}

type Link<A> = Rc<StateData<A>>;

#[derive(Clone)]
struct State<A> {
    head: Link<A>,
    hash: u32,
}

impl<A> State<A> {
    pub fn new(v: Vec<A>) -> Self {
        Self {
            head: Rc::new(StateData::Start(v)),
            hash: 0,
        }
    }

    pub fn count_nodes(&self) -> usize {
        let mut head = &self.head;
        fn aux<A>(head: &Link<A>) -> usize {
            match head.as_ref() {
                StateData::Start(_) => 0,
                StateData::Next(_, link) => 1 + aux(link),
            }
        }
        aux(head)
    }
}

impl<A> std::hash::Hash for State<A> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        state.write_u32(self.hash)
    }
}

impl<A: PartialEq> PartialEq for State<A> {
    fn eq(&self, other: &Self) -> bool {
        self.head == other.head
    }
}

impl<A: PartialEq> Eq for State<A> {}

struct StateIter<'a, A: 'a> {
    current: &'a StateData<A>,
    vector_idx: usize,
}

impl<'a, A> State<A> {
    fn iter(&'a self) -> StateIter<'a, A> {
        StateIter {
            current: &*self.head,
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
