use crate::abstract_plants::*;
use crate::solvers::greedy_base::min_generations;
use pathfinding::directed::astar::astar;
use std::rc::Rc;

pub fn breeding_program<A, B, K, S>(n_loci: usize, pop_0: Vec<A>, ideotype: A) -> Option<WGen<A, B>>
where
    A: Genotype<B> + SingleChrom + Diploid<B> + Clone,
    B: Gamete<A> + SingleChrom + Haploid,
    K: Crosspoint<A, B, usize>,
    S: Segment<A, B> + HaploidSegment<A, B> + Clone,
{
    // test for feasibility
    // create new population / lift genotypes
    let start: Vec<WGen<A, B>> = pop_0.iter().map(|x| x.lift_a()).collect();
    let state_0: State<WGen<A, B>> = State::new(start);
    // Initialise priority queue

    fn sucessors() {}

    let heuristic = |state: &State<WGen<A, B>>| min_generations::<A, B, K, S>(n_loci, &state.iter().map(|wx| wx.genotype.clone()).collect());

    fn succuess() {}

    //astar(&state_0)
    None
}

enum StateData<A> {
    Start(Vec<A>),
    Next(A, Rc<StateData<A>>),
}

type Link<A> = Rc<StateData<A>>;

struct State<A> {
    head: Link<A>,
}

impl<A> State<A> {
    pub fn new(v: Vec<A>) -> Self {
        Self {
            head: Rc::new(StateData::Start(v)),
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

struct StateIter<'a, A: 'a> {
    current: &'a StateData<A>,
    vector_idx: usize
}

impl<'a, A> State<A> {
    fn iter(&'a self) -> StateIter<'a, A> {
        StateIter {
            current: &*self.head,
            vector_idx: 0
        }
    }
}

impl<'a, A> IntoIterator for State<A> {
    type Item = A;
    type IntoIter = StateIter<'a, A>;
    fn into_iter(self) -> Self::IntoIter {
    }
}

impl<'a, A> Iterator for StateIter<'a, A> {
    type Item = &'a A;

    fn next(&mut self) -> Option<Self::Item> {
        match self.current {
            StateData::Start(ref v) => {
                if self.vector_idx >= v.len() { return None; }
                let x = &v[self.vector_idx];
                self.vector_idx += 1;
                Some(x)
            },
            StateData::Next(ref x, ref tail) => {
                self.current = tail.as_ref();
                Some(x)
            },
        }
    }
}
