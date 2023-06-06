use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solvers::enumerator_dominance::filter_non_dominating_key;
use crate::solvers::greedy_base::min_generations;
use pathfinding::directed::astar::astar;
use std::collections::{HashMap, HashSet};
use std::rc::Rc;

type WG = WGen<SingleChromGenotype, SingleChromGamete>;

pub fn breeding_program(
    n_loci: usize,
    pop_0: Vec<SingleChromGenotype>,
    ideotype: SingleChromGenotype,
) -> Option<WGen<SingleChromGenotype, SingleChromGamete>>
where
{
    // test for feasibility
    // create new population / lift genotypes
    let start: Vec<Rc<WG>> = pop_0.iter().map(|x| Rc::new(x.lift_a())).collect();
    let state_0: State<Rc<WG>> = State::new(start);
    // Initialise priority queue

    let successors = |state: &State<Rc<WG>>| {
        let mut out = vec![];
        // generate single node extensions
        out.extend(successors_single_node_extensios(state));
        // generate zigzag consolidating extensions
        out.extend(successors_zigzag_extensions(state));
        out
    };

    let heuristic = |state: &State<Rc<WG>>| {
        heuristic::<SingleChromGenotype, SingleChromGamete, SegmentBitVec>(n_loci, state)
    };

    let success = |state: &State<Rc<WG>>| match state.head.as_ref() {
        StateData::Start(ref v) => v.iter().any(|wx| wx.genotype == ideotype),
        StateData::Next(x, tail) => x.genotype == ideotype,
    };

    let (v, obj) = astar(&state_0, successors, heuristic, success)?;
    println!("obj: {}", obj);
    let x = v.last()?;
    clone_last_from_state(x.to_owned())
}

fn successors_single_node_extensios(
    state: &State<Rc<WG>>,
) -> impl IntoIterator<Item = (State<Rc<WG>>, usize)> {
    vec![]
}

fn successors_zigzag_extensions(
    state: &State<Rc<WG>>,
) -> impl IntoIterator<Item = (State<Rc<WG>>, usize)> {
    // for each genotype we want to identify the zigzags that are present and create new states
    // that consolidate those zigzags.
    // Question:
    // - is there any symmetry breaking or dominance that can be exploited on these structures
    // - given a maximal zigzag and its state after consolidation, are there subproblems that can
    // be exploited or
    // - True cost is in the number of odd zigzags
    vec![]
}

/// Computes a dp array of the min number of crossings to create a genotype x' that can produce a
/// gamete up to index j from index i where i=0 or i-1 is homozygous 0.
pub fn zigzag_costs_r(x: &SingleChromGenotype) -> Vec<usize> {
    let n_loci = x.get_n_loci(0);
    let column = |i| {
        (
            x.get(false, i).unwrap() as u8,
            x.get(true, i).unwrap() as u8,
        )
    };
    let mut out: Vec<usize> = vec![0; n_loci];

    let mut i = 0;
    while i < n_loci {
        // iterate over all contiguous (0, 0) columns
        while i < n_loci && column(i) == (0, 0) {
            out[i] = 0;
            i += 1;
        }
        if i == n_loci {
            break;
        }

        // iterate over all contiguous (1, 1) columns
        out[i] = 1;
        while i < n_loci && column(i) == (1, 1) {
            out[i] = 1;
            i += 1;
        }
        if i < n_loci && column(i) != (0, 0) {
            out[i] = 1;
        }

        // Now we are going to launch from a heterozygous loci into the rest of a zigzag.
        while i < n_loci && column(i) != (0, 0) {
            let mut j = i + 1;
            while j < n_loci && (column(j) == column(i) || column(j) == (1, 1)) {
                out[j] = out[i];
                j += 1;
            }
            if j == n_loci || column(j) == (0, 0) {
                i = j;
                break;
            }
            out[j] = out[i] + 1;
            i = j;
        }
    }
    for v in out.iter_mut() {
        *v = (*v + 1) >> 1;
    }
    out
}

fn zigzag_switch_points(x: &SingleChromGenotype) -> Vec<bool> {
    let n_loci = x.get_n_loci(0);
    let mut out: Vec<bool> = vec![false; n_loci]; // n_loci 0 would cause n_loci - 1 to overflow
    let f = |c: bool, i: usize| x.get(c, i).unwrap();
    let column = |i: usize| {
        (
            x.get(false, i).unwrap() as u8,
            x.get(true, i).unwrap() as u8,
        )
    };

    let mut i = 0;
    while i < n_loci {
        // iterate over all contiguous (0, 0) columns
        while i < n_loci && column(i) == (0, 0) {
            i += 1;
        }

        // iterate over all contiguous (1, 1) columns
        while i < n_loci && column(i) == (1, 1) {
            i += 1;
        }

        // Now we are going to launch from a heterozygous loci into the rest of a zigzag.
        while i < n_loci && column(i) != (0, 0) {
            let mut j = i + 1;
            while j < n_loci && (column(j) == column(i) || column(j) == (1, 1)) {
                j += 1;
            }
            if j == n_loci || column(j) == (0, 0) {
                i = j;
                break;
            } else {
                out[j] = true;
            }
            i = j;
        }
    }
    out
}

fn consodiate_zigzag(state: &State<Rc<WG>>, x: Rc<WG>, i: usize, j: usize) -> State<WG> {
    // assert that i..j is a zigzag
    let dp = zigzag_costs_r(&x.genotype);
    assert!((i..j).all(|k| dp[k] > 0));
    // assert! zigzag can't be consolidated immediately

    // extract segments and assume that segments are minimal
    // all crossings in the first phase have the same start chromosome
    type Val = Rc<WGam<SingleChromGenotype, SingleChromGamete>>;
    let mut hash_map: HashMap<SingleChromGamete, Val> = HashMap::new();

    let segments_s0: Vec<SegmentMC<Val>> =
        segment_from_range_genotype(&x.genotype, i, j)
            .iter()
            .map(|c| SegmentMC {
                s: c.s,
                e: c.e,
                g: hash_map
                    .entry(c.g.clone())
                    .or_insert_with(|| {
                        Rc::new({
                            let mut wg = WGam::new(c.g.clone());
                            wg.history.push(x.clone());
                            wg
                        })
                    })
                    .clone(),
            })
            .collect();
    // Extract min covering subset
    // Join segments
    assert!({ // all adjacent pairs of segments are joinable
        (0..segments_s0.len()-1).all(|i| {
            let SegmentMC { s: s_x, e: e_x, g: _ } = segments_s0[i];
            let SegmentMC { s: s_y, e: e_y, g: _ } = segments_s0[i+1];
            s_x < s_y && s_y <= e_x && e_x < e_y
        })
    });
    fn join_segments(segments_s0: &Vec<SegmentMC<Val>>, i: usize, j: usize) -> SegmentMC<Val> {
        // The base case can't return the raw segments as they are owned by the vector.
        // Although, we could clone the base segments.
        assert!(segments_s0.len() > 0);
        assert!(i < j && j < segments_s0.len());
        if j == i + 1 {
            return segments_s0[i].clone();
        } else {
            let mid = (i + j) >> 1;
            let res1 = join_segments(segments_s0, i, mid);
            let res2 = join_segments(segments_s0, mid, j);
            SegmentMC::<Val>::join(&res1, &res2)
        }
    }
    // extract genotype from the final segment
    unimplemented!()
}

fn non_dominated_gametes(x: Rc<WG>) -> Vec<Rc<WGam<SingleChromGenotype, SingleChromGamete>>> {
    let mut hashset: HashSet<SingleChromGamete> = HashSet::new();
    segment_from_range_genotype(&x.genotype, 0, x.genotype.get_n_loci(0))
        .iter()
        .map(|c| c.g.clone())
        .collect::<HashSet<SingleChromGamete>>()
        .iter()
        .map(|g| {
            Rc::new({
                let mut wg = WGam::new(g.clone());
                wg.history.push(x.clone());
                wg
            })
        })
        .collect()
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

fn naive_successors<A, B, K>(state: &State<WGen<A, B>>) -> Vec<(State<WGen<A, B>>, usize)> {
    let mut out = vec![];
    let xs = state.iter().map(|x| x).collect::<Vec<&WGen<A, B>>>();
    for i in 0..xs.len() {
        for j in i..xs.len() {
            for z in non_dominated_crossings::<A, B, K>(xs[i], xs[j]) {
                out.push((
                    State {
                        head: Rc::new(StateData::Next(z, state.head.clone())),
                        hash: todo!(),
                    },
                    1,
                ))
            }
        }
    }
    out
}

pub fn non_dominated_crossings<A, B, K>(x: &WGen<A, B>, y: &WGen<A, B>) -> Vec<WGen<A, B>> {
    unimplemented!()
}

fn heuristic<A, B, S>(n_loci: usize, state: &State<Rc<WGen<A, B>>>) -> usize
where
    A: Genotype<B> + SingleChrom + Diploid<B> + Clone,
    B: Gamete<A> + SingleChrom + Haploid,
    S: HaploidSegment<A, B> + Clone,
{
    let n_generations: usize = todo!();
    let n_crossings_div2: usize = todo!();
    n_generations.max(n_crossings_div2)
}

fn clone_last_from_state<A: Clone, B>(state: State<Rc<WGen<A, B>>>) -> Option<WGen<A, B>> {
    match state.head.as_ref() {
        StateData::Start(v) => v.last().map(|x| x.as_ref().clone()),
        StateData::Next(x, n) => Some(x.as_ref().clone()),
    }
}

#[derive(Clone, Copy)]
struct SegmentMC<T> {
    pub s: usize,
    pub e: usize,
    pub g: T,
}

impl SegmentMC<SingleChromGamete> {
    fn join(&self, other: &Self) -> Self {
        SegmentMC {
            s: self.s,
            e: other.e,
            g: {
                let z = SingleChromGenotype::from_gametes(&self.g, &other.g);
                let k_z = CrosspointBitVec::new(false, other.e);
                k_z.cross(&z)
            },
        }
    }
}

impl SegmentMC<Rc<SingleChromGamete>> {
    fn join(&self, other: &Self) -> Self {
        SegmentMC {
            s: self.s,
            e: other.e,
            g: {
                let z = Rc::new(SingleChromGenotype::from_gametes(&self.g, &other.g));
                let k_z = CrosspointBitVec::new(false, other.s);
                Rc::new(k_z.cross(z.as_ref()))
            },
        }
    }
}

impl SegmentMC<Rc<WGam<SingleChromGenotype, SingleChromGamete>>> {
    fn join(&self, other: &Self) -> Self {
        SegmentMC {
            s: self.s,
            e: other.e,
            g: {
                todo!()
            },
        }
    }
}

fn segment_from_range_genotype(
    x: &SingleChromGenotype,
    start: usize,
    end: usize,
) -> Vec<SegmentMC<SingleChromGamete>> {
    let q1_true = segment_from_range_gamete(&x.upper(), start, end);
    let q2_true = segment_from_range_gamete(&x.lower(), start, end);
    let mut q = (&q1_true, &q2_true);

    let mut used1_true = vec![false; q1_true.len()];
    let mut used2_true = vec![false; q2_true.len()];
    let mut used1 = &mut used1_true;
    let mut used2 = &mut used2_true;

    let mut i = 0;
    let mut j = 0;

    let mut out: Vec<SegmentMC<SingleChromGamete>> =
        Vec::with_capacity(q1_true.len() + q2_true.len());

    // oh dear god, does this even work?
    while i < q.0.len() && j < q.1.len() {
        let c1 = q.0[i].clone();
        let c2 = q.1[j].clone();
        let (s1, e1) = (c1.s, c1.e);
        let (s2, e2) = (c2.s, c2.e);

        if s1 > s2 || s1 == s2 && e1 < e2 {
            q = (q.1, q.0);
            (used1, used2) = (used2, used1);
            (i, j) = (j, i);
        } else if s1 <= s2 && e1 >= e2 {
            j += 1;
        } else if e1 + 1 < s2 {
            assert!(s1 < e1 + 1 && e1 + 1 < s2 && s2 <= e2);
            if !used1[i] {
                out.push(c1);
                used1[i] = true;
            }
            i += 1;
        } else {
            out.push(SegmentMC::<SingleChromGamete>::join(&c1, &c2));
            used1[i] = true;
            used2[j] = true;
            // break early if possible
            if e2 == end - 1 {
                i = q.0.len();
                j = q.1.len();
            } else {
                i += 1;
            }
        }
    }
    q.0[i..].iter().for_each(|c| out.push(c.clone()));
    q.1[j..].iter().for_each(|c| out.push(c.clone()));
    // TODO: resize to exact len //
    out
}

fn segment_from_range_gamete(
    gx: &SingleChromGamete,
    i: usize,
    j: usize,
) -> Vec<SegmentMC<SingleChromGamete>> {
    let mut out: Vec<SegmentMC<SingleChromGamete>> = Vec::with_capacity(gx.get_n_loci(0));
    let mut indexes = (i..j).filter(|k| gx.index(*k) == Allele::O);
    if let Some(s) = indexes.next() {
        let mut current_segment = SegmentMC {
            s,
            e: s + 1,
            g: gx.to_owned(),
        };
        let mut last = s;
        for k in indexes {
            if k > last + 1 {
                out.push(current_segment.clone());
                current_segment = SegmentMC {
                    s: k,
                    e: k + 1,
                    g: gx.to_owned(),
                };
            }
            current_segment.e = k + 1;
            last = k;
        }
        out.push(current_segment.clone());
    }
    // TODO: resize out to exact capacity //
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn breeding_program_test() {
        assert_ne!(
            breeding_program(
                3,
                vec![SingleChromGenotype::from_str("010", "101",)],
                SingleChromGenotype::ideotype(3)
            ),
            None
        );
    }

    #[test]
    fn zigzag_cost_test() {
        assert_eq!(
            zigzag_costs_r(&SingleChromGenotype::from_str(
                "0010110010001110010",
                "0101011100110011010"
            )),
            vec![0, 1, 1, 2, 2, 2, 3, 3, 3, 0, 1, 1, 1, 1, 1, 2, 0, 1, 0]
        );
    }

    #[test]
    fn heuristic_test() {
        panic!()
    }
}
