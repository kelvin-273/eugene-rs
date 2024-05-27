use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solution::*;
use crate::solvers::base_min_generations_segment;
use crate::solvers::base_min_generations_segment::Segment;
use pathfinding::directed::astar::astar;
use pyo3::prelude::*;
use std::collections::HashMap;
use std::rc::Rc;

type WGe = WGenS2<SingleChromGenotype, SingleChromGamete>;
type WGa = WGamS2<SingleChromGenotype, SingleChromGamete>;

/// Runs a breeding program given `n_loci` and `pop_0` where `pop_0` is a population of single
/// chromosome diploid genotypes with `n_loci` loci.
#[pyo3::pyfunction]
pub fn breeding_program_python(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<bool>>>,
) -> PyResult<
    Option<(
        Vec<Vec<Vec<i32>>>,
        Vec<&'static str>,
        Vec<usize>,
        Vec<usize>,
        usize,
    )>,
> {
    let ideotype = SingleChromGenotype::ideotype(n_loci);
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
    let res = breeding_program(n_loci, pop_0, ideotype);
    match res {
        None => Ok(None),
        Some(x_star) => {
            let sol = x_star.to_base_sol(n_loci, Objective::Crossings);
            Ok(Some((
                sol.tree_data,
                sol.tree_type,
                sol.tree_left,
                sol.tree_right,
                sol.objective,
            )))
        }
    }
}

pub fn breeding_program(
    n_loci: usize,
    pop_0: Vec<SingleChromGenotype>,
    ideotype: SingleChromGenotype,
) -> Option<WGe>
where
{
    // test for feasibility
    // create new population / lift genotypes
    let start: Vec<WGe> = pop_0.iter().map(|x| WGenS2::new(x.clone())).collect();
    let state_0: State<WGe> = State::new(start);
    // Initialise priority queue

    let successors = |state: &State<WGe>| successors_single_node_extensions(state);

    let heuristic = |state: &State<WGe>| heuristic(n_loci, state);

    let success = |state: &State<WGe>| match state.head.as_ref() {
        StateData::Start(ref v) => v.iter().any(|wx| wx.genotype() == &ideotype),
        StateData::Next(x, _tail) => x.genotype() == &ideotype,
    };

    let (v, _obj) = astar(&state_0, successors, heuristic, success)?;
    let x = v.last()?;
    clone_last_from_state::<SingleChromGenotype, SingleChromGenotype>(x)
    //match x.head.as_ref() {
    //StateData::Next(wx, _) => None,
    //StateData::Start(v) => None,
    //}
}

fn successors_single_node_extensions(
    state: &State<WGe>,
) -> impl IntoIterator<Item = (State<WGe>, usize)> {
    // vector of the genotype indices and rc gametes
    let gametes: Vec<(usize, WGa)> = state
        .iter()
        .enumerate()
        .flat_map(|(i, x)| non_dominated_gametes(x).into_iter().map(move |gx| (i, gx)))
        .collect();
    let first_j = gametes
        .iter()
        .enumerate()
        .find(|(_i_gamete, (i, _gx))| i >= &state.most_recent_parent)
        .map_or(0, |(i, _igx)| i);

    let wg_star = gametes
        .iter()
        .filter(|(_i, wg)| wg.gamete().alleles().iter().all(|a| *a == Allele::O))
        .next();
    if let Some((_, wg)) = wg_star {
        let wx_star = WGenS2::from_gametes(wg, wg);
        let new_state = state.push(wx_star);
        // TODO: halt the astar when this is found <11-06-23> //
        return vec![(new_state, 2)];
    }

    let n_gametes = gametes.len();
    let out = (0..gametes.len())
        .flat_map(|i| ((i + 1).max(first_j)..n_gametes).map(move |j| (i, j)))
        .map(|(i, j)| {
            let (_idx, wg1) = &gametes[i];
            let (idy, wg2) = &gametes[j];
            let wz = WGenS2::from_gametes(&wg1, &wg2);
            let new_state = state.push_with_most_recent_parent(wz, *idy);
            (new_state, 2)
        })
        .collect::<Vec<(State<WGe>, usize)>>();
    out
}

fn non_dominated_gametes(x: &WGe) -> Vec<WGa> {
    let res = segments_from_range_genotype(x.genotype(), 0, x.genotype().get_n_loci(0));
    res.iter()
        .map(|c| WGamS2::new_from_genotype(c.g.clone(), x.clone()))
        .collect()
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

    pub fn push_with_most_recent_parent(&self, x: A, most_recent_parent: usize) -> Self {
        Self {
            most_recent_parent,
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
                        most_recent_parent: j,
                        depth: state.depth + 1,
                    },
                    1,
                ))
            }
        }
    }
    out
}

pub fn non_dominated_crossings<A, B, K>(_x: &WGen<A, B>, _y: &WGen<A, B>) -> Vec<WGen<A, B>> {
    unimplemented!()
}

fn heuristic(n_loci: usize, state: &State<WGe>) -> usize {
    let n_segments = base_min_generations_segment::min_covering_segments(
        n_loci,
        &state.iter().map(|wgx| wgx.genotype().clone()).collect(),
    )
    .len();
    (n_segments + 1) >> 1
}

fn clone_last_from_state<A, B>(state: &State<WGe>) -> Option<WGe> {
    match state.head.as_ref() {
        StateData::Start(v) => v.last().map(|x| x.clone()),
        StateData::Next(x, _n) => Some(x.clone()),
    }
}

fn segments_from_range_genotype(
    x: &SingleChromGenotype,
    start: usize,
    end: usize,
) -> Vec<Segment<SingleChromGamete>> {
    let q1_true = segment_from_range_gamete(&x.upper(), start, end);
    let q2_true = segment_from_range_gamete(&x.lower(), start, end);
    let mut q = (&q1_true, &q2_true);
    let mut used1_true = vec![false; q1_true.len()];
    let mut used2_true = vec![false; q2_true.len()];
    let mut used1 = &mut used1_true;
    let mut used2 = &mut used2_true;

    let mut i = 0;
    let mut j = 0;

    let mut out: Vec<Segment<SingleChromGamete>> =
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
            out.push(c1.join(&c2));
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
) -> Vec<Segment<SingleChromGamete>> {
    let mut out: Vec<Segment<SingleChromGamete>> = Vec::with_capacity(gx.get_n_loci(0));
    let mut indexes = (i..j).filter(|k| gx.index(*k) == Allele::O);
    if let Some(s) = indexes.next() {
        let mut current_segment = Segment {
            s,
            e: s,
            g: gx.to_owned(),
        };
        let mut last = s;
        for k in indexes {
            if k > last + 1 {
                out.push(current_segment.clone());
                current_segment = Segment {
                    s: k,
                    e: k,
                    g: gx.to_owned(),
                };
            }
            current_segment.e = k;
            last = k;
        }
        out.push(current_segment.clone());
    }
    // TODO: resize out to exact capacity //
    out
}

/// Computes a dp array of the min number of crossings to create a genotype x' that can produce a
/// gamete up to index j from index i where i=0 or i-1 is homozygous 0.
fn _zigzag_costs_r(x: &SingleChromGenotype) -> Vec<usize> {
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

fn _zigzag_switch_points(x: &SingleChromGenotype) -> Vec<bool> {
    let n_loci = x.get_n_loci(0);
    let mut out: Vec<bool> = vec![false; n_loci]; // n_loci 0 would cause n_loci - 1 to overflow
    let _f = |c: bool, i: usize| x.get(c, i).unwrap();
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

fn _consodiate_zigzag(_state: &State<WGe>, x: WGe, i: usize, j: usize) -> State<WGe> {
    // assert that i..j is a zigzag
    let dp = _zigzag_costs_r(&x.genotype());
    assert!((i..j).all(|k| dp[k] > 0));
    // assert! zigzag can't be consolidated immediately

    // extract segments and assume that segments are minimal
    // all crossings in the first phase have the same start chromosome
    let mut hash_map: HashMap<SingleChromGamete, WGa> = HashMap::new();

    let segments_s0: Vec<Segment<WGa>> = segments_from_range_genotype(&x.genotype(), i, j)
        .iter()
        .map(|c| Segment {
            s: c.s,
            e: c.e,
            g: hash_map
                .entry(c.g.clone())
                .or_insert_with(|| WGa::new_from_genotype(c.g.clone(), x.clone()))
                .clone(),
        })
        .collect();
    // Extract min covering subset
    // Join segments
    assert!({
        // all adjacent pairs of segments are joinable
        (0..segments_s0.len() - 1).all(|i| {
            let Segment {
                s: s_x,
                e: e_x,
                g: _,
            } = segments_s0[i];
            let Segment {
                s: s_y,
                e: e_y,
                g: _,
            } = segments_s0[i + 1];
            s_x < s_y && s_y <= e_x && e_x < e_y
        })
    });
    fn join_segments(segments_s0: &Vec<Segment<WGa>>, i: usize, j: usize) -> Segment<WGa> {
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
            res1.join(&res2)
        }
    }
    // extract genotype from the final segment
    unimplemented!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn breeding_program_test() {
        let res = breeding_program(
            3,
            vec![SingleChromGenotype::from_str("010", "101")],
            SingleChromGenotype::ideotype(3),
        );
        assert_ne!(None, res);
        let res_val = res.unwrap();
        assert_eq!(2, res_val.generations());
        assert_eq!(2, res_val.crossings());
    }

    #[test]
    fn segment_gametes_test() {
        let gx = SingleChromGamete::from_str("000111011001110");
        let res = segment_from_range_gamete(&gx, 0, gx.get_n_loci(0));
        assert_eq!(
            vec![(3, 5), (7, 8), (11, 13)],
            res.iter()
                .map(|c| (c.s, c.e))
                .collect::<Vec<(usize, usize)>>()
        );
        assert!(res.iter().all(|c| c.g == gx));
    }

    #[test]
    fn segments_basic_test() {
        let x = SingleChromGenotype::from_str("010", "101");
        let res = segments_from_range_genotype(&x, 0, 3);
        dbg!(&res);
        assert_eq!(res.len(), 2);
    }

    #[test]
    fn breeding_program_zigzag_test() {
        let res = breeding_program(
            4,
            vec![
                SingleChromGenotype::from_str("1010", "1010"),
                SingleChromGenotype::from_str("0101", "0101"),
            ],
            SingleChromGenotype::ideotype(4),
        );
        assert_ne!(None, res);
        let sol = res
            .expect("feasible instance returned None")
            .to_base_sol(4, Objective::Crossings);
        assert_eq!(3, sol.objective);
    }

    #[test]
    fn successors_zigzag_test() {
        let pop_0 = vec![
            WGenS2::new(SingleChromGenotype::from_str("1010", "1010")),
            WGenS2::new(SingleChromGenotype::from_str("0101", "0101")),
        ];
        let state = State::new(pop_0);
        let new_states: Vec<(State<WGenS2<SingleChromGenotype, SingleChromGamete>>, usize)> =
            successors_single_node_extensions(&state)
                .into_iter()
                .collect();
        assert_eq!(
            vec![&SingleChromGenotype::from_str("1010", "0101")],
            new_states
                .iter()
                .map(|wx| wx.0.iter().next().unwrap().genotype())
                .collect::<Vec<_>>()
        );
        assert_eq!(1, new_states.len());
    }

    #[test]
    fn zigzag_cost_test() {
        assert_eq!(
            _zigzag_costs_r(&SingleChromGenotype::from_str(
                "0010110010001110010",
                "0101011100110011010"
            )),
            vec![0, 1, 1, 2, 2, 2, 3, 3, 3, 0, 1, 1, 1, 1, 1, 2, 0, 1, 0]
        );
    }
}
