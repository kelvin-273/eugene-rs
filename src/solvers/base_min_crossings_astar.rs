use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solvers::base_min_generations_segment;
use crate::solution::*;
use pathfinding::directed::astar::astar;
use pyo3::prelude::*;
use std::collections::hash_map::DefaultHasher;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::hash::Hasher;
use std::rc::Rc;

type WG = WGen<SingleChromGenotype, SingleChromGamete>;

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
            let sol = BaseSolution::min_cross_from_wgen(n_loci, &x_star);
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
) -> Option<WGen<SingleChromGenotype, SingleChromGamete>>
where
{
    // test for feasibility
    // create new population / lift genotypes
    let start: Vec<Rc<WG>> = pop_0.iter().map(|x| Rc::new(x.lift_a())).collect();
    let state_0: State<Rc<WG>> = State::new(start);
    // Initialise priority queue

    let successors = |state: &State<Rc<WG>>| {
        successors_single_node_extensions(state)
    };

    let heuristic = |state: &State<Rc<WG>>| heuristic(n_loci, state);

    let success = |state: &State<Rc<WG>>| {
        match state.head.as_ref() {
            StateData::Start(ref v) => v.iter().any(|wx| wx.genotype == ideotype),
            StateData::Next(x, _tail) => x.genotype == ideotype,
        }
    };

    let (v, _obj) = astar(&state_0, successors, heuristic, success)?;
    let x = v.last()?;
    clone_last_from_state(x.to_owned())
}

fn successors_single_node_extensions(
    state: &State<Rc<WG>>,
) -> impl IntoIterator<Item = (State<Rc<WG>>, usize)> {
    // vector of the genotype indices and rc gametes
    let gametes: Vec<(usize, Rc<WGam<SingleChromGenotype, SingleChromGamete>>)> = state
        .iter()
        .enumerate()
        .flat_map(|(i, x)| {
            non_dominated_gametes(x.clone())
                .into_iter()
                .map(move |gx| (i, gx))
        })
        .collect();
    let first_j = gametes
        .iter()
        .enumerate()
        .find(|(_i_gamete, (i, _gx))| i >= &state.most_recent_parent)
        .map_or(0, |(i, _igx)| i);

    let g_star = gametes
        .iter()
        .filter(|(_i, wg)| wg.gamete.alleles().iter().all(|a| *a == Allele::O))
        .next();
    if let Some((_, wg)) = g_star {
        let x_star = SingleChromGenotype::from_gametes(&wg.gamete, &wg.gamete);
        let mut wz = WGen::new(x_star);
        wz.history = Some((wg.clone(), wg.clone()));
        let new_state = state.push(Rc::new(wz));
        // TODO: halt the astar when this is found <11-06-23> //
        return vec![(new_state, 2)];
    }

    let n_gametes = gametes.len();
    let out = (0..gametes.len())
        .flat_map(|i| ((i + 1).max(first_j)..n_gametes).map(move |j| (i, j)))
        .map(|(i, j)| {
            let (_idx, wg1) = (&gametes[i]).clone();
            let (idy, wg2) = (&gametes[j]).clone();
            let z = SingleChromGenotype::from_gametes(&wg1.gamete, &wg2.gamete);
            let mut wz = WGen::new(z);
            wz.history = Some((wg1, wg2));
            let new_state = state.push_with_most_recent_parent(Rc::new(wz), idy);
            (new_state, 2)
        })
        .collect::<Vec<(State<Rc<WG>>, usize)>>();
    out
}

fn non_dominated_gametes(x: Rc<WG>) -> Vec<Rc<WGam<SingleChromGenotype, SingleChromGamete>>> {
    let _hashset: HashSet<SingleChromGamete> = HashSet::new();
    let res = segments_from_range_genotype(&x.genotype, 0, x.genotype.get_n_loci(0));
    res.iter()
        .map(|c| c.g.clone())
        //.collect::<HashSet<SingleChromGamete>>()
        //.iter()
        .map(|g| {
            Rc::new({
                let mut wg = WGam::new(g.clone());
                wg.history.push(x.clone());
                wg
            })
        })
        .collect()
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

fn _consodiate_zigzag(_state: &State<Rc<WG>>, x: Rc<WG>, i: usize, j: usize) -> State<WG> {
    // assert that i..j is a zigzag
    let dp = _zigzag_costs_r(&x.genotype);
    assert!((i..j).all(|k| dp[k] > 0));
    // assert! zigzag can't be consolidated immediately

    // extract segments and assume that segments are minimal
    // all crossings in the first phase have the same start chromosome
    type Val = Rc<WGam<SingleChromGenotype, SingleChromGamete>>;
    let mut hash_map: HashMap<SingleChromGamete, Val> = HashMap::new();

    let segments_s0: Vec<SegmentMC<Val>> = segments_from_range_genotype(&x.genotype, i, j)
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
    assert!({
        // all adjacent pairs of segments are joinable
        (0..segments_s0.len() - 1).all(|i| {
            let SegmentMC {
                s: s_x,
                e: e_x,
                g: _,
            } = segments_s0[i];
            let SegmentMC {
                s: s_y,
                e: e_y,
                g: _,
            } = segments_s0[i + 1];
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

impl<A: Clone, B> Clone for WGen<A, B> {
    fn clone(&self) -> Self {
        WGen {
            genotype: self.genotype.clone(),
            history: self.history.clone(),
        }
    }
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

    pub fn count_nodes(&self) -> usize {
        let head = &self.head;
        fn aux<A>(head: &Link<A>) -> usize {
            match head.as_ref() {
                StateData::Start(_) => 0,
                StateData::Next(_, link) => 1 + aux(link),
            }
        }
        aux(head)
    }

    pub fn hash(&self) -> u64
    where
        A: Hash,
    {
        let mut hasher1: DefaultHasher = Default::default();
        let mut intermediat_v: Vec<u64> = self
            .iter()
            .map(|x| {
                x.hash(&mut hasher1);
                hasher1.finish()
            })
            .collect();
        intermediat_v.sort();

        let mut hasher2: DefaultHasher = Default::default();
        intermediat_v.hash(&mut hasher2);
        hasher2.finish()
    }
}

impl<A> std::hash::Hash for State<A> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        state.write_usize(self.most_recent_parent)
        //state.write_u8(0)
    }
}

impl PartialEq for State<Rc<WG>> {
    fn eq(&self, other: &Self) -> bool {
        //match (self.head.as_ref(), other.head.as_ref()) {
        //    (StateData::Start(ref v1), StateData::Start(ref v2)) => todo!(),
        //    (StateData::Start(_), StateData::Next(_, _)) => false,
        //    (StateData::Next(_, _), StateData::Start(_)) => false,
        //    (StateData::Next(x1, _), StateData::Next(x2, _)) => x1.genotype == x2.genotype,
        //}
        let mut v1: Vec<_> = self.iter().map(|x| x.genotype.clone()).collect();
        let mut v2: Vec<_> = other.iter().map(|x| x.genotype.clone()).collect();
        v1.sort();
        v2.sort();
        v1 == v2
    }
}

//impl<A: PartialEq> Eq for State<A> {}

impl Eq for State<Rc<WG>> {}

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

fn heuristic(
    n_loci: usize,
    state: &State<Rc<WGen<SingleChromGenotype, SingleChromGamete>>>,
) -> usize {
    let n_segments = base_min_generations_segment::min_covering_segments::<
        SingleChromGenotype,
        SingleChromGamete,
        SegmentMC<Rc<WGam<SingleChromGenotype, SingleChromGamete>>>,
    >(
        n_loci,
        &state.iter().map(|wgx| wgx.genotype.clone()).collect(),
    )
    .len();
    ((n_segments + 1) >> 1).max((n_segments as f32).log2().ceil() as usize)
}

fn clone_last_from_state<A: Clone, B>(state: State<Rc<WGen<A, B>>>) -> Option<WGen<A, B>> {
    match state.head.as_ref() {
        StateData::Start(v) => v.last().map(|x| x.as_ref().clone()),
        StateData::Next(x, _n) => Some(x.as_ref().clone()),
    }
}

#[derive(Debug, Clone, Copy)]
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
                let k_z = CrosspointBitVec::new(false, other.s);
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

impl Segment<SingleChromGenotype, SingleChromGamete>
    for SegmentMC<Rc<WGam<SingleChromGenotype, SingleChromGamete>>>
{
}

impl HaploidSegment<SingleChromGenotype, SingleChromGamete>
    for SegmentMC<Rc<WGam<SingleChromGenotype, SingleChromGamete>>>
{
    fn start(&self) -> usize {
        self.s
    }
    fn end(&self) -> usize {
        self.e
    }
    fn gamete(&self) -> Rc<WGam<SingleChromGenotype, SingleChromGamete>> {
        self.g.clone()
    }
    fn from_start_end_gamete(
        s: usize,
        e: usize,
        g: Rc<WGam<SingleChromGenotype, SingleChromGamete>>,
    ) -> Self {
        Self { s, e, g }
    }
    fn join(&self, other: &Self) -> Self {
        assert!(self.s < other.s && other.s <= self.e + 1 && self.e < other.e);
        let z = SingleChromGenotype::from_gametes(&self.g.gamete, &other.g.gamete);
        let mut wz = z.lift_a();
        wz.history = Some((self.g.clone(), other.g.clone()));
        let gz = CrosspointBitVec::new(false, other.s).cross(&z);
        let wgz = WGam::new(gz);
        SegmentMC {
            s: self.s,
            e: other.e,
            g: Rc::new(wgz),
        }
    }
}

fn segments_from_range_genotype(
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
            e: s,
            g: gx.to_owned(),
        };
        let mut last = s;
        for k in indexes {
            if k > last + 1 {
                out.push(current_segment.clone());
                current_segment = SegmentMC {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::extra::analysis;

    #[test]
    fn breeding_program_test() {
        let res = breeding_program(
            3,
            vec![SingleChromGenotype::from_str("010", "101")],
            SingleChromGenotype::ideotype(3),
        );
        assert_ne!(None, res);
        let res_val = Rc::new(res.unwrap()).extract_first();
        assert_eq!(2, analysis::generations(res_val.as_ref()));
        assert_eq!(2, analysis::crossings(res_val.as_ref()));
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

    #[test]
    fn segments_basic_test() {
        let x = SingleChromGenotype::from_str("010", "101");
        let res = segments_from_range_genotype(&x, 0, 3);
        dbg!(&res);
        assert_eq!(res.len(), 2)
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
}
