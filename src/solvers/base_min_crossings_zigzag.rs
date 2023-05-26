#[derive(Debug, Copy, Clone)]
struct Segment<T> {
    start: usize,
    end: usize,
    data: T,
}

fn min_covering_segments<'a, T, I>(n_loci: usize, segments: I) -> Vec<&'a Segment<T>>
where
    T: Copy,
    I: IntoIterator<Item = &'a Segment<T>>,
{
    let mut segment_pigeonholes: Vec<Option<&'a Segment<T>>> = (0..n_loci).map(|_| None).collect();
    for c in segments {
        let s = c.start;
        let e = c.end;
        if segment_pigeonholes[s].is_none() {
            segment_pigeonholes[s] = Some(c);
        } else {
            segment_pigeonholes[s]
                .as_mut()
                .map(|c_old| match c_old.end < e {
                    true => *c_old = c,
                    false => {}
                });
        }
    }
    segment_pigeonholes.iter().filter_map(|c| *c).collect()
}

struct State<'a, T> {
    segments: &'a Vec<Vec<Segment<T>>>,
    selected_groups: Vec<bool>,
}

impl<'a, T> State<'a, T> {
    pub fn new(segment_groups: &'a Vec<Vec<Segment<T>>>) -> Self {
        Self {
            selected_groups: (0..segment_groups.len()).map(|_| true).collect(),
            segments: segment_groups,
        }
    }

    pub fn selected_segments(&self) -> impl IntoIterator<Item = &Segment<T>> {
        self.selected_groups
            .iter()
            .enumerate()
            .filter(|(_, &x)| x)
            .flat_map(|(i, _)| self.segments[i].iter())
    }

    pub fn deselect(&mut self, i: usize) {
        self.selected_groups[i] = false;
    }

    pub fn reselect(&mut self, i: usize) {
        self.selected_groups[i] = true;
    }

    pub fn remaining_groups(&self) -> Vec<usize> {
        self.selected_groups
            .iter()
            .enumerate()
            .filter(|(_, b)| **b)
            .map(|(i, _)| i)
            .collect()
    }

    pub fn n_groups(&self) -> usize {
        self.selected_groups.len()
    }
}

fn min_cover_segments_from_groups<T>(
    n_loci: usize,
    segment_groups: Vec<Vec<Segment<T>>>,
) -> Vec<Vec<Segment<T>>> {
    let n_groups = segment_groups.len();
    let mut state = State::<T> {
        segments: &segment_groups,
        selected_groups: (0..segment_groups.len()).map(|_| true).collect(),
    };

    fn f_obj<T>(state: &State<T>) -> usize {
        unimplemented!()
    };

    fn aux<T>(state: &mut State<T>, current_obj: usize) -> usize {
        let mut out = current_obj;
        for i in state.remaining_groups() {
            state.deselect(i);
            let new_obj = f_obj(&state);
            if new_obj <= current_obj {
                out = out.min(aux(state, new_obj));
            }
            state.reselect(i);
        }
        out
    }

    unimplemented!();
}

pub fn breeding_program<A, B>(n_loci: usize, pop_0: Vec<A>, ideotype: A) {
    unimplemented!();
}
