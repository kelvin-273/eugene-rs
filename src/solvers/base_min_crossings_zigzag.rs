#[derive(Debug, Copy, Clone)]
struct Segment<T> {
    start: usize,
    end: usize,
    data: T,
}

fn min_covering_segments<T, I>(n_loci: usize, segments: I) -> Vec<Segment<T>> where
    T: Copy,
    I: IntoIterator<Item = Segment<T>>
{
    let mut segment_pigeonholes: Vec<Option<Segment<T>>> = (0..n_loci).map(|_| None).collect();
    for c in segments {
        let s = c.start;
        let e = c.end;
        if segment_pigeonholes[s].is_none() {
            segment_pigeonholes[s] = Some(c);
        } else {
            segment_pigeonholes[s].as_mut().map(|c_old| {
                if c_old.end < e {
                    *c_old = c
                } else {
                }
            });
        }
    }
    segment_pigeonholes.iter().filter_map(|c| *c).collect()
}

fn min_cover_segments_from_groups<T>(
    n_loci: usize,
    semgent_groups: Vec<Vec<Segment<T>>>,
) -> Vec<Vec<Segment<T>>> {
    unimplemented!();
}

pub fn breeding_program<A, B>(n_loci: usize, pop_0: Vec<A>, ideotype: A) {
    unimplemented!();
}
