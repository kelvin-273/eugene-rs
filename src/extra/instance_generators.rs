use crate::plants::bit_array;

/// Generates a population of `SingleChromGenotype`s from an array
pub fn distarray_to_homo(v: &[usize]) -> Vec<bit_array::SingleChromGenotype> {
    dbg!(v);
    let n_loci = v.len();
    let n_pop0 = v.iter().max().map(|&length| length + 1).unwrap_or(0);
    let mut xs = vec![vec![(false, false); n_loci]; n_pop0];
    for (i, x) in v.iter().enumerate() {
        xs[*x][i] = (true, true);
    }
    dbg!(&xs);
    xs.iter()
        .map(|x| bit_array::SingleChromGenotype::new(x.clone()))
        .collect()
}

pub fn read_arrays() -> std::io::Result<Vec<(usize, Vec<bit_array::SingleChromGenotype>)>> {
    use std::io;
    let mut out = Vec::new();
    let mut lines = io::stdin().lines();
    while let Some(Ok(s)) = lines.next() {
        out.push(parse_homozygous(s.as_str()));
    }
    Ok(out)
}

pub fn parse_homozygous(s: &str) -> (usize, Vec<bit_array::SingleChromGenotype>) {
    let mut v = Vec::new();
    for w in s.split_whitespace() {
        let x: usize = w.parse().unwrap();
        v.push(x);
    }
    (v.len(), distarray_to_homo(&v))
}
