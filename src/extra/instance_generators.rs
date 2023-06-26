use crate::abstract_plants;
use crate::plants::bit_array;

/// A struct to iterate over Homozygous instances
///
/// We keep an internal representation of the sources of each allele, a stack frame to simulate
/// function calls and a flag to tell if we should stop iterating.
/// The stack should store the index that we are currently up to and the max value we have pushed
/// so far.
struct HomozygousInstances {
    n_loci: usize,
    key: Vec<usize>,
    stack: Vec<usize>,
    i: usize,
    complete: bool,
}

impl HomozygousInstances {
    pub fn new(n_loci: usize) -> Self {
        let mut stack = Vec::with_capacity(n_loci);
        let mut complete = false;

        if n_loci == 0 {
            complete = true;
        } else {
            stack.push(0)
        }

        Self {
            n_loci,
            key: vec![0; n_loci],
            stack,
            i: 1,
            complete,
        }
    }
}

impl Iterator for HomozygousInstances {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.complete {
            return None;
        }
        if self.stack.is_empty() {
            self.complete = true;
            return None;
        }

        // either its the first key,
        while self.i < self.n_loci {
            let c_max = *self.stack.last().unwrap();
            let prev = self.key[self.i - 1];
            let c = if prev == 0 { 1 } else { 0 };
            self.stack.push(c_max.max(c + 1));
            self.i += 1;
        }
        // somewhere in the middle,
        // or we're about to turn off
        None
    }
}

fn array_to_homo(v: &Vec<usize>) -> Vec<bit_array::SingleChromGenotype> {
    dbg!(v);
    let n_loci = v.len();
    let n_pop0 = v.iter().max().map(|&length| length + 1).unwrap_or(0);
    let mut xs = vec![vec![(false, false); n_loci]; n_pop0];
    dbg!(&xs);
    for (i, x) in v.iter().enumerate() {
        xs[*x][i] = (true, true);
    }
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
    (v.len(), array_to_homo(&v))
}
