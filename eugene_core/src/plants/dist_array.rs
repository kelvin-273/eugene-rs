use std::ops::{Deref, DerefMut};

// Description: A simple array-like structure that describes which gamete owns which locus.
// Each locus is owned by exactly one gamete, and each gamete can own multiple loci.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct DistArray {
    data: Vec<usize>,
}

impl DistArray {
    pub fn new(xs: &[usize]) -> Self {
        Self {
            data: xs.to_vec(),
        }
    }

    #[inline]
    pub fn n_loci(&self) -> usize {
        self.data.len()
    }

    pub fn n_pop(&self) -> usize {
        // assert that the values in data are in 0..n_pop
        assert!((0..self.n_loci()).all(|i| self[i] == 0 || self.contains(&(self[i] - 1))));
        self.data.iter().max().map_or(0, |&max| max + 1)
    }
}

impl Deref for DistArray {
    type Target = [usize];

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl DerefMut for DistArray {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

impl From<Vec<usize>> for DistArray {
    fn from(v: Vec<usize>) -> Self {
        Self { data: v }
    }
}

// Generate all possible redistributions of the DistArray xs
// Returns a vector of tuples (DistArray, gx, gy) where gx and gy are the two gametes that were
// merged to create the new DistArray. Assumes that xs has at least 2 gametes.
pub fn generate_redistributions(xs: &DistArray) -> Vec<(DistArray, usize, usize)> {
    let n_loci = xs.n_loci();
    let n_pop = xs.n_pop();
    let mut zs = xs.clone();
    if n_pop == 1 {
        return vec![];
    }
    let mut available_gz_values = (n_pop - 2..n_loci).collect::<Vec<_>>();
    let mut out = Vec::new();
    // Backtracking for raw redistributions
    // Assumes that gx and gy are already fixed
    fn bt(
        xs: &DistArray,
        out: &mut Vec<(DistArray, usize, usize)>,
        gx: usize,
        gy: usize,
        zs: &mut DistArray,
        i: usize,
        available_gz_values: &[usize],
        j_max: usize,
    ) {
        if i >= xs.len() {
            out.push((zs.clone(), gx, gy));
        } else if xs[i] != gx && xs[i] != gy {
            bt(xs, out, gx, gy, zs, i + 1, available_gz_values, j_max);
        } else {
            for (j, &gz) in available_gz_values.iter().take(j_max + 1).enumerate() {
                // set value
                zs[i] = gz;
                // check if set value would require more than 1 crossover
                if requires_multipoint(zs, xs, i) {
                    continue;
                }
                // TODO: check if set value would create a dominated gamete
                // Backtrack
                if j == j_max {
                    bt(xs, out, gx, gy, zs, i + 1, available_gz_values, j_max + 1);
                } else {
                    bt(xs, out, gx, gy, zs, i + 1, available_gz_values, j_max);
                }
            }
        }
    }
    for gx in 0..n_pop - 1 {
        for gy in gx + 1..n_pop {
            zs[..n_loci].copy_from_slice(&xs[..n_loci]);
            available_gz_values[0] = gx;
            available_gz_values[1] = gy;
            bt(xs, &mut out, gx, gy, &mut zs, 0, &available_gz_values, 0);
        }
    }
    if out.is_empty() {
        dbg!("No redistributions found for {:?}", xs);
    }
    out
}

pub fn requires_multipoint(zs: &DistArray, xs: &DistArray, i: usize) -> bool {
    let g = zs[i];
    let mut swaps = -1;
    let mut prev_value = None;
    for (&gz, &gx) in zs.iter().take(i + 1).zip(xs.iter()) {
        if gz == g && Some(gx) != prev_value {
            prev_value = Some(gx);
            swaps += 1;
        }
        if swaps > 1 {
            return true;
        }
    }
    false
}

pub fn simplify_dist_array(xs: &DistArray) -> DistArray {
    let mut mapping = vec![None; xs.n_pop()];
    let mut next = 0;
    let mut new_data = Vec::with_capacity(xs.n_loci());
    for &x in xs.iter() {
        if mapping[x].is_none() {
            mapping[x] = Some(next);
            next += 1;
        }
        new_data.push(mapping[x].unwrap());
    }
    DistArray::new(&new_data)
}

pub fn first_full_join(xs: &DistArray) -> Option<(DistArray, usize, usize)> {
    let n_loci = xs.n_loci();
    let n_pop = xs.n_pop();

    let mut start_points = vec![n_loci; n_pop];
    let mut end_points = vec![0; n_pop];
    for (i, x) in xs.iter().enumerate() {
        end_points[*x] = i;
    }
    for (i, x) in xs.iter().enumerate().rev() {
        start_points[*x] = i;
    }

    // find the (gx,gy) with the smallest start_points[gx] such that
    // a full join occurs between gx and gy
    for (gy, sy) in start_points.iter().enumerate() {
        if *sy > 0 && end_points[xs[sy - 1]] == sy - 1 {
            let gx = xs[sy - 1];
            let mut out = xs.clone();
            for i in *sy..end_points[gy] + 1 {
                if xs[i] == gy {
                    out[i] = gx;
                }
            }
            return Some((out, gx, gy));
        }
    }
    None
}

pub fn random_redistribution(xs: &DistArray, rng: &mut impl rand::Rng) -> (DistArray, usize, usize) {
    let n_pop = xs.n_pop();
    assert!(n_pop >= 2);
    let (gx, gy) = loop {
        let gx = rng.gen_range(0..n_pop);
        let gy = rng.gen_range(0..n_pop);
        if gx != gy {
            break (gx.min(gy), gx.max(gy));
        }
    };
    let mut zs = xs.clone();
    let mut available_gz_values = (n_pop - 2..xs.n_loci()).collect::<Vec<_>>();
    available_gz_values[0] = gx;
    available_gz_values[1] = gy;
    let mut j_max = 0;
    for i in 0..xs.n_loci() {
        if xs[i] != gx && xs[i] != gy {
            continue;
        }
        let mut j = rng.gen_range(0..=j_max);
        zs[i] = available_gz_values[j];
        while requires_multipoint(&zs, xs, i) {
            j = rng.gen_range(0..=j_max);
            zs[i] = available_gz_values[j];
        }
        if j == j_max {
            j_max += 1;
        }
    }
    (zs, gx, gy)
}

pub struct DistArrayGenerator {
    n_loci: usize,
    xs: DistArray,
    xs_max: Vec<usize>,
    initialised: bool,
}

impl DistArrayGenerator {
    pub fn new(n_loci: usize) -> Self {
        Self {
            n_loci,
            xs: DistArray::new(&vec![0; n_loci]),
            xs_max: vec![0; n_loci],
            initialised: false,
        }
    }
}

impl Iterator for DistArrayGenerator {
    type Item = DistArray;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.initialised {
            for i in 0..self.n_loci {
                self.xs[i] = i % 2;
                self.xs_max[i] = if i == 0 { 0 } else { 1 };
            }
            self.initialised = true;
            return Some(self.xs.clone())
        }
        for i in (0..self.n_loci).rev() {
            if i > 0 && self.xs[i] <= self.xs_max[i-1] {
                self.xs[i] += 1 + (self.xs[i] + 1 == self.xs[i-1]) as usize;
                self.xs_max[i] = self.xs_max[i-1].max(self.xs[i]);
                for j in i + 1..self.n_loci {
                    self.xs_max[j] = self.xs_max[i];
                    self.xs[j] = if self.xs[j-1] == 0 { 1 } else { 0 };
                }
                return Some(self.xs.clone());
            }
            else if i > 0 && self.xs_max[i] == self.xs_max[i-1] + 1 {
                continue;
            } else if i > 0 {
                panic!("Either xs_max[i] < xs_max[i-1] - 1 or xs_max[i] > xs_max[i-1] + 1")
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dist_array() {
        let d = DistArray::new(&[0, 1, 0, 2, 1]);
        assert_eq!(d.n_loci(), 5);
        assert_eq!(d.n_pop(), 3);
        assert_eq!(d.data, vec![0, 1, 0, 2, 1]);
    }

    #[test]
    fn test_generate_redistributions() {
        let d = DistArray::new(&[0, 1, 0, 2, 1]);
        let r = generate_redistributions(&d);
        for (dist, gx, gy) in r {
            println!("Dist: {:?}, gx: {}, gy: {}", dist, gx, gy);
        }
    }

    #[test]
    fn test_generate_distribute_arrays() {
        let n_loci = 4;
        let mut distarray_generator = DistArrayGenerator::new(n_loci);
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 0, 1])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 0, 2])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 2, 0])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 2, 1])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 2, 3])));
        assert_eq!(distarray_generator.next(), None);

        let n_loci = 5;
        let mut distarray_generator = DistArrayGenerator::new(n_loci);
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 0, 1, 0])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 0, 1, 2])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 0, 2, 0])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 0, 2, 1])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 0, 2, 3])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 2, 0, 1])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 2, 0, 2])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 2, 0, 3])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 2, 1, 0])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 2, 1, 2])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 2, 1, 3])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 2, 3, 0])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 2, 3, 1])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 2, 3, 2])));
        assert_eq!(distarray_generator.next(), Some(DistArray::new(&[0, 1, 2, 3, 4])));
        assert_eq!(distarray_generator.next(), None);
    }
}
