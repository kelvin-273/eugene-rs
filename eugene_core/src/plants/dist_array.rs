//! A structure to represent the distribution of loci among gametes.
//!
//! This module provides the `DistArray` struct, which is an array-like structure that
//! describes which gamete owns which locus. Each locus is owned by exactly one gamete,
//! and each gamete can own multiple loci. The module includes methods for manipulating
//! and analyzing `DistArray` instances, such as finding successors and predecessors
//! in lexicographic order, generating redistributions, and checking for canonical form.
//! It also includes utility functions for counting distinct distributions and computing
//! ordinal indices.
//!
//! # Examples
//! ```
//! use dist_array::dist_array;
//! let xs = dist_array![0, 1, 0, 1, 2];
//! let succ = xs.successor().unwrap();
//! assert_eq!(succ, dist_array![0, 1, 2, 0, 1]);
//! ```

use rand::Rng;
use serde::{Deserialize, Serialize};
use std::iter::FromIterator;
use std::ops::{Deref, DerefMut};

// Description: A simple array-like structure that describes which gamete owns which locus.
// Each locus is owned by exactly one gamete, and each gamete can own multiple loci.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Deserialize, Serialize)]
pub struct DistArray(Vec<usize>);

impl DistArray {
    /// Returns the number of loci. Is equivalent to `self.len()`.
    ///
    /// # Returns
    /// - `usize`: The number of loci in the DistArray.
    ///
    /// # Examples
    /// ```
    /// use dist_array::dist_array;
    /// let xs = dist_array![0, 1, 0, 1, 2];
    /// assert_eq!(xs.n_loci(), 5);
    /// ```
    #[inline]
    pub fn n_loci(&self) -> usize {
        self.len()
    }

    /// Returns the number of gametes in the population implied by the DistArray.
    pub fn n_pop(&self) -> usize {
        // assert that the values in data are in 0..n_pop
        assert!((0..self.n_loci()).all(|i| self[i] == 0 || self.contains(&(self[i] - 1))));
        self.0.iter().max().map_or(0, |&max| max + 1)
    }

    /// Checks if the current DistArray has `other` as its successor in lexicographic order.
    ///
    /// # Arguments
    /// - `other`: The DistArray to check as a successor.
    /// 
    /// # Returns
    /// - `bool`: True if `other` is the successor of the current DistArray, otherwise false.
    ///
    /// # Examples
    /// ```
    /// use dist_array::dist_array;
    /// let xs = dist_array![0, 1, 0, 1, 2];
    /// let succ = dist_array![0, 1, 2, 0, 1];
    /// assert!(xs.has_successor(&succ));
    /// ```
    pub fn has_successor(&self, other: &Self) -> bool {
        let n = self.len();
        if other.len() < n || other.len() > n + 1 {
            false
        } else if other.len() == n + 1 {
            *self == (0..n).collect::<DistArray>().into()
                && *other == (0..n + 1).map(|x| x % 2).collect::<DistArray>().into()
        } else {
            let mut i = 0;
            while i < n && self[i] == other[i] {
                i += 1;
            }
            if i == n || other[i] != self[i] + 1 {
                return false;
            }
            let mut a = 0;
            i += 1;
            while i < n && other[i] == a {
                a = (a + 1) % 2;
                i += 1;
            }
            return i == n;
        }
    }

    /// Returns the successor DistArray in lexicographic order.
    ///
    /// # Returns
    /// - `Option<DistArray>`: The successor DistArray if it exists, otherwise None.
    ///
    /// # Panics
    /// - Panics if the DistArray is not in canonical form.
    ///
    /// # Examples
    /// ```
    /// use dist_array::dist_array;
    /// let xs = dist_array![0, 1, 0, 1, 2];
    /// let succ = xs.successor().unwrap();
    /// assert_eq!(succ, dist_array![0, 1, 2, 0, 1]);
    /// ```
    pub fn successor(&self) -> Option<Self>
    where
        Self: Sized,
    {
        assert!(is_canonical_dist_array(self));
        let n = self.len();
        let mut xs = self.0.clone();

        let mut xs_max = xs.clone();
        for i in 1..n {
            xs_max[i] = xs_max[i - 1].max(xs[i]);
        }

        for i in (0..n).rev() {
            if i > 0 && xs[i] <= xs_max[i - 1] {
                xs[i] += 1 + (xs[i] + 1 == xs[i - 1]) as usize;
                for j in i + 1..n {
                    xs[j] = if xs[j - 1] == 0 { 1 } else { 0 };
                }
                return Some(DistArray::from(xs));
            } else if i > 0 && xs_max[i] == xs_max[i - 1] + 1 {
                continue;
            } else if i > 0 {
                panic!("Either xs_max[i] < xs_max[i-1] - 1 or xs_max[i] > xs_max[i-1] + 1")
            }
        }
        Some((0..n + 1).map(|x| x % 2).collect::<DistArray>().into())
    }

    /// Returns the successor DistArray in lexicographic order.
    ///
    /// # Returns
    /// - `Option<DistArray>`: The successor DistArray if it exists, otherwise None.
    ///
    /// # Panics
    /// - Panics if the DistArray is not in canonical form.
    ///
    /// # Examples
    /// ```
    /// use dist_array::dist_array;
    /// let xs = dist_array![0, 1, 2, 0, 1];
    /// let pred = xs.predecessor();
    /// assert_eq!(pred, Some(dist_array![0, 1, 0, 1, 2]));
    /// let ys = dist_array![0];
    /// assert_eq!(ys.predecessor(), None)
    /// ```
    pub fn predecessor(&self) -> Option<Self>
    where
        Self: Sized,
    {
        assert!(is_canonical_dist_array(self));
        let n = self.len();
        if n == 1 {
            return None;
        }

        let mut xs = self.clone();

        for i in (2..n).rev() {
            if xs[i] > 0 && xs[i] != xs[i-1] + 1 {
                xs[i] -= 1;
                for j in 0..n-i {
                   xs[i + 1 + j] = j & 1;
                }
                return Some(xs);
            } else {
                todo!("other cases")
            }
        }
        Some((0..n-1).collect::<DistArray>().into())
    }
}

/// Checks if a DistArray is in canonical form.
///
/// A DistArray is in canonical form if:
/// - The first element is 0.
/// - Each subsequent element is at most 1 greater than the maximum of all previous elements.
/// - No two consecutive elements are equal.
///
/// # Arguments
/// - `xs`: The DistArray to check.
///
/// # Returns
/// - `bool`: True if the DistArray is in canonical form, otherwise false.
///
/// # Examples
/// ```
/// use dist_array::{dist_array, is_canonical_dist_array};
/// let xs = dist_array![0, 1, 0, 2, 1];
/// assert!(is_canonical_dist_array(&xs));
/// let ys = dist_array![0, 2, 1, 0];
/// assert!(!is_canonical_dist_array(&ys));
/// ```
pub fn is_canonical_dist_array(xs: &DistArray) -> bool {
    let n = xs.len();
    if n == 0  || xs[0] != 0 {
        return false;
    }
    let mut x_max = 0;
    for i in 1..n {
        if xs[i] > x_max + 1 || xs[i] == xs[i-1] {
            return false
        }
        x_max = x_max.max(xs[i])
    }
    true
}

impl Deref for DistArray {
    type Target = [usize];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for DistArray {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl From<Vec<usize>> for DistArray {
    fn from(v: Vec<usize>) -> Self {
        Self(v)
    }
}

impl Into<Vec<usize>> for DistArray {
    fn into(self) -> Vec<usize> {
        self.0
    }
}

impl FromIterator<usize> for DistArray {
    fn from_iter<I: IntoIterator<Item = usize>>(iter: I) -> Self {
        let data: Vec<usize> = iter.into_iter().collect();
        Self(data)
    }
}

impl<'a> FromIterator<&'a usize> for DistArray {
    fn from_iter<I: IntoIterator<Item = &'a usize>>(iter: I) -> Self {
        let data: Vec<usize> = iter.into_iter().cloned().collect();
        Self(data)
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
    DistArray(new_data)
}

pub fn canonical_dist_array(xs: &DistArray) -> DistArray {
    let mut zs = xs.clone();
    let n_loci = zs.len();
    let n_pop = xs.0.iter().max().unwrap() + 1;
    let mut mapping = vec![None; n_pop];
    mapping[zs[0]] = Some(0);
    let mut x_max = 0;
    for i in 0..n_loci {
        zs[i] = *mapping[zs[i]].get_or_insert_with(|| {
            x_max += 1;
            x_max
        });
    }
    zs
}

#[derive(Clone, Debug)]
pub struct OrderedDistArray(pub DistArray);

impl PartialEq for OrderedDistArray {
    fn eq(&self, other: &Self) -> bool {
        self.0.0 == other.0.0
    }
}

impl Eq for OrderedDistArray {}

impl std::hash::Hash for OrderedDistArray {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        for &x in &self.0 .0 {
            x.hash(state);
        }
    }
}

impl std::ops::Deref for OrderedDistArray {
    type Target = DistArray;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<DistArray> for OrderedDistArray {
    fn from(d: DistArray) -> Self {
        OrderedDistArray(d)
    }
}

impl From<OrderedDistArray> for DistArray {
    fn from(od: OrderedDistArray) -> Self {
        od.0
    }
}

impl PartialOrd for OrderedDistArray {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OrderedDistArray {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.n_loci() < other.n_loci() {
            return std::cmp::Ordering::Less;
        } else if self.n_loci() > other.n_loci() {
            return std::cmp::Ordering::Greater;
        }
        let n = self.n_loci().min(other.n_loci());
        for i in 0..n {
            if self[i] < other[i] {
                return std::cmp::Ordering::Less;
            } else if self[i] > other[i] {
                return std::cmp::Ordering::Greater;
            }
        }
        std::cmp::Ordering::Equal
    }
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

pub fn random_redistribution(xs: &DistArray, rng: &mut impl Rng) -> (DistArray, usize, usize) {
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
            xs: DistArray(vec![0; n_loci]),
            xs_max: vec![0; n_loci],
            initialised: false,
        }
    }

    pub fn from_dist_array(xs: &DistArray) -> Self {
        let n_loci = xs.n_loci();
        Self {
            n_loci,
            xs: xs.clone(),
            xs_max: vec![0; n_loci],
            initialised: false,
        }
    }
}

impl Default for DistArrayGenerator {
    fn default() -> Self {
        Self::new(0)
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
            return Some(self.xs.clone());
        }
        for i in (0..self.n_loci).rev() {
            if i > 0 && self.xs[i] <= self.xs_max[i - 1] {
                self.xs[i] += 1 + (self.xs[i] + 1 == self.xs[i - 1]) as usize;
                self.xs_max[i] = self.xs_max[i - 1].max(self.xs[i]);
                for j in i + 1..self.n_loci {
                    self.xs_max[j] = self.xs_max[i];
                    self.xs[j] = if self.xs[j - 1] == 0 { 1 } else { 0 };
                }
                return Some(self.xs.clone());
            } else if i > 0 && self.xs_max[i] == self.xs_max[i - 1] + 1 {
                continue;
            } else if i > 0 {
                panic!("Either xs_max[i] < xs_max[i-1] - 1 or xs_max[i] > xs_max[i-1] + 1")
            }
        }
        None
    }
}

/// Counts the number of distinct DistArrays with a given number of loci.
///
/// # Arguments
/// - `n_loci`: The number of loci.
///
/// # Returns
/// - `usize`: The number of distinct DistArrays with the given number of loci.
///
/// # Panics
/// - Panics if `n_loci` is greater than or equal to 26, as this would cause overflow.
///
/// # Examples
/// ```
/// use dist_array::count_distribute_instances;
/// let count = dist_array::count_distribute_instances(1);
/// assert_eq!(count, 1);
/// let count = dist_array::count_distribute_instances(3);
/// assert_eq!(count, 2);
/// let count = dist_array::count_distribute_instances(5);
/// assert_eq!(count, 15);
/// ```
pub fn count_distribute_instances(n_loci: usize) -> usize {
    if n_loci < 26 {
        H_LOOKUP[n_loci - 1][0]
    } else {
        panic!("n_loci too large, will cause overflow");
    }
}

/// Computes the ordinal index of a DistArray in the lexicographic ordering of all possible
/// DistArrays with the same number of loci.
///
/// # Arguments
/// - `xs`: The DistArray to compute the ordinal index for.
///
/// # Returns
/// - `usize`: The ordinal index of the DistArray.
///
/// # Examples
/// ```
/// use dist_array::dist_array;
/// let xs = dist_array![0, 1, 0, 2, 1];
/// let index = dist_array::ordinal_index(&xs);
/// assert_eq!(index, 4);
/// ```
pub fn ordinal_index(xs: &DistArray) -> usize {
    let n = xs.n_loci();
    let mut total = 1;
    assert!(xs[0] == 0);
    let mut x_max = 0;
    for i in 1..n {
        let x_next = x_max.max(xs[i]);
        let m = xs[i] - (xs[i - 1] < xs[i]) as usize;
        total += m * H_LOOKUP[n - i - 1][x_max];
        x_max = x_next;
    }
    total
}

const H_LOOKUP: [[usize; 26]; 26] = [
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 0], [2, 5, 10, 17, 26, 37, 50, 65, 82, 101, 122, 145, 170, 197, 226, 257, 290, 325, 362, 401, 442, 485, 530, 577, 0, 0], [5, 15, 37, 77, 141, 235, 365, 537, 757, 1031, 1365, 1765, 2237, 2787, 3421, 4145, 4965, 5887, 6917, 8061, 9325, 10715, 12237, 0, 0, 0], [15, 52, 151, 372, 799, 1540, 2727, 4516, 7087, 10644, 15415, 21652, 29631, 39652, 52039, 67140, 85327, 106996, 132567, 162484, 197215, 237252, 0, 0, 0, 0], [52, 203, 674, 1915, 4736, 10427, 20878, 38699, 67340, 111211, 175802, 267803, 395224, 567515, 795686, 1092427, 1472228, 1951499, 2548690, 3284411, 4181552, 0, 0, 0, 0, 0], [203, 877, 3263, 10481, 29371, 73013, 163967, 338233, 649931, 1176701, 2025823, 3341057, 5310203, 8173381, 12232031, 17858633, 25507147, 35724173, 49160831, 66585361, 0, 0, 0, 0, 0, 0], [877, 4140, 17007, 60814, 190497, 529032, 1322035, 3017562, 6376149, 12616132, 23599287, 42061830, 71895817, 118485984, 189107067, 293386642, 443838525, 656471772, 951480319, 0, 0, 0, 0, 0, 0, 0], [4140, 21147, 94828, 372939, 1291020, 3967195, 10949772, 27499083, 63625324, 137144475, 278054700, 534575947, 981235788, 1729424859, 2940885580, 4844638155, 7757888172, 12111500443, 0, 0, 0, 0, 0, 0, 0, 0], [21147, 115975, 562595, 2409837, 9131275, 30785747, 93197715, 256118905, 646147067, 1512354975, 3315122947, 6861571205, 13504254315, 25423408747, 46017036275, 80427460497, 136237711195, 0, 0, 0, 0, 0, 0, 0, 0, 0], [115975, 678570, 3535027, 16360786, 67310847, 247126450, 815305195, 2438979402, 6681531511, 16926317722, 40012800675, 88981537570, 187474460527, 376521349986, 724665968347, 1342649618650, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [678570, 4213597, 23430840, 116393205, 516369838, 2050937445, 7330810572, 23754387325, 70378569810, 192349660173, 489109544320, 1166271373797, 2626214876310, 5619443518165, 11487973175508, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [4213597, 27644437, 163254885, 865549453, 4116416797, 17585497797, 67739250757, 236659281085, 755378218653, 2220256485877, 6057366816997, 15455199988077, 37134022033885, 84540738911653, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [27644437, 190899322, 1192059223, 6713065156, 34051164985, 155666739742, 643094785627, 2411993186248, 8263282235101, 26039675189890, 76028868158047, 207141221902732, 530149003318273, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [190899322, 1382958545, 9097183602, 54190360453, 291871399682, 1421428484337, 6270561900010, 25147234538837, 92145933070698, 310385944867057, 967429903483202, 2808702444248325, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1382958545, 10480142147, 72384727657, 454442481041, 2588914083065, 13377704321695, 62770605938897, 268176574842557, 1047553409432641, 3760903407286715, 12483001479080345, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [10480142147, 82864869804, 599211936355, 3952241526188, 23733360653955, 129659127547372, 644800210475939, 2924789433330540, 12141330682747843, 46331132144660780, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [82864869804, 682076806159, 5150665398898, 35590085232519, 224592570163192, 1293095848212799, 6793590696186174, 32614856716061623, 143461777606643524, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [682076806159, 5832742205057, 45891416030315, 331362825860749, 2191466128865567, 13259069937250169, 73376400893178667, 371765774619074885, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [5832742205057, 51724158235372, 423145657921379, 3185554606447814, 22024934452712437, 139671750579429512, 812024179978146887, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [51724158235372, 474869816156751, 4031845922290572, 31581598272055879, 227771488390279260, 1510382932875294447, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [474869816156751, 4506715738447323, 39645290116637023, 322516283206446897, 2421468886436411487, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [4506715738447323, 44152005855084346, 401806863439720943, 3389017736055752178, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [44152005855084346, 445958869294805289, 4192631462935194064, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [445958869294805289, 4638590332229999353, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [4638590332229999353, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
];

/// Macro to create a DistArray from a list of usize values.
///
/// # Examples
/// ```
/// use dist_array::dist_array;
/// let d = dist_array![0, 1, 0, 2, 1];
/// assert_eq!(d.n_loci(), 5);
/// assert_eq!(d.n_pop(), 3);
/// assert_eq!(d.0, vec![0, 1, 0, 2, 1]);
/// ```
#[macro_export]
macro_rules! dist_array {
    ($($x:expr),*) => {
        DistArray::from(vec![$($x),*])
    };
}
pub use dist_array;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dist_array() {
        let d = DistArray(vec![0, 1, 0, 2, 1]);
        assert_eq!(d.n_loci(), 5);
        assert_eq!(d.n_pop(), 3);
        assert_eq!(d.0, vec![0, 1, 0, 2, 1]);
    }

    #[test]
    fn test_dist_array_macro() {
        let d = dist_array![0, 1, 0, 2, 1];
        assert_eq!(d.n_loci(), 5);
        assert_eq!(d.n_pop(), 3);
        assert_eq!(d.0, vec![0, 1, 0, 2, 1]);
    }

    #[test]
    fn test_generate_redistributions() {
        let d = DistArray(vec![0, 1, 0, 2, 1]);
        let r = generate_redistributions(&d);
        for (dist, gx, gy) in r {
            println!("Dist: {:?}, gx: {}, gy: {}", dist, gx, gy);
        }
    }

    #[test]
    fn test_generate_distribute_arrays() {
        let n_loci = 4;
        let mut xs_gen = DistArrayGenerator::new(n_loci);
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 0, 1])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 0, 2])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 2, 0])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 2, 1])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 2, 3])));
        assert_eq!(xs_gen.next(), None);

        let n_loci = 5;
        let mut xs_gen = DistArrayGenerator::new(n_loci);
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 0, 1, 0])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 0, 1, 2])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 0, 2, 0])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 0, 2, 1])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 0, 2, 3])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 2, 0, 1])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 2, 0, 2])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 2, 0, 3])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 2, 1, 0])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 2, 1, 2])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 2, 1, 3])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 2, 3, 0])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 2, 3, 1])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 2, 3, 2])));
        assert_eq!(xs_gen.next(), Some(DistArray(vec![0, 1, 2, 3, 4])));
        assert_eq!(xs_gen.next(), None);
    }

    #[test]
    fn canonical_dist_array_test() {
        assert_eq!(
            canonical_dist_array(&DistArray(vec![5, 2, 6, 2, 5, 4])),
            DistArray(vec![0, 1, 2, 1, 0, 3])
        );
    }

    #[test]
    fn ordinal_index_test() {
        for n_loci in 1..11 {
            for (i, xs) in DistArrayGenerator::new(n_loci).enumerate() {
                let idx = ordinal_index(&xs);
                assert_eq!(idx, i + 1, "Failed at xs = {:?}", xs);
            }
        }
    }
}
