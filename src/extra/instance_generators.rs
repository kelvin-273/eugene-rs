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
