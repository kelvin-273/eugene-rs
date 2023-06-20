use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use rand::prelude::*;
use std::rc::Rc;

type A = SingleChromGenotype;
type B = SingleChromGamete;
type S = (usize, usize, Rc<WGam<A, B>>);

pub struct WSegS<S> {
    segment: S,
    history: Option<(Rc<WSegS<S>>, Rc<WSegS<S>>)>,
}

pub fn segment_tree(wx: WGenS<A, B>) -> Option<WSegS<S>> where
{
    unimplemented!();
}

pub fn random_crossing_schedule_binary_tree_depth(depth: usize) -> WGenS<A, B> {
    unimplemented!();
}

pub fn random_crossing_schedule_binary_tree_nodes(n_loci: usize, nodes: usize) -> WGenS<A, B> {
    // if the number of nodes is makes a valid binary tree attached to the ideotype
    if todo!() {
        return todo!();
    }
    let mut v: Vec<A> = Vec::with_capacity(nodes + 1);
    let mut p: Vec<Option<(usize, usize)>> = vec![None; nodes + 1];

    v[0] = A::ideotype(n_loci);
    if nodes == 0 {
        return todo!();
    }
    p[0] = Some((1, 1));

    fn aux(
        v: &mut Vec<SingleChromGenotype>,
        p: &mut Vec<Option<(usize, usize)>>,
        i: usize,
        j: usize,
    ) {
        let n = j - i;
        if n == 1 {}
    }

    let mut rng = thread_rng();
    // make a node for the ideotype
    // split the remaining nodes between the left and right
}

fn random_parent(gx: &B) -> A {
    unimplemented!();
}
