use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solvers::greedy_base;
use crate::extra::analysis;
use std::rc::Rc;

fn breeding_program(
    n_loci: usize,
    pop_0: Vec<SingleChromGenotype>,
) -> Option<WGen<SingleChromGenotype, SingleChromGamete>> {
    let initial_solution = greedy_base::breeding_program::<
        SingleChromGenotype,
        SingleChromGamete,
        CrosspointBitVec,
        SegmentBitVec,
    >(n_loci, pop_0.clone(), SingleChromGenotype::ideotype(n_loci))?;
    
    let base_objective = analysis::crossings(Rc::new(initial_solution.clone()).extract_first().as_ref());
    Some(initial_solution)
}
