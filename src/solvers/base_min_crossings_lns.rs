use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solvers::base_min_generations_segment;

pub fn breeding_program(
    n_loci: usize,
    pop_0: &Vec<SingleChromGenotype>,
) -> Option<WGen<SingleChromGenotype, SingleChromGamete>> {
    let _initial_solution = base_min_generations_segment::breeding_program(n_loci, pop_0)?;
    unimplemented!();
}
