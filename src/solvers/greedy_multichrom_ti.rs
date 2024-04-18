use crate::abstract_plants::*;

fn breeding_program<A, B, K>(_pop_0: Vec<A>, _ideotype: A) -> Option<WGen<A, B>>
where
    A: Genotype<B>,
    B: Gamete<A>,
{
    // Create a vector of vectors of segments
    // Each tree is going to be created individually bottom up before being joined together
    None
}
