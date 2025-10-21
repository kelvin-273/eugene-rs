//! A submodule for the different solvers.
//!
//! Each solver has a function `breeding_program` which solves solves instances of various crossing
//! schedule instances.

/// Finds the min crossings for base problem using A*
pub mod base_min_crossings_astar;
/// Finds the min crossings for base problem distribute instances using A*
pub mod base_min_crossings_distribute_astar;
///// Finds the min crossings for base problem using MIP
////pub mod base_min_crossings_mip;
/// Solves by breeding all possible progeny in each generation
pub mod base_min_generations_enumerator;
/// Solves by breeding all non-dominating progeny in each generation
pub mod base_min_generations_enumerator_dominance;
/// Finds the min generations for base problem
pub mod base_min_generations_segment;
/// Finds the min generations for multi-chromosome diploid plants
pub mod greedy_multichrom_ti;

/// Finds the min resources for base problem using Beam Search
pub mod base_min_resources_beam_search;

/// ## Heuristic methods

/// Find a crossing schedule by repeated performing random crossings until the ideotype is created
pub mod base_heuristic_random_selection;
