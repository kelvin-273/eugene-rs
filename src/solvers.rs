//! A submodule for the different solvers.
//!
//! Each solver has a function `breeding_program` which solves solves instances of various crossing
//! schedule instances.

/// Solves by breeding all possible progeny in each generation
pub mod base_min_generations_enumerator;
pub mod base_min_generations_enumerator2;
/// Solves by breeding all non-dominating progeny in each generation
pub mod base_min_generations_enumerator_dominance;
/// Finds the min generations for base problem
pub mod base_min_generations_segment;
/// Finds the min generations for multi-chromosome diploid plants
pub mod greedy_multichrom_ti;
/// Finds the min crossings for base problem using A*
pub mod base_min_crossings_astar;
/// Finds the min crossings for base problem using MIP
//pub mod base_min_crossings_mip;
/// Finds the min crossings for base problem using zigzag decomposition
pub mod base_min_crossings_zigzag;
/// Finds the min crossings for base problem using LNS approaches
pub mod base_min_crossings_lns;
/// Finds the min crossings for base problem using backtracking
pub mod base_min_crossings_backtracking;
/// Finds the min crossings for base problem distribute instances using A*
pub mod base_min_crossings_distribute_astar;
