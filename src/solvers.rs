//! A submodule for the different solvers.
//!
//! Each solver has a function `breeding_program` which solves solves instances of various crossing
//! schedule instances.

/// Finds the min crossings for base problem using A*
pub mod base_min_crossings_astar;
/// Finds the min crossings for base problem using zigzag decomposition
pub mod base_min_crossings_zigzag;
/// Solves by breeding all possible progeny in each generation
pub mod enumerator;
/// Solves by breeding all non-dominating progeny in each generation
pub mod enumerator_dominance;
/// Finds the min generations for base problem
pub mod greedy_base;
/// Finds the min generations for multi-chromosome diploid plants
pub mod greedy_multichrom_ti;
