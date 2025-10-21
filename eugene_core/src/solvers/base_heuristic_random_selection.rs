use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solution::{BaseSolution, Objective};
use crate::solvers::base_min_generations_enumerator_dominance::filter_non_dominating_key;

/// Finds a crossing schedule by repeatedly performing random crossings until the ideotype is
/// created.
/// (This function is a shell for `breeding_program_random` with a default number of crossings.)
///
/// # Arguments
/// * `n_loci` - The number of loci in the genotype.
/// * `pop_0` - The initial population of single chromosome genotypes.
/// # Returns
/// An `Option<BaseSolution>` which contains the solution if the ideotype is found, or `None` if
/// not.
///
/// # Panics
/// Panics if `n_loci` is zero or if `pop_0` is empty.
/// # Example
/// ```
/// use crate::plants::bit_array::*;
/// use crate::solvers::base_heuristic_random_selection::breeding_program;
/// let n_loci = 10;
/// let pop_0 = vec![SingleChromGenotype::from_str("11001100", "00110011")];
/// let solution = breeding_program(n_loci, &pop_0);
/// assert!(solution.is_some());
/// ```
pub fn breeding_program(n_loci: usize, pop_0: &[SingleChromGenotype]) -> Option<BaseSolution> {
    let mut rng = rand::thread_rng();
    breeding_program_random(n_loci, pop_0, Some(20), &mut rng)
}

/// Finds a crossing schedule by repeatedly performing random crossings until the ideotype is
/// created.
///
/// # Arguments
/// * `n_loci` - The number of loci in the genotype.
/// * `pop_0` - The initial population of single chromosome genotypes.
/// # Returns
/// An `Option<BaseSolution>` which contains the solution if the ideotype is found, or `None` if
/// not.
///
/// # Panics
/// Panics if `n_loci` is zero or if `pop_0` is empty.
/// # Example
/// ```
/// use crate::plants::bit_array::*;
/// use crate::solvers::base_heuristic_random_selection::breeding_program;
/// let n_loci = 10;
/// let pop_0 = vec![SingleChromGenotype::from_str("11001100", "00110011")];
/// let solution = breeding_program_random(n_loci, &pop_0);
/// assert!(solution.is_some());
/// ```
pub fn breeding_program_random(
    n_loci: usize,
    pop_0: &[SingleChromGenotype],
    k_cross: Option<usize>,
    rng: impl rand::Rng,
) -> Option<BaseSolution> {
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    // TODO: check for infeasibility <18-07-24> //
    repeated_breeding_random(n_loci, pop_0, k_cross, rng)
        .iter()
        .find(|wx| wx.genotype() == &ideotype)
        .map(|wx| wx.to_base_sol(n_loci, Objective::Crossings))
}

/// Repeatedly performs random crossings until the ideotype is created.
/// # Arguments
/// * `n_loci` - The number of loci in the genotype.
/// * `pop_0` - The initial population of single chromosome genotypes.
/// * `k_cross` - An optional maximum number of crossings to perform.
/// # Returns
/// A vector of `WGen<SingleChromGenotype, SingleChromGamete>` which contains the population after
/// the crossings.
pub fn repeated_breeding_random(
    n_loci: usize,
    pop_0: &[SingleChromGenotype],
    k_cross: Option<usize>,
    mut rng: impl rand::Rng,
) -> Vec<WGen<SingleChromGenotype, SingleChromGamete>> {
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    let mut pop: Vec<WGen<SingleChromGenotype, SingleChromGamete>> =
        pop_0.iter().cloned().map(WGen::new).collect();
    let mut n_crossings = 0;
    while !pop.iter().any(|wx| wx.genotype() == &ideotype)
        && k_cross.map_or(true, |k_cross| n_crossings < k_cross)
    {
        let n_pop = pop.len();
        let i = rng.gen_range(0..n_pop);
        let j = rng.gen_range(i..n_pop);
        let wx = &pop[i];
        let wy = &pop[j];
        let wz = WGen::from_gametes(
            &wx.cross(CrosspointBitVec::new(
                Chrom::from(rng.gen::<bool>()),
                rng.gen_range(0..n_loci),
            )),
            &wy.cross(CrosspointBitVec::new(
                Chrom::from(rng.gen::<bool>()),
                rng.gen_range(0..n_loci),
            )),
        );
        pop.push(wz);
        n_crossings += 1;
    }
    pop
}

/// Finds a crossing schedule by repeatedly performing random crossings until the ideotype is
/// created. Crossings are performed only between non-dominating gametes.
/// # Arguments
/// * `n_loci` - The number of loci in the genotype.
/// * `pop_0` - The initial population of single chromosome genotypes.
/// * `k_cross` - An optional maximum number of crossings to perform.
/// # Returns
/// An `Option<BaseSolution>` which contains the solution if the ideotype is found, or `None` if
/// not.
/// # Panics
/// Panics if `n_loci` is zero or if `pop_0` is empty.
/// # Example
/// ```
/// use crate::plants::bit_array::*;
/// use crate::solvers::base_heuristic_random_selection::breeding_program_random_dominance;
/// let n_loci = 10;
/// let pop_0 = vec![SingleChromGenotype::from_str("11001100", "00110011")];
/// let solution = breeding_program_random_dominance(n_loci, &pop_0, Some(20), &mut
/// rand::thread_rng());
/// assert!(solution.is_some());
/// ```
pub fn breeding_program_random_dominance(
    n_loci: usize,
    pop_0: &[SingleChromGenotype],
    k_cross: Option<usize>,
    rng: impl rand::Rng,
) -> Option<BaseSolution> {
    // TODO: check for infeasibility <18-07-24> //
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    repeated_breeding_random_dominance(n_loci, pop_0, k_cross, rng)
        .iter()
        .find(|wx| wx.genotype() == &ideotype)
        .map(|wx| wx.to_base_sol(n_loci, Objective::Crossings))
}

/// Repeatedly performs random crossings until the ideotype is created.
/// Crossings are performed only between non-dominating gametes.
/// # Arguments
/// * `n_loci` - The number of loci in the genotype.
/// * `pop_0` - The initial population of single chromosome genotypes.
/// * `k_cross` - An optional maximum number of crossings to perform.
/// # Returns
/// A vector of `WGen<SingleChromGenotype, SingleChromGamete>` which contains the population after
/// the crossings.
/// # Panics
/// Panics if `n_loci` is zero or if `pop_0` is empty.
/// # Example
/// ```
/// use crate::plants::bit_array::*;
/// use crate::solvers::base_heuristic_random_selection::repeated_breeding_random_dominance;
/// let n_loci = 10;
/// let pop_0 = vec![SingleChromGenotype::from_str("11001100", "00110011")];
/// let pop = repeated_breeding_random_dominance(n_loci, &pop_0, Some(20), &mut
/// rand::thread_rng());
/// assert!(pop.iter().any(|wx| wx.genotype() == &SingleChromGenotype::ideotype(n_loci)));
/// ```
pub fn repeated_breeding_random_dominance(
    n_loci: usize,
    pop_0: &[SingleChromGenotype],
    k_cross: Option<usize>,
    mut rng: impl rand::Rng,
) -> Vec<WGen<SingleChromGenotype, SingleChromGamete>> {
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    let mut pop: Vec<WGen<SingleChromGenotype, SingleChromGamete>> =
        pop_0.iter().cloned().map(WGen::new).collect();
    let mut n_crossings = 0;
    while !pop.iter().any(|wx| wx.genotype() == &ideotype)
        && k_cross.map_or(true, |k_cross| n_crossings < k_cross)
    {
        let non_dominating_wgametes = filter_non_dominating_key::<_, _, _, DomGamete>(
            pop.iter()
                .flat_map(|wx| CrosspointBitVec::crosspoints(&n_loci).map(|k| wx.cross(k))),
            |wg| wg.gamete(),
        );
        let i = rng.gen_range(0..non_dominating_wgametes.len());
        let j = rng.gen_range(i..non_dominating_wgametes.len());
        pop.push(WGen::from_gametes(
            &non_dominating_wgametes[i],
            &non_dominating_wgametes[j],
        ));
        n_crossings += 1;
    }
    pop
}
