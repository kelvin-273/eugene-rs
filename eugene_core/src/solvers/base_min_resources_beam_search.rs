use crate::extra::resources::{cost_of_crossing, RecRate};
use crate::plants::bit_array::SingleChromGenotype;
use crate::plants::dist_array::DistArray;
use crate::solution::BaseSolution;

#[derive(Debug, Clone, Default)]
pub struct Output {
    //expansions: usize,
    //pushed_nodes: usize,
    //children_created_by_branching: usize,
    //objective: Option<usize>,
}

pub fn breeding_program(
    n_loci: usize,
    pop_0: &[SingleChromGenotype],
    recombination_rates: &[f64],
    gamma: f64,
) -> Option<BaseSolution>
where
{
    let mut _best_solution: Option<BaseSolution> = None;
    let mut best_objective = f64::INFINITY;
    let mut _output = Output::default();

    let subproblem_generator = SubproblemGenerator::new(n_loci, pop_0, recombination_rates);

    for (xs, rec_rates) in subproblem_generator {
        let lower_bound = compute_lower_bound(&xs, &rec_rates, gamma);
        if lower_bound >= best_objective {
            continue;
        }
        // Solve the subproblem
        let (_res, objective) = match solve_subproblem(&xs, rec_rates, gamma) {
            None => continue,
            Some(res) => res,
        };
        if objective < best_objective {
            best_objective = objective;
            // Update best solution
        }
    }

    // Placeholder implementation
    Some(BaseSolution {
        tree_data: vec![],
        tree_type: vec![],
        tree_left: vec![],
        tree_right: vec![],
        objective: 0,
    })
}

struct SubproblemGenerator {
    // Define the structure of the subproblem generator
}

impl SubproblemGenerator {
    pub fn new(
        _n_loci: usize,
        _pop_0: &[SingleChromGenotype],
        _recombination_rates: &[f64],
    ) -> Self {
        // Initialize the subproblem generator
        SubproblemGenerator {}
    }
}

impl Iterator for SubproblemGenerator {
    type Item = (DistArray, RecRate);

    fn next(&mut self) -> Option<Self::Item> {
        // Generate the next subproblem
        todo!()
    }
}

fn compute_lower_bound(_xs: &DistArray, _recombination_rates: &RecRate, _gamma: f64) -> f64 {
    todo!()
}

fn solve_subproblem(_xs: &DistArray, _recombination_rates: RecRate, _gamma: f64) -> Option<(Vec<usize>, f64)> {
    unimplemented!()
}

struct _Node {
    // Define the structure of a node in the search tree
}

struct _SearchTree {
    // Define the structure of the search tree
}

struct _State {
    xs: DistArray,
    pr_gx: Vec<f64>,
    resourses_used: f64,
}

impl _State {
    pub fn _next_state_given_redistribution(
        &self,
        gamma: f64,
        _zs: DistArray,
        gx: usize,
        gy: usize,
    ) -> Self {
        let px = self.pr_gx[gx];
        let py = self.pr_gx[gy];
        let _next_cost = self.resourses_used + cost_of_crossing(gamma, px, py);
        let mut _next_pr_gx = self.pr_gx.clone();
        // Update based on new distribution
        todo!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plants::bit_array::SingleChromGenotype;

    #[test]
    fn test_breeding_program_2_loci() {
        let n_loci = 2;
        let pop_0 = vec![
            SingleChromGenotype::from_str("01", "01"),
            SingleChromGenotype::from_str("10", "10"),
        ];
        let recombination_rates = vec![0.1];
        let gamma = 0.5;
        let result = breeding_program(n_loci, &pop_0, &recombination_rates, gamma);
        assert!(result.is_some());
        assert_eq!(result.unwrap().objective, 1);
    }

    #[test]
    fn test_breeding_program_3_loci() {
        let n_loci = 3;
        let pop_0 = vec![
            SingleChromGenotype::from_str("000", "111"),
            SingleChromGenotype::from_str("011", "000"),
        ];
        let recombination_rates = vec![0.1, 0.2];
        let gamma = 0.5;
        let result = breeding_program(n_loci, &pop_0, &recombination_rates, gamma);
        assert!(result.is_some());
    }
}
