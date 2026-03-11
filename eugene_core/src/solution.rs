use std::ops::Index;

use bit_vec::BitVec;

use crate::{
    abstract_plants::{Allele, Diploid, Genotype, Haploid, IndexAllele, WGen},
    extra::resources::RecRate,
};

#[derive(Debug, PartialEq, Eq)]
pub struct BaseSolution {
    pub tree_data: Vec<Vec<Vec<i32>>>,
    pub tree_type: Vec<&'static str>,
    pub tree_left: Vec<usize>,
    pub tree_right: Vec<usize>,
    pub objective: usize,
}

#[derive(Debug, PartialEq, Eq)]
pub enum Objective {
    Generations,
    Crossings,
}

impl BaseSolution {}

#[derive(Debug, PartialEq, Eq, Clone)]
enum TreeType {
    Node,
    Leaf,
}

#[derive(Debug, PartialEq, Eq)]
pub struct CrossingSchedule {
    n_loci: usize,
    tree_data: Vec<[BitVec; 2]>,
    tree_type: Vec<TreeType>,
    tree_left: Vec<usize>,
    tree_right: Vec<usize>,
}

/// Compact representation of a crossing schedule.
///
/// Maintains the following invariants:
/// - crossing schedule is non-empty
/// - elements are stored `Node`s first and `Leave`s second
/// - the target if constructed is stored in index 0
/// - crossing schedule is an acyclic graph
/// - parents of each node are stored to the right of the node
impl CrossingSchedule {
    pub fn len(&self) -> usize {
        self.tree_type.len()
    }

    #[inline]
    pub fn n_loci(&self) -> usize {
        self.n_loci
    }

    fn genotype_view(&self, i: usize) -> Option<GenotypeView> {
        if i >= self.len() {
            return None;
        }
        Some(GenotypeView {
            crossing_schedule: self,
            idx: i,
        })
    }

    /// Computes the number of generations required to create the first node in the crossing
    /// schedule. Assumes the elements in the crossing schedule are ordered with the target in
    /// index 0 and leaves in the end.
    pub fn generations(&self) -> usize {
        let mut dp = vec![0; self.len()];
        for i in (0..self.len()).rev() {
            dp[i] = match self.tree_type[i] {
                TreeType::Leaf => 0,
                TreeType::Node => 1 + dp[self.tree_left[i]].max(dp[self.tree_right[i]]),
            }
        }
        dp[0]
    }

    pub fn crossings(&self) -> usize {
        self.tree_type
            .iter()
            .filter(|typ| **typ == TreeType::Node)
            .count()
    }

    pub fn resources(&self, rec_rate: &RecRate, gamma: f64) -> usize {
        let mut out = 0;
        let mut i = 0;
        while self.tree_type[i] == TreeType::Node {
            out += rec_rate
                .crossing_resources(
                    gamma,
                    &self.genotype_view(self.tree_left[i]).unwrap(),
                    &self.genotype_view(self.tree_right[i]).unwrap(),
                    &self.genotype_view(i).unwrap(),
                )
                .ceil() as usize;
            i += 1;
        }
        out
    }
}

struct GenotypeView<'a> {
    crossing_schedule: &'a CrossingSchedule,
    idx: usize,
}

struct GameteView<'a> {
    crossing_schedule: &'a CrossingSchedule,
    idx: usize,
    chrom: bool,
}

impl<'a> Diploid<GameteView<'a>> for GenotypeView<'a> {
    fn upper(&self) -> GameteView<'a> {
        GameteView {
            crossing_schedule: self.crossing_schedule,
            idx: self.idx,
            chrom: false,
        }
    }

    fn lower(&self) -> GameteView<'a> {
        GameteView {
            crossing_schedule: self.crossing_schedule,
            idx: self.idx,
            chrom: true,
        }
    }
}

impl<'a> Haploid for GameteView<'a> {
    fn alleles(&self) -> Vec<Allele> {
        self.crossing_schedule.tree_data[self.idx][self.chrom as usize]
            .iter()
            .map(Into::into)
            .collect()
    }
}

impl<'a> IndexAllele<usize> for GameteView<'a> {
    fn index(&self, idx: usize) -> Allele {
        self.crossing_schedule.tree_data[self.idx][self.chrom as usize]
            .get(idx)
            .expect("idx out of bounds")
            .into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abstract_plants::*;
    use crate::plants::bit_array::*;

    #[test]
    fn wgens_to_base_sol_ideotype_test() {
        let x = SingleChromGenotype::ideotype(4);
        let wx: WGen<SingleChromGenotype, SingleChromGamete> = WGen::new(x);
        let sol = wx.to_base_sol(4, Objective::Crossings);
        assert_eq!(0, sol.objective);
        assert_eq!(vec!["Leaf"], sol.tree_type);
        assert_eq!(vec![0], sol.tree_left);
        assert_eq!(vec![0], sol.tree_right);
    }

    #[test]
    fn wgen_to_base_sol_pair_test() {
        let n_loci = 2;
        let xl = SingleChromGenotype::from_str("10", "10");
        let xr = SingleChromGenotype::from_str("01", "01");

        let wxl = WGen::new(xl);
        let wxr = WGen::new(xr);

        let gl = SingleChromGamete::from_str("10");
        let gr = SingleChromGamete::from_str("01");

        let wgl = WGam::new_from_genotype(gl, wxl);
        let wgr = WGam::new_from_genotype(gr, wxr);

        let wxz = WGen::from_gametes(&wgl, &wgr);
        let wgz = WGam::new_from_genotype(SingleChromGamete::from_str("11"), wxz);

        let wx_star = WGen::from_gametes(&wgz, &wgz);

        let sol_c = wx_star.to_base_sol(n_loci, Objective::Crossings);
        assert_eq!(2, sol_c.objective);
        assert_eq!(vec!["Node", "Node", "Leaf", "Leaf"], sol_c.tree_type);
        assert_eq!(vec![2, 3, 0, 0], sol_c.tree_left);
        assert_eq!(vec![2, 4, 0, 0], sol_c.tree_right);

        let sol_g = wx_star.to_base_sol(n_loci, Objective::Generations);
        assert_eq!(2, sol_g.objective);
        assert_eq!(vec!["Node", "Node", "Leaf", "Leaf"], sol_g.tree_type);
    }

    #[test]
    fn wgens_to_base_sol_3loci_test() {
        let n_loci = 3;
        let x = SingleChromGenotype::from_str("011", "111");
        let g = SingleChromGamete::from_str("111");
        let wx = WGen::new(x);
        let wg = WGam::new_from_genotype(g, wx);
        let wz = WGen::from_gametes(&wg, &wg);

        let sol_c = wz.to_base_sol(n_loci, Objective::Crossings);
        assert_eq!(1, sol_c.objective);
        assert_eq!(vec!["Node", "Leaf"], sol_c.tree_type);
        assert_eq!(vec![2, 0], sol_c.tree_left);
        assert_eq!(vec![2, 0], sol_c.tree_right);

        let sol_g = wz.to_base_sol(n_loci, Objective::Generations);
        assert_eq!(1, sol_g.objective);
        assert_eq!(vec!["Node", "Leaf"], sol_g.tree_type);
        assert_eq!(vec![2, 0], sol_g.tree_left);
        assert_eq!(vec![2, 0], sol_g.tree_right);
    }
}
