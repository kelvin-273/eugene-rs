use std::{
    collections::{HashMap, HashSet, VecDeque},
    ops::Index,
};

use bit_vec::BitVec;

use crate::{
    abstract_plants::{Allele, Diploid, Genotype, Haploid, IndexAllele, WGen},
    extra::resources::RecRate,
    plants::bit_array::{SingleChromGamete, SingleChromGenotype},
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
pub enum TreeType {
    Node,
    Leaf,
}

#[derive(Debug, PartialEq, Eq, Clone)]
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
    pub fn new(
        n_loci: usize,
        tree_data: Vec<[BitVec; 2]>,
        tree_type: Vec<TreeType>,
        tree_left: Vec<usize>,
        tree_right: Vec<usize>,
    ) -> Self {
        assert!(!tree_type.is_empty());
        assert_eq!(tree_data.len(), tree_type.len());
        assert_eq!(tree_left.len(), tree_type.len());
        assert_eq!(tree_right.len(), tree_type.len());
        Self {
            n_loci,
            tree_data,
            tree_type,
            tree_left,
            tree_right,
        }
    }

    pub fn len(&self) -> usize {
        self.tree_type.len()
    }

    #[inline]
    pub fn n_loci(&self) -> usize {
        self.n_loci
    }

    pub fn genotype_view(&'_ self, i: usize) -> Option<GenotypeView<'_>> {
        if i >= self.len() {
            return None;
        }
        Some(GenotypeView {
            crossing_schedule: self,
            idx: i,
        })
    }

    pub fn left_parent(&'_ self, i: usize) -> Option<GenotypeView<'_>> {
        if i >= self.len() || self.tree_type[i] == TreeType::Leaf {
            return None;
        }
        Some(GenotypeView {
            crossing_schedule: self,
            idx: self.tree_left[i],
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
            out += rec_rate.crossing_resources(
                gamma,
                &self.genotype_view(self.tree_left[i]).unwrap(),
                &self.genotype_view(self.tree_right[i]).unwrap(),
                &self.genotype_view(i).unwrap(),
            );
            i += 1;
        }
        out
    }
}

impl From<WGen<SingleChromGenotype, SingleChromGamete>> for CrossingSchedule {
    fn from(wtarget: WGen<SingleChromGenotype, SingleChromGamete>) -> Self {
        let n_loci = wtarget.genotype().n_loci();

        let mut tree_refs = Vec::new();
        let mut tree_type = Vec::new();
        let mut tree_left = Vec::new();
        let mut tree_right = Vec::new();

        let mut hash_map = HashMap::new();

        // DFS to construct the crossing schedule tree from the target WGen. Maintains a hash set
        // to avoid visiting the same WGen twice and to assign indices to WGen in the crossing
        // schedule.
        fn visit(
            wx: WGen<SingleChromGenotype, SingleChromGamete>,
            hash_map: &mut HashMap<SingleChromGenotype, usize>,
            tree_refs: &mut Vec<WGen<SingleChromGenotype, SingleChromGamete>>,
            tree_type: &mut Vec<TreeType>,
            tree_left: &mut Vec<usize>,
            tree_right: &mut Vec<usize>,
        ) -> usize {
            if let Some(idx) = hash_map.get(wx.genotype()) {
                return *idx;
            }
            match wx.history() {
                None => {
                    tree_type.push(TreeType::Leaf);
                    tree_left.push(0);
                    tree_right.push(0);
                }
                Some((wgx, wgy)) => {
                    let left_idx = visit(
                        wgx.history().unwrap(),
                        hash_map,
                        tree_refs,
                        tree_type,
                        tree_left,
                        tree_right,
                    );
                    let right_idx = visit(
                        wgy.history().unwrap(),
                        hash_map,
                        tree_refs,
                        tree_type,
                        tree_left,
                        tree_right,
                    );
                    tree_type.push(TreeType::Node);
                    tree_left.push(left_idx);
                    tree_right.push(right_idx);
                }
            };
            let idx = tree_refs.len();
            tree_refs.push(wx.clone());
            hash_map.insert(wx.genotype().clone(), idx);
            idx
        }

        visit(
            wtarget,
            &mut hash_map,
            &mut tree_refs,
            &mut tree_type,
            &mut tree_left,
            &mut tree_right,
        );

        tree_refs.reverse();
        tree_type.reverse();
        tree_left.reverse();
        tree_right.reverse();

        let tree_data = tree_refs
            .into_iter()
            .map(|wx| wx.genotype().clone().into())
            .collect();
        tree_left
            .iter_mut()
            .zip(tree_type.iter())
            .for_each(|(left, typ)| {
                *left = match typ {
                    TreeType::Node => tree_type.len() - *left - 1,
                    TreeType::Leaf => 0,
                }
            });
        tree_right
            .iter_mut()
            .zip(tree_type.iter())
            .for_each(|(right, typ)| {
                *right = match typ {
                    TreeType::Node => tree_type.len() - *right - 1,
                    TreeType::Leaf => 0,
                }
            });

        Self {
            n_loci,
            tree_data,
            tree_type,
            tree_left,
            tree_right,
        }
    }
}

#[derive(Debug)]
pub struct GenotypeView<'a> {
    crossing_schedule: &'a CrossingSchedule,
    idx: usize,
}

#[derive(Debug)]
pub struct GameteView<'a> {
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

impl From<&CrossingSchedule>
    for (
        Vec<Vec<Vec<i32>>>,
        Vec<&'static str>,
        Vec<usize>,
        Vec<usize>,
    )
{
    fn from(cs: &CrossingSchedule) -> Self {
        let tree_data = cs
            .tree_data
            .iter()
            .map(|x| {
                x.iter()
                    .map(|xi| xi.iter().map(|b| b as i32).collect())
                    .collect()
            })
            .collect();
        let tree_type = cs
            .tree_type
            .iter()
            .map(|ty| match ty {
                TreeType::Node => "Node",
                TreeType::Leaf => "Leaf",
            })
            .collect();
        let tree_left = cs.tree_left.clone();
        let tree_right = cs.tree_right.clone();

        (tree_data, tree_type, tree_left, tree_right)
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

    #[test]
    fn crossing_schedule_from_wgen_test() {
        let x1 = SingleChromGenotype::from_str("10101", "01010");
        let wx1 = WGen::new(x1);

        let g1 = SingleChromGamete::from_str("01101");
        let g2 = SingleChromGamete::from_str("01011");
        let wg1 = WGam::new_from_genotype(g1, wx1.clone());
        let wg2 = WGam::new_from_genotype(g2, wx1.clone());
        let wx2 = WGen::from_gametes(&wg1, &wg2);

        let g3 = SingleChromGamete::from_str("10101");
        let g4 = SingleChromGamete::from_str("01111");

        let wg3 = WGam::new_from_genotype(g3, wx1.clone());
        let wg4 = WGam::new_from_genotype(g4, wx2.clone());
        let wx3 = WGen::from_gametes(&wg3, &wg4);

        let schedule = CrossingSchedule::from(wx3);
        assert_eq!(schedule.n_loci(), 5);
        assert_eq!(schedule.len(), 3);
        assert_eq!(
            schedule.tree_type,
            vec![TreeType::Node, TreeType::Node, TreeType::Leaf]
        );
        assert_eq!(schedule.tree_left, vec![2, 2, 0]);
        assert_eq!(schedule.tree_right, vec![1, 2, 0]);
    }

    #[test]
    fn crossing_schedule_generations_test() {
        let x1 = SingleChromGenotype::from_str("10101", "01010");
        let wx1 = WGen::new(x1);

        let g1 = SingleChromGamete::from_str("01101");
        let g2 = SingleChromGamete::from_str("01011");
        let wg1 = WGam::new_from_genotype(g1, wx1.clone());
        let wg2 = WGam::new_from_genotype(g2, wx1.clone());
        let wx2 = WGen::from_gametes(&wg1, &wg2);

        let g3 = SingleChromGamete::from_str("10101");
        let g4 = SingleChromGamete::from_str("01111");

        let wg3 = WGam::new_from_genotype(g3, wx1.clone());
        let wg4 = WGam::new_from_genotype(g4, wx2.clone());
        let wx3 = WGen::from_gametes(&wg3, &wg4);

        let schedule = CrossingSchedule::from(wx3);
        assert_eq!(schedule.generations(), 2);
    }

    #[test]
    fn crossing_schedule_crossings_test() {
        let x1 = SingleChromGenotype::from_str("10101", "01010");
        let wx1 = WGen::new(x1);

        let g1 = SingleChromGamete::from_str("01101");
        let g2 = SingleChromGamete::from_str("01011");
        let wg1 = WGam::new_from_genotype(g1, wx1.clone());
        let wg2 = WGam::new_from_genotype(g2, wx1.clone());
        let wx2 = WGen::from_gametes(&wg1, &wg2);

        let g3 = SingleChromGamete::from_str("10101");
        let g4 = SingleChromGamete::from_str("01111");

        let wg3 = WGam::new_from_genotype(g3, wx1.clone());
        let wg4 = WGam::new_from_genotype(g4, wx2.clone());
        let wx3 = WGen::from_gametes(&wg3, &wg4);

        let schedule = CrossingSchedule::from(wx3);
        assert_eq!(schedule.crossings(), 2);
    }

    #[test]
    fn crossing_schedule_resources_test() {
        let x1 = SingleChromGenotype::from_str("10", "10");
        let x2 = SingleChromGenotype::from_str("01", "01");

        let g1 = SingleChromGamete::from_str("10");
        let g2 = SingleChromGamete::from_str("01");
        let g3 = SingleChromGamete::from_str("11");

        let wx1 = WGen::new(x1);
        let wx2 = WGen::new(x2);
        let wg1 = WGam::new_from_genotype(g1, wx1.clone());
        let wg2 = WGam::new_from_genotype(g2, wx2.clone());
        let wx3 = WGen::from_gametes(&wg1, &wg2);
        let wg3 = WGam::new_from_genotype(g3.clone(), wx3.clone());
        let wx4 = WGen::from_gametes(&wg3, &wg2);
        let wg4 = WGam::new_from_genotype(g3, wx4.clone());
        let wx5 = WGen::from_gametes(&wg4, &wg4);
        // Head of a shorter crossing schedule with the same target as wx5.
        let wx6 = WGen::from_gametes(&wg3, &wg3);

        let rec_rate = RecRate::from(vec![0.1]);
        let gamma = 0.99;
        let crossing_schedule0 = CrossingSchedule::from(wx3);
        let crossing_schedule1 = CrossingSchedule::from(wx5);
        let crossing_schedule2 = CrossingSchedule::from(wx6);
        assert_eq!(crossing_schedule0.resources(&rec_rate, gamma), 1);
        assert_eq!(crossing_schedule1.resources(&rec_rate, gamma), 108);
        assert_eq!(crossing_schedule2.resources(&rec_rate, gamma), 1841);
    }

    #[test]
    fn crossing_schedule_view_test() {
        use Allele::*;
        let x1 = SingleChromGenotype::from_str("10101", "01010");
        let wx1 = WGen::new(x1);

        let g1 = SingleChromGamete::from_str("01101");
        let g2 = SingleChromGamete::from_str("01011");
        let wg1 = WGam::new_from_genotype(g1, wx1.clone());
        let wg2 = WGam::new_from_genotype(g2, wx1.clone());
        let wx2 = WGen::from_gametes(&wg1, &wg2);

        let g3 = SingleChromGamete::from_str("10101");
        let g4 = SingleChromGamete::from_str("01111");

        let wg3 = WGam::new_from_genotype(g3, wx1.clone());
        let wg4 = WGam::new_from_genotype(g4, wx2.clone());
        let wx3 = WGen::from_gametes(&wg3, &wg4);

        let schedule = CrossingSchedule::from(wx3);

        let vz = schedule
            .genotype_view(1)
            .expect("One'th element yields None");

        assert_eq!(vz.upper().index(0), Z);
        assert_eq!(vz.upper().index(1), O);
        assert_eq!(vz.upper().index(2), O);
        assert_eq!(vz.upper().index(3), Z);
        assert_eq!(vz.upper().index(4), O);

        assert_eq!(vz.lower().index(0), Z);
        assert_eq!(vz.lower().index(1), O);
        assert_eq!(vz.lower().index(2), Z);
        assert_eq!(vz.lower().index(3), O);
        assert_eq!(vz.lower().index(4), O);
    }
}
