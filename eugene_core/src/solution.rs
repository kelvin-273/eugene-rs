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
