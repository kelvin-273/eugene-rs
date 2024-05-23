use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use std::collections::{HashMap, VecDeque};

pub struct BaseSolution {
    pub tree_data: Vec<Vec<Vec<i32>>>,
    pub tree_type: Vec<&'static str>,
    pub tree_left: Vec<usize>,
    pub tree_right: Vec<usize>,
    pub objective: usize,
}

enum Objective {
    Generations,
    Crossings,
}

impl BaseSolution {
    #[deprecated]
    pub fn from_wgen(
        n_loci: usize,
        x_star: &WGen<SingleChromGenotype, SingleChromGamete>,
    ) -> BaseSolution {
        wgen_to_base_sol(n_loci, &x_star, Objective::Crossings)
    }

    pub fn min_gen_from_wgen(
        n_loci: usize,
        x_star: &WGen<SingleChromGenotype, SingleChromGamete>,
    ) -> BaseSolution {
        wgen_to_base_sol(n_loci, &x_star, Objective::Generations)
    }

    pub fn min_cross_from_wgen(
        n_loci: usize,
        x_star: &WGen<SingleChromGenotype, SingleChromGamete>,
    ) -> BaseSolution {
        wgen_to_base_sol(n_loci, &x_star, Objective::Crossings)
    }

    pub fn min_gen_from_wgens(
        n_loci: usize,
        x_star: &WGenS<SingleChromGenotype, SingleChromGamete>,
    ) -> BaseSolution {
        wgens_to_base_sol(n_loci, &x_star, Objective::Generations)
    }

    pub fn min_cross_from_wgens(
        n_loci: usize,
        x_star: &WGenS<SingleChromGenotype, SingleChromGamete>,
    ) -> BaseSolution {
        wgens_to_base_sol(n_loci, &x_star, Objective::Crossings)
    }
}

fn wgen_to_base_sol(
    n_loci: usize,
    x_star: &WGen<SingleChromGenotype, SingleChromGamete>,
    obj: Objective,
) -> BaseSolution {
    // get the number of crossings and cells required
    let mut q_node = VecDeque::new();
    let mut n_cells = 1;
    let mut n_cross = 0;

    /*
     * map every genotype to a unique index in the order they were seen starting from 0
     */

    // create index_map
    let mut index_map: HashMap<SingleChromGenotype, usize> = HashMap::new();
    index_map.insert(x_star.genotype.clone(), 0);

    q_node.push_back(x_star);
    if x_star.history.is_some() {
        n_cross += 1;
    }

    // using q_node to process all genotypes for the counting
    while let Some(z) = q_node.pop_front() {
        if let Some((wgx, wgy)) = z.history.as_ref() {
            // QUESTION: how is it that a single & here
            // makes all the difference?
            let x = wgx.history.first().unwrap().as_ref();
            if !index_map.contains_key(&x.genotype) {
                index_map.insert(x.genotype.clone(), n_cells);
                n_cells += 1;
                q_node.push_back(x);
                if x.history.is_some() {
                    n_cross += 1;
                }
            }

            let y = wgy.history.first().unwrap().as_ref();
            if !index_map.contains_key(&y.genotype) {
                index_map.insert(y.genotype.clone(), n_cells);
                n_cells += 1;
                q_node.push_back(y);
                if y.history.is_some() {
                    n_cross += 1;
                }
            }
        }
    }
    index_map.drain();

    // initialise arrays
    let mut tree_data = vec![vec![vec![0; n_loci]; 2]; n_cells];
    let mut tree_type = vec!["Null"; n_cells];
    let mut tree_left = vec![0; n_cells];
    let mut tree_right = vec![0; n_cells];

    // flush index_map because nodes and leaves have mixed order
    index_map.drain();
    index_map.insert(x_star.genotype.clone(), 0);

    q_node.push_back(&x_star);
    let mut i_node = 0;
    let mut i_leaf = n_cross;

    // using q_node to process all genotypes for the counting
    while let Some(z) = q_node.pop_front() {
        let &i_z = index_map.get(&z.genotype).unwrap();

        // fill in data
        for (j, a) in z.genotype.upper().alleles().iter().enumerate() {
            tree_data[i_z][0][j] = match a {
                Allele::Z => 0,
                Allele::O => 1,
            };
        }
        for (j, a) in z.genotype.lower().alleles().iter().enumerate() {
            tree_data[i_z][1][j] = match a {
                Allele::Z => 0,
                Allele::O => 1,
            };
        }

        if let Some((wgx, wgy)) = z.history.as_ref() {
            tree_type[i_z] = "Node";

            let x = wgx.history.first().unwrap().as_ref();
            if !index_map.contains_key(&x.genotype) {
                if x.history.is_some() {
                    i_node += 1;
                    index_map.insert(x.genotype.clone(), i_node);
                } else {
                    index_map.insert(x.genotype.clone(), i_leaf);
                    i_leaf += 1;
                }
                q_node.push_back(x);
            }
            let &i_x = index_map.get(&x.genotype).unwrap();
            tree_left[i_z] = i_x + 1;

            let y = wgy.history.first().unwrap().as_ref();
            if !index_map.contains_key(&y.genotype) {
                if y.history.is_some() {
                    i_node += 1;
                    index_map.insert(y.genotype.clone(), i_node);
                } else {
                    index_map.insert(y.genotype.clone(), i_leaf);
                    i_leaf += 1;
                }
                q_node.push_back(y);
            }
            let &i_y = index_map.get(&y.genotype).unwrap();
            tree_right[i_z] = i_y + 1;
        } else {
            tree_type[i_z] = "Leaf";
        }
    }

    // compute generations if required
    let objective = match obj {
        Objective::Crossings => n_cross,
        Objective::Generations => {
            let mut generations_dp = vec![-1; n_cells];
            fn aux(
                tree_type: &Vec<&str>,
                tree_left: &Vec<usize>,
                tree_right: &Vec<usize>,
                generations_dp: &mut Vec<i32>,
                i: usize,
            ) -> i32 {
                if generations_dp[i] < 0 {
                    generations_dp[i] = match tree_type[i] {
                        "Leaf" => 0,
                        "Node" => {
                            let gen_l = aux(
                                tree_type,
                                tree_left,
                                tree_right,
                                generations_dp,
                                tree_left[i] - 1,
                            );
                            let gen_r = aux(
                                tree_type,
                                tree_left,
                                tree_right,
                                generations_dp,
                                tree_right[i] - 1,
                            );
                            gen_l.max(gen_r) + 1
                        }
                        _ => panic!("tree type is neither Leaf nor Node"),
                    }
                }
                generations_dp[i]
            }
            aux(&tree_type, &tree_left, &tree_right, &mut generations_dp, 0) as usize
        }
    };

    BaseSolution {
        tree_data,
        tree_type,
        tree_left,
        tree_right,
        objective,
    }
}

fn wgens_to_base_sol(
    n_loci: usize,
    x_star: &WGenS<SingleChromGenotype, SingleChromGamete>,
    obj: Objective,
) -> BaseSolution {
    // get the number of crossings and cells required
    let mut q_node = VecDeque::new();
    let mut n_cells = 1;
    let mut n_cross = 0;

    /*
     * map every genotype to a unique index in the order they were seen starting from 0
     */

    // create index_map
    let mut index_map: HashMap<SingleChromGenotype, usize> = HashMap::new();
    index_map.insert(x_star.genotype.clone(), 0);

    q_node.push_back(x_star);
    if x_star.history.is_some() {
        n_cross += 1;
    }

    // using q_node to process all genotypes for the counting
    while let Some(z) = q_node.pop_front() {
        if let Some((wgx, wgy)) = z.history.as_ref() {
            // QUESTION: how is it that a single & here
            // makes all the difference?
            let x = wgx
                .history
                .as_ref()
                .expect("gamete wgx should have a history");
            if !index_map.contains_key(&x.genotype) {
                index_map.insert(x.genotype.clone(), n_cells);
                n_cells += 1;
                q_node.push_back(x);
                if x.history.is_some() {
                    n_cross += 1;
                }
            }

            let y = wgy
                .history
                .as_ref()
                .expect("gamete wgy should have a history");
            if !index_map.contains_key(&y.genotype) {
                index_map.insert(y.genotype.clone(), n_cells);
                n_cells += 1;
                q_node.push_back(y);
                if y.history.is_some() {
                    n_cross += 1;
                }
            }
        }
    }
    index_map.drain();

    // initialise arrays
    let mut tree_data = vec![vec![vec![0; n_loci]; 2]; n_cells];
    let mut tree_type = vec!["Null"; n_cells];
    let mut tree_left = vec![0; n_cells];
    let mut tree_right = vec![0; n_cells];

    // flush index_map because nodes and leaves have mixed order
    index_map.drain();
    index_map.insert(x_star.genotype.clone(), 0);

    q_node.push_back(&x_star);
    let mut i_node = 0;
    let mut i_leaf = n_cross;

    // using q_node to process all genotypes for the counting
    while let Some(z) = q_node.pop_front() {
        let &i_z = index_map
            .get(&z.genotype)
            .expect("z.genotype is not in the index_map");

        // fill in data
        for (j, a) in z.genotype.upper().alleles().iter().enumerate() {
            tree_data[i_z][0][j] = match a {
                Allele::Z => 0,
                Allele::O => 1,
            };
        }
        for (j, a) in z.genotype.lower().alleles().iter().enumerate() {
            tree_data[i_z][1][j] = match a {
                Allele::Z => 0,
                Allele::O => 1,
            };
        }

        if let Some((wgx, wgy)) = z.history.as_ref() {
            tree_type[i_z] = "Node";

            let x = wgx
                .history
                .as_ref()
                .expect("wgx should have a history (in the queue)");
            if !index_map.contains_key(&x.genotype) {
                if x.history.is_some() {
                    i_node += 1;
                    index_map.insert(x.genotype.clone(), i_node);
                } else {
                    index_map.insert(x.genotype.clone(), i_leaf);
                    i_leaf += 1;
                }
                q_node.push_back(x);
            }
            let &i_x = index_map
                .get(&x.genotype)
                .expect("x.genotype should be in the index_map");
            tree_left[i_z] = i_x + 1;

            let y = wgy
                .history
                .as_ref()
                .expect("wgy should have a history (in the queue)");
            if !index_map.contains_key(&y.genotype) {
                if y.history.is_some() {
                    i_node += 1;
                    index_map.insert(y.genotype.clone(), i_node);
                } else {
                    index_map.insert(y.genotype.clone(), i_leaf);
                    i_leaf += 1;
                }
                q_node.push_back(y);
            }
            let &i_y = index_map
                .get(&y.genotype)
                .expect("y.genotype should be in the index_map");
            tree_right[i_z] = i_y + 1;
        } else {
            tree_type[i_z] = "Leaf";
        }
    }

    // compute generations if required
    let objective = match obj {
        Objective::Crossings => n_cross,
        Objective::Generations => {
            let mut generations_dp = vec![-1; n_cells];
            fn aux(
                tree_type: &Vec<&str>,
                tree_left: &Vec<usize>,
                tree_right: &Vec<usize>,
                generations_dp: &mut Vec<i32>,
                i: usize,
            ) -> i32 {
                if generations_dp[i] < 0 {
                    generations_dp[i] = match tree_type[i] {
                        "Leaf" => 0,
                        "Node" => {
                            let gen_l = aux(
                                tree_type,
                                tree_left,
                                tree_right,
                                generations_dp,
                                tree_left[i] - 1,
                            );
                            let gen_r = aux(
                                tree_type,
                                tree_left,
                                tree_right,
                                generations_dp,
                                tree_right[i] - 1,
                            );
                            gen_l.max(gen_r) + 1
                        }
                        _ => panic!("tree type is neither Leaf nor Node"),
                    }
                }
                generations_dp[i]
            }
            aux(&tree_type, &tree_left, &tree_right, &mut generations_dp, 0) as usize
        }
    };

    BaseSolution {
        tree_data,
        tree_type,
        tree_left,
        tree_right,
        objective,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    use std::rc::Rc;

    #[test]
    fn wgens_to_base_sol_ideotype_test() {
        let x = SingleChromGenotype::ideotype(4);
        let wx: WGenS<SingleChromGenotype, SingleChromGamete> = WGenS::new(x);
        let sol = BaseSolution::min_cross_from_wgens(4, &wx);
        assert_eq!(0, sol.objective);
        assert_eq!(vec!["Leaf"], sol.tree_type);
        assert_eq!(vec![0], sol.tree_left);
        assert_eq!(vec![0], sol.tree_right);
    }

    #[test]
    fn wgens_to_base_sol_pair_test() {
        let xl = SingleChromGenotype::from_str("10", "10");
        let xr = SingleChromGenotype::from_str("01", "01");

        let wxl = Rc::new(WGenS::new(xl));
        let wxr = Rc::new(WGenS::new(xr));

        let gl = SingleChromGamete::from_str("10");
        let gr = SingleChromGamete::from_str("01");

        let wgl = Rc::new(WGamS {
            gamete: gl,
            history: Some(wxl),
        });
        let wgr = Rc::new(WGamS {
            gamete: gr,
            history: Some(wxr),
        });

        let wxz = Rc::new(WGenS {
            genotype: SingleChromGenotype::from_gametes(&wgl.gamete, &wgr.gamete),
            history: Some((wgl, wgr)),
        });

        let wgz = Rc::new(WGamS {
            gamete: SingleChromGamete::from_str("11"),
            history: Some(wxz)
        });

        let wx_star = WGenS {
            genotype: SingleChromGenotype::from_gametes(&wgz.gamete, &wgz.gamete),
            history: Some((wgz.clone(), wgz)),
        };

        let sol_c = BaseSolution::min_cross_from_wgens(2, &wx_star);
        assert_eq!(2, sol_c.objective);
        assert_eq!(vec!["Node", "Node", "Leaf", "Leaf"], sol_c.tree_type);
        assert_eq!(vec![2, 3, 0, 0], sol_c.tree_left);
        assert_eq!(vec![2, 4, 0, 0], sol_c.tree_right);

        let sol_g = BaseSolution::min_gen_from_wgens(2, &wx_star);
        assert_eq!(2, sol_g.objective);
        assert_eq!(vec!["Node", "Node", "Leaf", "Leaf"], sol_g.tree_type);
    }

    #[test]
    fn wgens_to_base_sol_3loci_test() {

        let x = SingleChromGenotype::from_str("011", "111");
        let z = SingleChromGenotype::from_str("111", "111");
        let g = SingleChromGamete::from_str("111");
        let wx = Rc::new(WGenS::new(x));
        let wg = Rc::new(WGamS { gamete: g, history: Some(wx) });
        let wz = Rc::new(WGenS { genotype: z, history: Some((wg.clone(), wg)) });

        let sol_c = BaseSolution::min_cross_from_wgens(3, &wz);
        assert_eq!(1, sol_c.objective);
        assert_eq!(vec!["Node", "Leaf"], sol_c.tree_type);
        assert_eq!(vec![2, 0], sol_c.tree_left);
        assert_eq!(vec![2, 0], sol_c.tree_right);

        let sol_c = BaseSolution::min_gen_from_wgens(3, &wz);
        assert_eq!(1, sol_c.objective);
        assert_eq!(vec!["Node", "Leaf"], sol_c.tree_type);
        assert_eq!(vec![2, 0], sol_c.tree_left);
        assert_eq!(vec![2, 0], sol_c.tree_right);
    }
}
