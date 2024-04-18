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
}

fn wgen_to_base_sol(
    n_loci: usize,
    x_star: &WGen<SingleChromGenotype, SingleChromGamete>,
    obj: Objective,
) -> BaseSolution {
    println!("x_star: {:#?}", x_star);
    // get the number of crossings and cells required
    let mut q_node = VecDeque::new();
    let mut q_leaf = VecDeque::new();
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

    //println!("n_cells: {:#?}", n_cells);
    //println!("n_cross: {:#?}", n_cross);

    // initialise arrays
    let mut tree_data = vec![vec![vec![0; n_loci]; 2]; n_cells];
    let mut tree_type = vec!["Null"; n_cells];
    let mut tree_left = vec![0; n_cells];
    let mut tree_right = vec![0; n_cells];

    // flush index_map because nodes and leaves have mixed order
    index_map.drain();
    index_map.insert(x_star.genotype.clone(), 0);

    if x_star.history.is_none() {
        q_leaf.push_back(&x_star);
    } else {
        q_node.push_back(&x_star);
    }
    let mut i_node = 0;
    let mut i_leaf = n_cross;

    // using q_node to process all genotypes for the counting
    while let Some(z) = q_node.pop_front() {
        //println!("z: {:#?}", z.genotype);
        //println!("z: {:#?}", z);
        let &i_z = index_map.get(&z.genotype).unwrap();
        //println!("i_z: {:#?}", i_z);

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
            //println!("x: {:#?}", x.genotype);
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
            tree_left[i_z] = i_x;

            let y = wgy.history.first().unwrap().as_ref();
            //println!("y: {:#?}", y.genotype);
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
            tree_right[i_z] = i_y;
        } else {
            tree_type[i_z] = "Leaf";
        }
    }

    println!("tree_type:\t{:?}", tree_type);
    println!("tree_left:\t{:?}", tree_left);
    println!("tree_right:\t{:?}", tree_right);
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
                //println!("i: {:#?}", i);
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
