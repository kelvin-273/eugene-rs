use eugene_core::plants::dist_array::{DistArray, OrderedDistArray};
use eugene_core::solvers::base_min_crossings_distribute_astar::{
    breeding_program_distribute_general, Config,
};
use eugene_core::{extra::instance_generators, plants::dist_array};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
struct CompletedInstances {}

impl CompletedInstances {
    pub fn new() -> Self {
        Self {}
    }

    pub fn get(&self, _instance: &Vec<usize>) -> Option<bool> {
        unimplemented!()
    }

    pub fn _insert(&mut self, _instance: &Vec<usize>, _value: bool) {
        unimplemented!()
    }
}

pub fn main() {
    let mut db = CompletedInstances::new();
    let conf_opt = Config::new(false, false, false, None);
    let conf_heu = Config::new(true, false, false, None);

    for n_loci in 2.. {
        for instance in dist_array::DistArrayGenerator::new(n_loci) {
            let xs: Vec<_> = instance.into();

            print!("Solving instance: {:?}", xs);

            // Reverse check
            let rev_xs: DistArray = dist_array::canonical_dist_array(&xs.iter().rev().collect());
            if OrderedDistArray(xs.clone().into()) <= OrderedDistArray(rev_xs) {
                print!(" Already completed in reverse. Skipping.\n");
                continue;
            }
            if let Some(completed) = db.get(&xs) {
                print!(
                    " Already completed with result: {}. Skipping.\n",
                    if completed { "HEU=OPT" } else { "HEU!=OPT" }
                );
                continue;
            }

            let res_heu = breeding_program_distribute_general(&xs, &conf_heu);
            print!(".");
            if res_heu.is_none() {
                println!(" HEU returned None");
                break;
            }
            print!(".");
            let res_opt = breeding_program_distribute_general(&xs, &conf_opt);
            if res_opt.is_none() {
                println!(" OPT returned None");
                break;
            }
            println!("");
        }
    }
}
