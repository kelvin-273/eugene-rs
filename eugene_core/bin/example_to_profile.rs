use eugene_core::solvers::base_min_crossings_distribute_astar;

pub fn main() {
    let v = vec![0, 1, 0, 2, 3, 1, 2, 0, 3, 1, 4, 2, 0, 3, 1, 2, 4];
    println!("Testing {:?}", &v);
    println!("Length: {}", v.len());
    let res = base_min_crossings_distribute_astar::breeding_program_distribute(&v);
    println!("{:?}", res);
}
