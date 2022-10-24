use eugene::*;
use bitvec::prelude::*;

fn main() {
    println!("Hello, world!");
    //let x = bits![1, 0, 0, 1, 1, 1];
    let x = bits![1,0,1,0,1,1];
    println!("{:?}", x);
    //println!("{:?}", x[3]);
    //println!("{:?}", x[1]);
    //println!("{:?}", x[2]);
    x.iter().for_each(|i| println!("{}", i))
}
