#![allow(unused)]

pub mod abstract_plants;
pub mod enumerator;
pub mod enumerator_dominance;
pub mod greedy_base;
pub mod greedy_multichrom_ti;
pub mod visualisation;
pub mod plant_implementations;

mod play {
    #[cfg(test)]
    mod tests {

        #[test]
        fn masking_test() {
            use bitvec::prelude::*;
            let x = bitvec![3, 5, 1, 3, 1, 3, 4, 6, 1, 3, 44];
            println!("{:?}", x);
            let y: u8 = 0b0110_1011;
            let _z: BitArray<u8, Lsb0> = BitArray::from(y);
        }
    }
}
