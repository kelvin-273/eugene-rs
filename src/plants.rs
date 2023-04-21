pub mod u64_and_u32;

pub mod vec_of_bools;

pub mod bigint;

pub mod bit_array;

#[cfg(test)]
mod tests {
    use super::*;
    use rand::distributions::{Distribution, Uniform};
    use rand::prelude::*;

    #[test]
    fn vec_of_bools_random_plant_generates_gamete_with_alleles_from_self() {
        use vec_of_bools::*;
        use crate::abstract_plants::Crosspoint;
        let mut rng = thread_rng();
        let n_loci_distribution = Uniform::from(1..10);
        for _ in (0..1000) {
            //let n_loci = random.randrange(1, 10)
            let n_loci = n_loci_distribution.sample(&mut rng);
            let n_loci_2 = n_loci_distribution.sample(&mut rng);
            //x = plant.generate_random_plant(n_loci)
            let x: SingleChromGenotype;
            //let g: SingleChromGamete = ;
            //cs = format(x.create_gamete(), f"0{n_loci}b")
            //for i in range(n_loci):
                //assert cs[i] == su[i] or cs[i] == sl[i]
        }
    }
}
