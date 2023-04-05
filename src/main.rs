use std::io;
use std::rc::Rc;
use eugene::enumerator_dominance;
use eugene::visualisation;

pub fn main() -> io::Result<()> {
    _main2()
}

pub fn _main2() -> io::Result<()> {
    use eugene::plant_implementations::bit_array::*;
    let res = enumerator_dominance::breeding_program_timeout_gametes(
        SingleChromGenotype::init_pop_random(13, 6)
        );
    Ok(())
}

pub fn _main1() -> io::Result<()> {
    use eugene::plant_implementations::vec_of_bools::*;
    let n_loci = 8;
    let pop0: Vec<SingleChromGenotype> = vec![
        SingleChromGenotype::from_str(
            "01010101",
            "10101010",
        )
    ];
    if let Some(res) = enumerator_dominance::breeding_program::<
        SingleChromGenotype,
        SingleChromGamete,
        CrosspointSingleVob,
        WeakDomSingle,
        usize,
    >(pop0, SingleChromGenotype::ideotype(n_loci), n_loci) {
        print!("breeding_program successful");

        visualisation::draw_tree_genotypes(Rc::new(res).extract_first())?;
    }
    Ok(())
}

fn _main0() -> io::Result<()> {
    use eugene::plant_implementations::vec_of_bools::*;
    let n_loci = 10;
    let pop0: Vec<SingleChromGenotype> = vec![
        SingleChromGenotype::from_str(
            "0000011000",
            "0000000000",
            ),
        SingleChromGenotype::from_str(
            "0000000100",
            "0000000000",
            ),
        SingleChromGenotype::from_str(
            "0000000101",
            "0000000010",
            ),
        SingleChromGenotype::from_str(
            "1100000000",
            "0000000000",
            ),
        SingleChromGenotype::from_str(
            "0010000000",
            "0000000000",
            ),
        SingleChromGenotype::from_str(
            "0010100000",
            "0001000000",
            ),
    ];
    if let Some(res) = enumerator_dominance::breeding_program::<
        SingleChromGenotype,
        SingleChromGamete,
        CrosspointSingleVob,
        WeakDomSingle,
        usize,
    >(pop0, SingleChromGenotype::ideotype(n_loci), n_loci) {
        print!("breeding_program successful");

        visualisation::draw_tree_genotypes(Rc::new(res).extract_first())?;
    }
    Ok(())
}
