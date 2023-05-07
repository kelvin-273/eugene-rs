use eugene::solvers::*;
use eugene::visualisation;
use rand::prelude::*;
use std::io;
use std::rc::Rc;
use std::env;

pub fn main() -> io::Result<()> {
    env::set_var("RUST_BACKTRACE", "1");
    _main4()
}

pub fn _main4() -> io::Result<()> {
    use eugene::plants::bit_array::*;
    let n_loci = 20;
    let n_pop = 1024;
    let _res = greedy_base::breeding_program::<
        SingleChromGenotype,
        SingleChromGamete,
        CrosspointBitVec,
        SegmentBitVec,
        >(
        n_loci,
        SingleChromGenotype::init_pop_random(&mut thread_rng(), n_loci, n_pop),
        SingleChromGenotype::ideotype(n_loci),
        );
    Ok(())
}

pub fn _main3() -> io::Result<()> {
    use eugene::plants::bit_array::*;
    let n_loci = 11;
    let n_pop = 6;
    let _res = enumerator::breeding_program::<
        SingleChromGenotype,
        SingleChromGamete,
        CrosspointBitVec,
        usize,
        >(
        SingleChromGenotype::init_pop_random(&mut thread_rng(), n_loci, n_pop),
        SingleChromGenotype::ideotype(n_loci),
        n_loci
        );
    Ok(())
}

pub fn _main2() -> io::Result<()> {
    use eugene::plants::bit_array::*;
    let n_loci = 16;
    let n_pop = 6;
    let _res = enumerator_dominance::breeding_program_timeout_gametes::<
        SingleChromGenotype,
        SingleChromGamete,
        CrosspointBitVec,
        DomGamete,
        usize,
    >(
        SingleChromGenotype::init_pop_random(&mut thread_rng(), n_loci, n_pop),
        SingleChromGenotype::ideotype(n_loci),
        None,
        n_loci,
    );
    Ok(())
}

pub fn _main1() -> io::Result<()> {
    use eugene::plants::vec_of_bools::*;
    let n_loci = 8;
    let pop0: Vec<SingleChromGenotype> =
        vec![SingleChromGenotype::from_str("01010101", "10101010")];
    if let Some(res) = enumerator_dominance::breeding_program::<
        SingleChromGenotype,
        SingleChromGamete,
        CrosspointSingleVob,
        WeakDomSingle,
        usize,
    >(pop0, SingleChromGenotype::ideotype(n_loci), n_loci)
    {
        print!("breeding_program successful");

        visualisation::draw_tree_genotypes(Rc::new(res).extract_first())?;
    }
    Ok(())
}

fn _main0() -> io::Result<()> {
    use eugene::plants::vec_of_bools::*;
    let n_loci = 10;
    let pop0: Vec<SingleChromGenotype> = vec![
        SingleChromGenotype::from_str("0000011000", "0000000000"),
        SingleChromGenotype::from_str("0000000100", "0000000000"),
        SingleChromGenotype::from_str("0000000101", "0000000010"),
        SingleChromGenotype::from_str("1100000000", "0000000000"),
        SingleChromGenotype::from_str("0010000000", "0000000000"),
        SingleChromGenotype::from_str("0010100000", "0001000000"),
    ];
    if let Some(res) = enumerator_dominance::breeding_program::<
        SingleChromGenotype,
        SingleChromGamete,
        CrosspointSingleVob,
        WeakDomSingle,
        usize,
    >(pop0, SingleChromGenotype::ideotype(n_loci), n_loci)
    {
        print!("breeding_program successful");

        visualisation::draw_tree_genotypes(Rc::new(res).extract_first())?;
    }
    Ok(())
}
