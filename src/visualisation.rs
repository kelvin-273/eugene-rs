use std::rc::Rc;
use std::collections::HashMap;
use std::os::unix::prelude::*;
use std::fs;
use std::io;
use svg::Document;
use svg::node::element::Rectangle;
use svg::node::element::Path;
use svg::node::element::Group;
use svg::node::element::path::Data;

use crate::abstract_plants::*;

fn allele_colour(allele: Allele) -> String {
    match allele {
        Allele::O => "yellow".to_owned(),
        Allele::Z => "blue".to_owned()
    }
}

pub fn draw_gen_base<A, B>(x: &A, offset: (i32, i32)) -> Group where
    A: Diploid<B>,
    B: BioSize + Haploid
{
    Group::new()
        .add(draw_gam_base::<A, B>(&x.upper(), (0, 0)))
        .add(draw_gam_base::<A, B>(&x.lower(), (0, 10)))
}

pub fn draw_gam_base<A, B>(x: &B, offset: (i32, i32)) -> Group where
    B: BioSize + Haploid
{

    let g = Group::new();
    let size = x.get_sizes();
    x.alleles().iter().enumerate().map(|(i, a)| {
        Rectangle::new()
            .set("width", 10)
            .set("height", 10)
            .set("x", (i as i32)*10 + offset.0)
            .set("y", (i as i32)*0 + offset.1)
            .set("stroke", "black")
            .set("fill", allele_colour(*a))
    }).fold(g, |d, square| {
        d.add(square)
    })
}

pub fn draw_base<A, B>(filename: &str, x: &A) where
    A: BioSize + Diploid<B>,
    B: BioSize + Haploid
{
    let size = x.get_sizes();
    let doc = Document::new()
        // TODO: replace this a viewBox that covers everything //
        .set("viewBox", (0, 0, size[0].1 * 10 + 1, 2*10 + 1));

    let out = doc.add(draw_gen_base(x, (0, 0)));
    svg::save(filename, &out);
}

pub fn draw_tree_genotypes_base<A, B>(tree: Rc<WGenS<A, B>>) -> io::Result<()> where
    A: BioSize + Diploid<B>,
    B: BioSize + Haploid
{
    // collect all the unique genotypes a hashmap of all the unique genotypes
    let mut h: HashMap<A, usize> = HashMap::new();
    fn aux<A, B>(node: Rc<WGenS<A, B>>, h: &mut HashMap<A, usize>) {
        unimplemented!();
    }
    aux(tree.clone(), &mut h);

    // create a folder in tmp
    fs::create_dir("/tmp/tree-image")?;
    //
    // for every genotype create an svg for the genotype
    for (key, val) in h.iter() {
        draw_base(format!("/tmp/tree-image/node{}.svg", val).as_str(), key);
    }

    // traverse the tree collecting edges from genotype to genotype
    // create dotfile with nodes containing links to the svg files
    // run dot
    unimplemented!();
}
