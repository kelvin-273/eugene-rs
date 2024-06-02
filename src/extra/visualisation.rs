use std::collections::HashMap;
use std::fs;
use std::hash::Hash;
use std::io;
use std::io::prelude::*;
use svg::node::element::Group;
use svg::node::element::Rectangle;
use svg::Document;

use crate::abstract_plants::*;
use crate::solution::BaseSolution;

fn allele_colour(allele: Allele) -> String {
    match allele {
        Allele::O => "yellow".to_owned(),
        Allele::Z => "blue".to_owned(),
    }
}

pub fn draw_gen_base<A, B>(x: &A, offset: (i32, i32)) -> Group
where
    A: Diploid<B>,
    B: BioSize + Haploid,
{
    Group::new()
        .add(draw_gam_base::<B>(&x.upper(), (offset.0, offset.1)))
        .add(draw_gam_base::<B>(&x.lower(), (offset.0, offset.1 + 10)))
}

pub fn draw_gam_base<B>(x: &B, offset: (i32, i32)) -> Group
where
    B: BioSize + Haploid,
{
    let g = Group::new();
    //let size = x.get_sizes();
    x.alleles()
        .iter()
        .enumerate()
        .map(|(i, a)| {
            Rectangle::new()
                .set("width", 10)
                .set("height", 10)
                .set("x", (i as i32) * 10 + offset.0)
                .set("y", (i as i32) * 0 + offset.1)
                .set("stroke", "black")
                .set("fill", allele_colour(*a))
        })
        .fold(g, |d, square| d.add(square))
}

pub fn draw_base<A, B>(filename: &str, x: &A) -> io::Result<()>
where
    A: BioSize + Diploid<B>,
    B: BioSize + Haploid,
{
    //let doc = Document::new()
    //    // TODO: replace this a viewBox that covers everything //
    //    .set("viewBox", (0, 0, size[0].1 * 10 + 1, 2*10 + 1));
    let doc = Document::new();

    let out = doc.add(draw_gen_base(x, (0, 0)));
    svg::save(filename, &out)
}

pub const BLOCKSIZE: usize = 10;
pub trait Draw {
    /// Returns (width, height) of the viewBox needed to contain self.
    fn view_box_size(&self) -> Option<(usize, usize)>;

    fn draw(&self) -> Group;

    fn draw_to_file(&self, s: &str) -> io::Result<()> {
        let out = match self.view_box_size() {
            None => Document::new().add(self.draw()),
            Some((w, h)) => Document::new()
                .set("viewBox", (0, 0, w, h))
                .add(self.draw()),
        };
        svg::save(s, &out)
    }
}

pub fn draw_tree_genotypes<A, B>(tree: &BaseSolution) -> io::Result<()>
where
    A: Draw + Eq + Hash + Clone,
    B: Draw,
{
    unimplemented!()
}

pub fn draw_tree_gametes<A, B>(tree: &BaseSolution) -> io::Result<()>
where
    B: Draw + Eq + Hash + Clone,
{
    unimplemented!()
}
