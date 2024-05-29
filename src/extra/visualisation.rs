use std::collections::HashMap;
use std::fs;
use std::hash::Hash;
use std::io;
use std::io::prelude::*;
use std::rc::Rc;
use svg::node::element::Group;
use svg::node::element::Rectangle;
use svg::Document;

use crate::abstract_plants::*;

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

pub fn draw_tree_genotypes<A, B>(tree: Rc<WGenS<A, B>>) -> io::Result<()>
where
    A: Draw + Eq + Hash + Clone,
    B: Draw,
{
    // collect all the unique genotypes a hashmap of all the unique genotypes
    let mut h: HashMap<A, usize> = HashMap::new();
    fn aux<A, B>(node: Rc<WGenS<A, B>>, h: &mut HashMap<A, usize>)
    where
        A: Eq + Hash + Clone,
    {
        let n = h.len();
        h.entry(node.genotype.clone()).or_insert(n);
        if let Some((lg, rg)) = node.history.as_ref() {
            if let Some(g) = lg.history.clone() {
                aux(g, h)
            }
            if let Some(g) = rg.history.clone() {
                aux(g, h)
            }
        }
    }
    aux(tree.clone(), &mut h);

    // create a folder in tmp
    fs::create_dir_all("/tmp/tree-image")?;
    //
    // for every genotype create an svg for the genotype
    for (key, val) in h.iter() {
        key.draw_to_file(format!("/tmp/tree-image/node{}.svg", val).as_str())?;
    }

    // traverse the tree collecting edges from genotype to genotype
    // create dotfile with nodes containing links to the svg files
    // run dot
    //let mut v: Vec<(usize, usize)> = Vec::new();
    fn collect_edges<A, B>(
        node: Rc<WGenS<A, B>>,
        h: &HashMap<A, usize>,
        v: &mut Vec<(usize, usize)>,
    ) where
        A: Hash + Eq,
    {
        let z = h.get(&node.genotype).unwrap();
        node.history.clone().map(|(lg, rg)| {
            let g = &lg.history.clone()?.genotype;
            let x = h.get(g).unwrap();
            v.push((*x, *z));
            collect_edges(lg.history.clone()?, h, v);

            let g = &rg.history.clone()?.genotype;
            let y = h.get(g).unwrap();
            v.push((*y, *z));
            collect_edges(rg.history.clone()?, h, v);
            Some(v)
        });
    }
    let mut v = Vec::new();
    collect_edges(tree.clone(), &h, &mut v);

    // create dotviz file
    let mut buffer_dot_file = fs::File::create("/tmp/tree-image/tree.dot")?;
    buffer_dot_file.write(b"digraph G {\n")?;
    for val in h.values() {
        buffer_dot_file.write(
            format!(
                "\tx{} [label=\"\", image=\"/tmp/tree-image/node{}.svg\"]\n",
                *val, *val
            )
            .as_bytes(),
        )?;
    }
    for (a, b) in v {
        buffer_dot_file.write(format!("\tx{} -> x{}\n", a, b).as_bytes())?;
    }

    buffer_dot_file.write(b"}")?;
    Ok(())
}

pub fn draw_tree_gametes<A, B>(tree: Rc<WGenS<A, B>>) -> io::Result<()>
where
    B: Draw + Eq + Hash + Clone,
{
    // Create a hash map from gamete
    let mut h: HashMap<B, usize> = HashMap::new();
    fn aux<A, B>(node: Rc<WGamS<A, B>>, h: &mut HashMap<B, usize>)
    where
        B: Eq + Hash + Clone,
    {
        let n = h.len();
        h.entry(node.gamete.clone()).or_insert(n);
        if let Some(wx) = node.history.clone() {
            if let Some((lg, rg)) = wx.history.clone() {
                aux(lg, h);
                aux(rg, h);
            }
        }
    }
    if let Some((lg, rg)) = tree.history.as_ref() {
        aux(lg.clone(), &mut h);
        if lg.gamete != rg.gamete {
            aux(rg.clone(), &mut h);
        }
    };

    // create a folder in tmp
    fs::create_dir("/tmp/tree-image")?;
    //
    // for every genotype create an svg for the genotype
    for (key, val) in h.iter() {
        key.draw_to_file(format!("/tmp/tree-image/node{}.svg", val).as_str())?;
    }

    // traverse the tree collecting edges from genotype to genotype
    // create dotfile with nodes containing links to the svg files
    // run dot
    let mut v: Vec<(usize, usize)> = Vec::new();
    fn collect_edges<A, B>(
        node: Rc<WGamS<A, B>>,
        h: &HashMap<B, usize>,
        v: &mut Vec<(usize, usize)>,
    ) where
        B: Hash + Eq,
    {
        let gz = h.get(&node.gamete).unwrap();
        match &node.history {
            None => None,
            Some(wx) => wx.history.as_ref().map(|(lg, rg)| {
                let x = h.get(&lg.gamete).unwrap();
                v.push((*x, *gz));
                collect_edges(lg.clone(), h, v);

                let y = h.get(&rg.gamete).unwrap();
                v.push((*y, *gz));
                collect_edges(rg.clone(), h, v);
            }),
        };
    }
    match tree.history.as_ref() {
        None => {}
        Some((lg, rg)) => {
            collect_edges(lg.clone(), &h, &mut v);
            if lg.gamete != rg.gamete {
                collect_edges(rg.clone(), &h, &mut v);
            }
        }
    };

    // create dotviz file
    Ok(())
}
