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
                .set("width", BLOCKSIZE)
                .set("height", BLOCKSIZE)
                .set("x", (i as i32) * BLOCKSIZE as i32 + offset.0)
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

pub fn draw_individual_genotypes<A, B>(pop: &Vec<WGen<A, B>>) -> Group
where
    A: Genotype<B> + Diploid<B>,
    B: Gamete<A> + Haploid,
{
    pop.iter()
        .enumerate()
        .fold((0, Group::new()), |(total_rows, g_curr), (i, x)| {
            let x_group = [
                x.genotype().upper().alleles(),
                x.genotype().lower().alleles(),
            ]
            .iter()
            .enumerate()
            .fold(Group::new(), |x_group, (j, row)| {
                x_group.add(
                    row.iter()
                        .enumerate()
                        .fold(Group::new(), |row_group, (p, a)| {
                            row_group.add(
                                Rectangle::new()
                                    .set("width", BLOCKSIZE)
                                    .set("height", BLOCKSIZE)
                                    .set("x", BLOCKSIZE * p as usize)
                                    .set("y", BLOCKSIZE * (total_rows + i))
                                    .set("stroke", "black")
                                    .set("fill", allele_colour(*a)),
                            )
                        }),
                )
            });
            (total_rows + 2 + 1, g_curr.add(x_group))
        })
        .1
}

use crate::plants::bit_array::SingleChromGenotype;

pub fn create_random_selection(
    n_loci: usize,
    pop_0: Vec<SingleChromGenotype>,
    filename: Option<&str>,
    k_cross: Option<usize>,
) {
    use crate::solvers::base_heuristic_random_selection as ehue;
    let filename = filename.unwrap_or("/tmp/output.svg");
    let mut rng = rand::thread_rng();
    let output_pop = ehue::repeated_breeding_random(n_loci, &pop_0, k_cross, &mut rng);
    let output_svg = draw_individual_genotypes(&output_pop);
    let _ = svg::save(filename, &output_svg);
}

pub fn create_random_selection_dominance(
    n_loci: usize,
    pop_0: Vec<SingleChromGenotype>,
    filename: Option<&str>,
    k_cross: Option<usize>,
) {
    use crate::solvers::base_heuristic_random_selection as ehue;
    let filename = filename.unwrap_or("/tmp/output.svg");
    let mut rng = rand::thread_rng();
    let output_pop = ehue::repeated_breeding_random_dominance(n_loci, &pop_0, k_cross, &mut rng);
    let output_svg = draw_individual_genotypes(&output_pop);
    let _ = svg::save(filename, &output_svg);
}

#[pyo3::pyfunction]
pub fn create_random_selection_py(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<bool>>>,
    filename: Option<&str>,
    k_cross: Option<usize>,
) {
    todo!();
}

mod bezier {
    #![allow(unused)]

    #[derive(Debug, Clone, Copy, PartialEq)]
    struct Point {
        x: f64,
        y: f64,
    }

    impl Point {
        pub fn new(x: f64, y: f64) -> Self {
            Self { x, y }
        }
    }

    impl std::ops::Add for Point {
        type Output = Point;
        fn add(self, rhs: Self) -> Self::Output {
            Self {
                x: self.x + rhs.x,
                y: self.y + rhs.y,
            }
        }
    }

    impl std::ops::Sub for Point {
        type Output = Point;
        fn sub(self, rhs: Self) -> Self::Output {
            Self {
                x: self.x - rhs.x,
                y: self.y - rhs.y,
            }
        }
    }

    impl std::ops::Mul<f64> for Point {
        type Output = Point;
        fn mul(self, rhs: f64) -> Self::Output {
            Self {
                x: self.x * rhs,
                y: self.y * rhs,
            }
        }
    }

    struct Bezier {
        p: (Point, Point, Point, Point),
    }

    impl Bezier {
        pub fn new(p0: Point, p1: Point, p2: Point, p3: Point) -> Self {
            Self {
                p: (p0, p1, p2, p3),
            }
        }

        fn evaluate_y(&self, t: f64) -> Point {
            self.p.0 * (1.0 - t).powi(3)
                + self.p.1 * (1.0 - t).powi(2) * t
                + self.p.2 * (1.0 - t).powi(1) * t.powi(2)
                + self.p.3 * t.powi(3)
        }

        fn evaluate_dy(self, t: f64) -> Point {
            (self.p.1 - self.p.0) * (1.0 - t).powi(2) * 3.0
                + (self.p.2 - self.p.1) * (1.0 - t) * t * 6.0
                + (self.p.3 - self.p.2) * t.powi(2) * 3.0
        }

        fn evaluate_ddy(self, t: f64) -> Point {
            (self.p.2 - self.p.1 * 2.0 + self.p.0) * (1.0 - t) * 6.0
                + (self.p.3 - self.p.2 * 2.0 + self.p.1) * t * 6.0
        }

        /// Returns true if the two Bézier curves intersect
        fn intersects_with(&self, other: &Self) -> bool {
            fn bounding_box_int(lhs: &Bezier, rhs: &Bezier) -> bool {
                unimplemented!()
            }
            unimplemented!()
        }

        fn bisect(&self) -> (Self, Self) {
            unimplemented!()
        }
    }
}
