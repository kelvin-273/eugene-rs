use crate::abstract_plants::*;
use crate::plants::bit_array::*;
use crate::solution::BaseSolution;
use grb::prelude::*;
use grb::Error;
use pyo3::prelude::*;

/// Runs a breeding program given `n_loci` and `pop_0` where `pop_0` is a population of single
/// chromosome diploid genotypes with `n_loci` loci.
#[pyo3::pyfunction]
pub fn breeding_program_python_mip(
    n_loci: usize,
    pop_0: Vec<Vec<Vec<bool>>>,
) -> PyResult<
    Option<(
        Vec<Vec<Vec<i32>>>,
        Vec<&'static str>,
        Vec<usize>,
        Vec<usize>,
        usize,
    )>,
> {
    let ideotype = SingleChromGenotype::ideotype(n_loci);
    let pop_0 = pop_0
        .iter()
        .map(|x| {
            SingleChromGenotype::new(
                x[0].iter()
                    .zip(x[1].iter())
                    .map(|(a, b)| (*a, *b))
                    .collect(),
            )
        })
        .collect();
    let res = breeding_program(n_loci, pop_0, ideotype);
    Ok(res.ok().map(|x| {
        (
            x.tree_data,
            x.tree_type,
            x.tree_left,
            x.tree_right,
            x.objective,
        )
    }))
}

pub fn breeding_program(
    n_loci: usize,
    pop_0: Vec<SingleChromGenotype>,
    _ideotype: SingleChromGenotype,
) -> Result<BaseSolution, Error>
where
{
    let mut m = Model::new("MinCrossings")?;

    let t_max = upper_bound(n_loci, &pop_0);
    let n_max = 1 << n_loci;

    let u: Vec<Vec<Vec<Var>>> = (0..n_max)
        .map(|_| {
            (0..n_max)
                .map(|_| {
                    (0..=t_max)
                        .map(|_| add_binvar!(m))
                        .collect::<Result<Vec<_>, Error>>()
                })
                .collect::<Result<Vec<_>, Error>>()
        })
        .collect::<Result<Vec<_>, Error>>()?;

    let v: Vec<Vec<Var>> = (0..n_max)
        .map(|_| {
            (0..t_max)
                .map(|_| add_binvar!(m))
                .collect::<Result<Vec<_>, Error>>()
        })
        .collect::<Result<Vec<_>, Error>>()?;

    let w: Vec<Vec<Vec<Var>>> = (0..n_max)
        .map(|_| {
            (0..n_max)
                .map(|_| {
                    (0..t_max)
                        .map(|_| add_binvar!(m))
                        .collect::<Result<Vec<_>, Error>>()
                })
                .collect::<Result<Vec<_>, Error>>()
        })
        .collect::<Result<Vec<_>, Error>>()?;

    // Initial population
    let mut gen_0: Vec<Vec<bool>> = vec![vec![false; n_max]; n_max];
    let f = |xc: SingleChromGamete| {
        xc.alleles().iter().fold(0, |acc, b| {
            (acc << 1)
                | (match b {
                    Allele::Z => 0,
                    Allele::O => 1,
                })
        })
    };
    for x in pop_0 {
        gen_0[f(x.upper())][f(x.lower())] = true;
    }
    for gx in 0..n_max {
        for gy in 0..n_max {
            m.add_constr("P_0", c!(u[gx][gy][0] == gen_0[gx][gy] as i16))?;
        }
    }

    // target
    m.add_constr(
        &"target",
        c!(u[(1 << n_loci) - 1][(1 << n_loci) - 1][t_max] == 1),
    )?;

    m.set_objective(
        w.iter().flatten().into_iter().flatten().grb_sum(),
        ModelSense::Minimize,
    )?;

    // A genotype can only be available if it is created
    for gx in 0..n_max {
        for gy in 0..n_max {
            for t in 0..t_max {
                m.add_constr(
                    "creation",
                    c!(w[gx][gy][t] == u[gx][gy][t + 1] - u[gx][gy][t]),
                )?;
            }
        }
    }

    // gx is available if and only if there is some x that can produce gx
    // v[gx][t] <-> exists(u[x] : gx in x)
    // ≡ -v[gx][t] <-> forall(-x : gx in x)

    // -v[gx][t] -> forall(-u[x] : gx in x)
    // ≡ v[gx][t] \/ forall(-u[x] : gx in x)
    // ≡ forall(v[gx][t] \/ -u[x] : gx in x)
    for gx in 0..n_max {
        for gy in 0..n_max {
            for gz in gen_recombined_gametes(n_loci, gx, gy) {
                for t in 0..t_max {
                    m.add_constr("existence implies gamete", c!(v[gz][t] >= u[gx][gy][t]))?;
                }
            }
        }
    }

    // v[gx][t] -> exists(u[x] : gx in x)
    // ≡ -v[gx][t] \/ exists(u[x] : gx in x)
    for gz in 0..n_max {
        let parents = gen_parent_gametes(n_loci, gz);
        for t in 0..t_max {
            m.add_constr(
                "gamete implies existence",
                c!((1 - v[gz][t]) + parents.iter().map(|(gx, gy)| u[*gx][*gy][t]).grb_sum() >= 1),
            )?;
        }
    }

    // The gametes required to create x must be available
    // in order for x to be created
    for gx in 0..n_max {
        for gy in 0..n_max {
            for t in 0..t_max {
                m.add_constr("gamete required for creation", c!(w[gx][gy][t] <= v[gx][t]))?;
                m.add_constr("gamete required for creation", c!(w[gx][gy][t] <= v[gy][t]))?;
            }
        }
    }

    m.optimize()?;

    Ok(BaseSolution {
        tree_data: vec![],
        tree_type: vec![],
        tree_left: vec![],
        tree_right: vec![],
        objective: m.get_attr(attr::ObjVal)? as usize,
    })
}

fn upper_bound(n_loci: usize, _pop_0: &Vec<SingleChromGenotype>) -> usize {
    n_loci
}

fn is_child_gamete(n_loci: usize, gx: usize, gy: usize, gz: usize) -> bool {
    let mask_0 = (1 << n_loci) - 1;
    for k in 0..n_loci {
        let mask_l = (mask_0 >> k) << k;
        let mask_r = mask_0 ^ mask_l;
        if ((gx & mask_l) | (gy & mask_r)) == gz {
            return true;
        }
        if ((gy & mask_l) | (gx & mask_r)) == gz {
            return true;
        }
    }
    false
}

fn gen_parent_gametes(n_loci: usize, gz: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let n = 1 << n_loci;
    for gx in 0..n {
        for gy in 0..n {
            if is_child_gamete(n_loci, gx, gy, gz) {
                out.push((gx, gy));
            }
        }
    }
    out
}

fn gen_recombined_gametes(n_loci: usize, gx: usize, gy: usize) -> Vec<usize> {
    let mut out = Vec::with_capacity(2 * n_loci);
    let mask_0 = (1 << n_loci) - 1;
    for k in 0..n_loci {
        let mask_l = (mask_0 >> k) << k;
        let mask_r = mask_0 ^ mask_l;
        out.push((gx & mask_l) | (gy & mask_r));
        out.push((gy & mask_l) | (gx & mask_r));
    }
    out
}
