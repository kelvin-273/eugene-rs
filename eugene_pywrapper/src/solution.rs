use eugene_core::extra::resources::RecRate;
use eugene_core::plants::bit_array::SingleChromGenotype;
use eugene_core::solution::CrossingSchedule;
use num_bigint::BigUint;
use pyo3::prelude::*;

pub type PyBaseSolution = PyResult<
    Option<(
        Vec<Vec<Vec<i32>>>,
        Vec<&'static str>,
        Vec<usize>,
        Vec<usize>,
        usize,
    )>,
>;

pub type PyBaseSolutionNoObjective = PyResult<
    Option<(
        Vec<Vec<Vec<i32>>>,
        Vec<&'static str>,
        Vec<usize>,
        Vec<usize>,
    )>,
>;

#[pyclass]
pub struct PyRecRate {
    rec_rate: RecRate,
}

#[pymethods]
impl PyRecRate {
    #[new]
    fn new(rec_rate: Vec<f64>) -> Self {
        Self {
            rec_rate: RecRate::new(rec_rate),
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.rec_rate)
    }

    fn n_loci(&self) -> usize {
        self.rec_rate.n_loci()
    }

    fn cost_of_crossing(&self, gamma: f64, x: BigUint, y: BigUint, z: BigUint) -> usize {
        let n_loci = self.rec_rate.n_loci();
        let x = genotype_from_biguint(x, n_loci);
        let y = genotype_from_biguint(y, n_loci);
        let z = genotype_from_biguint(z, n_loci);
        self.rec_rate.crossing_resources(gamma, &x, &y, &z)
    }
}

fn genotype_from_biguint(x: BigUint, n_loci: usize) -> SingleChromGenotype {
    let mask = (BigUint::from(1u64) << n_loci) - 1u64;
    let x1 = (&x >> n_loci) & mask.clone();
    let x2 = x & mask;

    let mut v = vec![(false, false); n_loci];
    for (i, a) in x1.iter_u32_digits().enumerate() {
        for j in 0..32 {
            if i * 32 + j >= n_loci {
                break;
            }
            v[i * 32 + j].0 = (a >> j) & 1 == 1;
        }
    }
    for (i, b) in x2.iter_u32_digits().enumerate() {
        for j in 0..32 {
            if i * 32 + j >= n_loci {
                break;
            }
            v[i * 32 + j].1 = (b >> j) & 1 == 1;
        }
    }
    v.reverse();
    SingleChromGenotype::new(v)
}

fn gamete_from_biguint(x: BigUint, n_loci: usize) -> SingleChromGenotype {
    let mask = (BigUint::from(1u64) << n_loci) - 1u64;
    let x = x & mask;
    let mut v: Vec<bool> = x
        .iter_u32_digits()
        .flat_map(|a| (0..32).map(move |i| (a >> i) & 1 == 1))
        .take(n_loci)
        .collect();
    v.reverse();
    SingleChromGenotype::new(v.into_iter().map(|b| (b, b)).collect())
}

#[pyclass]
pub struct PyCrossingSchedule {
    crossing_schedule: CrossingSchedule,
}

impl PyCrossingSchedule {
    pub fn new(crossing_schedule: CrossingSchedule) -> Self {
        Self { crossing_schedule }
    }
}

#[pymethods]
impl PyCrossingSchedule {
    fn to_base_solution_no_obj(&self) -> PyBaseSolutionNoObjective {
        Ok(Some((&self.crossing_schedule).into()))
    }

    fn generations(&self) -> usize {
        self.crossing_schedule.generations()
    }

    fn crossings(&self) -> usize {
        self.crossing_schedule.crossings()
    }

    fn resources(&self, rec_rate: &PyRecRate, gamma: f64) -> usize {
        self.crossing_schedule.resources(&rec_rate.rec_rate, gamma)
    }
}
