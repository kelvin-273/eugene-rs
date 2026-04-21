use bit_vec::BitVec;
use eugene_core::extra::resources::{RecRate, SinglePointRecProb};
use eugene_core::plants::bit_array::{SingleChromGamete, SingleChromGenotype};
use eugene_core::solution::{CrossingSchedule, TreeType};
use num_bigint::BigUint;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

type PyCrossingScheduleState = (
    usize,
    Vec<Vec<Vec<i32>>>,
    Vec<String>,
    Vec<usize>,
    Vec<usize>,
);

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

#[pyclass]
pub struct PySinglePointRecProb {
    rec_prob: SinglePointRecProb,
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

#[pymethods]
impl PySinglePointRecProb {
    #[new]
    fn new(rec_prob: Vec<f64>) -> Self {
        Self {
            rec_prob: SinglePointRecProb::new(rec_prob),
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.rec_prob)
    }

    fn n_loci(&self) -> usize {
        self.rec_prob.n_loci()
    }

    fn cost_of_crossing(&self, gamma: f64, x: BigUint, y: BigUint, z: BigUint) -> usize {
        let n_loci = self.rec_prob.n_loci();
        let x = genotype_from_biguint(x, n_loci);
        let y = genotype_from_biguint(y, n_loci);
        let z = genotype_from_biguint(z, n_loci);
        self.rec_prob.crossing_resources(gamma, &x, &y, &z)
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

fn gamete_from_biguint(x: BigUint, n_loci: usize) -> SingleChromGamete {
    let mask = (BigUint::from(1u64) << n_loci) - 1u64;
    let x = x & mask;
    let mut v: Vec<bool> = x
        .iter_u32_digits()
        .flat_map(|a| (0..32).map(move |i| (a >> i) & 1 == 1))
        .take(n_loci)
        .collect();
    v.reverse();
    SingleChromGamete::new(&v)
}

#[pyclass(module = "eugene_pywrapper")]
#[derive(Debug)]
pub struct PyCrossingSchedule {
    crossing_schedule: CrossingSchedule,
}

impl PyCrossingSchedule {
    pub fn new(crossing_schedule: CrossingSchedule) -> Self {
        Self { crossing_schedule }
    }

    fn state(&self) -> PyCrossingScheduleState {
        let (tree_data, tree_type, tree_left, tree_right) = (&self.crossing_schedule).into();
        (
            self.crossing_schedule.n_loci(),
            tree_data,
            tree_type.into_iter().map(str::to_owned).collect(),
            tree_left,
            tree_right,
        )
    }

    fn from_state(state: PyCrossingScheduleState) -> PyResult<Self> {
        let (n_loci, tree_data, tree_type, tree_left, tree_right) = state;
        let tree_data = tree_data
            .into_iter()
            .enumerate()
            .map(|(genotype_idx, genotype)| {
                let genotype = genotype
                    .into_iter()
                    .enumerate()
                    .map(|(chromosome_idx, chromosome)| {
                        if chromosome.len() != n_loci {
                            return Err(PyValueError::new_err(format!(
                                "chromosome {chromosome_idx} in genotype {genotype_idx} has length {}, expected {n_loci}",
                                chromosome.len()
                            )));
                        }
                        chromosome
                            .into_iter()
                            .map(|allele| match allele {
                                0 => Ok(false),
                                1 => Ok(true),
                                _ => Err(PyValueError::new_err(format!(
                                    "allele values must be 0 or 1, got {allele}"
                                ))),
                            })
                            .collect::<PyResult<BitVec>>()
                    })
                    .collect::<PyResult<Vec<_>>>()?;
                let genotype: [BitVec; 2] = genotype.try_into().map_err(|_| {
                    PyValueError::new_err(format!(
                        "genotype {genotype_idx} must contain exactly two chromosomes"
                    ))
                })?;
                Ok(genotype)
            })
            .collect::<PyResult<Vec<_>>>()?;
        let tree_type = tree_type
            .into_iter()
            .map(|typ| match typ.as_str() {
                "Node" => Ok(TreeType::Node),
                "Leaf" => Ok(TreeType::Leaf),
                _ => Err(PyValueError::new_err(format!(
                    "tree type must be 'Node' or 'Leaf', got {typ:?}"
                ))),
            })
            .collect::<PyResult<Vec<_>>>()?;

        Ok(Self {
            crossing_schedule: CrossingSchedule::new(
                n_loci, tree_data, tree_type, tree_left, tree_right,
            ),
        })
    }
}

#[pymethods]
impl PyCrossingSchedule {
    #[new]
    fn py_new(state: PyCrossingScheduleState) -> PyResult<Self> {
        Self::from_state(state)
    }

    fn to_base_solution_no_obj(&self) -> PyBaseSolutionNoObjective {
        Ok(Some((&self.crossing_schedule).into()))
    }

    fn generations(&self) -> usize {
        self.crossing_schedule.generations()
    }

    fn crossings(&self) -> usize {
        self.crossing_schedule.crossings()
    }

    fn resources(&self, rec_rate: Vec<f64>, gamma: f64) -> usize {
        self.crossing_schedule
            .resources(&RecRate::new(rec_rate), gamma)
    }

    fn resources_single_point(&self, rec_prob: Vec<f64>, gamma: f64) -> usize {
        self.crossing_schedule
            .resources(&SinglePointRecProb::new(rec_prob), gamma)
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.crossing_schedule)
    }

    fn __getnewargs__(&self) -> (PyCrossingScheduleState,) {
        (self.state(),)
    }

    fn __getstate__(&self) -> PyCrossingScheduleState {
        self.state()
    }

    fn __setstate__(&mut self, state: PyCrossingScheduleState) -> PyResult<()> {
        *self = Self::from_state(state)?;
        Ok(())
    }
}
