use crate::solution::{PyBaseSolution, PyCrossingSchedule};
use eugene_core::solvers::base_min_crossings_distribute_astar;
use eugene_core::solvers::base_min_crossings_distribute_astar::Config;
use pyo3::PyResult;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

#[derive(pyo3::IntoPyObject, Debug)]
pub struct Output {
    expansions: usize,
    pushed_nodes: usize,
    children_created_by_branching: usize,
    crossing_schedule: Option<PyCrossingSchedule>,
    objective: Option<usize>,
}

/// Computes an optimal crossing schedule given a distribute array `xs`.
///
/// For now, only the objective is computed.
#[pyo3::pyfunction]
pub fn breeding_program_distribute_minres_python(
    xs: Vec<usize>,
    timeout: Option<u64>,
) -> PyResult<Option<PyCrossingSchedule>> {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = base_min_crossings_distribute_astar::breeding_program_distribute_general(
            &xs.into(),
            &Config::new(true, true, false, None),
        )
        .map(|res| PyCrossingSchedule::new(res.crossing_schedule().unwrap().clone()));
        let _ = tx.send(res);
    });
    let res = rx
        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
        .ok()
        .flatten();
    Ok(res)
}

/// Returns the number of crossings required for distribute array `xs` along with several other
/// runtime statistics.
#[pyo3::pyfunction]
#[pyo3(signature = (xs, timeout=None, diving=false, full_join=false, dominance=false, debug_trace_file=None))]
pub fn breeding_program_distribute_general_python(
    xs: Vec<usize>,
    timeout: Option<u64>,
    diving: bool,
    full_join: bool,
    dominance: bool,
    debug_trace_file: Option<String>,
) -> pyo3::PyResult<Option<Output>> {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = base_min_crossings_distribute_astar::breeding_program_distribute_general(
            &xs.into(),
            &Config::new(full_join, dominance, diving, debug_trace_file),
        );
        let _ = tx.send(res);
    });
    let res = rx
        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
        .ok()
        .flatten();
    Ok(res.map(|res| Output {
        expansions: res.expansions(),
        pushed_nodes: res.pushed_nodes(),
        children_created_by_branching: res.children_created_by_branching(),
        crossing_schedule: res
            .crossing_schedule()
            .map(|cs| PyCrossingSchedule::new(cs.clone())),
        objective: res.objective(),
    }))
}

/// Computes an optimal crossing schedule given a distribute array `xs`.
///
/// For now, only the objective is computed.
#[pyo3::pyfunction]
pub fn breeding_program_distribute_python(xs: Vec<usize>, timeout: Option<u64>) -> PyBaseSolution {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = base_min_crossings_distribute_astar::breeding_program_distribute(&xs.into());
        let _ = tx.send(res);
    });
    let res = rx
        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
        .ok()
        .flatten();
    match res {
        None => Ok(None),
        Some(sol) => Ok(Some((
            sol.tree_data,
            sol.tree_type,
            sol.tree_left,
            sol.tree_right,
            sol.objective,
        ))),
    }
}

#[pyo3::pyfunction]
pub fn breeding_program_distribute_no_full_join_python(
    xs: Vec<usize>,
    timeout: Option<u64>,
) -> PyBaseSolution {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res = base_min_crossings_distribute_astar::breeding_program_distribute_no_full_join(
            &xs.into(),
        );
        let _ = tx.send(res);
    });
    let res = rx
        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
        .ok()
        .flatten();
    match res {
        None => Ok(None),
        Some(sol) => Ok(Some((
            sol.tree_data,
            sol.tree_type,
            sol.tree_left,
            sol.tree_right,
            sol.objective,
        ))),
    }
}

#[pyo3::pyfunction]
pub fn breeding_program_distribute_dominance_python(
    xs: Vec<usize>,
    timeout: Option<u64>,
) -> PyBaseSolution {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res =
            base_min_crossings_distribute_astar::breeding_program_distribute_dominance(&xs.into());
        let _ = tx.send(res);
    });
    let res = rx
        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
        .ok()
        .flatten();
    match res {
        None => Ok(None),
        Some(sol) => Ok(Some((
            sol.tree_data,
            sol.tree_type,
            sol.tree_left,
            sol.tree_right,
            sol.objective,
        ))),
    }
}

#[pyo3::pyfunction]
pub fn breeding_program_distribute_dominance_no_full_join_python(
    xs: Vec<usize>,
    timeout: Option<u64>,
) -> PyBaseSolution {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res =
            base_min_crossings_distribute_astar::breeding_program_distribute_dominance_no_full_join(
                &xs.into(),
            );
        tx.send(res)
    });
    let res = rx
        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
        .ok()
        .flatten();
    match res {
        None => Ok(None),
        Some(sol) => Ok(Some((
            sol.tree_data,
            sol.tree_type,
            sol.tree_left,
            sol.tree_right,
            sol.objective,
        ))),
    }
}

#[pyo3::pyfunction]
pub fn breeding_program_distribute_diving_python(
    xs: Vec<usize>,
    timeout: Option<u64>,
) -> PyBaseSolution {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let res =
            base_min_crossings_distribute_astar::breeding_program_distribute_diving(&xs.into());
        tx.send(res)
    });
    let res = rx
        .recv_timeout(Duration::new(timeout.unwrap_or(u64::MAX), 0))
        .ok()
        .flatten();
    match res {
        None => Ok(None),
        Some(sol) => Ok(Some((
            sol.tree_data,
            sol.tree_type,
            sol.tree_left,
            sol.tree_right,
            sol.objective,
        ))),
    }
}
