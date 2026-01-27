use std::ops::Range;

use eugene_core::plants::dist_array;
use eugene_core::plants::dist_array::{dist_array, DistArray, OrderedDistArray};
use eugene_core::solvers::base_min_crossings_distribute_astar::{
    breeding_program_distribute_general, Config,
};
use btree_range_map::{AsRange, IntoRange, RangeMap, Measure};
use range_traits::MaybeBounded;
use serde::{Serialize, Deserialize};

/// A wrapper around DistArray to implement Ord and Measure traits.
///
//// This struct provides a way to compare DistArray instances and calculate distances
/// between them based on their ordinal indices.
///
/// # Examples
/// ```
/// let da1 = dist_array![0, 1, 0, 2, 1];
/// let da2 = dist_array![0, 1, 2, 0];
/// let wrapper1 = DistArrayWrapper::new(da1);
/// let wrapper2 = DistArrayWrapper::new(da2);
///
/// assert!(wrapper1 > wrapper2);
/// assert_eq!(wrapper1.distance(&wrapper2), 5);
/// ```
#[derive(Debug, PartialEq, Eq, Hash, Clone, Serialize, Deserialize)]
struct DistArrayWrapper(DistArray);

impl DistArrayWrapper {
    fn new(da: DistArray) -> Self {
        DistArrayWrapper(da)
    }

    pub fn to_ordinal(&self) -> usize {
        dist_array::ordinal_index(&self.0)
    }
}

impl std::ops::Deref for DistArrayWrapper {
    type Target = DistArray;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for DistArrayWrapper {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl PartialOrd for DistArrayWrapper {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        // Assumes that dist arrays are in canonical form
        if self.0.len() != other.0.len() {
            Some(self.0.len().cmp(&other.0.len()))
        } else {
            Some(self.0.cmp(&other.0))
        }
    }
}

impl Ord for DistArrayWrapper {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

impl Measure for DistArrayWrapper {
    type Len = usize;

    fn len(&self) -> usize {
        1
    }

    fn distance(&self, other: &Self) -> usize {
        if self.n_loci() > other.n_loci() || (self.n_loci() == other.n_loci() && self >= other) {
            other.distance(self)
        } else {
            let ordinal_self = self.to_ordinal();
            let ordinal_other = other.to_ordinal();
            let mut dist = 0;
            for n_loci in self.n_loci()..other.n_loci() {
                dist += dist_array::count_distribute_instances(n_loci);
            }
            dist + (ordinal_other - ordinal_self)
        }
    }
}

enum DistArrayRange {
    Empty,
    UpperBounded(DistArrayWrapper),
    LowerBounded(DistArrayWrapper),
    DoubleBounded(DistArrayWrapper, DistArrayWrapper),
}

fn f() -> Range<DistArrayWrapper> {
    let start = DistArrayWrapper::new(dist_array![0, 1, 0, 1]);
    let end = DistArrayWrapper::new(dist_array![0, 1, 2, 3]);
    start..end
}

impl range_traits::Enum for DistArrayWrapper {
    fn succ(&self) -> Option<Self> {
        self.successor().map(DistArrayWrapper)
    }

    fn pred(&self) -> Option<Self> {
        self.predecessor().map(DistArrayWrapper)
    }
}


/// Represents the result of a completed instance comparison between HEU and OPT.
///
/// # Variants
/// - `HeuEqualsOpt`: Indicates that the heuristic solution equals the optimal solution.
/// - `HeuNotEqualsOpt`: Indicates that the heuristic solution does not equal the optimal
/// solution.
/// - `Timeout`: Indicates that the comparison resulted in a timeout.
///
/// # Examples
/// ```
/// let result = CompletedInstanceResult::HeuEqualsOpt;
/// assert_eq!(result.to_bool(), Some(true));
///
/// let result = CompletedInstanceResult::HeuNotEqualsOpt;
/// assert_eq!(result.to_bool(), Some(false));
///
/// let result = CompletedInstanceResult::Timeout;
/// assert_eq!(result.to_bool(), None);
///
/// let result = CompletedInstanceResult::from_bool(true);
/// assert_eq!(result, CompletedInstanceResult::HeuEqualsOpt);
///
/// let result = CompletedInstanceResult::from_bool(false);
/// assert_eq!(result, CompletedInstanceResult::HeuNotEqualsOpt);
/// ```
#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash, Serialize, Deserialize)]
enum CompletedInstanceResult {
    HeuEqualsOpt,
    HeuNotEqualsOpt,
    Timeout,
}

impl CompletedInstanceResult {
    /// Creates a CompletedInstanceResult from a `bool`.
    ///
    /// # Arguments
    /// - `b`: A boolean indicating whether HEU equals OPT.
    ///
    /// # Returns
    /// - `CompletedInstanceResult`: The corresponding CompletedInstanceResult.
    ///
    /// # Examples
    /// ```
    /// let result = CompletedInstanceResult::from_bool(true);
    /// assert_eq!(result, CompletedInstanceResult::HeuEqualsOpt);
    ///
    /// let result = CompletedInstanceResult::from_bool(false);
    /// assert_eq!(result, CompletedInstanceResult::HeuNotEqualsOpt);
    /// ```
    pub fn from_bool(b: bool) -> Self {
        if b {
            CompletedInstanceResult::HeuEqualsOpt
        } else {
            CompletedInstanceResult::HeuNotEqualsOpt
        }
    }

    /// Converts the result to an `Option<bool>`.
    ///
    /// # Returns
    /// - `Some(true)`: if HEU equals OPT.
    /// - `Some(false)`: if HEU not equals OPT.
    /// - `None`: if the result is a timeout.
    ///
    /// # Examples
    /// ```
    /// let result = CompletedInstanceResult::HeuEqualsOpt;
    /// assert_eq!(result.to_bool(), Some(true));
    ///
    /// let result = CompletedInstanceResult::HeuNotEqualsOpt;
    /// assert_eq!(result.to_bool(), Some(false));
    ///
    /// let result = CompletedInstanceResult::Timeout;
    /// assert_eq!(result.to_bool(), None);
    /// ```
    pub fn to_bool(&self) -> Option<bool> {
        match self {
            CompletedInstanceResult::HeuEqualsOpt => Some(true),
            CompletedInstanceResult::HeuNotEqualsOpt => Some(false),
            CompletedInstanceResult::Timeout => None,
        }
    }

    /// Returns true if the result indicates a timeout.
    ///
    /// # Returns
    /// - `bool`: True if the result is a timeout, false otherwise.
    ///
    /// # Examples
    /// ```
    /// let result = CompletedInstanceResult::Timeout;
    /// assert!(result.is_timeout());
    /// ```
    pub fn is_timeout(&self) -> bool {
        matches!(self, CompletedInstanceResult::Timeout)
    }

    /// Returns true if the result indicates HEU equals OPT.
    ///
    /// # Returns
    /// - `bool`: True if the result is HEU equals OPT, false otherwise.
    ///
    /// # Examples
    /// ```
    /// let result = CompletedInstanceResult::HeuEqualsOpt;
    /// assert!(result.is_heu_equals_opt());
    /// ```
    pub fn is_heu_equals_opt(&self) -> bool {
        matches!(self, CompletedInstanceResult::HeuEqualsOpt)
    }

    /// Returns true if the result indicates HEU not equals OPT.
    ///
    /// # Returns
    /// - `bool`: True if the result is HEU not equals OPT, false otherwise.
    ///
    /// # Examples
    /// ```
    /// let result = CompletedInstanceResult::HeuNotEqualsOpt;
    /// assert!(result.is_heu_not_equals_opt());
    /// ```
    pub fn is_heu_not_equals_opt(&self) -> bool {
        matches!(self, CompletedInstanceResult::HeuNotEqualsOpt)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CompletedInstances {
    btree: RangeMap<DistArrayWrapper, CompletedInstanceResult>,
}

impl CompletedInstances {
    pub fn new() -> Self {
        CompletedInstances {
            btree: RangeMap::new(),
        }
    }

    pub fn get(&self, da: &DistArrayWrapper) -> Option<CompletedInstanceResult> {
        self.btree.get(da.clone()).copied()
    }

    pub fn insert(&mut self, da: DistArrayWrapper, result: bool) {
        let res_enum = CompletedInstanceResult::from_bool(result);
        self.btree.insert(da, res_enum);
    }
}

pub fn main() {
    let mut db = CompletedInstances::new();
    let conf_opt = Config::new(false, false, false, None);
    let conf_heu = Config::new(true, false, false, None);

    for n_loci in 2.. {
        for instance in dist_array::DistArrayGenerator::new(n_loci) {
            let xs: Vec<_> = instance.into();

            print!("Solving instance: {:?}", xs);

            // Reverse check
            let rev_xs: OrderedDistArray =
                dist_array::canonical_dist_array(&xs.iter().rev().collect()).into();
            if OrderedDistArray(xs.clone().into()) <= rev_xs {
                print!(" Already completed in reverse. Skipping.\n");
                continue;
            }
            if let Some(result) = db.get(&DistArrayWrapper(
                OrderedDistArray(xs.clone().into()).into(),
            )) {
                print!(
                    " Already completed with result: {}. Skipping.\n",
                    match result {
                        CompletedInstanceResult::HeuEqualsOpt => "HEU=OPT",
                        CompletedInstanceResult::HeuNotEqualsOpt => "HEU!=OPT",
                        CompletedInstanceResult::Timeout => "TIMEOUT",
                    }
                );
                continue;
            }

            let Some(res_heu) = breeding_program_distribute_general(&xs, &conf_heu) else {
                println!(" HEU returned None");
                break;
            };
            let Some(res_opt) = breeding_program_distribute_general(&xs, &conf_opt) else {
                println!(" OPT returned None");
                break;
            };

            if !(res_heu.objective() == res_opt.objective()) {
                println!(" HEU!=OPT");
            } else {
                println!(" HEU=OPT");
            }
            let res = db.get(&DistArrayWrapper({
                let rev = xs.iter().rev().cloned().collect::<Vec<usize>>();
                dist_array::canonical_dist_array(&DistArray::from(rev))
            }));
            db.insert(
                DistArrayWrapper(OrderedDistArray(xs.into()).into()),
                res_heu.objective() == res_opt.objective()
            );
        }
    }
}
