use crate::abstract_plants::*;
use std::collections::HashSet;
use std::hash::Hash;

pub fn generations<A, B>(wz: &WGenS<A, B>) -> usize {
    // This is cursed.
    (|| {
        let (gx, gy) = wz.history?;
        let wx = gx.history?.as_ref();
        let wy = gy.history?.as_ref();
        Some(generations(wx).max(generations(wy)) + 1)
    })().unwrap_or(0)
}

#[derive(Hash, PartialEq, Eq)]
struct Triple<A> {
    x: A,
    y: A,
    z: A,
}

/// Computes the number of crossings in a crossing schedule.
///
/// Assumes that each crossing can only be created once, i.e. only one (x, y, z) ∀x,y,z.
pub fn crossings<A, B>(wz: &WGenS<A, B>) -> usize
where
    A: Clone + Hash + Eq,
{
    fn aux<A, B>(s: &mut HashSet<Triple<A>>, wz: &WGenS<A, B>) -> usize
    where
        A: Clone + Hash + Eq,
    {
        (|| {
            let (wgx, wgy) = wz.history?;
            let triple = Triple {
                x: wgx.history?.genotype.clone(),
                y: wgy.history?.genotype.clone(),
                z: wz.genotype.clone(),
            };
            if s.contains(&triple) {
                Some(0)
            } else {
                s.insert(triple);
                let resx = aux(s, wgx.history?.as_ref());
                let resy = aux(s, wgy.history?.as_ref());
                Some(1 + resx + resy)
            }
        })().unwrap_or(0)
    }
    let mut s: HashSet<Triple<A>> = HashSet::new();
    aux(&mut s, &wz)
}
