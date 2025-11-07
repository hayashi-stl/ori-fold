use std::hash::Hash;

use indexmap::{IndexMap, indexmap};

/// A "permutation" resulting from various operations
/// performed on vectors, stored sparsely.
/// 
/// Unlike an actual permutation, elements can be repeated
/// and generated from nowhere.
pub struct Permutation<I>(IndexMap<I, Option<I>>);

impl<I: Eq + Clone + Hash> Permutation<I> {
    /// Creates a permutation mapping new indices to old indices
    pub fn new(map: IndexMap<I, Option<I>>) -> Self {
        Self(map)
    }

    /// Returns the identity permutation
    pub fn identity() -> Self {
        Self::new(indexmap!{})
    }

    /// Returns the permutation that results from
    /// applying `self` and then applying `then`.
    pub fn combine(&self, mut then: Self) -> Self {
        for v in then.0.values_mut() {
            *v = v.clone().and_then(|v| self.0.get(&v).cloned().unwrap_or(Some(v)))
        }
        then
    }
}