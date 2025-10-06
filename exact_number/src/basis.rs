use std::{cell::RefCell, cmp::Ordering, fmt::{Debug, Display}, ops::Deref, sync::{Arc, LazyLock, Mutex}};

use malachite::{base::num::basic::traits::{One, Zero}, Integer};
use malachite::base::num::arithmetic::traits::{Sign, Abs};

use crate::{interval::{Fixed, Interval}, pslq::{checked_integer_relation, checked_integer_relation_product, IntegerRelationError}, rat::Rat};

pub type SqrtExprSum = Vec<(Integer, SqrtExpr)>;
#[derive(Clone, PartialEq, Eq)]
pub enum SqrtExpr {
    Int(Integer),
    Sum(SqrtExprSum)
}

impl SqrtExpr {
    pub const ONE: SqrtExpr = SqrtExpr::Int(Integer::ONE);
}

impl Display for SqrtExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Int(a) => write!(f, "√{}", a)?,
            Self::Sum(terms) => {
                write!(f, "√(")?;
                write!(f, "{}{}", terms[0].0, terms[0].1)?;
                for (a, sqrt) in terms.iter().skip(1) {
                    write!(f, " {} {}{}", if a.sign() == Ordering::Less {"-"} else {"+"}, a.abs(), sqrt)?
                }
                write!(f, ")")?
            }
        };
        Ok(())
    }
}

impl Debug for SqrtExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

pub(crate) static BASES: LazyLock<Mutex<RefCell<Vec<Arc<Basis>>>>> = LazyLock::new(|| Mutex::new(RefCell::new(vec![])));

/// A basis of with-rational-coefficients linearly independent numbers. The first element must be 1.
#[derive(Debug)]
pub struct Basis {
    pub exprs: Vec<SqrtExpr>,
    /// Matrices stored sparsely, used for multiplication
    /// Specifically, for a component i, when multiplying a by b,
    /// we have (a*b)[i] = a * matrices[i] * b.
    /// The matrices are stored sparsely as (value, component to multiply in a, component to multiply in b)
    pub matrices: Vec<Vec<(Rat, usize, usize)>>,
    /// Approximations for the basis elements.
    pub float_approx: Vec<Interval<f64>>,
    pub fixed_approx: Mutex<RefCell<Vec<Vec<Interval<Fixed>>>>>,
}

/// An error from trying to check a basis
#[derive(Clone, Debug)]
pub enum BasisError {
    /// Linear combination that proves linear dependence
    LinearlyDependent(Vec<Integer>),
    /// (a, b, error), where a * b is linearly independent with the "basis"
    NotClosed(SqrtExpr, SqrtExpr, IntegerRelationError),
    /// The first element is not 1.
    FirstElementNotOne
}

impl Basis {
    fn calc_matrices(basis: &[SqrtExpr]) -> Result<Vec<Vec<(Rat, usize, usize)>>, BasisError> {
        let mut result = vec![vec![]; basis.len()];
        for i in 0..basis.len() {
            for j in 0..basis.len() {
                if i == 0 {
                    // First basis term is 1.
                    result[j].push((Rat::ONE, i, j));
                    continue;
                } else if j == 0 {
                    // Second basis term is 1.
                    result[i].push((Rat::ONE, i, j));
                    continue;
                } else if i == j {
                    // Basis terms equal; easy to square a square root
                    match &basis[i] {
                        SqrtExpr::Int(int) => {
                            // Square sqrt(integer)? Easy!
                            result[0].push((int.clone().into(), i, j));
                            continue;
                        },
                        SqrtExpr::Sum(terms) => {
                            let indexes = terms.iter()
                                .map(|(coeff, sqrt)| {
                                    if let Some(i) = basis.iter().position(|b| b == sqrt) {Some((i, coeff))} else {None}
                                })
                                .collect::<Option<Vec<_>>>();
                            if let Some(indexes) = indexes {
                                // Easily found in basis; assign and dip
                                for (index, coeff) in indexes {
                                    result[index].push((coeff.clone().into(), i, j));
                                }
                                continue;
                            }
                        }
                    }
                }
                // Welp; couldn't easily find the expression in terms of the basis.
                // Time to bust out PSLQ
                let coeffs = checked_integer_relation_product(basis, &basis[i], &basis[j])
                    .map_err(|err| BasisError::NotClosed(basis[i].clone(), basis[j].clone(), err))?;
                for (index, coeff) in coeffs.into_iter().enumerate() {
                    if coeff != Rat::ZERO {
                        result[index].push((coeff.into(), i, j));
                    }
                }
            }
        }
        //for (i, row) in result.iter().enumerate() {
        //    println!("{}: [", i);
        //    for entry in row.iter() {
        //        println!("    ({}, {}, {})", entry.0, entry.1, entry.2);
        //    }
        //    println!("]");
        //}
        Ok(result)
    }

    /// Checks that this is actually a closed basis.
    /// Guarantees that the basis is closed, but doesn't guarantee
    /// that it's linearly independent, even though it tries its best.
    fn new_checked(basis: Vec<SqrtExpr>) -> Result<Self, BasisError> {
        if basis.len() == 0 || basis[0] != SqrtExpr::ONE {
            Err(BasisError::FirstElementNotOne)?;
        }
        if let Ok(coeffs) = checked_integer_relation(&basis) {
            Err(BasisError::LinearlyDependent(coeffs))?;
        }

        let matrices = Self::calc_matrices(&basis)?;

        // Use the exact interval for the first one instead of the silly fat interval that sqrt() will give
        let approximations = std::iter::once(Interval::<f64>::from_integer(Integer::ONE, &0.0))
            .chain(basis.iter().skip(1).cloned().map(|b| Interval::<f64>::from_sqrt_expr(b, &0.0)))
            .collect::<Vec<_>>();
        
        Ok(Self { exprs: basis, matrices, float_approx: approximations, fixed_approx: Mutex::new(RefCell::new(vec![])) })
    }

    pub(crate) fn precision(level: usize) -> u64 {
        64u64 << level
    }

    pub(crate) fn fixed_approx<'a>(exprs: &[SqrtExpr], cell: &'a mut Vec<Vec<Interval<Fixed>>>, level: usize) -> &'a [Interval<Fixed>] {
        let mut curr = cell.len();
        if level >= curr {
            cell.resize_with(level + 1, || {
                let approx = exprs.iter().cloned()
                    .map(|b| Interval::<Fixed>::from_sqrt_expr(b, &Fixed(Integer::ZERO, Self::precision(curr))))
                    .collect::<Vec<_>>();
                curr += 1;
                approx
            });
        }
        &cell[level]
    }

    /// Panics if proven to not be a basis.
    pub(crate) fn new(basis: Vec<SqrtExpr>) -> Self {
        Self::new_checked(basis).unwrap()
    }

    pub(crate) fn new_arc_checked(basis: Vec<SqrtExpr>) -> Result<ArcBasis, BasisError> {
        //println!("Basis: {:?}", basis);
        let basis_lock = BASES.lock().expect("BASES mutex is poisoned");
        let mut basis_vec = basis_lock.borrow_mut();
        let existing = basis_vec.iter().find(|b| &b.as_ref().exprs == &basis);
        let arc = if let Some(existing) = existing.cloned() {
            existing
        } else {
            let arc = Arc::new(Basis::new_checked(basis)?);
            basis_vec.push(arc.clone());
            arc
        };
        Ok(ArcBasis(arc))
    }

    fn new_arc(basis: Vec<SqrtExpr>) -> ArcBasis {
        Self::new_arc_checked(basis).unwrap()
    }
}

impl PartialEq for Basis {
    fn eq(&self, other: &Self) -> bool {
        self.exprs == other.exprs
    }
}

impl Eq for Basis {}

#[derive(Clone, Debug)]
pub struct ArcBasis(Arc<Basis>);

impl ArcBasis {
    pub(crate) fn into_unified(a: Option<Self>, b: Option<Self>) -> Option<Self> {
        match (a, b) {
            (None, None) => None,
            (None, Some(basis)) => Some(basis),
            (Some(basis), None) => Some(basis),
            (Some(a), Some(b)) => 
                if a == b {
                    Some(a)
                } else {
                    panic!("Basis {:?} doesn't match basis {:?}", a, b)
                }
        }
    }

    pub(crate) fn into_unified_vr(a: Option<Self>, b: Option<&Self>) -> Option<Self> {
        match (a, b) {
            (None, None) => None,
            (None, Some(basis)) => Some(basis.clone()),
            (Some(basis), None) => Some(basis),
            (Some(a), Some(b)) => 
                if &a == b {
                    Some(a)
                } else {
                    panic!("Basis {:?} doesn't match basis {:?}", a, b)
                }
        }
    }

    pub(crate) fn into_unified_rv(a: Option<&Self>, b: Option<Self>) -> Option<Self> {
        match (a, b) {
            (None, None) => None,
            (None, Some(basis)) => Some(basis),
            (Some(basis), None) => Some(basis.clone()),
            (Some(a), Some(b)) => 
                if a == &b {
                    Some(b)
                } else {
                    panic!("Basis {:?} doesn't match basis {:?}", a, b)
                }
        }
    }

    pub(crate) fn to_unified<'a>(a: Option<&'a Self>, b: Option<&'a Self>) -> Option<&'a Self> {
        match (a, b) {
            (None, None) => None,
            (None, Some(basis)) => Some(basis),
            (Some(basis), None) => Some(basis),
            (Some(a), Some(b)) => 
                if a == b {
                    Some(a)
                } else {
                    panic!("Basis {:?} doesn't match basis {:?}", a, b)
                }
        }
    }

    pub(crate) fn assert_compatible(&self, other: &ArcBasis) -> () {
        assert_eq!(self, other)
    }
}

impl PartialEq for ArcBasis {
    fn eq(&self, other: &Self) -> bool {
        Arc::as_ptr(self) == Arc::as_ptr(other)
    }
}

impl Eq for ArcBasis {}

impl Deref for ArcBasis {
    type Target = Arc<Basis>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

//impl Clone for ArcBasis {
//    fn clone(&self) -> Self {
//        let _basis_lock = BASES.lock().expect("BASES mutex is poisoned");
//        Self(self.0.clone())
//    }
//}

impl Drop for ArcBasis {
    fn drop(&mut self) {
        let basis_lock = BASES.lock().expect("BASES mutex is poisoned");
        // If `strong_count == 2`:
        //     Only references left are this one and the one in BASES.
        //     This means that no other threads have a clone of this `ArcBasis` and thus
        //     can't change the reference count by calling `clone()`.
        //     They can try to change the reference count by calling `Basis::new_arc`, but
        //     that's protected by the same lock.
        //
        // If `strong_count > 2`:
        //     The only way another thread can change the reference count while this is locked
        //     is by calling `clone()`, which can only increase the reference count.
        if Arc::strong_count(&self.0) == 2 {
            let mut basis_vec = basis_lock.borrow_mut();
            if let Some((i, _)) = basis_vec.iter().enumerate().find(|(_, b)| &b.as_ref().exprs == &self.exprs) {
                basis_vec.swap_remove(i);
            }
        }
    }
}
