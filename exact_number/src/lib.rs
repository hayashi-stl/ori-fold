use std::cell::{RefCell, RefMut};
use std::cmp::Ordering;
use std::fmt::{Debug, Display};
use std::ops::{Add, AddAssign, Deref, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::{Arc, LazyLock, Mutex, MutexGuard};

use malachite::base::num::basic::traits::{One, Zero as _};
use malachite::Integer;
use malachite::base::num::arithmetic::traits::{Abs, NegAssign, Sign, Square};
use nalgebra::{DMatrix, DVector, DVectorView, RealField};
use num::{Signed, Zero};

use crate::interval::{Fixed, Interval, IntervalContent};
use crate::pslq::{checked_integer_relation, checked_integer_relation_product, IntegerRelationError};
use crate::rat::Rat;

mod pslq;
mod algebraic;
mod rat;
mod interval;
//mod matrix;

#[macro_export]
macro_rules! sqrt_expr_helper {
    ($( { $as:literal $an:literal sqrt $($sqrt:tt)* } )*) => {
        crate::SqrtExpr::Sum(vec![$(
            (
                malachite::Integer::from_sign_and_abs($as, ($an as u64).into()),
                crate::sqrt_expr_parse!($($sqrt)*)
            ),
        )*])
    };
}

#[macro_export]
macro_rules! sqrt_expr_parse {
    ($int:literal) => { crate::SqrtExpr::Int(malachite::Integer::from($int as i64)) };
    ($( { $($init:tt)* } )* + sqrt $($tt:tt)*) => { crate::sqrt_expr_parse!($( { $($init)* } )* [+] 1 sqrt $($tt)*) };
    ($( { $($init:tt)* } )* - sqrt $($tt:tt)*) => { crate::sqrt_expr_parse!($( { $($init)* } )* [-] 1 sqrt $($tt)*) };
    ($( { $($init:tt)* } )* + $($tt:tt)*) => { crate::sqrt_expr_parse!($( { $($init)* } )* [+] $($tt)*) };
    ($( { $($init:tt)* } )* - $($tt:tt)*) => { crate::sqrt_expr_parse!($( { $($init)* } )* [-] $($tt)*) };
    ($( { $($init:tt)* } )* sqrt $($tt:tt)*) => { crate::sqrt_expr_parse!($( { $($init)* } )* [+] 1 sqrt $($tt)*) };
    ($an:literal $($tt:tt)*) => { crate::sqrt_expr_parse!(+ $an $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal) => { crate::sqrt_expr_parse!($( { $($init)* } )* [+] $an sqrt(1)) };
    ($( { $($init:tt)* } )* [-] $an:literal) => { crate::sqrt_expr_parse!($( { $($init)* } )* [-] $an sqrt(1)) };
    ($( { $($init:tt)* } )* [+] $an:literal + $($tt:tt)*) => { crate::sqrt_expr_parse!($( { $($init)* } )* [+] $an sqrt(1) + $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal + $($tt:tt)*) => { crate::sqrt_expr_parse!($( { $($init)* } )* [-] $an sqrt(1) + $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal - $($tt:tt)*) => { crate::sqrt_expr_parse!($( { $($init)* } )* [+] $an sqrt(1) - $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal - $($tt:tt)*) => { crate::sqrt_expr_parse!($( { $($init)* } )* [-] $an sqrt(1) - $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal sqrt $sqrt:literal $($tt:tt)*) => { crate::sqrt_expr_parse!($( { $($init)* } )* [+] $an sqrt ($sqrt) $($tt)* ) };
    ($( { $($init:tt)* } )* [-] $an:literal sqrt $sqrt:literal $($tt:tt)*) => { crate::sqrt_expr_parse!($( { $($init)* } )* [-] $an sqrt ($sqrt) $($tt)* ) };
    ($( { $($init:tt)* } )* [+] $an:literal sqrt ($($sqrt:tt)*) $($tt:tt)*) => { crate::sqrt_expr_parse!($( { $($init)* } )* { true $an sqrt $($sqrt)* } $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal sqrt ($($sqrt:tt)*) $($tt:tt)*) => { crate::sqrt_expr_parse!($( { $($init)* } )* { false $an sqrt $($sqrt)* } $($tt)*) };
    ($( { $($terms:tt)* } )*) => { crate::sqrt_expr_helper!($( { $($terms)* } )*)}
}

#[macro_export]
macro_rules! sqrt_expr {
    ($($tt:tt)*) => {
        crate::sqrt_expr_parse!($($tt)*)
    };
}

#[macro_export]
macro_rules! expr_helper {
    ($( { $as:literal $an:literal / $ad:literal sqrt $($sqrt:tt)* })*) => {
        crate::BasedExpr::from_terms(vec![$(
            (
                crate::rat::Rat::from_sign_and_naturals($as, malachite::Natural::from($an as u64), malachite::Natural::from($ad as u64)),
                crate::sqrt_expr_parse!($($sqrt)*)
            ),
        )*])
    };
}

#[macro_export]
macro_rules! expr_parse {
    ($( { $($init:tt)* } )* + sqrt $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [+] 1 / 1 sqrt $($tt)*) };
    ($( { $($init:tt)* } )* - sqrt $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [-] 1 / 1 sqrt $($tt)*) };
    ($( { $($init:tt)* } )* + $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [+] $($tt)*) };
    ($( { $($init:tt)* } )* - $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [-] $($tt)*) };
    ($( { $($init:tt)* } )* sqrt $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [+] 1 / 1 sqrt $($tt)*) };
    ($an:literal $($tt:tt)*) => { crate::expr_parse!(+ $an $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal) => { crate::expr_parse!($( { $($init)* } )* [+] $an / 1) };
    ($( { $($init:tt)* } )* [-] $an:literal) => { crate::expr_parse!($( { $($init)* } )* [-] $an / 1) };
    ($( { $($init:tt)* } )* [+] $an:literal + $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [+] $an / 1 + $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal + $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [-] $an / 1 + $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal - $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [+] $an / 1 - $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal - $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [-] $an / 1 - $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal sqrt $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [+] $an / 1 sqrt $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal sqrt $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [-] $an / 1 sqrt $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal / $ad:literal) => { crate::expr_parse!($( { $($init)* } )* [+] $an / $ad sqrt(1)) };
    ($( { $($init:tt)* } )* [-] $an:literal / $ad:literal) => { crate::expr_parse!($( { $($init)* } )* [-] $an / $ad sqrt(1)) };
    ($( { $($init:tt)* } )* [+] $an:literal / $ad:literal + $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [+] $an / $ad sqrt(1) + $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal / $ad:literal + $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [-] $an / $ad sqrt(1) + $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal / $ad:literal - $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [+] $an / $ad sqrt(1) - $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal / $ad:literal - $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [-] $an / $ad sqrt(1) - $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal / $ad:literal sqrt $sqrt:literal $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [+] $an / $ad sqrt ($sqrt) $($tt)* ) };
    ($( { $($init:tt)* } )* [-] $an:literal / $ad:literal sqrt $sqrt:literal $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* [-] $an / $ad sqrt ($sqrt) $($tt)* ) };
    ($( { $($init:tt)* } )* [+] $an:literal / $ad:literal sqrt ($($sqrt:tt)*) $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* { true $an / $ad sqrt $($sqrt)* } $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal / $ad:literal sqrt ($($sqrt:tt)*) $($tt:tt)*) => { crate::expr_parse!($( { $($init)* } )* { false $an / $ad sqrt $($sqrt)* } $($tt)*) };
    ($( { $($terms:tt)* } )*) => { crate::expr_helper!($( { $($terms)* } )*)}
}

#[macro_export]
macro_rules! based_expr {
    ($($tt:tt)*) => {
        crate::expr_parse!($($tt)*)
    };
}

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

static BASES: LazyLock<Mutex<RefCell<Vec<Arc<Basis>>>>> = LazyLock::new(|| Mutex::new(RefCell::new(vec![])));

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

    fn precision(level: usize) -> u64 {
        64u64 << level
    }

    fn fixed_approx<'a>(exprs: &[SqrtExpr], cell: &'a mut Vec<Vec<Interval<Fixed>>>, level: usize) -> &'a [Interval<Fixed>] {
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
    fn new(basis: Vec<SqrtExpr>) -> Self {
        Self::new_checked(basis).unwrap()
    }

    fn new_arc_checked(basis: Vec<SqrtExpr>) -> Result<ArcBasis, BasisError> {
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
    fn into_unified(a: Option<Self>, b: Option<Self>) -> Option<Self> {
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

    fn into_unified_vr(a: Option<Self>, b: Option<&Self>) -> Option<Self> {
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

    fn into_unified_rv(a: Option<&Self>, b: Option<Self>) -> Option<Self> {
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

    fn to_unified<'a>(a: Option<&'a Self>, b: Option<&'a Self>) -> Option<&'a Self> {
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

    fn assert_compatible(&self, other: &ArcBasis) -> () {
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

/// An expression of a number, with a basis.
/// Two values with different bases are incompatible and mixing them results in a panic.
#[derive(Clone)]
pub enum BasedExpr {
    /// The basis is unknown. Used when a value is needed and the basis can't be passed in
    Undefined(Rat),
    /// The basis is the list of `SqrtExpr`. The first element of the basis must be 1.
    Based(DVector<Rat>, ArcBasis)
}

impl Default for BasedExpr {
    fn default() -> Self {
        Self::Undefined(Rat::ZERO)
    }
}

fn longer_first<T>(lhs: Vec<T>, rhs: Vec<T>) -> (Vec<T>, Vec<T>) {
    if lhs.len() < rhs.len() { (rhs, lhs) } else { (lhs, rhs) }
}

fn longer_first_ref<'a, T>(lhs: &'a [T], rhs: &'a [T]) -> (&'a [T], &'a [T]) {
    if lhs.len() < rhs.len() { (rhs, lhs) } else { (lhs, rhs) }
}

impl BasedExpr {
    /// Constructs a BasedExpr out of a rational without defining the basis yet
    pub fn new_undefined(q: Rat) -> Self {
        Self::Undefined(q)
    }

    /// Constructs a BasedExpr from its terms. The basis is defined by the second values of each entry.
    pub fn from_terms_checked(terms: impl IntoIterator<Item = (Rat, SqrtExpr)>) -> Result<Self, BasisError> {
        let (mut coeffs, mut basis): (Vec<_>, Vec<_>) = terms.into_iter().unzip();
        if basis[0] != SqrtExpr::ONE {
            coeffs.insert(0, 0u64.into());
            basis.insert(0, SqrtExpr::ONE);
        }
        Ok(BasedExpr::Based(coeffs.into(), Basis::new_arc_checked(basis)?))
    }

    /// Constructs a BasedExpr from its terms. The basis is defined by the second values of each entry.
    /// Panics if the basis is proven to not be a basis.
    pub fn from_terms(terms: impl IntoIterator<Item = (Rat, SqrtExpr)>) -> Self {
        Self::from_terms_checked(terms).unwrap()
    }

    /// Gets the sign of this expression.
    /// Returns `Ordering::Greater` if greater than 0, `Ordering::Equal` if equal to 0,
    /// and `Ordering::Less` if less than 0.
    pub fn sign(&self) -> Ordering {
        match self {
            Self::Undefined(q) => q.sign(),
            Self::Based(coeffs, basis) => {
                if coeffs.iter().all(|c| c.is_zero()) {
                    return Ordering::Equal;
                }
                unimplemented!()
            }
        }
    }

    fn interval_f64(coeffs: &DVector<Rat>, basis: &Basis) -> Interval<f64> {
        coeffs.iter().zip(basis.float_approx.iter())
            .map(|(coeff, elem)| Interval::from_rational(coeff.clone(), &0.0) * elem)
            .sum()
    }

    /// The precision is 64 << `precision_level`
    fn interval_fixed(coeffs: &DVector<Rat>, basis: &Basis, precision_level: usize) -> Interval<Fixed> {
        let precision = Basis::precision(precision_level);
        let lock = basis.fixed_approx.lock().unwrap();
        let mut borrow = lock.borrow_mut();
        coeffs.iter().zip(Basis::fixed_approx(&basis.exprs, &mut borrow, precision_level).iter())
            .map(|(coeff, elem)| Interval::from_rational(coeff.clone(), &Fixed(Integer::ZERO, precision)) * elem)
            .sum()
    }

    fn based_rational(q: Rat, basis: ArcBasis) -> Self {
        let mut coeffs = vec![Rat::ZERO; basis.exprs.len()];
        coeffs[0] = q;
        Self::Based(coeffs.into(), basis)
    }

    fn has_basis(&self) -> bool {
        match self {
            Self::Undefined(_) => false,
            Self::Based(_, basis) => true
        }
    }

    fn basis(&self) -> Option<&ArcBasis> {
        match self {
            Self::Undefined(_) => None,
            Self::Based(_, basis) => Some(basis)
        }
    }

    fn assert_compatible(&self, other: &BasedExpr) -> () {
        if let (Self::Based(_, basis_a), Self::Based(_, basis_b)) = (self, other) {
            basis_a.assert_compatible(basis_b);
        }
    }

    fn into_coeffs_basis(self) -> (DVector<Rat>, Option<ArcBasis>) {
        match self {
            Self::Undefined(q) => (vec![q].into(), None),
            Self::Based(coeffs, basis) => (coeffs, Some(basis))
        }
    }

    fn to_coeffs_basis(&'_ self) -> (DVectorView<'_, Rat>, Option<&'_ ArcBasis>) {
        match self {
            Self::Undefined(q) => (DVectorView::from_slice(std::slice::from_ref(q), 1), None),
            Self::Based(coeffs, basis) => (coeffs.as_view(), Some(basis))
        }
    }

    fn into_coeffs_basis_2(self, other: BasedExpr) -> (DVector<Rat>, DVector<Rat>, Option<ArcBasis>) {
        let (coeffs_a, basis_a) = self.into_coeffs_basis();
        let (coeffs_b, basis_b) = other.into_coeffs_basis();
        (coeffs_a, coeffs_b, ArcBasis::into_unified(basis_a, basis_b))
    }

    fn into_coeffs_basis_2_vr(self, other: &'_ BasedExpr) -> (DVector<Rat>, DVectorView<'_, Rat>, Option<ArcBasis>) {
        let (coeffs_a, basis_a) = self.into_coeffs_basis();
        let (coeffs_b, basis_b) = other.to_coeffs_basis();
        (coeffs_a, coeffs_b, ArcBasis::into_unified_vr(basis_a, basis_b))
    }

    fn into_coeffs_basis_2_rv(&'_ self, other: BasedExpr) -> (DVectorView<'_, Rat>, DVector<Rat>, Option<ArcBasis>) {
        let (coeffs_a, basis_a) = self.to_coeffs_basis();
        let (coeffs_b, basis_b) = other.into_coeffs_basis();
        (coeffs_a, coeffs_b, ArcBasis::into_unified_rv(basis_a, basis_b))
    }

    fn into_coeffs_basis_2_rr<'a>(&'a self, other: &'a BasedExpr) -> (DVectorView<'a, Rat>, DVectorView<'a, Rat>, Option<ArcBasis>) {
        let (coeffs_a, basis_a) = self.to_coeffs_basis();
        let (coeffs_b, basis_b) = other.to_coeffs_basis();
        (coeffs_a, coeffs_b, ArcBasis::to_unified(basis_a, basis_b).cloned())
    }

    fn to_coeffs_basis_2<'a>(&'a self, other: &'a BasedExpr) -> (DVectorView<'a, Rat>, DVectorView<'a, Rat>, Option<&'a ArcBasis>) {
        let (coeffs_a, basis_a) = self.to_coeffs_basis();
        let (coeffs_b, basis_b) = other.to_coeffs_basis();
        (coeffs_a, coeffs_b, ArcBasis::to_unified(basis_a, basis_b))
    }
}

impl Zero for BasedExpr {
    fn zero() -> Self {
        Self::new_undefined(Rat::ZERO)
    }

    fn is_zero(&self) -> bool {
        match self {
            Self::Undefined(q) => q == &Rat::ZERO,
            Self::Based(coeffs, _) => coeffs.iter().all(|c| c == &Rat::ZERO)
        }
    }
}

impl Display for BasedExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Undefined(q) => write!(f, "{}", q),
            Self::Based(coeffs, basis) => {
                let basis = basis.as_ref();
                write!(f, "{}{}", coeffs[0], basis.exprs[0])?;
                for (a, sqrt) in coeffs.iter().zip(basis.exprs.iter()).skip(1) {
                    write!(f, " {} {}{}", if a.sign() == Ordering::Less {"-"} else {"+"}, Abs::abs(a), sqrt)?
                }
                Ok(())
            }
        }?;
        Ok(())
    }
}

impl Debug for BasedExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

impl PartialEq for BasedExpr {
    fn eq(&self, other: &Self) -> bool {
        if !self.has_basis() && other.has_basis() { return other == self }
        match (self, other) {
            (BasedExpr::Undefined(a), BasedExpr::Undefined(b)) => a == b,
            (BasedExpr::Based(a, _), BasedExpr::Undefined(b)) =>
                &a[0] == b && a.iter().skip(1).all(|q| q == &Rat::ZERO),
            (BasedExpr::Undefined(_), BasedExpr::Based(_, _)) => unreachable!(),
            (BasedExpr::Based(a, basis_a), BasedExpr::Based(b, basis_b)) => {
                basis_a.assert_compatible(&basis_b);
                a.iter().zip(b.iter()).all(|(a, b)| a == b)
            }
        }
    }
}

impl Eq for BasedExpr {}

impl AddAssign<BasedExpr> for BasedExpr {
    fn add_assign(&mut self, mut rhs: BasedExpr) {
        if !self.has_basis() && rhs.has_basis() {
            std::mem::swap(self, &mut rhs);
        }
        match (self, rhs) {
            (BasedExpr::Undefined(a), BasedExpr::Undefined(b)) => *a += b,
            (BasedExpr::Based(a, _), BasedExpr::Undefined(b)) => a[0] += b,
            (BasedExpr::Undefined(_), BasedExpr::Based(_, _)) => unreachable!(),
            (BasedExpr::Based(a, basis_a), BasedExpr::Based(b, basis_b)) => {
                basis_a.assert_compatible(&basis_b);
                *a += b;
            }
        }
    }
}

impl AddAssign<&BasedExpr> for BasedExpr {
    fn add_assign(&mut self, rhs: &BasedExpr) {
        if !self.has_basis() && rhs.has_basis() {
            let new_rhs = std::mem::replace(self, rhs.clone());
            *self += new_rhs;
            return;
        }
        match (self, rhs) {
            (BasedExpr::Undefined(a), BasedExpr::Undefined(b)) => *a += b,
            (BasedExpr::Based(a, _), BasedExpr::Undefined(b)) => a[0] += b,
            (BasedExpr::Undefined(_), BasedExpr::Based(_, _)) => unreachable!(),
            (BasedExpr::Based(a, basis_a), BasedExpr::Based(b, basis_b)) => {
                basis_a.assert_compatible(&basis_b);
                *a += b;
            }
        }
    }
}

impl Add<BasedExpr> for BasedExpr {
    type Output = BasedExpr;

    fn add(mut self, rhs: BasedExpr) -> Self::Output {
        self += rhs;
        self
    }
}

impl Add<&BasedExpr> for BasedExpr {
    type Output = BasedExpr;

    fn add(mut self, rhs: &BasedExpr) -> Self::Output {
        self += rhs;
        self
    }
}

impl<'a> Add<BasedExpr> for &'a BasedExpr {
    type Output = BasedExpr;

    fn add(self, rhs: BasedExpr) -> Self::Output {
        rhs + self
    }
}

impl<'a> Add<&BasedExpr> for &'a BasedExpr {
    type Output = BasedExpr;

    fn add(self, rhs: &BasedExpr) -> Self::Output {
        if !self.has_basis() && rhs.has_basis() { return rhs + self }
        self.clone() + rhs
    }
}

impl NegAssign for BasedExpr {
    fn neg_assign(&mut self) {
        match self {
            BasedExpr::Undefined(q) => q.neg_assign(),
            BasedExpr::Based(coeffs, _) => coeffs.iter_mut().for_each(|q| q.neg_assign()),
        }
    }
}

impl Neg for BasedExpr {
    type Output = BasedExpr;

    fn neg(mut self) -> Self::Output {
        self.neg_assign();
        self
    }
}

impl<'a> Neg for &'a BasedExpr {
    type Output = BasedExpr;

    fn neg(self) -> Self::Output {
        -self.clone()
    }
}

impl SubAssign<BasedExpr> for BasedExpr {
    fn sub_assign(&mut self, mut rhs: BasedExpr) {
        if !self.has_basis() && rhs.has_basis() {
            std::mem::swap(self, &mut rhs);
            self.neg_assign();
            *self += rhs;
            return;
        }
        match (self, rhs) {
            (BasedExpr::Undefined(a), BasedExpr::Undefined(b)) => *a -= b,
            (BasedExpr::Based(a, _), BasedExpr::Undefined(b)) => a[0] -= b,
            (BasedExpr::Undefined(_), BasedExpr::Based(_, _)) => unreachable!(),
            (BasedExpr::Based(a, basis_a), BasedExpr::Based(b, basis_b)) => {
                basis_a.assert_compatible(&basis_b);
                *a -= b;
            }
        }
    }
}

impl SubAssign<&BasedExpr> for BasedExpr {
    fn sub_assign(&mut self, rhs: &BasedExpr) {
        if !self.has_basis() && rhs.has_basis() {
            let new_rhs = std::mem::replace(self, -rhs);
            *self += new_rhs;
            return;
        }
        match (self, rhs) {
            (BasedExpr::Undefined(a), BasedExpr::Undefined(b)) => *a -= b,
            (BasedExpr::Based(a, _), BasedExpr::Undefined(b)) => a[0] -= b,
            (BasedExpr::Undefined(_), BasedExpr::Based(_, _)) => unreachable!(),
            (BasedExpr::Based(a, basis_a), BasedExpr::Based(b, basis_b)) => {
                basis_a.assert_compatible(&basis_b);
                *a -= b;
            }
        }
    }
}

impl Sub<BasedExpr> for BasedExpr {
    type Output = BasedExpr;

    fn sub(mut self, rhs: BasedExpr) -> Self::Output {
        self -= rhs;
        self
    }
}

impl Sub<&BasedExpr> for BasedExpr {
    type Output = BasedExpr;

    fn sub(mut self, rhs: &BasedExpr) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<'a> Sub<BasedExpr> for &'a BasedExpr {
    type Output = BasedExpr;

    fn sub(self, rhs: BasedExpr) -> Self::Output {
        -rhs + self
    }
}

impl<'a> Sub<&BasedExpr> for &'a BasedExpr {
    type Output = BasedExpr;

    fn sub(self, rhs: &BasedExpr) -> Self::Output {
        if !self.has_basis() && rhs.has_basis() { return -rhs + self }
        self.clone() - rhs
    }
}

impl MulAssign<BasedExpr> for BasedExpr {
    fn mul_assign(&mut self, mut rhs: BasedExpr) {
        if !self.has_basis() && rhs.has_basis() {
            std::mem::swap(self, &mut rhs);
        }
        match (self, rhs) {
            (BasedExpr::Undefined(a), BasedExpr::Undefined(b)) => *a *= b,
            (BasedExpr::Based(a, _), BasedExpr::Undefined(b)) => a.iter_mut().for_each(|a| *a *= &b),
            (BasedExpr::Undefined(_), BasedExpr::Based(_, _)) => unreachable!(),
            (BasedExpr::Based(a, basis_a), BasedExpr::Based(b, basis_b)) => {
                basis_a.assert_compatible(&basis_b);
                *a = DVector::from_iterator(a.len(), basis_a.matrices.iter().map(|mtx| {
                    mtx.iter().map(|(coeff, ai, bi)| {
                        coeff * &a[*ai] * &b[*bi]
                    }).sum()
                }));
            }
        }
    }
}

impl MulAssign<&BasedExpr> for BasedExpr {
    fn mul_assign(&mut self, rhs: &BasedExpr) {
        if !self.has_basis() && rhs.has_basis() {
            let new_rhs = std::mem::replace(self, rhs.clone());
            *self *= new_rhs;
            return;
        }
        match (self, rhs) {
            (BasedExpr::Undefined(a), BasedExpr::Undefined(b)) => *a *= b,
            (BasedExpr::Based(a, _), BasedExpr::Undefined(b)) => a.iter_mut().for_each(|a| *a *= b),
            (BasedExpr::Undefined(_), BasedExpr::Based(_, _)) => unreachable!(),
            (BasedExpr::Based(a, basis_a), BasedExpr::Based(b, basis_b)) => {
                basis_a.assert_compatible(&basis_b);
                *a = DVector::from_iterator(a.len(), basis_a.matrices.iter().map(|mtx| {
                    mtx.iter().map(|(coeff, ai, bi)| {
                        coeff * &a[*ai] * &b[*bi]
                    }).sum()
                }));
            }
        }
    }
}

impl Mul<BasedExpr> for BasedExpr {
    type Output = BasedExpr;

    fn mul(mut self, rhs: BasedExpr) -> Self::Output {
        self *= rhs;
        self
    }
}

impl Mul<&BasedExpr> for BasedExpr {
    type Output = BasedExpr;

    fn mul(mut self, rhs: &BasedExpr) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<'a> Mul<BasedExpr> for &'a BasedExpr {
    type Output = BasedExpr;

    fn mul(self, rhs: BasedExpr) -> Self::Output {
        rhs * self
    }
}

impl<'a> Mul<&BasedExpr> for &'a BasedExpr {
    type Output = BasedExpr;

    fn mul(self, rhs: &BasedExpr) -> Self::Output {
        if !self.has_basis() && rhs.has_basis() { return rhs * self }
        self.clone() * rhs
    }
}

impl DivAssign<BasedExpr> for BasedExpr {
    fn div_assign(&mut self, rhs: BasedExpr) {
        if let (None, Some(basis)) = (self.basis(), rhs.basis()) {
            *self = BasedExpr::based_rational(if let BasedExpr::Undefined(q) = std::mem::take(self) {q} else {unreachable!()}, basis.clone());
        }
        match (self, rhs) {
            (BasedExpr::Undefined(a), BasedExpr::Undefined(b)) => *a /= b,
            (BasedExpr::Based(a, _), BasedExpr::Undefined(b)) => a.iter_mut().for_each(|a| *a /= &b),
            (BasedExpr::Undefined(_), BasedExpr::Based(_, _)) => unreachable!(),
            (BasedExpr::Based(a, basis_a), BasedExpr::Based(b, basis_b)) => {
                // The dreaded division
                // TODO: Make this faster if necessary.
                basis_a.assert_compatible(&basis_b);
                let mut mtx = DMatrix::repeat(a.len(), a.len(), Rat::ZERO);
                for (row_i, row) in basis_a.matrices.iter().enumerate() {
                    for (coeff, a_i, b_i) in row {
                        mtx[(row_i, *b_i)] += coeff * &b[*a_i];
                    }
                }
                let result = mtx.try_inverse_mut();
                if !result { panic!("Tried to divide by 0"); }
                *a = mtx * std::mem::take(a);
            }
        }
    }
}

impl DivAssign<&BasedExpr> for BasedExpr {
    fn div_assign(&mut self, rhs: &BasedExpr) {
        if let (None, Some(basis)) = (self.basis(), rhs.basis()) {
            *self = BasedExpr::based_rational(if let BasedExpr::Undefined(q) = std::mem::take(self) {q} else {unreachable!()}, basis.clone());
        }
        match (self, rhs) {
            (BasedExpr::Undefined(a), BasedExpr::Undefined(b)) => *a /= b,
            (BasedExpr::Based(a, _), BasedExpr::Undefined(b)) => a.iter_mut().for_each(|a| *a /= b),
            (BasedExpr::Undefined(_), BasedExpr::Based(_, _)) => unreachable!(),
            (BasedExpr::Based(a, basis_a), BasedExpr::Based(b, basis_b)) => {
                // The dreaded division
                // TODO: Make this faster if necessary.
                basis_a.assert_compatible(&basis_b);
                let mut mtx = DMatrix::repeat(a.len(), a.len(), Rat::ZERO);
                for (row_i, row) in basis_a.matrices.iter().enumerate() {
                    for (coeff, a_i, b_i) in row {
                        mtx[(row_i, *b_i)] += coeff * &b[*a_i];
                    }
                }
                let result = mtx.try_inverse_mut();
                if !result { panic!("Tried to divide by 0"); }
                *a = mtx * std::mem::take(a);
            }
        }
    }
}

impl Div<BasedExpr> for BasedExpr {
    type Output = BasedExpr;

    fn div(mut self, rhs: BasedExpr) -> Self::Output {
        self /= rhs;
        self
    }
}

impl Div<&BasedExpr> for BasedExpr {
    type Output = BasedExpr;

    fn div(mut self, rhs: &BasedExpr) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<'a> Div<BasedExpr> for &'a BasedExpr {
    type Output = BasedExpr;

    fn div(self, mut rhs: BasedExpr) -> Self::Output {
        let denom = std::mem::replace(&mut rhs, self.clone());
        rhs / denom
    }
}

impl<'a> Div<&BasedExpr> for &'a BasedExpr {
    type Output = BasedExpr;

    fn div(self, rhs: &BasedExpr) -> Self::Output {
        if let (None, Some(basis)) = (self.basis(), rhs.basis()) {
            let numer = BasedExpr::based_rational(if let BasedExpr::Undefined(q) = self {q.clone()} else {unreachable!()}, basis.clone());
            return numer / rhs
        }
        self.clone() / rhs
    }
}

impl PartialOrd for BasedExpr {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let diff = self - other;
        match diff {
            BasedExpr::Undefined(a) => a.partial_cmp(&Rat::ZERO),
            BasedExpr::Based(a, basis) => {
                // Easy cases
                if a.iter().all(|a| a == &Rat::ZERO) {
                    return Some(Ordering::Equal);
                } if a.len() == 1 {
                    return a[0].partial_cmp(&Rat::ZERO);
                } else if a.len() == 2 {
                    // The basis has to be [1, sqrt(c)] for some square-free integer c >= 2
                    // However, it could be represented as, e.g. [1, sqrt(sqrt(4))].
                    // Why you would do that, I don't know.
                    if let SqrtExpr::Int(c) = &basis.exprs[1] {
                        return (a[0].signum() * (&a[0]).square() + a[1].signum() * (&a[1]).square() * Rat::from(c)).partial_cmp(&Rat::ZERO);
                    }
                }

                let float_approx = Self::interval_f64(&a, basis.as_ref());
                let ordering = float_approx.cmp_zero();
                if ordering.is_some() { return ordering; }

                let mut level = 0;
                // Increase the interval; we need separation!
                loop {
                    let approx = Self::interval_fixed(&a, basis.as_ref(), level);
                    let ordering = approx.cmp_zero();
                    if ordering.is_some() { return ordering; }
                    level += 1;
                }
            }
        }
    }
}

impl Ord for BasedExpr {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

//impl Signed for BasedExpr {
//    fn abs(&self) -> Self {
//        todo!()
//    }
//
//    fn abs_sub(&self, other: &Self) -> Self {
//        todo!()
//    }
//
//    fn signum(&self) -> Self {
//        todo!()
//    }
//
//    fn is_positive(&self) -> bool {
//        todo!()
//    }
//
//    fn is_negative(&self) -> bool {
//        todo!()
//    }
//}

#[cfg(test)]
mod tests {
    use std::sync::RwLock;

    use super::*;

    static LOCK: LazyLock<RwLock<()>> = LazyLock::new(|| RwLock::new(()));

    #[test]
    fn test_macro_int() {
        let _lock = LOCK.read();
        let result = based_expr!(1);
        let expected = BasedExpr::from_terms(vec![(Rat::from(1i64), SqrtExpr::Int(Integer::from(1i64)))]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_negative_int() {
        let _lock = LOCK.read();
        let result = based_expr!(-3);
        let expected = BasedExpr::from_terms(vec![(Rat::from(-3i64), SqrtExpr::Int(Integer::from(1i64)))]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_rational() {
        let _lock = LOCK.read();
        let result = based_expr!(3/5);
        let expected = BasedExpr::from_terms(vec![(Rat::from_integers(Integer::from(3i64), Integer::from(5i64)), SqrtExpr::Int(Integer::from(1i64)))]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_simple_sqrt() {
        let _lock = LOCK.read();
        let result = based_expr!(sqrt 3);
        let expected = BasedExpr::from_terms(vec![(Rat::from(1i64), SqrtExpr::Int(Integer::from(3i64)))]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_plus_sqrt() {
        let _lock = LOCK.read();
        let result = based_expr!(-4 + 2 sqrt 2);
        let expected = BasedExpr::from_terms(vec![
            (Rat::from(-4i64), SqrtExpr::Int(Integer::from(1i64))),
            (Rat::from(2i64), SqrtExpr::Int(Integer::from(2i64))),
            ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_rational_sqrt() {
        let _lock = LOCK.read();
        let result = based_expr!(1/2 + 1/2 sqrt 5);
        let expected = BasedExpr::from_terms(vec![
            (Rat::from_integers(1i64.into(), 2i64.into()), SqrtExpr::Int(Integer::from(1i64))),
            (Rat::from_integers(1i64.into(), 2i64.into()), SqrtExpr::Int(Integer::from(5i64))),
            ]);
        assert_eq!(result, expected);
    }

    //#[test]
    //fn test_macro_frrt() {
    //    let _lock = LOCK.read();
    //    let result = based_expr!(sqrt (sqrt 5));
    //    let expected = BasedExpr::from_terms(vec![
    //        (Rat::from(1i64), SqrtExpr::Sum(vec![
    //            (Integer::from(1i64), SqrtExpr::Int(Integer::from(5i64))),
    //        ])),
    //    ]);
    //    assert_eq!(result, expected);
    //}

    #[test]
    fn test_macro_complicated_sqrt() {
        let _lock = LOCK.read();
        // Real-life example! tan 11.25°
        let result = based_expr!(-1 - sqrt 2 + 2 sqrt(2 - sqrt 2) + sqrt(4 - 2 sqrt 2));
        let expected = BasedExpr::from_terms(vec![
            (Rat::from(-1i64), SqrtExpr::Int(Integer::from(1i64))),
            (Rat::from(-1i64), SqrtExpr::Int(Integer::from(2i64))),
            (Rat::from(2i64), SqrtExpr::Sum(vec![
                (Integer::from(2i64), SqrtExpr::Int(Integer::from(1i64))),
                (Integer::from(-1i64), SqrtExpr::Int(Integer::from(2i64))),
            ])),
            (Rat::from(1i64), SqrtExpr::Sum(vec![
                (Integer::from(4i64), SqrtExpr::Int(Integer::from(1i64))),
                (Integer::from(-2i64), SqrtExpr::Int(Integer::from(2i64))),
            ])),
        ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_bases_list() {
        let _lock = LOCK.write();
        let basis = Basis::new(vec![sqrt_expr_parse!(1), sqrt_expr_parse!(2)]);
        let expr = based_expr!(1 + sqrt 2);

        {
            let basis_lock = BASES.lock().unwrap();
            let basis_vec = basis_lock.borrow();
            assert_eq!(basis_vec.len(), 1);
            assert_eq!(basis_vec[0].exprs, basis.exprs);
        }

        let expr2 = expr.clone();
        {
            let basis_lock = BASES.lock().unwrap();
            let basis_vec = basis_lock.borrow();
            assert_eq!(basis_vec.len(), 1);
        }

        let expr3 = based_expr!(3 + 2 sqrt 2); // same basis
        {
            let basis_lock = BASES.lock().unwrap();
            let basis_vec = basis_lock.borrow();
            assert_eq!(basis_vec.len(), 1);
        }

        let expr4 = based_expr!(5 + sqrt 5); // different basis
        {
            let basis_lock = BASES.lock().unwrap();
            let basis_vec = basis_lock.borrow();
            assert_eq!(basis_vec.len(), 2);
        }

        drop(expr);
        drop(expr2);
        drop(expr3);
        drop(expr4);
        {
            let basis_lock = BASES.lock().unwrap();
            let basis_vec = basis_lock.borrow();
            assert_eq!(basis_vec.len(), 0);
        }
    }

    fn add_test(a: BasedExpr, b: BasedExpr, expected: BasedExpr) {
        assert_eq!(a.clone() + b.clone(), expected);
        assert_eq!(a.clone() + &b, expected);
        assert_eq!(&a + b.clone(), expected);
        assert_eq!(&a + &b, expected);
    }

    fn radd_test(a: BasedExpr, b: BasedExpr, expected: BasedExpr) {
        assert_eq!(b.clone() + a.clone(), expected);
        assert_eq!(b.clone() + &a, expected);
        assert_eq!(&b + a.clone(), expected);
        assert_eq!(&b + &a, expected);
    }

    fn sub_test(a: BasedExpr, b: BasedExpr, expected: BasedExpr) {
        assert_eq!(a.clone() - b.clone(), expected);
        assert_eq!(a.clone() - &b, expected);
        assert_eq!(&a - b.clone(), expected);
        assert_eq!(&a - &b, expected);
    }

    fn rsub_test(a: BasedExpr, b: BasedExpr, expected: BasedExpr) {
        assert_eq!(b.clone() - a.clone(), expected);
        assert_eq!(b.clone() - &a, expected);
        assert_eq!(&b - a.clone(), expected);
        assert_eq!(&b - &a, expected);
    }

    fn mul_test(a: BasedExpr, b: BasedExpr, expected: BasedExpr) {
        assert_eq!(a.clone() * b.clone(), expected);
        assert_eq!(a.clone() * &b, expected);
        assert_eq!(&a * b.clone(), expected);
        assert_eq!(&a * &b, expected);
    }

    fn rmul_test(a: BasedExpr, b: BasedExpr, expected: BasedExpr) {
        assert_eq!(b.clone() * a.clone(), expected);
        assert_eq!(b.clone() * &a, expected);
        assert_eq!(&b * a.clone(), expected);
        assert_eq!(&b * &a, expected);
    }

    fn div_test(a: BasedExpr, b: BasedExpr, expected: BasedExpr) {
        assert_eq!(a.clone() / b.clone(), expected);
        assert_eq!(a.clone() / &b, expected);
        assert_eq!(&a / b.clone(), expected);
        assert_eq!(&a / &b, expected);
    }

    fn rdiv_test(a: BasedExpr, b: BasedExpr, expected: BasedExpr) {
        assert_eq!(b.clone() / a.clone(), expected);
        assert_eq!(b.clone() / &a, expected);
        assert_eq!(&b / a.clone(), expected);
        assert_eq!(&b / &a, expected);
    }

    fn cmp_test(a: BasedExpr, b: BasedExpr, expected: Ordering) {
        assert_eq!(a.cmp(&b), expected);
    }

    fn rcmp_test(a: BasedExpr, b: BasedExpr, expected: Ordering) {
        assert_eq!(b.cmp(&a), expected);
    }

    // Assuming a perfectly working Rational library and just testing BasedExpr-specific behavior
    #[test]
    fn test_based_expr_add() {
        let _lock = LOCK.read();
        add_test(BasedExpr::new_undefined((-7i64).into()), BasedExpr::new_undefined((-5i64).into()), BasedExpr::new_undefined((-12i64).into()));
        add_test(BasedExpr::new_undefined(2.into()), BasedExpr::new_undefined(2.into()), BasedExpr::new_undefined(4.into()));

        add_test(BasedExpr::new_undefined((-7i64).into()), based_expr!(-5), based_expr!(-12));
        add_test(BasedExpr::new_undefined((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(-4 + sqrt 2));
        add_test(BasedExpr::new_undefined(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(121 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6));

        radd_test(BasedExpr::new_undefined((-7i64).into()), based_expr!(-5), based_expr!(-12));
        radd_test(BasedExpr::new_undefined((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(-4 + sqrt 2));
        radd_test(BasedExpr::new_undefined(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(121 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6));

        add_test(based_expr!(4), based_expr!(8), based_expr!(12));
        add_test(based_expr!(1 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(4 + 3 sqrt 2));
        add_test(based_expr!(1/2 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(7/2 + 3 sqrt 2));
        add_test(based_expr!(1 + sqrt 2 - 4 sqrt 5 + 8 sqrt 10), based_expr!(-1 + 100 sqrt 2 + 50 sqrt 5 + 7 sqrt 10),
            based_expr!(0 + 101 sqrt 2 + 46 sqrt 5 + 15 sqrt 10));
    }

    #[test]
    fn test_based_expr_sub() {
        let _lock = LOCK.read();
        sub_test(BasedExpr::new_undefined((-7i64).into()), BasedExpr::new_undefined((-5i64).into()), BasedExpr::new_undefined((-2i64).into()));
        sub_test(BasedExpr::new_undefined(2.into()), BasedExpr::new_undefined(2.into()), BasedExpr::new_undefined(0.into()));

        sub_test(BasedExpr::new_undefined((-7i64).into()), based_expr!(-5), based_expr!(-2));
        sub_test(BasedExpr::new_undefined((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(-6 - sqrt 2));
        sub_test(BasedExpr::new_undefined(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(119 - 2 sqrt 2 - 3 sqrt 3 - 6 sqrt 6));

        rsub_test(BasedExpr::new_undefined((-7i64).into()), based_expr!(-5), based_expr!(2));
        rsub_test(BasedExpr::new_undefined((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(6 + sqrt 2));
        rsub_test(BasedExpr::new_undefined(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(-119 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6));

        sub_test(based_expr!(4), based_expr!(8), based_expr!(-4));
        sub_test(based_expr!(1 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(-2 - sqrt 2));
        sub_test(based_expr!(1/2 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(-5/2 - sqrt 2));
        sub_test(based_expr!(1 + sqrt 2 - 4 sqrt 5 + 8 sqrt 10), based_expr!(-1 + 100 sqrt 2 + 50 sqrt 5 + 7 sqrt 10),
            based_expr!(2 - 99 sqrt 2 - 54 sqrt 5 + sqrt 10));
    }

    #[test]
    fn test_based_expr_mul() {
        let _lock = LOCK.read();
        mul_test(BasedExpr::new_undefined((-7i64).into()), BasedExpr::new_undefined((-5i64).into()), BasedExpr::new_undefined(35.into()));
        mul_test(BasedExpr::new_undefined(2.into()), BasedExpr::new_undefined(2.into()), BasedExpr::new_undefined(4.into()));

        mul_test(BasedExpr::new_undefined((-7i64).into()), based_expr!(-5), based_expr!(35));
        mul_test(BasedExpr::new_undefined((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(-5 - 5 sqrt 2));
        mul_test(BasedExpr::new_undefined(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(120 + 240 sqrt 2 + 360 sqrt 3 + 720 sqrt 6));

        rmul_test(BasedExpr::new_undefined((-7i64).into()), based_expr!(-5), based_expr!(35));
        rmul_test(BasedExpr::new_undefined((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(-5 - 5 sqrt 2));
        rmul_test(BasedExpr::new_undefined(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(120 + 240 sqrt 2 + 360 sqrt 3 + 720 sqrt 6));

        mul_test(based_expr!(4), based_expr!(8), based_expr!(32));
        mul_test(based_expr!(1 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(7 + 5 sqrt 2));
        mul_test(based_expr!(1/2 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(11/2 + 4 sqrt 2));
        mul_test(based_expr!(1 + sqrt 2 - 4 sqrt 5 + 8 sqrt 10), based_expr!(-1 + 100 sqrt 2 + 50 sqrt 5 + 7 sqrt 10),
            based_expr!(-241 + 1959 sqrt 2 + 1668 sqrt 5 - 351 sqrt 10));
    }

    #[test]
    fn test_based_expr_div() {
        let _lock = LOCK.read();
        div_test(BasedExpr::new_undefined((-7i64).into()), BasedExpr::new_undefined((-5i64).into()),
            BasedExpr::new_undefined(Rat::from_signeds(7, 5)));
        div_test(BasedExpr::new_undefined(2.into()), BasedExpr::new_undefined(2.into()), BasedExpr::new_undefined(1.into()));

        div_test(BasedExpr::new_undefined((-7i64).into()), based_expr!(-5), based_expr!(7/5));
        div_test(BasedExpr::new_undefined((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(5 - 5 sqrt 2));
        div_test(BasedExpr::new_undefined(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(60/91 - 120/91 sqrt 2 - 180/91 sqrt 3 + 360/91 sqrt 6));

        rdiv_test(BasedExpr::new_undefined((-7i64).into()), based_expr!(-5), based_expr!(5/7));
        rdiv_test(BasedExpr::new_undefined((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(-1/5 - 1/5 sqrt 2));
        rdiv_test(BasedExpr::new_undefined(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(1/120 + 1/60 sqrt 2 + 1/40 sqrt 3 + 1/20 sqrt 6));

        div_test(based_expr!(4), based_expr!(8), based_expr!(1/2));
        div_test(based_expr!(1 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(-1 + sqrt 2));
        div_test(based_expr!(1/2 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(-5/2 + 2 sqrt 2));
        div_test(based_expr!(1 + sqrt 2 - 4 sqrt 5 + 8 sqrt 10), based_expr!(-1 + 100 sqrt 2 + 50 sqrt 5 + 7 sqrt 10),
        // Kaboom!
            based_expr!(-8551371/21774121 - 9982071/21774121 sqrt 2 + 7355940/21774121 sqrt 5 + 2437885/21774121 sqrt 10));
    }

    #[test]
    fn test_based_expr_div_even_more_complicated_basis_simpler() {
        let _lock = LOCK.read();
        // The dreaded 5.625°. Is it even possible?!?!
        div_test(
            based_expr!(1 + 0 sqrt 2 + 0 sqrt(2 + sqrt 2) + 0 sqrt(2 - sqrt 2)
                + 0 sqrt(2 + sqrt(2 + sqrt 2)) + 0 sqrt(2 + sqrt(2 - sqrt 2)) + 0 sqrt(2 - sqrt(2 + sqrt 2)) + 0 sqrt(2 - sqrt(2 - sqrt 2))),
            based_expr!(0 + 0 sqrt 2 + 0 sqrt(2 + sqrt 2) + 0 sqrt(2 - sqrt 2)
                + 1 sqrt(2 + sqrt(2 + sqrt 2)) + 0 sqrt(2 + sqrt(2 - sqrt 2)) + 0 sqrt(2 - sqrt(2 + sqrt 2)) + 0 sqrt(2 - sqrt(2 - sqrt 2))),
        based_expr!(
                + 0
                + 0 sqrt 2
                + 0 sqrt(2 + sqrt 2)
                + 0 sqrt(2 - sqrt 2)
                + 1/2 sqrt(2 + sqrt(2 + sqrt 2))
                - 1/2 sqrt(2 + sqrt(2 - sqrt 2))
                - 1/2 sqrt(2 - sqrt(2 + sqrt 2))
                + 1/2 sqrt(2 - sqrt(2 - sqrt 2))
            )
        );
    }

    #[test]
    fn test_based_expr_div_even_more_complicated_basis() {
        let _lock = LOCK.read();
        // The dreaded 5.625°. Is it even possible?!?!
        div_test(
            based_expr!(4 + 0 sqrt 2 + 0 sqrt(2 + sqrt 2) + 0 sqrt(2 - sqrt 2)
                + 0 sqrt(2 + sqrt(2 + sqrt 2)) + 0 sqrt(2 + sqrt(2 - sqrt 2)) + 0 sqrt(2 - sqrt(2 + sqrt 2)) + 0 sqrt(2 - sqrt(2 - sqrt 2))),
            based_expr!(0 + 1 sqrt 2 + 2 sqrt(2 + sqrt 2) + 3 sqrt(2 - sqrt 2)
                + 4 sqrt(2 + sqrt(2 + sqrt 2)) + 5 sqrt(2 + sqrt(2 - sqrt 2)) + 6 sqrt(2 - sqrt(2 + sqrt 2)) + 7 sqrt(2 - sqrt(2 - sqrt 2))),
            based_expr!(
                - 7719012/6424543
                - 5185348/6424543 sqrt 2
                + 7011758/6424543 sqrt(2 + sqrt 2)
                + 2347666/6424543 sqrt(2 - sqrt 2)
                + 1738972/6424543 sqrt(2 + sqrt(2 + sqrt 2))
                - 4398496/6424543 sqrt(2 + sqrt(2 - sqrt 2))
                - 5442656/6424543 sqrt(2 - sqrt(2 + sqrt 2))
                + 6380064/6424543 sqrt(2 - sqrt(2 - sqrt 2))
            )
        );
    }

    #[test]
    fn test_based_expr_cmp() {
        let _lock = LOCK.read();

        cmp_test(BasedExpr::new_undefined((-7i64).into()), BasedExpr::new_undefined((-5i64).into()), Ordering::Less);
        cmp_test(BasedExpr::new_undefined(4.into()), BasedExpr::new_undefined(0.into()), Ordering::Greater);
        cmp_test(BasedExpr::new_undefined(1.into()), BasedExpr::new_undefined(1.into()), Ordering::Equal);

        cmp_test(BasedExpr::new_undefined(7.into()), based_expr!(0 + 5 sqrt 2), Ordering::Less);
        cmp_test(BasedExpr::new_undefined(6.into()), based_expr!(6), Ordering::Equal);
        cmp_test(BasedExpr::new_undefined(6.into()), based_expr!(6 + 0 sqrt 2), Ordering::Equal);

        rcmp_test(BasedExpr::new_undefined(7.into()), based_expr!(0 + 5 sqrt 2), Ordering::Greater);
        rcmp_test(BasedExpr::new_undefined(6.into()), based_expr!(6), Ordering::Equal);
        rcmp_test(BasedExpr::new_undefined(6.into()), based_expr!(6 + 0 sqrt 2), Ordering::Equal);

        cmp_test(based_expr!(1 + 2 sqrt 3 - 5 sqrt 5 + 0 sqrt 15), based_expr!(1 + 2 sqrt 3 - 5 sqrt 5 + 0 sqrt 15), Ordering::Equal);
        cmp_test(based_expr!(3), based_expr!(5/2), Ordering::Greater);
        cmp_test(based_expr!(5/2), based_expr!(3), Ordering::Less);
        cmp_test(based_expr!(4), based_expr!(4), Ordering::Equal);
        cmp_test(based_expr!(2 + sqrt 2), based_expr!(0 + 0 sqrt 2), Ordering::Greater);
        cmp_test(based_expr!(2 + 2 sqrt 2), based_expr!(0 + 0 sqrt 2), Ordering::Greater);
        cmp_test(based_expr!(2 - sqrt 2), based_expr!(0 + 0 sqrt 2), Ordering::Greater);
        cmp_test(based_expr!(2 - 2 sqrt 2), based_expr!(0 + 0 sqrt 2), Ordering::Less);
        cmp_test(based_expr!(-2 + sqrt 2), based_expr!(0 + 0 sqrt 2), Ordering::Less);
        cmp_test(based_expr!(-2 + 2 sqrt 2), based_expr!(0 + 0 sqrt 2), Ordering::Greater);
        cmp_test(based_expr!(-2 - sqrt 2), based_expr!(0 + 0 sqrt 2), Ordering::Less);
        cmp_test(based_expr!(-2 - 2 sqrt 2), based_expr!(0 + 0 sqrt 2), Ordering::Less);
        // Nice easily separable cases
        cmp_test(based_expr!( 16 + 11 sqrt 2 + 7 sqrt 5 + 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Greater);
        cmp_test(based_expr!( 16 + 11 sqrt 2 + 7 sqrt 5 - 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Greater);
        cmp_test(based_expr!( 16 + 11 sqrt 2 - 7 sqrt 5 + 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Greater);
        cmp_test(based_expr!( 16 + 11 sqrt 2 - 7 sqrt 5 - 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Greater);
        cmp_test(based_expr!( 16 - 11 sqrt 2 + 7 sqrt 5 + 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Greater);
        cmp_test(based_expr!( 16 - 11 sqrt 2 + 7 sqrt 5 - 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Greater);
        cmp_test(based_expr!( 16 - 11 sqrt 2 - 7 sqrt 5 + 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Greater);
        cmp_test(based_expr!( 16 - 11 sqrt 2 - 7 sqrt 5 - 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Less);
        cmp_test(based_expr!(-16 + 11 sqrt 2 + 7 sqrt 5 + 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Greater);
        cmp_test(based_expr!(-16 + 11 sqrt 2 + 7 sqrt 5 - 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Less);
        cmp_test(based_expr!(-16 + 11 sqrt 2 - 7 sqrt 5 + 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Less);
        cmp_test(based_expr!(-16 + 11 sqrt 2 - 7 sqrt 5 - 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Less);
        cmp_test(based_expr!(-16 - 11 sqrt 2 + 7 sqrt 5 + 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Less);
        cmp_test(based_expr!(-16 - 11 sqrt 2 + 7 sqrt 5 - 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Less);
        cmp_test(based_expr!(-16 - 11 sqrt 2 - 7 sqrt 5 + 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Less);
        cmp_test(based_expr!(-16 - 11 sqrt 2 - 7 sqrt 5 - 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Less);
    }

    #[test]
    fn test_based_expr_cmp_hard_to_separate() {
        let _lock = LOCK.read();
        // Requires level 0 (64 bits) precision
        let mut expr = based_expr!(1 + 0 sqrt 2 + 0 sqrt 5 - 364585791794594742/1152921504606846976 sqrt 10);
        cmp_test(expr.clone(), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Greater);
        cmp_test(-expr, based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Less);
    }

    #[test]
    fn test_based_expr_cmp_tiny_rational() {
        let _lock = LOCK.read();
        // A devilish case: requires quite a lot of precision with fixed point
        let mut expr = based_expr!(1 / 18446744073709551615 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6);
        expr = &expr * &expr; // 128
        expr = &expr * &expr; // 256
        expr = &expr * &expr; // 512
        expr = &expr * &expr; // 1024
        expr = &expr * &expr; // 2048. Now it's too small for a f64.
        cmp_test(expr.clone(), based_expr!(0 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6), Ordering::Greater);
        cmp_test(-expr, based_expr!(0 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6), Ordering::Less);
    }
}