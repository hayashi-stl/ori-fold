use std::cell::RefCell;
use std::cmp::Ordering;
use std::fmt::{Debug, Display};
use std::ops::{Add, AddAssign, Deref, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::{Arc, LazyLock, Mutex};

use malachite::base::num::basic::traits::{One, Zero};
use malachite::{rational::Rational, Integer};
use malachite::base::num::arithmetic::traits::{Abs, NegAssign, Sign};
use nalgebra::{DMatrix, DVector, DVectorSlice, DVectorView};

use crate::pslq::{checked_integer_relation, checked_integer_relation_product, IntegerRelationError};

mod pslq;
mod algebraic;

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
                malachite::rational::Rational::from_sign_and_naturals($as, malachite::Natural::from($an as u64), malachite::Natural::from($ad as u64)),
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
#[derive(Clone, Debug)]
pub struct Basis {
    pub exprs: Vec<SqrtExpr>,
    /// Matrices stored sparsely, used for multiplication
    /// Specifically, for a component i, when multiplying a by b,
    /// we have (a*b)[i] = a * matrices[i] * b.
    /// The matrices are stored sparsely as (value, component to multiply in a, component to multiply in b)
    pub matrices: Vec<Vec<(Rational, usize, usize)>>,
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
    fn calc_matrices(basis: &[SqrtExpr]) -> Result<Vec<Vec<(Rational, usize, usize)>>, BasisError> {
        let mut result = vec![vec![]; basis.len()];
        for i in 0..basis.len() {
            for j in 0..basis.len() {
                if i == 0 {
                    // First basis term is 1.
                    result[j].push((Rational::ONE, i, j));
                    continue;
                } else if j == 0 {
                    // Second basis term is 1.
                    result[i].push((Rational::ONE, i, j));
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
                    if coeff != Rational::ZERO {
                        result[index].push((coeff.into(), i, j));
                    }
                }
            }
        }
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
        Ok(Self { exprs: basis, matrices })
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
    Undefined(Rational),
    /// The basis is the list of `SqrtExpr`. The first element of the basis must be 1.
    Based(DVector<Rational>, ArcBasis)
}

impl Default for BasedExpr {
    fn default() -> Self {
        Self::Undefined(Rational::ZERO)
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
    pub fn new_undefined(q: Rational) -> Self {
        Self::Undefined(q)
    }

    /// Constructs a BasedExpr from its terms. The basis is defined by the second values of each entry.
    pub fn from_terms_checked(terms: impl IntoIterator<Item = (Rational, SqrtExpr)>) -> Result<Self, BasisError> {
        let (mut coeffs, mut basis): (Vec<_>, Vec<_>) = terms.into_iter().unzip();
        if basis[0] != SqrtExpr::ONE {
            coeffs.insert(0, 0u64.into());
            basis.insert(0, SqrtExpr::ONE);
        }
        Ok(BasedExpr::Based(coeffs.into(), Basis::new_arc_checked(basis)?))
    }

    /// Constructs a BasedExpr from its terms. The basis is defined by the second values of each entry.
    /// Panics if the basis is proven to not be a basis.
    pub fn from_terms(terms: impl IntoIterator<Item = (Rational, SqrtExpr)>) -> Self {
        Self::from_terms_checked(terms).unwrap()
    }

    fn based_rational(q: Rational, basis: ArcBasis) -> Self {
        let mut coeffs = vec![Rational::ZERO; basis.exprs.len()];
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

    fn into_coeffs_basis(self) -> (DVector<Rational>, Option<ArcBasis>) {
        match self {
            Self::Undefined(q) => (vec![q].into(), None),
            Self::Based(coeffs, basis) => (coeffs, Some(basis))
        }
    }

    fn to_coeffs_basis(&'_ self) -> (DVectorView<'_, Rational>, Option<&'_ ArcBasis>) {
        match self {
            Self::Undefined(q) => (DVectorView::from_slice(std::slice::from_ref(q), 1), None),
            Self::Based(coeffs, basis) => (coeffs.as_view(), Some(basis))
        }
    }

    fn into_coeffs_basis_2(self, other: BasedExpr) -> (DVector<Rational>, DVector<Rational>, Option<ArcBasis>) {
        let (coeffs_a, basis_a) = self.into_coeffs_basis();
        let (coeffs_b, basis_b) = other.into_coeffs_basis();
        (coeffs_a, coeffs_b, ArcBasis::into_unified(basis_a, basis_b))
    }

    fn into_coeffs_basis_2_vr(self, other: &'_ BasedExpr) -> (DVector<Rational>, DVectorView<'_, Rational>, Option<ArcBasis>) {
        let (coeffs_a, basis_a) = self.into_coeffs_basis();
        let (coeffs_b, basis_b) = other.to_coeffs_basis();
        (coeffs_a, coeffs_b, ArcBasis::into_unified_vr(basis_a, basis_b))
    }

    fn into_coeffs_basis_2_rv(&'_ self, other: BasedExpr) -> (DVectorView<'_, Rational>, DVector<Rational>, Option<ArcBasis>) {
        let (coeffs_a, basis_a) = self.to_coeffs_basis();
        let (coeffs_b, basis_b) = other.into_coeffs_basis();
        (coeffs_a, coeffs_b, ArcBasis::into_unified_rv(basis_a, basis_b))
    }

    fn into_coeffs_basis_2_rr<'a>(&'a self, other: &'a BasedExpr) -> (DVectorView<'a, Rational>, DVectorView<'a, Rational>, Option<ArcBasis>) {
        let (coeffs_a, basis_a) = self.to_coeffs_basis();
        let (coeffs_b, basis_b) = other.to_coeffs_basis();
        (coeffs_a, coeffs_b, ArcBasis::to_unified(basis_a, basis_b).cloned())
    }

    fn to_coeffs_basis_2<'a>(&'a self, other: &'a BasedExpr) -> (DVectorView<'a, Rational>, DVectorView<'a, Rational>, Option<&'a ArcBasis>) {
        let (coeffs_a, basis_a) = self.to_coeffs_basis();
        let (coeffs_b, basis_b) = other.to_coeffs_basis();
        (coeffs_a, coeffs_b, ArcBasis::to_unified(basis_a, basis_b))
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
                    write!(f, " {} {}{}", if a.sign() == Ordering::Less {"-"} else {"+"}, a.abs(), sqrt)?
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
                &a[0] == b && a.iter().skip(1).all(|q| q == &Rational::ZERO),
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

//impl DivAssign<BasedExpr> for BasedExpr {
//    fn div_assign(&mut self, rhs: BasedExpr) {
//        if let (None, Some(basis)) = (self.basis(), rhs.basis()) {
//            *self = BasedExpr::based_rational(if let BasedExpr::Undefined(q) = std::mem::take(self) {q} else {unreachable!()}, basis.clone());
//        }
//        match (self, rhs) {
//            (BasedExpr::Undefined(a), BasedExpr::Undefined(b)) => *a /= b,
//            (BasedExpr::Based(a, _), BasedExpr::Undefined(b)) => a.iter_mut().for_each(|a| *a /= &b),
//            (BasedExpr::Undefined(_), BasedExpr::Based(_, _)) => unreachable!(),
//            (BasedExpr::Based(a, basis_a), BasedExpr::Based(b, basis_b)) => {
//                // The dreaded division
//                // TODO: Make this faster if necessary.
//                basis_a.assert_compatible(&basis_b);
//                let mut mtx = DMatrix::repeat(a.len(), a.len(), Rational::ZERO);
//                for (row_i, row) in basis_a.matrices.iter().enumerate() {
//                    for (coeff, a_i, b_i) in row {
//                        mtx[(row_i, *b_i)] = coeff * &a[*a_i];
//                    }
//                }
//                let solution = mtx.solve(DVector::from_iterator(a.len(), b));
//            }
//        }
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
        let expected = BasedExpr::from_terms(vec![(Rational::from(1i64), SqrtExpr::Int(Integer::from(1i64)))]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_negative_int() {
        let _lock = LOCK.read();
        let result = based_expr!(-3);
        let expected = BasedExpr::from_terms(vec![(Rational::from(-3i64), SqrtExpr::Int(Integer::from(1i64)))]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_rational() {
        let _lock = LOCK.read();
        let result = based_expr!(3/5);
        let expected = BasedExpr::from_terms(vec![(Rational::from_integers(Integer::from(3i64), Integer::from(5i64)), SqrtExpr::Int(Integer::from(1i64)))]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_simple_sqrt() {
        let _lock = LOCK.read();
        let result = based_expr!(sqrt 3);
        let expected = BasedExpr::from_terms(vec![(Rational::from(1i64), SqrtExpr::Int(Integer::from(3i64)))]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_plus_sqrt() {
        let _lock = LOCK.read();
        let result = based_expr!(-4 + 2 sqrt 2);
        let expected = BasedExpr::from_terms(vec![
            (Rational::from(-4i64), SqrtExpr::Int(Integer::from(1i64))),
            (Rational::from(2i64), SqrtExpr::Int(Integer::from(2i64))),
            ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_rational_sqrt() {
        let _lock = LOCK.read();
        let result = based_expr!(1/2 + 1/2 sqrt 5);
        let expected = BasedExpr::from_terms(vec![
            (Rational::from_integers(1i64.into(), 2i64.into()), SqrtExpr::Int(Integer::from(1i64))),
            (Rational::from_integers(1i64.into(), 2i64.into()), SqrtExpr::Int(Integer::from(5i64))),
            ]);
        assert_eq!(result, expected);
    }

    //#[test]
    //fn test_macro_frrt() {
    //    let _lock = LOCK.read();
    //    let result = based_expr!(sqrt (sqrt 5));
    //    let expected = BasedExpr::from_terms(vec![
    //        (Rational::from(1i64), SqrtExpr::Sum(vec![
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
            (Rational::from(-1i64), SqrtExpr::Int(Integer::from(1i64))),
            (Rational::from(-1i64), SqrtExpr::Int(Integer::from(2i64))),
            (Rational::from(2i64), SqrtExpr::Sum(vec![
                (Integer::from(2i64), SqrtExpr::Int(Integer::from(1i64))),
                (Integer::from(-1i64), SqrtExpr::Int(Integer::from(2i64))),
            ])),
            (Rational::from(1i64), SqrtExpr::Sum(vec![
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

    #[test]
    fn test_add_based_based() {
        let _lock = LOCK.read();
        let a = based_expr!(1 + sqrt 2);
        let b = based_expr!(3 + 2 sqrt 2);
        let expected = based_expr!(4 + 3 sqrt 2);
        
        let result = a.clone() + b.clone();
        assert_eq!(result, expected);
        let result = a.clone() + &b;
        assert_eq!(result, expected);
        let result = &a + b.clone();
        assert_eq!(result, expected);
        let result = &a + &b;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_add_based_undefined() {
        let _lock = LOCK.read();
        let a = based_expr!(1 + sqrt 2);
        let b = BasedExpr::new_undefined((-5i64).into());
        let expected = based_expr!(-4 + sqrt 2);
        
        let result = a.clone() + b.clone();
        assert_eq!(result, expected);
        let result = a.clone() + &b;
        assert_eq!(result, expected);
        let result = &a + b.clone();
        assert_eq!(result, expected);
        let result = &a + &b;
        assert_eq!(result, expected);
        let result = b.clone() + a.clone();
        assert_eq!(result, expected);
        let result = b.clone() + &a;
        assert_eq!(result, expected);
        let result = &b + a.clone();
        assert_eq!(result, expected);
        let result = &b + &a;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_add_undefined_undefined() {
        let _lock = LOCK.read();
        let a = BasedExpr::new_undefined((-7i64).into());
        let b = BasedExpr::new_undefined((-5i64).into());
        let expected = BasedExpr::new_undefined((-12i64).into());
        
        let result = a.clone() + b.clone();
        assert_eq!(result, expected);
        let result = a.clone() + &b;
        assert_eq!(result, expected);
        let result = &a + b.clone();
        assert_eq!(result, expected);
        let result = &a + &b;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_sub_based_based() {
        let _lock = LOCK.read();
        let a = based_expr!(1 - sqrt 5);
        let b = based_expr!(6 - 2 sqrt 5);
        let expected = based_expr!(-5 + sqrt 5);
        
        let result = a.clone() - b.clone();
        assert_eq!(result, expected);
        let result = a.clone() - &b;
        assert_eq!(result, expected);
        let result = &a - b.clone();
        assert_eq!(result, expected);
        let result = &a - &b;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_sub_based_undefined() {
        let _lock = LOCK.read();
        let a = based_expr!(1 - sqrt 17);
        let b = BasedExpr::new_undefined((-5i64).into());
        let expected = based_expr!(6 - sqrt 17);
        
        let result = a.clone() - b.clone();
        assert_eq!(result, expected);
        let result = a.clone() - &b;
        assert_eq!(result, expected);
        let result = &a - b.clone();
        assert_eq!(result, expected);
        let result = &a - &b;
        assert_eq!(result, expected);

        let expected = based_expr!(-6 + sqrt 17);
        let result = b.clone() - a.clone();
        assert_eq!(result, expected);
        let result = b.clone() - &a;
        assert_eq!(result, expected);
        let result = &b - a.clone();
        assert_eq!(result, expected);
        let result = &b - &a;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_sub_undefined_undefined() {
        let _lock = LOCK.read();
        let a = BasedExpr::new_undefined((-7i64).into());
        let b = BasedExpr::new_undefined((-5i64).into());
        let expected = BasedExpr::new_undefined((-2i64).into());
        
        let result = a.clone() - b.clone();
        assert_eq!(result, expected);
        let result = a.clone() - &b;
        assert_eq!(result, expected);
        let result = &a - b.clone();
        assert_eq!(result, expected);
        let result = &a - &b;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_mul_based_based() {
        let _lock = LOCK.read();
        let a = based_expr!(1 - sqrt 5);
        let b = based_expr!(6 - 2 sqrt 5);
        let expected = based_expr!(16 - 8 sqrt 5);
        
        let result = a.clone() * b.clone();
        assert_eq!(result, expected);
        let result = a.clone() * &b;
        assert_eq!(result, expected);
        let result = &a * b.clone();
        assert_eq!(result, expected);
        let result = &a * &b;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_mul_based_based_by_0() {
        let _lock = LOCK.read();
        let a = based_expr!(1 - sqrt 5);
        let b = based_expr!(0 + 0 sqrt 5);
        let expected = based_expr!(0 + 0 sqrt 5);
        
        let result = a.clone() * b.clone();
        assert_eq!(result, expected);
        let result = a.clone() * &b;
        assert_eq!(result, expected);
        let result = &a * b.clone();
        assert_eq!(result, expected);
        let result = &a * &b;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_mul_based_based_more_complicated_basis() {
        let _lock = LOCK.read();
        // Requires PSLQ to prove closure, at least for now.
        let a = based_expr!(0 + 0 sqrt 2 + 1 sqrt(2 + sqrt 2) + 0 sqrt(2 - sqrt 2));
        let b = based_expr!(0 + 1 sqrt 2 + 0 sqrt(2 + sqrt 2) + 0 sqrt(2 - sqrt 2));
        let expected = based_expr!(0 + 0 sqrt 2 + 1 sqrt(2 + sqrt 2) + 1 sqrt(2 - sqrt 2));
        
        let result = a.clone() * b.clone();
        assert_eq!(result, expected);
        let result = a.clone() * &b;
        assert_eq!(result, expected);
        let result = &a * b.clone();
        assert_eq!(result, expected);
        let result = &a * &b;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_mul_based_undefined() {
        let _lock = LOCK.read();
        let a = based_expr!(1 - sqrt 17);
        let b = BasedExpr::new_undefined((-5i64).into());
        let expected = based_expr!(-5 + 5 sqrt 17);
        
        let result = a.clone() * b.clone();
        assert_eq!(result, expected);
        let result = a.clone() * &b;
        assert_eq!(result, expected);
        let result = &a * b.clone();
        assert_eq!(result, expected);
        let result = &a * &b;
        assert_eq!(result, expected);
        let result = b.clone() * a.clone();
        assert_eq!(result, expected);
        let result = b.clone() * &a;
        assert_eq!(result, expected);
        let result = &b * a.clone();
        assert_eq!(result, expected);
        let result = &b * &a;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_mul_undefined_undefined() {
        let _lock = LOCK.read();
        let a = BasedExpr::new_undefined((-121i64).into());
        let b = BasedExpr::new_undefined((11i64).into());
        let expected = BasedExpr::new_undefined((-1331i64).into());
        
        let result = a.clone() * b.clone();
        assert_eq!(result, expected);
        let result = a.clone() * &b;
        assert_eq!(result, expected);
        let result = &a * b.clone();
        assert_eq!(result, expected);
        let result = &a * &b;
        assert_eq!(result, expected);
    }
}