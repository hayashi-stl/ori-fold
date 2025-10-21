use std::cmp::Ordering;
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::sync::{LazyLock, RwLock};

pub use malachite;
use malachite::base::num::basic::traits::{Zero as _};
use malachite::Integer;
use malachite::base::num::arithmetic::traits::{Abs, Sign, Square};
use nalgebra::{DVector, DVectorView};
use num::{Num, One, Signed, Zero};

use crate::basis::{ArcBasis, Basis, BasisError, SqrtExpr};
use crate::interval::{Fixed, Interval};
use crate::rat::Rat;

mod pslq;
mod algebraic;
pub mod rat;
mod interval;
pub mod basis;
pub mod math;
pub mod conversion;
mod parse;
//mod matrix;

#[macro_export]
macro_rules! sqrt_expr_helper {
    ($( { $as:literal $an:literal sqrt $($sqrt:tt)* } )*) => {
        $crate::basis::SqrtExpr::Sum(vec![$(
            (
                $crate::malachite::Integer::from_sign_and_abs($as, ($an as u64).into()),
                $crate::sqrt_expr_parse!($($sqrt)*)
            ),
        )*])
    };
}

#[macro_export]
macro_rules! sqrt_expr_parse {
    ($int:literal) => { $crate::basis::SqrtExpr::Int($crate::malachite::Integer::from($int as i64)) };
    ($( { $($init:tt)* } )* + sqrt $($tt:tt)*) => { $crate::sqrt_expr_parse!($( { $($init)* } )* [+] 1 sqrt $($tt)*) };
    ($( { $($init:tt)* } )* - sqrt $($tt:tt)*) => { $crate::sqrt_expr_parse!($( { $($init)* } )* [-] 1 sqrt $($tt)*) };
    ($( { $($init:tt)* } )* + $($tt:tt)*) => { $crate::sqrt_expr_parse!($( { $($init)* } )* [+] $($tt)*) };
    ($( { $($init:tt)* } )* - $($tt:tt)*) => { $crate::sqrt_expr_parse!($( { $($init)* } )* [-] $($tt)*) };
    ($( { $($init:tt)* } )* sqrt $($tt:tt)*) => { $crate::sqrt_expr_parse!($( { $($init)* } )* [+] 1 sqrt $($tt)*) };
    ($an:literal $($tt:tt)*) => { $crate::sqrt_expr_parse!(+ $an $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal) => { $crate::sqrt_expr_parse!($( { $($init)* } )* [+] $an sqrt(1)) };
    ($( { $($init:tt)* } )* [-] $an:literal) => { $crate::sqrt_expr_parse!($( { $($init)* } )* [-] $an sqrt(1)) };
    ($( { $($init:tt)* } )* [+] $an:literal + $($tt:tt)*) => { $crate::sqrt_expr_parse!($( { $($init)* } )* [+] $an sqrt(1) + $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal + $($tt:tt)*) => { $crate::sqrt_expr_parse!($( { $($init)* } )* [-] $an sqrt(1) + $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal - $($tt:tt)*) => { $crate::sqrt_expr_parse!($( { $($init)* } )* [+] $an sqrt(1) - $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal - $($tt:tt)*) => { $crate::sqrt_expr_parse!($( { $($init)* } )* [-] $an sqrt(1) - $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal sqrt $sqrt:literal $($tt:tt)*) => { $crate::sqrt_expr_parse!($( { $($init)* } )* [+] $an sqrt ($sqrt) $($tt)* ) };
    ($( { $($init:tt)* } )* [-] $an:literal sqrt $sqrt:literal $($tt:tt)*) => { $crate::sqrt_expr_parse!($( { $($init)* } )* [-] $an sqrt ($sqrt) $($tt)* ) };
    ($( { $($init:tt)* } )* [+] $an:literal sqrt ($($sqrt:tt)*) $($tt:tt)*) => { $crate::sqrt_expr_parse!($( { $($init)* } )* { true $an sqrt $($sqrt)* } $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal sqrt ($($sqrt:tt)*) $($tt:tt)*) => { $crate::sqrt_expr_parse!($( { $($init)* } )* { false $an sqrt $($sqrt)* } $($tt)*) };
    ($( { $($terms:tt)* } )*) => { $crate::sqrt_expr_helper!($( { $($terms)* } )*)}
}

#[macro_export]
macro_rules! sqrt_expr {
    ($($tt:tt)*) => {
        $crate::sqrt_expr_parse!($($tt)*)
    };
}

#[macro_export]
macro_rules! expr_helper {
    ($( { $as:literal $an:literal / $ad:literal sqrt $($sqrt:tt)* })*) => {
        $crate::BasedExpr::from_terms(vec![$(
            (
                $crate::rat::Rat::from_sign_and_naturals($as, $crate::malachite::Natural::from($an as u64), $crate::malachite::Natural::from($ad as u64)),
                $crate::sqrt_expr_parse!($($sqrt)*)
            ),
        )*])
    };
}

#[macro_export]
macro_rules! expr_parse {
    ($( { $($init:tt)* } )* + sqrt $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [+] 1 / 1 sqrt $($tt)*) };
    ($( { $($init:tt)* } )* - sqrt $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [-] 1 / 1 sqrt $($tt)*) };
    ($( { $($init:tt)* } )* + $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [+] $($tt)*) };
    ($( { $($init:tt)* } )* - $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [-] $($tt)*) };
    ($( { $($init:tt)* } )* sqrt $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [+] 1 / 1 sqrt $($tt)*) };
    ($an:literal $($tt:tt)*) => { $crate::expr_parse!(+ $an $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal) => { $crate::expr_parse!($( { $($init)* } )* [+] $an / 1) };
    ($( { $($init:tt)* } )* [-] $an:literal) => { $crate::expr_parse!($( { $($init)* } )* [-] $an / 1) };
    ($( { $($init:tt)* } )* [+] $an:literal + $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [+] $an / 1 + $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal + $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [-] $an / 1 + $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal - $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [+] $an / 1 - $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal - $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [-] $an / 1 - $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal sqrt $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [+] $an / 1 sqrt $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal sqrt $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [-] $an / 1 sqrt $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal / $ad:literal) => { $crate::expr_parse!($( { $($init)* } )* [+] $an / $ad sqrt(1)) };
    ($( { $($init:tt)* } )* [-] $an:literal / $ad:literal) => { $crate::expr_parse!($( { $($init)* } )* [-] $an / $ad sqrt(1)) };
    ($( { $($init:tt)* } )* [+] $an:literal / $ad:literal + $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [+] $an / $ad sqrt(1) + $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal / $ad:literal + $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [-] $an / $ad sqrt(1) + $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal / $ad:literal - $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [+] $an / $ad sqrt(1) - $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal / $ad:literal - $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [-] $an / $ad sqrt(1) - $($tt)*) };
    ($( { $($init:tt)* } )* [+] $an:literal / $ad:literal sqrt $sqrt:literal $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [+] $an / $ad sqrt ($sqrt) $($tt)* ) };
    ($( { $($init:tt)* } )* [-] $an:literal / $ad:literal sqrt $sqrt:literal $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* [-] $an / $ad sqrt ($sqrt) $($tt)* ) };
    ($( { $($init:tt)* } )* [+] $an:literal / $ad:literal sqrt ($($sqrt:tt)*) $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* { true $an / $ad sqrt $($sqrt)* } $($tt)*) };
    ($( { $($init:tt)* } )* [-] $an:literal / $ad:literal sqrt ($($sqrt:tt)*) $($tt:tt)*) => { $crate::expr_parse!($( { $($init)* } )* { false $an / $ad sqrt $($sqrt)* } $($tt)*) };
    ($( { $($terms:tt)* } )*) => { $crate::expr_helper!($( { $($terms)* } )*)}
}


/// A macro for generating expression "literals".
/// 
/// Some examples:
/// ```rust
/// based_expr!(1)
/// based_expr!(1/2)
/// based_expr!(1 + sqrt 2)
/// based_expr!(-3 + 2/5 sqrt 6)
/// based_expr!(0 + 0 sqrt 2 + 0 sqrt(2 + sqrt 2) + 1 sqrt(2 - sqrt 2))
/// ```
/// At least for now, a full basis must be specified, so you're not allowed to do,
/// ```rust
/// based_expr!(sqrt 17)
/// ```
/// for example, because $[\sqrt{17}]$ does not form a closed basis. You must do
/// ```rust
/// based_expr!(0 + sqrt 17)
/// ```
/// In addition, the first element of the basis must be 1, so
/// ```rust
/// based_expr!(sqrt 2 + 1)
/// ```
/// is not allowed either.
#[macro_export]
macro_rules! based_expr {
    ($($tt:tt)*) => {
        $crate::expr_parse!($($tt)*)
    };
}

/// An expression of an algebraic number with a basis. Technically a member of a
/// field extension of the real numbers.
/// 
/// For example, $3 + 2\sqrt{2}$ has the basis $\left\langle 1, \sqrt{2}\right\rangle$
/// with coefficients $\\left\langle 3, 2\right\rangle$.
/// 
/// It is assumed that the first element of the basis is 1.
/// 
/// Two values with different bases are incompatible and mixing them results in a panic.
/// Treat them like different types.
/// 
/// Sometimes, a value is needed in a context where no basis can be provided. In such
/// cases, a *baseless number* (`BasedExpr::Baseless`) is created. This is
/// compatible with any basis, but any operation done where a basis is involved
/// will result in a *based number* (`BasedExpr::Based`).
/// 
/// For this type, the `%` operator *floors*, which is not the behavior
/// for good ol' primitive integers.
#[derive(Clone)]
pub enum BasedExpr {
    /// The basis is unknown. Used when a value is needed and there's nowhere to get a basis from
    Baseless(Rat),
    /// The basis is the list of `SqrtExpr`. The first element of the basis is 1.
    Based(DVector<Rat>, ArcBasis)
}

impl Default for BasedExpr {
    /// The default value, a baseless 0.
    fn default() -> Self {
        Self::BASELESS_ZERO
    }
}

fn longer_first<T>(lhs: Vec<T>, rhs: Vec<T>) -> (Vec<T>, Vec<T>) {
    if lhs.len() < rhs.len() { (rhs, lhs) } else { (lhs, rhs) }
}

fn longer_first_ref<'a, T>(lhs: &'a [T], rhs: &'a [T]) -> (&'a [T], &'a [T]) {
    if lhs.len() < rhs.len() { (rhs, lhs) } else { (lhs, rhs) }
}

impl BasedExpr {
    pub const BASELESS_ZERO: Self = Self::Baseless(Rat::ZERO);
    pub const BASELESS_ONE: Self = Self::Baseless(Rat::ONE);

    /// Constructs a baseless BasedExpr (one without a defined basis)
    pub fn new_baseless(q: Rat) -> Self {
        Self::Baseless(q)
    }

    /// Constructs a BasedExpr from its terms. The basis is defined by the second values of each entry.
    /// 
    /// Remember, the first term of the basis must be 1.
    /// 
    /// Returns an `Err` if the basis:
    /// * has something other than 1 as the first term, or
    /// * is proven to have a linear dependency, or
    /// * is not proven to be closed under multiplication
    pub fn try_from_terms(terms: impl IntoIterator<Item = (Rat, SqrtExpr)>) -> Result<Self, BasisError> {
        let (mut coeffs, mut basis): (Vec<_>, Vec<_>) = terms.into_iter().unzip();
        if basis[0] != SqrtExpr::ONE {
            coeffs.insert(0, 0u64.into());
            basis.insert(0, SqrtExpr::ONE);
        }
        Ok(BasedExpr::Based(coeffs.into(), Basis::new_arc_checked(basis)?))
    }

    /// Constructs a BasedExpr from its terms. The basis is defined by the second values of each entry.
    /// 
    /// Remember, the first term of the basis must be 1.
    /// 
    /// # Panics
    /// Panics if the basis is:
    /// * has something other than 1 as the first term, or
    /// * is proven to have a linear dependency, or
    /// * is not proven to be closed under multiplication
    pub fn from_terms(terms: impl IntoIterator<Item = (Rat, SqrtExpr)>) -> Self {
        Self::try_from_terms(terms).unwrap()
    }

    pub fn try_from_coeffs_and_basis(mut coeffs: Vec<Rat>, basis: ArcBasis) -> Option<Self> {
        if coeffs.len() > basis.len() { None? }
        coeffs.resize(basis.len(), Rat::ZERO);
        Some(BasedExpr::Based(DVector::from_vec(coeffs), basis))
    }

    pub fn based_zero(basis: ArcBasis) -> Self {
        Self::based_rational(Rat::ZERO, basis)
    }

    pub fn based_one(basis: ArcBasis) -> Self {
        Self::based_rational(Rat::ONE, basis)
    }

    /// Gets the sign of this expression.
    /// Returns `Ordering::Greater` if greater than 0, `Ordering::Equal` if equal to 0,
    /// and `Ordering::Less` if less than 0.
    pub fn sign(&self) -> Ordering {
        match self {
            Self::Baseless(q) => q.sign(),
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

    // impl for<T: IntervalContent> Fn(Interval<T>) -> Option<U> is not allowed
    fn narrow_down<U>(coeffs: &DVector<Rat>, basis: &Basis,
        result_fn_f64: Option<impl FnOnce(Interval<f64>) -> Option<U>>, mut result_fn_fixed: impl FnMut(Interval<Fixed>) -> Option<U>) -> U
    {
        if let Some(result_fn_f64) = result_fn_f64 {
            let float_approx = Self::interval_f64(&coeffs, basis);
            let result = result_fn_f64(float_approx);
            if let Some(result) = result { return result };
        }

        let mut level = 0;
        // Increase the interval; we need separation!
        loop {
            let approx = Self::interval_fixed(&coeffs, basis, level);
            let result = result_fn_fixed(approx);
            if let Some(result) = result { return result };
            level += 1;
        }
    }

    pub fn based_rational(q: Rat, basis: ArcBasis) -> Self {
        let mut coeffs = vec![Rat::ZERO; basis.exprs.len()];
        coeffs[0] = q;
        Self::Based(coeffs.into(), basis)
    }

    fn has_basis(&self) -> bool {
        match self {
            Self::Baseless(_) => false,
            Self::Based(_, basis) => true
        }
    }

    fn basis(&self) -> Option<&ArcBasis> {
        match self {
            Self::Baseless(_) => None,
            Self::Based(_, basis) => Some(basis)
        }
    }

    fn with_basis(self, basis: Option<ArcBasis>) -> Self {
        if let Some(basis) = basis {
            match self {
                Self::Baseless(q) => {
                    let mut coeffs = DVector::repeat(basis.exprs.len(), Rat::ZERO);
                    coeffs[0] = q;
                    Self::Based(coeffs, basis)
                },
                Self::Based(coeffs, b) => {
                    b.assert_compatible(&basis);
                    Self::Based(coeffs, b)
                }
            }
        } else { self }
    }

    /// Generates a 0 with `self`'s basis (or a baseless 0 if `self` is baseless)
    pub fn into_zero(mut self) -> Self {
        self.set_zero();
        self
    }

    /// Generates a 0 with `self`'s basis (or a baseless 0 if `self` is baseless)
    pub fn to_zero(&self) -> Self {
        Self::BASELESS_ZERO.with_basis(self.basis().cloned())
    }

    /// Generates a 1 with `self`'s basis (or a baseless 1 if `self` is baseless)
    pub fn into_one(mut self) -> Self {
        self.set_one();
        self
    }

    /// Generates a 1 with `self`'s basis (or a baseless 1 if `self` is baseless)
    pub fn to_one(&self) -> Self {
        Self::BASELESS_ONE.with_basis(self.basis().cloned())
    }

    /// Compares this to 0
    pub fn cmp_zero(&self) -> Ordering {
        match self {
            BasedExpr::Baseless(a) => a.cmp(&Rat::ZERO),
            BasedExpr::Based(a, basis) => {
                // Easy cases
                if a.iter().all(|a| a == &Rat::ZERO) {
                    return Ordering::Equal;
                } if a.len() == 1 {
                    return a[0].cmp(&Rat::ZERO);
                } else if a.len() == 2 {
                    // The basis has to be [1, sqrt(c)] for some square-free integer c >= 2
                    // However, it could be represented as, e.g. [1, sqrt(sqrt(4))].
                    // Why you would do that, I don't know.
                    if let SqrtExpr::Int(c) = &basis.exprs[1] {
                        return ((&a[0]).signum() * (&a[0]).square() + (&a[1]).signum() * (&a[1]).square() * Rat::from(c)).cmp(&Rat::ZERO);
                    }
                }

                Self::narrow_down(&a, &basis, Some(|k: Interval<f64>| k.cmp_zero()), |k| k.cmp_zero())
            }
        }
    }

    fn assert_compatible(&self, other: &BasedExpr) -> () {
        if let (Self::Based(_, basis_a), Self::Based(_, basis_b)) = (self, other) {
            basis_a.assert_compatible(basis_b);
        }
    }

    fn into_coeffs_basis(self) -> (DVector<Rat>, Option<ArcBasis>) {
        match self {
            Self::Baseless(q) => (vec![q].into(), None),
            Self::Based(coeffs, basis) => (coeffs, Some(basis))
        }
    }

    fn to_coeffs_basis(&'_ self) -> (DVectorView<'_, Rat>, Option<&'_ ArcBasis>) {
        match self {
            Self::Baseless(q) => (DVectorView::from_slice(std::slice::from_ref(q), 1), None),
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
    /// Generates a baseless 0
    fn zero() -> Self {
        Self::new_baseless(Rat::ZERO)
    }

    fn is_zero(&self) -> bool {
        match self {
            Self::Baseless(q) => q == &Rat::ZERO,
            Self::Based(coeffs, _) => coeffs.iter().all(|c| c == &Rat::ZERO)
        }
    }

    /// Sets this to 0, preserving the basis if there is one
    fn set_zero(&mut self) {
        match self {
            Self::Baseless(q) => *q = Rat::ZERO,
            Self::Based(coeffs, _) => coeffs.fill(Rat::ZERO)
        }
    }
}

impl One for BasedExpr {
    /// Generates a baseless 1
    fn one() -> Self {
        Self::new_baseless(Rat::ONE)
    }
    
    /// Sets this to 1, preserving the basis if there is one
    fn set_one(&mut self) {
        match self {
            Self::Baseless(q) => *q = Rat::ONE,
            Self::Based(coeffs, _) => {
                coeffs.fill(Rat::ZERO);
                coeffs[0] = Rat::ONE
            }
        }
    }
}

impl Display for BasedExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Baseless(q) => write!(f, "{}", q),
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
            (BasedExpr::Baseless(a), BasedExpr::Baseless(b)) => a == b,
            (BasedExpr::Based(a, _), BasedExpr::Baseless(b)) =>
                &a[0] == b && a.iter().skip(1).all(|q| q == &Rat::ZERO),
            (BasedExpr::Baseless(_), BasedExpr::Based(_, _)) => unreachable!(),
            (BasedExpr::Based(a, basis_a), BasedExpr::Based(b, basis_b)) => {
                basis_a.assert_compatible(&basis_b);
                a.iter().zip(b.iter()).all(|(a, b)| a == b)
            }
        }
    }
}

impl Eq for BasedExpr {}

impl Hash for BasedExpr {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        core::mem::discriminant(self).hash(state);
        match self {
            Self::Baseless(q) => q.hash(state),
            Self::Based(coeffs, _) => coeffs.hash(state)
        }
    }
}

impl PartialOrd for BasedExpr {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let diff = self - other;
        Some(diff.cmp_zero())
    }
}

impl Ord for BasedExpr {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

#[used]
static TEST_LOCK: LazyLock<RwLock<()>> = LazyLock::new(|| RwLock::new(()));

#[cfg(test)]
mod tests {
    use std::sync::RwLock;

    use crate::basis::BASES;

    use super::*;

    #[test]
    fn test_macro_int() {
        let _lock = TEST_LOCK.read();
        let result = based_expr!(1);
        let expected = BasedExpr::from_terms(vec![(Rat::from(1i64), SqrtExpr::Int(Integer::from(1i64)))]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_negative_int() {
        let _lock = TEST_LOCK.read();
        let result = based_expr!(-3);
        let expected = BasedExpr::from_terms(vec![(Rat::from(-3i64), SqrtExpr::Int(Integer::from(1i64)))]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_rational() {
        let _lock = TEST_LOCK.read();
        let result = based_expr!(3/5);
        let expected = BasedExpr::from_terms(vec![(Rat::from_integers(Integer::from(3i64), Integer::from(5i64)), SqrtExpr::Int(Integer::from(1i64)))]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_simple_sqrt() {
        let _lock = TEST_LOCK.read();
        let result = based_expr!(sqrt 3);
        let expected = BasedExpr::from_terms(vec![(Rat::from(1i64), SqrtExpr::Int(Integer::from(3i64)))]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_plus_sqrt() {
        let _lock = TEST_LOCK.read();
        let result = based_expr!(-4 + 2 sqrt 2);
        let expected = BasedExpr::from_terms(vec![
            (Rat::from(-4i64), SqrtExpr::Int(Integer::from(1i64))),
            (Rat::from(2i64), SqrtExpr::Int(Integer::from(2i64))),
            ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_macro_rational_sqrt() {
        let _lock = TEST_LOCK.read();
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
        let _lock = TEST_LOCK.read();
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
        let _lock = TEST_LOCK.write();
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
}