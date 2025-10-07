use std::{cmp::Ordering, iter::{Product, Sum}, ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign}};

use approx::{AbsDiffEq, RelativeEq, UlpsEq};
use malachite::{base::{num::{arithmetic::traits::{Abs, AbsDiff, Ceiling, CheckedDiv, Floor, NegAssign, Pow}, conversion::traits::RoundingFrom, logic::traits::SignificantBits}, rounding_modes::RoundingMode}, Integer};
use nalgebra::{ComplexField, DMatrix, DVector, Field, RealField, SimdValue};
use num::{FromPrimitive, Signed, Zero};
use simba::scalar::SubsetOf;

use crate::{conversion::{IntoBaseless, TryBaselessFrom}, rat::Rat, BasedExpr};

impl AddAssign<BasedExpr> for BasedExpr {
    fn add_assign(&mut self, mut rhs: BasedExpr) {
        if !self.has_basis() && rhs.has_basis() {
            std::mem::swap(self, &mut rhs);
        }
        match (self, rhs) {
            (BasedExpr::Baseless(a), BasedExpr::Baseless(b)) => *a += b,
            (BasedExpr::Based(a, _), BasedExpr::Baseless(b)) => a[0] += b,
            (BasedExpr::Baseless(_), BasedExpr::Based(_, _)) => unreachable!(),
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
            (BasedExpr::Baseless(a), BasedExpr::Baseless(b)) => *a += b,
            (BasedExpr::Based(a, _), BasedExpr::Baseless(b)) => a[0] += b,
            (BasedExpr::Baseless(_), BasedExpr::Based(_, _)) => unreachable!(),
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

impl Sum<BasedExpr> for BasedExpr {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(BasedExpr::BASELESS_ZERO, |a, b| a + b)
    }
}

impl<'a> Sum<&'a BasedExpr> for BasedExpr {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(BasedExpr::BASELESS_ZERO, |a, b| a + b)
    }
}

impl NegAssign for BasedExpr {
    fn neg_assign(&mut self) {
        match self {
            BasedExpr::Baseless(q) => q.neg_assign(),
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
            (BasedExpr::Baseless(a), BasedExpr::Baseless(b)) => *a -= b,
            (BasedExpr::Based(a, _), BasedExpr::Baseless(b)) => a[0] -= b,
            (BasedExpr::Baseless(_), BasedExpr::Based(_, _)) => unreachable!(),
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
            (BasedExpr::Baseless(a), BasedExpr::Baseless(b)) => *a -= b,
            (BasedExpr::Based(a, _), BasedExpr::Baseless(b)) => a[0] -= b,
            (BasedExpr::Baseless(_), BasedExpr::Based(_, _)) => unreachable!(),
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
            (BasedExpr::Baseless(a), BasedExpr::Baseless(b)) => *a *= b,
            (BasedExpr::Based(a, _), BasedExpr::Baseless(b)) => a.iter_mut().for_each(|a| *a *= &b),
            (BasedExpr::Baseless(_), BasedExpr::Based(_, _)) => unreachable!(),
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
            (BasedExpr::Baseless(a), BasedExpr::Baseless(b)) => *a *= b,
            (BasedExpr::Based(a, _), BasedExpr::Baseless(b)) => a.iter_mut().for_each(|a| *a *= b),
            (BasedExpr::Baseless(_), BasedExpr::Based(_, _)) => unreachable!(),
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

impl Product<BasedExpr> for BasedExpr {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(BasedExpr::BASELESS_ONE, |a, b| a * b)
    }
}

impl<'a> Product<&'a BasedExpr> for BasedExpr {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(BasedExpr::BASELESS_ONE, |a, b| a * b)
    }
}

impl DivAssign<BasedExpr> for BasedExpr {
    fn div_assign(&mut self, rhs: BasedExpr) {
        *self = std::mem::take(self) / rhs;
    }
}

impl DivAssign<&BasedExpr> for BasedExpr {
    fn div_assign(&mut self, rhs: &BasedExpr) {
        *self = std::mem::take(self) / rhs;
    }
}

impl Div<BasedExpr> for BasedExpr {
    type Output = BasedExpr;

    fn div(self, rhs: BasedExpr) -> Self::Output {
        self.checked_div(rhs).expect("division by 0")
    }
}

impl Div<&BasedExpr> for BasedExpr {
    type Output = BasedExpr;

    fn div(self, rhs: &BasedExpr) -> Self::Output {
        self.checked_div(rhs).expect("division by 0")
    }
}

impl<'a> Div<BasedExpr> for &'a BasedExpr {
    type Output = BasedExpr;

    fn div(self, rhs: BasedExpr) -> Self::Output {
        self.checked_div(rhs).expect("division by 0")
    }
}

impl<'a> Div<&BasedExpr> for &'a BasedExpr {
    type Output = BasedExpr;

    fn div(self, rhs: &BasedExpr) -> Self::Output {
        self.checked_div(rhs).expect("division by 0")
    }
}

impl CheckedDiv<BasedExpr> for BasedExpr {
    type Output = BasedExpr;

    fn checked_div(mut self, rhs: BasedExpr) -> Option<Self::Output> {
        if let (None, Some(basis)) = (self.basis(), rhs.basis()) {
            self = BasedExpr::based_rational(if let BasedExpr::Baseless(q) = self {q} else {unreachable!()}, basis.clone());
        }
        match (self, rhs) {
            (BasedExpr::Baseless(a), BasedExpr::Baseless(b)) => a.checked_div(b).map(BasedExpr::Baseless),
            (BasedExpr::Based(a, basis), BasedExpr::Baseless(b)) =>
                (!b.is_zero()).then(|| BasedExpr::Based(a / b, basis)),
            (BasedExpr::Baseless(_), BasedExpr::Based(_, _)) => unreachable!(),
            (BasedExpr::Based(a, basis_a), BasedExpr::Based(b, basis_b)) => {
                // The dreaded checked_division
                // TODO: Make this faster if necessary.
                basis_a.assert_compatible(&basis_b);
                let mut mtx = DMatrix::repeat(a.len(), a.len(), Rat::ZERO);
                for (row_i, row) in basis_a.matrices.iter().enumerate() {
                    for (coeff, a_i, b_i) in row {
                        mtx[(row_i, *b_i)] += coeff * &b[*a_i];
                    }
                }
                let result = mtx.try_inverse_mut();
                if !result { None? }
                Some(BasedExpr::Based(mtx * a, basis_a))
            }
        }
    }
}

impl CheckedDiv<&BasedExpr> for BasedExpr {
    type Output = BasedExpr;

    fn checked_div(mut self, rhs: &BasedExpr) -> Option<Self::Output> {
        if let (None, Some(basis)) = (self.basis(), rhs.basis()) {
            self = BasedExpr::based_rational(if let BasedExpr::Baseless(q) = self {q} else {unreachable!()}, basis.clone());
        }
        match (self, rhs) {
            (BasedExpr::Baseless(a), BasedExpr::Baseless(b)) => a.checked_div(b).map(BasedExpr::Baseless),
            (BasedExpr::Based(a, basis), BasedExpr::Baseless(b)) =>
                (!b.is_zero()).then(|| BasedExpr::Based(a / b.clone(), basis)),
            (BasedExpr::Baseless(_), BasedExpr::Based(_, _)) => unreachable!(),
            (BasedExpr::Based(a, basis_a), BasedExpr::Based(b, basis_b)) => {
                // The dreaded checked_division
                // TODO: Make this faster if necessary.
                basis_a.assert_compatible(&basis_b);
                let mut mtx = DMatrix::repeat(a.len(), a.len(), Rat::ZERO);
                for (row_i, row) in basis_a.matrices.iter().enumerate() {
                    for (coeff, a_i, b_i) in row {
                        mtx[(row_i, *b_i)] += coeff * &b[*a_i];
                    }
                }
                let result = mtx.try_inverse_mut();
                if !result { None? }
                Some(BasedExpr::Based(mtx * a, basis_a))
            }
        }
    }
}

impl<'a> CheckedDiv<BasedExpr> for &'a BasedExpr {
    type Output = BasedExpr;

    fn checked_div(self, mut rhs: BasedExpr) -> Option<Self::Output> {
        let denom = std::mem::replace(&mut rhs, self.clone());
        rhs.checked_div(denom)
    }
}

impl<'a> CheckedDiv<&BasedExpr> for &'a BasedExpr {
    type Output = BasedExpr;

    fn checked_div(self, rhs: &BasedExpr) -> Option<Self::Output> {
        if let (None, Some(basis)) = (self.basis(), rhs.basis()) {
            let numer = BasedExpr::based_rational(if let BasedExpr::Baseless(q) = self {q.clone()} else {unreachable!()}, basis.clone());
            return numer.checked_div(rhs)
        }
        self.clone().checked_div(rhs)
    }
}

impl Pow<i64> for BasedExpr {
    type Output = BasedExpr;

    fn pow(self, exp: i64) -> Self::Output {
        (&self).pow(exp)
    }
}

impl Pow<i64> for &BasedExpr {
    type Output = BasedExpr;

    fn pow(self, exp: i64) -> Self::Output {
        if exp < 0 {
            (self.to_one() / self).pow((-exp) as u64)
        } else {
            self.pow(exp as u64)
        }
    }
}

impl Pow<u64> for BasedExpr {
    type Output = BasedExpr;

    fn pow(self, exp: u64) -> Self::Output {
        (&self).pow(exp)
    }
}

impl Pow<u64> for &BasedExpr {
    type Output = BasedExpr;

    fn pow(self, exp: u64) -> Self::Output {
        let mut acc = self.to_one();
        for i in (0..exp.significant_bits()).rev() {
            if (exp >> i & 1) == 0 {
                acc = &acc * &acc;
            } else {
                acc = &acc * &acc * self;
            }
        }
        acc
    }
}

impl RemAssign<BasedExpr> for BasedExpr {
    /// Performs the `%=` operation.
    ///
    /// This *floors*, unlike with primitive integers.
    fn rem_assign(&mut self, rhs: Self) {
        let multiple = Floor::floor(&*self / &rhs) * rhs;
        *self -= multiple;
    }
}

impl RemAssign<&BasedExpr> for BasedExpr {
    /// Performs the `%=` operation.
    ///
    /// This *floors*, unlike with primitive integers.
    fn rem_assign(&mut self, rhs: &Self) {
        let multiple = Floor::floor(&*self / rhs) * rhs;
        *self -= multiple;
    }
}

impl Rem<BasedExpr> for BasedExpr {
    type Output = BasedExpr;

    /// Performs the `%` operation.
    ///
    /// This *floors*, unlike with primitive integers.
    fn rem(mut self, rhs: BasedExpr) -> Self::Output {
        self %= rhs;
        self
    }
}

impl Rem<&BasedExpr> for BasedExpr {
    type Output = BasedExpr;

    /// Performs the `%` operation.
    ///
    /// This *floors*, unlike with primitive integers.
    fn rem(mut self, rhs: &BasedExpr) -> Self::Output {
        self %= rhs;
        self
    }
}

impl<'a> Rem<BasedExpr> for &'a BasedExpr {
    type Output = BasedExpr;

    /// Performs the `%` operation.
    ///
    /// This *floors*, unlike with primitive integers.
    fn rem(self, rhs: BasedExpr) -> Self::Output {
        self - Floor::floor(self / &rhs) * rhs
    }
}

impl<'a> Rem<&BasedExpr> for &'a BasedExpr {
    type Output = BasedExpr;

    /// Performs the `%` operation.
    ///
    /// This *floors*, unlike with primitive integers.
    fn rem(self, rhs: &BasedExpr) -> Self::Output {
        self - Floor::floor(self / rhs) * rhs
    }
}

impl Abs for BasedExpr {
    type Output = BasedExpr;
    fn abs(self) -> Self::Output { if self.cmp_zero() == Ordering::Less { -self } else { self } }
}

impl Abs for &BasedExpr {
    type Output = BasedExpr;
    fn abs(self) -> Self::Output { if self.cmp_zero() == Ordering::Less { -self } else { self.clone() } }
}

impl AbsDiff<BasedExpr> for BasedExpr {
    type Output = BasedExpr;
    fn abs_diff(self, other: BasedExpr) -> Self::Output { Abs::abs(self - other) }
}

impl AbsDiff<&BasedExpr> for BasedExpr {
    type Output = BasedExpr;
    fn abs_diff(self, other: &BasedExpr) -> Self::Output { Abs::abs(self - other) }
}

impl AbsDiff<BasedExpr> for &BasedExpr {
    type Output = BasedExpr;
    fn abs_diff(self, other: BasedExpr) -> Self::Output { Abs::abs(self - other) }
}

impl AbsDiff<&BasedExpr> for &BasedExpr {
    type Output = BasedExpr;
    fn abs_diff(self, other: &BasedExpr) -> Self::Output { Abs::abs(self - other) }
}

impl Ceiling for BasedExpr {
    type Output = BasedExpr; // because we need to preserve the basis
    fn ceiling(self) -> Self::Output { (&self).round(RoundingMode::Ceiling) }
}

impl Ceiling for &BasedExpr {
    type Output = BasedExpr; // because we need to preserve the basis
    fn ceiling(self) -> Self::Output { self.round(RoundingMode::Ceiling) }
}

impl Floor for BasedExpr {
    type Output = BasedExpr; // because we need to preserve the basis
    fn floor(self) -> Self::Output { (&self).round(RoundingMode::Floor) }
}

impl Floor for &BasedExpr {
    type Output = BasedExpr; // because we need to preserve the basis
    fn floor(self) -> Self::Output { self.round(RoundingMode::Floor) }
}

// ComplexField traits (for better nalgebra support. Note that not all functions are valid.)

impl Signed for BasedExpr {
    fn abs(&self) -> Self {
        Abs::abs(self)
    }

    fn abs_sub(&self, other: &Self) -> Self {
        self.abs_diff(other)
    }

    fn signum(&self) -> Self {
        match self {
            BasedExpr::Baseless(a) => BasedExpr::Baseless(a.signum()),
            BasedExpr::Based(a, basis) => {
                let mut coeffs = DVector::repeat(a.len(), Rat::ZERO);
                coeffs[0] = match self.cmp_zero() {
                    Ordering::Less => Rat::MINUS_ONE,
                    Ordering::Equal => Rat::ZERO,
                    Ordering::Greater => Rat::ONE,
                };
                BasedExpr::Based(coeffs, basis.clone())
            }
        }
    }

    fn is_positive(&self) -> bool {
        self.cmp_zero() == Ordering::Greater
    }

    fn is_negative(&self) -> bool {
        self.cmp_zero() == Ordering::Less
    }
}

impl SimdValue for BasedExpr {
    const LANES: usize = 1;
    type Element = BasedExpr;
    type SimdBool = bool;

    fn splat(val: Self::Element) -> Self {
        val
    }

    fn extract(&self, _: usize) -> Self::Element {
        self.clone()
    }

    unsafe fn extract_unchecked(&self, _: usize) -> Self::Element {
        self.clone()
    }

    fn replace(&mut self, _: usize, val: Self::Element) {
        *self = val;
    }

    unsafe fn replace_unchecked(&mut self, _: usize, val: Self::Element) {
        *self = val;
    }

    fn select(self, cond: Self::SimdBool, other: Self) -> Self {
        if cond { self } else { other }
    }
}

impl Field for BasedExpr {}

impl AbsDiffEq for BasedExpr {
    type Epsilon = BasedExpr;

    /// Returns a baseless 0. These values are exact; no need for epsilon comparison.
    fn default_epsilon() -> Self::Epsilon {
        BasedExpr::BASELESS_ZERO
    }

    /// Equivalent to checking equality
    fn abs_diff_eq(&self, other: &Self, _: Self::Epsilon) -> bool {
        // Exact values; no need for fudgy comparisons
        self == other
    }
}

impl UlpsEq for BasedExpr {
    /// Returns 0. These values are exact; no need for epsilon comparison.
    fn default_max_ulps() -> u32 {
        0
    }

    /// Equivalent to checking equality
    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, _: u32) -> bool {
        // Exact values; no need for fudgy comparisons
        self.abs_diff_eq(other, epsilon)
    }
}

impl RelativeEq for BasedExpr {
    /// Returns a baseless 1. These values are exact; no need for epsilon comparison.
    fn default_max_relative() -> Self::Epsilon {
        Self::BASELESS_ONE
    }

    /// Equivalent to checking equality
    fn relative_eq(&self, other: &Self, _: Self::Epsilon, _: Self::Epsilon)
        -> bool {
        // Exact values; no need for fudgy comparisons
        self == other
    }
}

impl FromPrimitive for BasedExpr {
    /// Converts into a baseless number.
    fn from_i64(n: i64) -> Option<Self> {
        Some(n.into_baseless())
    }

    /// Converts into a baseless number.
    fn from_u64(n: u64) -> Option<Self> {
        Some(n.into_baseless())
    }
}

impl SubsetOf<BasedExpr> for BasedExpr {
    fn to_superset(&self) -> BasedExpr {
        self.clone()
    }

    fn from_superset_unchecked(element: &BasedExpr) -> Self {
        element.clone()
    }

    fn is_in_subset(_: &BasedExpr) -> bool {
        true
    }
}

impl SubsetOf<BasedExpr> for f64 {
    /// `f64` is technically not a subset of `BasedExpr`. This panics if
    /// `self` is not finite.
    fn to_superset(&self) -> BasedExpr {
        // A necessary evil to make this compatible with nalgebra. Some methods don't *really* need ComplexField.
        BasedExpr::try_baseless_from(*self).unwrap_or_else(|e| panic!("{e:?}: cannot convert {self} to BasedExpr"))
    }

    /// Uses nearest rounding, ties to even.
    fn from_superset_unchecked(element: &BasedExpr) -> Self {
        element.round_to_f64(RoundingMode::Nearest)
    }

    fn is_in_subset(element: &BasedExpr) -> bool {
        true
    }
}

impl SubsetOf<BasedExpr> for f32 {
    /// `f32` is technically not a subset of `BasedExpr`. This panics if
    /// `self` is not finite.
    fn to_superset(&self) -> BasedExpr {
        // A necessary evil to make this compatible with nalgebra. Some methods don't *really* need ComplexField.
        BasedExpr::try_baseless_from(*self).unwrap_or_else(|e| panic!("{e:?}: cannot convert {self} to BasedExpr"))
    }

    /// Uses nearest rounding, ties to even.
    fn from_superset_unchecked(element: &BasedExpr) -> Self {
        element.round_to_f64(RoundingMode::Nearest) as f32
    }

    fn is_in_subset(element: &BasedExpr) -> bool {
        true
    }
}

impl RealField for BasedExpr {
    fn is_sign_positive(&self) -> bool {
        self.is_positive()
    }

    fn is_sign_negative(&self) -> bool {
        self.is_negative()
    }

    fn copysign(self, sign: Self) -> Self {
        match sign.cmp_zero() {
            Ordering::Less => (-self).with_basis(sign.basis().cloned()),
            Ordering::Equal => sign.with_basis(self.basis().cloned()),
            Ordering::Greater => self.with_basis(sign.basis().cloned()),
        }
    }

    fn max(self, other: Self) -> Self {
        Ord::max(self, other)
    }

    fn min(self, other: Self) -> Self {
        Ord::min(self, other)
    }

    fn clamp(self, min: Self, max: Self) -> Self {
        Ord::clamp(self, min, max)
    }

    /// Fake function; do not call.
    fn atan2(self, other: Self) -> Self {
        panic!("cannot take atan2 of field extension elements {self}, {other}")
    }

    fn min_value() -> Option<Self> {
        None
    }

    fn max_value() -> Option<Self> {
        None
    }

    // And here comes another necessary evil
    // A bunch of irrational constants that can't be defined

    /// Fake function; do not call.
    fn pi() -> Self {
        panic!("π is not algebraic")
    }

    /// Fake function; do not call.
    fn two_pi() -> Self {
        panic!("2π is not algebraic")
    }

    /// Fake function; do not call.
    fn frac_pi_2() -> Self {
        panic!("π/2 is not algebraic")
    }

    /// Fake function; do not call.
    fn frac_pi_3() -> Self {
        panic!("π/3 is not algebraic")
    }

    /// Fake function; do not call.
    fn frac_pi_4() -> Self {
        panic!("π/4 is not algebraic")
    }

    /// Fake function; do not call.
    fn frac_pi_6() -> Self {
        panic!("π/6 is not algebraic")
    }

    /// Fake function; do not call.
    fn frac_pi_8() -> Self {
        panic!("π/8 is not algebraic")
    }

    /// Fake function; do not call.
    fn frac_1_pi() -> Self {
        panic!("1/π is not algebraic")
    }

    /// Fake function; do not call.
    fn frac_2_pi() -> Self {
        panic!("2/π is not algebraic")
    }

    /// Fake function; do not call.
    fn frac_2_sqrt_pi() -> Self {
        panic!("2/√π is not algebraic")
    }

    /// Fake function; do not call.
    fn e() -> Self {
        panic!("e is not algebraic")
    }

    /// Fake function; do not call.
    fn log2_e() -> Self {
        panic!("log_2(e) is not algebraic")
    }

    /// Fake function; do not call.
    fn log10_e() -> Self {
        panic!("log_10(e) is not algebraic")
    }

    /// Fake function; do not call.
    fn ln_2() -> Self {
        panic!("ln(2) is not algebraic")
    }

    /// Fake function; do not call.
    fn ln_10() -> Self {
        panic!("ln(10) is not algebraic")
    }
}

// For use in matrices
impl ComplexField for BasedExpr {
    type RealField = BasedExpr;

    #[doc = r" Builds a pure-real complex number from the given value."]
    fn from_real(re: Self::RealField) -> Self {
        re
    }

    #[doc = r" The real part of this complex number."]
    fn real(self) -> Self::RealField {
        self
    }

    #[doc = r" The imaginary part of this complex number."]
    fn imaginary(self) -> Self::RealField {
        self.into_zero()
    }

    #[doc = r" The modulus of this complex number."]
    fn modulus(self) -> Self::RealField {
        Abs::abs(self)
    }

    #[doc = r" The squared modulus of this complex number."]
    fn modulus_squared(self) -> Self::RealField {
        &self * &self
    }

    #[doc = r" The argument of this complex number."]
    fn argument(self) -> Self::RealField {
        self.into_zero()
    }

    #[doc = r" The sum of the absolute value of this complex number's real and imaginary part."]
    fn norm1(self) -> Self::RealField {
        Abs::abs(self)
    }

    #[doc = r" Multiplies this complex number by `factor`."]
    fn scale(self, factor: Self::RealField) -> Self {
        self * factor
    }

    #[doc = r" Divides this complex number by `factor`."]
    fn unscale(self, factor: Self::RealField) -> Self {
        self / factor
    }

    fn floor(self) -> Self {
        Floor::floor(self)
    }

    fn ceil(self) -> Self {
        self.ceiling()
    }

    fn round(self) -> Self {
        let basis = self.basis().cloned();
        BasedExpr::new_baseless(self.round_to_integer(RoundingMode::Nearest).into()).with_basis(basis)
    }

    fn trunc(self) -> Self {
        let basis = self.basis().cloned();
        BasedExpr::new_baseless(self.round_to_integer(RoundingMode::Down).into()).with_basis(basis)
    }

    fn fract(self) -> Self {
        self.clone() - self.trunc()
    }

    fn mul_add(self, a: Self, b: Self) -> Self {
        self * a + b
    }

    #[doc = r" The absolute value of this complex number: `self / self.signum()`."]
    #[doc = r""]
    #[doc = r" This is equivalent to `self.modulus()`."]
    fn abs(self) -> Self::RealField {
        Abs::abs(self)
    }

    #[doc = r" Computes (self.conjugate() * self + other.conjugate() * other).sqrt()"]
    /// Currently fake function; do not call.
    fn hypot(self, other: Self) -> Self::RealField {
        unimplemented!("cannot take hypothenuse of field extension yet. Some cases will stay undefined.")
    }

    fn recip(self) -> Self {
        Self::BASELESS_ONE / self
    }

    fn conjugate(self) -> Self {
        self
    }

    /// Fake function; do not call.
    fn sin(self) -> Self {
        panic!("cannot take sin of field extension {self}")
    }

    /// Fake function; do not call.
    fn cos(self) -> Self {
        panic!("cannot take cos of field extension {self}")
    }

    /// Fake function; do not call.
    fn sin_cos(self) -> (Self,Self) {
        panic!("cannot take sin_cos of field extension {self}")
    }

    /// Fake function; do not call.
    fn tan(self) -> Self {
        panic!("cannot take tan of field extension {self}")
    }

    /// Fake function; do not call.
    fn asin(self) -> Self {
        panic!("cannot take asin of field extension {self}")
    }

    /// Fake function; do not call.
    fn acos(self) -> Self {
        panic!("cannot take acos of field extension {self}")
    }

    /// Fake function; do not call.
    fn atan(self) -> Self {
        panic!("cannot take atan of field extension {self}")
    }

    /// Fake function; do not call.
    fn sinh(self) -> Self {
        panic!("cannot take sinh of field extension {self}")
    }

    /// Fake function; do not call.
    fn cosh(self) -> Self {
        panic!("cannot take cosh of field extension {self}")
    }

    /// Fake function; do not call.
    fn tanh(self) -> Self {
        panic!("cannot take tanh of field extension {self}")
    }

    /// Fake function; do not call.
    fn asinh(self) -> Self {
        panic!("cannot take asinh of field extension {self}")
    }

    /// Fake function; do not call.
    fn acosh(self) -> Self {
        panic!("cannot take acosh of field extension {self}")
    }

    /// Fake function; do not call.
    fn atanh(self) -> Self {
        panic!("cannot take atanh of field extension {self}")
    }

    /// Currently fake function; do not call.
    fn log(self, base: Self::RealField) -> Self {
        unimplemented!("cannot take log_{base} of field extension yet. Some cases will stay undefined.")
    }

    /// Currently fake function; do not call.
    fn log2(self) -> Self {
        unimplemented!("cannot take log_2 of field extension yet. Some cases will stay undefined.")
    }

    /// Currently fake function; do not call.
    fn log10(self) -> Self {
        unimplemented!("cannot take log_10 of field extension yet. Some cases will stay undefined.")
    }

    /// Fake function; do not call.
    fn ln(self) -> Self {
        panic!("cannot take ln_e of field extension {self}")
    }

    /// Fake function; do not call.
    fn ln_1p(self) -> Self {
        panic!("cannot take ln_1p of field extension {self}")
    }

    /// Currently fake function; do not call.
    fn sqrt(self) -> Self {
        panic!("cannot take sqrt of field extension yet. Some cases will stay undefined.")
    }

    /// Currently fake function; do not call.
    fn exp(self) -> Self {
        panic!("cannot take exp of field extension {self}")
    }

    /// Currently fake function; do not call.
    fn exp2(self) -> Self {
        panic!("cannot take exp2 of this field extension yet. Some cases will stay undefined.")
    }

    /// Fake function; do not call.
    fn exp_m1(self) -> Self {
        panic!("cannot take exp_m1 of field extension {self}")
    }

    fn powi(self, n: i32) -> Self {
        self.pow(n as i64)
    }

    /// Currently fake function; do not call.
    fn powf(self, n: Self::RealField) -> Self {
        panic!("cannot raise field extension to the {n}-th power yet. Some cases will stay undefined")
    }

    fn powc(self, n: Self) -> Self {
        self.powf(n)
    }

    /// Currently fake function; do not call.
    fn cbrt(self) -> Self {
        panic!("cannot take cbrt of field extension yet. Some cases will stay undefined")
    }

    fn is_finite(&self) -> bool {
        true
    }

    /// Currently fake function; always returns `None`
    fn try_sqrt(self) -> Option<Self> {
        // TODO: Define if ready
        None
    }
}

#[cfg(test)]
mod test {
    use std::cmp::Ordering;

    use crate::{based_expr, rat::Rat, BasedExpr};
    use malachite::base::num::arithmetic::traits::{Abs, AbsDiff, Ceiling, CheckedDiv, Floor, Pow};

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

    fn checked_div_test(a: BasedExpr, b: BasedExpr, expected: Option<BasedExpr>) {
        assert_eq!(a.clone().checked_div(b.clone()), expected);
        assert_eq!(a.clone().checked_div(&b), expected);
        assert_eq!((&a).checked_div(b.clone()), expected);
        assert_eq!((&a).checked_div(&b), expected);
    }

    fn rchecked_div_test(a: BasedExpr, b: BasedExpr, expected: Option<BasedExpr>) {
        assert_eq!(b.clone().checked_div(a.clone()), expected);
        assert_eq!(b.clone().checked_div(&a), expected);
        assert_eq!((&b).checked_div(a.clone()), expected);
        assert_eq!((&b).checked_div(&a), expected);
    }

    fn pow_test(a: BasedExpr, exp: i64, expected: BasedExpr) {
        assert_eq!(a.clone().pow(exp), expected);
        assert_eq!((&a).pow(exp), expected);
    }

    fn abs_test(a: BasedExpr, expected: BasedExpr) {
        assert_eq!(Abs::abs(a.clone()), expected);
        assert_eq!(Abs::abs(&a), expected);
    }

    fn abs_diff_test(a: BasedExpr, b: BasedExpr, expected: BasedExpr) {
        assert_eq!(a.clone().abs_diff(b.clone()), expected);
        assert_eq!(a.clone().abs_diff(&b), expected);
        assert_eq!((&a).abs_diff(b.clone()), expected);
        assert_eq!((&a).abs_diff(&b), expected);
    }

    fn rabs_diff_test(a: BasedExpr, b: BasedExpr, expected: BasedExpr) {
        assert_eq!(b.clone().abs_diff(a.clone()), expected);
        assert_eq!(b.clone().abs_diff(&a), expected);
        assert_eq!((&b).abs_diff(a.clone()), expected);
        assert_eq!((&b).abs_diff(&a), expected);
    }

    fn cmp_test(a: BasedExpr, b: BasedExpr, expected: Ordering) {
        assert_eq!(a.cmp(&b), expected);
    }

    fn rcmp_test(a: BasedExpr, b: BasedExpr, expected: Ordering) {
        assert_eq!(b.cmp(&a), expected);
    }

    fn floor_test(a: BasedExpr, expected: BasedExpr) {
        assert_eq!(a.floor(), expected);
    }

    fn ceiling_test(a: BasedExpr, expected: BasedExpr) {
        assert_eq!(a.ceiling(), expected);
    }

    // Assuming a perfectly working Rational library and just testing BasedExpr-specific behavior
    #[test]
    fn test_based_expr_add() {
        let _lock = crate::TEST_LOCK.read();
        add_test(BasedExpr::new_baseless((-7i64).into()), BasedExpr::new_baseless((-5i64).into()), BasedExpr::new_baseless((-12i64).into()));
        add_test(BasedExpr::new_baseless(2.into()), BasedExpr::new_baseless(2.into()), BasedExpr::new_baseless(4.into()));

        add_test(BasedExpr::new_baseless((-7i64).into()), based_expr!(-5), based_expr!(-12));
        add_test(BasedExpr::new_baseless((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(-4 + sqrt 2));
        add_test(BasedExpr::new_baseless(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(121 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6));

        radd_test(BasedExpr::new_baseless((-7i64).into()), based_expr!(-5), based_expr!(-12));
        radd_test(BasedExpr::new_baseless((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(-4 + sqrt 2));
        radd_test(BasedExpr::new_baseless(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(121 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6));

        add_test(based_expr!(4), based_expr!(8), based_expr!(12));
        add_test(based_expr!(1 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(4 + 3 sqrt 2));
        add_test(based_expr!(1/2 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(7/2 + 3 sqrt 2));
        add_test(based_expr!(1 + sqrt 2 - 4 sqrt 5 + 8 sqrt 10), based_expr!(-1 + 100 sqrt 2 + 50 sqrt 5 + 7 sqrt 10),
            based_expr!(0 + 101 sqrt 2 + 46 sqrt 5 + 15 sqrt 10));
    }

    #[test]
    fn test_based_expr_sub() {
        let _lock = crate::TEST_LOCK.read();
        sub_test(BasedExpr::new_baseless((-7i64).into()), BasedExpr::new_baseless((-5i64).into()), BasedExpr::new_baseless((-2i64).into()));
        sub_test(BasedExpr::new_baseless(2.into()), BasedExpr::new_baseless(2.into()), BasedExpr::new_baseless(0.into()));

        sub_test(BasedExpr::new_baseless((-7i64).into()), based_expr!(-5), based_expr!(-2));
        sub_test(BasedExpr::new_baseless((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(-6 - sqrt 2));
        sub_test(BasedExpr::new_baseless(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(119 - 2 sqrt 2 - 3 sqrt 3 - 6 sqrt 6));

        rsub_test(BasedExpr::new_baseless((-7i64).into()), based_expr!(-5), based_expr!(2));
        rsub_test(BasedExpr::new_baseless((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(6 + sqrt 2));
        rsub_test(BasedExpr::new_baseless(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(-119 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6));

        sub_test(based_expr!(4), based_expr!(8), based_expr!(-4));
        sub_test(based_expr!(1 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(-2 - sqrt 2));
        sub_test(based_expr!(1/2 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(-5/2 - sqrt 2));
        sub_test(based_expr!(1 + sqrt 2 - 4 sqrt 5 + 8 sqrt 10), based_expr!(-1 + 100 sqrt 2 + 50 sqrt 5 + 7 sqrt 10),
            based_expr!(2 - 99 sqrt 2 - 54 sqrt 5 + sqrt 10));
    }

    #[test]
    fn test_based_expr_abs_diff() {
        // Might as well borrow the subtraction test case list since it's there.
        let _lock = crate::TEST_LOCK.read();
        abs_diff_test(BasedExpr::new_baseless((-7i64).into()), BasedExpr::new_baseless((-5i64).into()), BasedExpr::new_baseless(2.into()));
        abs_diff_test(BasedExpr::new_baseless(2.into()), BasedExpr::new_baseless(2.into()), BasedExpr::new_baseless(0.into()));

        abs_diff_test(BasedExpr::new_baseless((-7i64).into()), based_expr!(-5), based_expr!(2));
        abs_diff_test(BasedExpr::new_baseless((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(6 + sqrt 2));
        abs_diff_test(BasedExpr::new_baseless(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(119 - 2 sqrt 2 - 3 sqrt 3 - 6 sqrt 6));

        rabs_diff_test(BasedExpr::new_baseless((-7i64).into()), based_expr!(-5), based_expr!(2));
        rabs_diff_test(BasedExpr::new_baseless((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(6 + sqrt 2));
        rabs_diff_test(BasedExpr::new_baseless(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(119 - 2 sqrt 2 - 3 sqrt 3 - 6 sqrt 6));

        abs_diff_test(based_expr!(4), based_expr!(8), based_expr!(4));
        abs_diff_test(based_expr!(1 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(2 + sqrt 2));
        abs_diff_test(based_expr!(1/2 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(5/2 + sqrt 2));
        abs_diff_test(based_expr!(1 + sqrt 2 - 4 sqrt 5 + 8 sqrt 10), based_expr!(-1 + 100 sqrt 2 + 50 sqrt 5 + 7 sqrt 10),
            based_expr!(-2 + 99 sqrt 2 + 54 sqrt 5 - sqrt 10));
    }

    #[test]
    fn test_based_expr_mul() {
        let _lock = crate::TEST_LOCK.read();
        mul_test(BasedExpr::new_baseless((-7i64).into()), BasedExpr::new_baseless((-5i64).into()), BasedExpr::new_baseless(35.into()));
        mul_test(BasedExpr::new_baseless(2.into()), BasedExpr::new_baseless(2.into()), BasedExpr::new_baseless(4.into()));

        mul_test(BasedExpr::new_baseless((-7i64).into()), based_expr!(-5), based_expr!(35));
        mul_test(BasedExpr::new_baseless((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(-5 - 5 sqrt 2));
        mul_test(BasedExpr::new_baseless(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(120 + 240 sqrt 2 + 360 sqrt 3 + 720 sqrt 6));

        rmul_test(BasedExpr::new_baseless((-7i64).into()), based_expr!(-5), based_expr!(35));
        rmul_test(BasedExpr::new_baseless((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(-5 - 5 sqrt 2));
        rmul_test(BasedExpr::new_baseless(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(120 + 240 sqrt 2 + 360 sqrt 3 + 720 sqrt 6));

        mul_test(based_expr!(4), based_expr!(8), based_expr!(32));
        mul_test(based_expr!(1 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(7 + 5 sqrt 2));
        mul_test(based_expr!(1/2 + sqrt 2), based_expr!(3 + 2 sqrt 2), based_expr!(11/2 + 4 sqrt 2));
        mul_test(based_expr!(1 + sqrt 2 - 4 sqrt 5 + 8 sqrt 10), based_expr!(-1 + 100 sqrt 2 + 50 sqrt 5 + 7 sqrt 10),
            based_expr!(-241 + 1959 sqrt 2 + 1668 sqrt 5 - 351 sqrt 10));
    }

    #[test]
    fn test_based_expr_checked_div() {
        let _lock = crate::TEST_LOCK.read();
        checked_div_test(BasedExpr::new_baseless((-7i64).into()), BasedExpr::new_baseless((-5i64).into()),
            Some(BasedExpr::new_baseless(Rat::from_signeds(7, 5))));
        checked_div_test(BasedExpr::new_baseless(2.into()), BasedExpr::new_baseless(2.into()), Some(BasedExpr::new_baseless(1.into())));
        checked_div_test(BasedExpr::new_baseless(3.into()), BasedExpr::new_baseless(0.into()), None);

        checked_div_test(BasedExpr::new_baseless((-7i64).into()), based_expr!(-5), Some(based_expr!(7/5)));
        checked_div_test(BasedExpr::new_baseless((-5i64).into()), based_expr!(1 + sqrt 2), Some(based_expr!(5 - 5 sqrt 2)));
        checked_div_test(BasedExpr::new_baseless(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            Some(based_expr!(60/91 - 120/91 sqrt 2 - 180/91 sqrt 3 + 360/91 sqrt 6)));
        checked_div_test(BasedExpr::new_baseless((-5i64).into()), based_expr!(0 + 0 sqrt 2), None);

        rchecked_div_test(BasedExpr::new_baseless((-7i64).into()), based_expr!(-5), Some(based_expr!(5/7)));
        rchecked_div_test(BasedExpr::new_baseless((-5i64).into()), based_expr!(1 + sqrt 2), Some(based_expr!(-1/5 - 1/5 sqrt 2)));
        rchecked_div_test(BasedExpr::new_baseless(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            Some(based_expr!(1/120 + 1/60 sqrt 2 + 1/40 sqrt 3 + 1/20 sqrt 6)));
        rchecked_div_test(BasedExpr::new_baseless(0.into()), based_expr!(3 + 1/2 sqrt 2), None);

        checked_div_test(based_expr!(4), based_expr!(8), Some(based_expr!(1/2)));
        checked_div_test(based_expr!(1 + sqrt 2), based_expr!(3 + 2 sqrt 2), Some(based_expr!(-1 + sqrt 2)));
        checked_div_test(based_expr!(1/2 + sqrt 2), based_expr!(3 + 2 sqrt 2), Some(based_expr!(-5/2 + 2 sqrt 2)));
        checked_div_test(based_expr!(1 + sqrt 2 - 4 sqrt 5 + 8 sqrt 10), based_expr!(-1 + 100 sqrt 2 + 50 sqrt 5 + 7 sqrt 10),
        // Kaboom!
            Some(based_expr!(-8551371/21774121 - 9982071/21774121 sqrt 2 + 7355940/21774121 sqrt 5 + 2437885/21774121 sqrt 10)));
        checked_div_test(based_expr!(1 + sqrt 2 - 4 sqrt 5 + 8 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), None);
    }

    #[test]
    fn test_based_expr_div() {
        let _lock = crate::TEST_LOCK.read();
        div_test(BasedExpr::new_baseless((-7i64).into()), BasedExpr::new_baseless((-5i64).into()),
            BasedExpr::new_baseless(Rat::from_signeds(7, 5)));
        div_test(BasedExpr::new_baseless(2.into()), BasedExpr::new_baseless(2.into()), BasedExpr::new_baseless(1.into()));

        div_test(BasedExpr::new_baseless((-7i64).into()), based_expr!(-5), based_expr!(7/5));
        div_test(BasedExpr::new_baseless((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(5 - 5 sqrt 2));
        div_test(BasedExpr::new_baseless(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
            based_expr!(60/91 - 120/91 sqrt 2 - 180/91 sqrt 3 + 360/91 sqrt 6));

        rdiv_test(BasedExpr::new_baseless((-7i64).into()), based_expr!(-5), based_expr!(5/7));
        rdiv_test(BasedExpr::new_baseless((-5i64).into()), based_expr!(1 + sqrt 2), based_expr!(-1/5 - 1/5 sqrt 2));
        rdiv_test(BasedExpr::new_baseless(120.into()), based_expr!(1 + 2 sqrt 2 + 3 sqrt 3 + 6 sqrt 6),
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
        let _lock = crate::TEST_LOCK.read();
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
        let _lock = crate::TEST_LOCK.read();
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
        let _lock = crate::TEST_LOCK.read();

        cmp_test(BasedExpr::new_baseless((-7i64).into()), BasedExpr::new_baseless((-5i64).into()), Ordering::Less);
        cmp_test(BasedExpr::new_baseless(4.into()), BasedExpr::new_baseless(0.into()), Ordering::Greater);
        cmp_test(BasedExpr::new_baseless(1.into()), BasedExpr::new_baseless(1.into()), Ordering::Equal);

        cmp_test(BasedExpr::new_baseless(7.into()), based_expr!(0 + 5 sqrt 2), Ordering::Less);
        cmp_test(BasedExpr::new_baseless(6.into()), based_expr!(6), Ordering::Equal);
        cmp_test(BasedExpr::new_baseless(6.into()), based_expr!(6 + 0 sqrt 2), Ordering::Equal);

        rcmp_test(BasedExpr::new_baseless(7.into()), based_expr!(0 + 5 sqrt 2), Ordering::Greater);
        rcmp_test(BasedExpr::new_baseless(6.into()), based_expr!(6), Ordering::Equal);
        rcmp_test(BasedExpr::new_baseless(6.into()), based_expr!(6 + 0 sqrt 2), Ordering::Equal);

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
        cmp_test(based_expr!( 16 + 11 sqrt 2 + 7 sqrt 5 + 5 sqrt 10), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Greater);
    }

    #[test]
    fn test_based_expr_cmp_hard_to_separate() {
        let _lock = crate::TEST_LOCK.read();
        // Requires level 0 (64 bits) precision
        let mut expr = based_expr!(1 + 0 sqrt 2 + 0 sqrt 5 - 364585791794594742/1152921504606846976 sqrt 10);
        cmp_test(expr.clone(), based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Greater);
        cmp_test(-expr, based_expr!(0 + 0 sqrt 2 + 0 sqrt 5 + 0 sqrt 10), Ordering::Less);
    }

    #[test]
    fn test_based_expr_cmp_huge() {
        let _lock = crate::TEST_LOCK.read();
        let mut expr = based_expr!(18446744073709551615 + 1 sqrt 2 + 1 sqrt 3 + 1 sqrt 6);
        expr = &expr * &expr; // 128
        expr = &expr * &expr; // 256
        expr = &expr * &expr; // 512
        expr = &expr * &expr; // 1024
        expr = &expr * &expr; // 2048. Now it's too big for a f64.
        cmp_test(expr.clone(), based_expr!(0 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6), Ordering::Greater);
        cmp_test(-expr, based_expr!(0 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6), Ordering::Less);
    }

    #[test]
    fn test_based_expr_cmp_tiny_rational() {
        let _lock = crate::TEST_LOCK.read();
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

    #[test]
    fn test_based_expr_floor() {
        floor_test(BasedExpr::new_baseless(4.into()), BasedExpr::new_baseless(4.into()));
        floor_test(BasedExpr::new_baseless(Rat::from_signeds(-7, 2)), BasedExpr::new_baseless((-4i64).into()));
        floor_test(BasedExpr::new_baseless(Rat::from_signeds(7, 2)), BasedExpr::new_baseless(3.into()));

        floor_test(based_expr!(1), based_expr!(1));
        floor_test(based_expr!(-16/3), based_expr!(-6));
        floor_test(based_expr!(16/3), based_expr!(5));

        floor_test(based_expr!(1 + 0 sqrt 19), based_expr!(1 + 0 sqrt 19));
        floor_test(based_expr!(-16/3 + 0 sqrt 19), based_expr!(-6 + 0 sqrt 19));
        floor_test(based_expr!(16/3 + 0 sqrt 19), based_expr!(5 + 0 sqrt 19));

        floor_test(based_expr!( 10 + 5 sqrt 2), based_expr!(17 + 0 sqrt 2));
        floor_test(based_expr!( 10 - 5 sqrt 2), based_expr!(2 + 0 sqrt 2));
        floor_test(based_expr!(-10 + 5 sqrt 2), based_expr!(-3 + 0 sqrt 2));
        floor_test(based_expr!(-10 - 5 sqrt 2), based_expr!(-18 + 0 sqrt 2));

        // Interesting: calculating sqrt 2 to arbitrary precision
        floor_test(based_expr!(0 - 10_000_000_000_000_000_000 sqrt 2), based_expr!(-14_142_135_623_730_950_489 + 0 sqrt 2));
        floor_test(based_expr!(0 + 10_000_000_000_000_000_000 sqrt 2), based_expr!(14_142_135_623_730_950_488 + 0 sqrt 2));
    }

    //#[test]
    //fn test_based_expr_floor_huge() {
    //    let _lock = LOCK.read();
    //    let mut expr = based_expr!(18446744073709551615 / 7 + 1 sqrt 2 + 2 sqrt 3 + 3 sqrt 6);
    //    expr = &expr * &expr; // 128
    //    expr = &expr * &expr; // 256
    //    expr = &expr * &expr; // 512
    //    expr = &expr * &expr; // 1024
    //    expr = &expr * &expr; // 2048. Now it's too big for a f64.
    //    floor_test(expr.clone(), based_expr!(29261315000490985866034265910655352362012703219793199536888319019821371131933806760325547444646635954444043342139809928318139696798544675925888188338221339885089369629447696498655821861600601720768456739971717834746537757852568694229208325968271188918862240764563977541245923183255772770994335182783063858980663169108890025394237468362301773429446826331614785772463900774182847877589130333071186724787988209609668792681925755067232598618056520190463661900957315389139194309000770415381977960484373066441133179907554455493079694231835246112771216149123492919052440327771892024748423837590563 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6));
    //}
    // Unfortunately, it's also too big for a literal. Let's just hope.

    #[test]
    fn test_based_expr_floor_extremely_hard_to_separate_rational() {
        let _lock = crate::TEST_LOCK.read();
        // A devilish case: requires quite a lot of precision with fixed point
        let mut expr = based_expr!(1 / 18446744073709551615 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6);
        expr = &expr * &expr; // 128
        expr = &expr * &expr; // 256
        expr = &expr * &expr; // 512
        expr = &expr * &expr; // 1024
        expr = &expr * &expr; // 2048. Now it's too small for a f64.
        expr += based_expr!(25 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6);
        floor_test(expr.clone(), based_expr!(25 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6));
        floor_test(-expr, based_expr!(-26 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6));
    }

    #[test]
    fn test_based_expr_ceiling() {
        ceiling_test(BasedExpr::new_baseless(4.into()), BasedExpr::new_baseless(4.into()));
        ceiling_test(BasedExpr::new_baseless(Rat::from_signeds(-7, 2)), BasedExpr::new_baseless((-3i64).into()));
        ceiling_test(BasedExpr::new_baseless(Rat::from_signeds(7, 2)), BasedExpr::new_baseless(4.into()));

        ceiling_test(based_expr!(1), based_expr!(1));
        ceiling_test(based_expr!(-16/3), based_expr!(-5));
        ceiling_test(based_expr!(16/3), based_expr!(6));

        ceiling_test(based_expr!(1 + 0 sqrt 19), based_expr!(1 + 0 sqrt 19));
        ceiling_test(based_expr!(-16/3 + 0 sqrt 19), based_expr!(-5 + 0 sqrt 19));
        ceiling_test(based_expr!(16/3 + 0 sqrt 19), based_expr!(6 + 0 sqrt 19));

        ceiling_test(based_expr!( 10 + 5 sqrt 2), based_expr!(18 + 0 sqrt 2));
        ceiling_test(based_expr!( 10 - 5 sqrt 2), based_expr!(3 + 0 sqrt 2));
        ceiling_test(based_expr!(-10 + 5 sqrt 2), based_expr!(-2 + 0 sqrt 2));
        ceiling_test(based_expr!(-10 - 5 sqrt 2), based_expr!(-17 + 0 sqrt 2));

        // Interesting: calculating sqrt 2 to arbitrary precision
        ceiling_test(based_expr!(0 - 10_000_000_000_000_000_000 sqrt 2), based_expr!(-14_142_135_623_730_950_488 + 0 sqrt 2));
        ceiling_test(based_expr!(0 + 10_000_000_000_000_000_000 sqrt 2), based_expr!(14_142_135_623_730_950_489 + 0 sqrt 2));
    }

    //#[test]
    //fn test_based_expr_ceiling_huge() {
    //    let _lock = LOCK.read();
    //    let mut expr = based_expr!(18446744073709551615 / 7 + 1 sqrt 2 + 2 sqrt 3 + 3 sqrt 6);
    //    expr = &expr * &expr; // 128
    //    expr = &expr * &expr; // 256
    //    expr = &expr * &expr; // 512
    //    expr = &expr * &expr; // 1024
    //    expr = &expr * &expr; // 2048. Now it's too big for a f64.
    //    ceiling_test(expr.clone(), based_expr!(29261315000490985866034265910655352362012703219793199536888319019821371131933806760325547444646635954444043342139809928318139696798544675925888188338221339885089369629447696498655821861600601720768456739971717834746537757852568694229208325968271188918862240764563977541245923183255772770994335182783063858980663169108890025394237468362301773429446826331614785772463900774182847877589130333071186724787988209609668792681925755067232598618056520190463661900957315389139194309000770415381977960484373066441133179907554455493079694231835246112771216149123492919052440327771892024748423837590564 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6));
    //}
    // Unfortunately, it's also too big for a literal. Let's just hope.

    #[test]
    fn test_based_expr_ceiling_extremely_hard_to_separate_rational() {
        let _lock = crate::TEST_LOCK.read();
        // A devilish case: requires quite a lot of precision with fixed point
        let mut expr = based_expr!(1 / 18446744073709551615 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6);
        expr = &expr * &expr; // 128
        expr = &expr * &expr; // 256
        expr = &expr * &expr; // 512
        expr = &expr * &expr; // 1024
        expr = &expr * &expr; // 2048. Now it's too small for a f64.
        expr += based_expr!(25 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6);
        ceiling_test(expr.clone(), based_expr!(26 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6));
        ceiling_test(-expr, based_expr!(-25 + 0 sqrt 2 + 0 sqrt 3 + 0 sqrt 6));
    }

    #[test]
    fn test_based_expr_abs() {
        // Based on cmp, so just check each case
        abs_test(based_expr!(1 - sqrt 2), based_expr!(-1 + sqrt 2));
        abs_test(based_expr!(-1 + sqrt 2), based_expr!(-1 + sqrt 2));
        abs_test(based_expr!(0 + 0 sqrt 3), based_expr!(0 + 0 sqrt 3));
    }

    #[test]
    fn test_based_expr_pow() {
        pow_test(BasedExpr::new_baseless(0.into()), 0, BasedExpr::new_baseless(1.into()));
        pow_test(BasedExpr::new_baseless(0.into()), 1, BasedExpr::new_baseless(0.into()));
        pow_test(BasedExpr::new_baseless(5.into()), 0, BasedExpr::new_baseless(1.into()));
        pow_test(BasedExpr::new_baseless(5.into()), 1, BasedExpr::new_baseless(5.into()));
        pow_test(BasedExpr::new_baseless(5.into()), 2, BasedExpr::new_baseless(25.into()));
        pow_test(BasedExpr::new_baseless(5.into()), -6, BasedExpr::new_baseless(Rat::from_signeds(1, 15625)));
        pow_test(BasedExpr::new_baseless((-5i64).into()), 1, BasedExpr::new_baseless((-5i64).into()));
        pow_test(BasedExpr::new_baseless((-5i64).into()), 2, BasedExpr::new_baseless(25.into()));

        pow_test(based_expr!(0), 0, based_expr!(1));
        pow_test(based_expr!(0), 1, based_expr!(0));
        pow_test(based_expr!(5), 0, based_expr!(1));
        pow_test(based_expr!(5), 1, based_expr!(5));
        pow_test(based_expr!(5), 2, based_expr!(25));
        pow_test(based_expr!(5), -6, based_expr!(1/15625));
        pow_test(based_expr!(-5), 1, based_expr!(-5));
        pow_test(based_expr!(-5), 2, based_expr!(25));

        pow_test(based_expr!(1 + sqrt 2), 2, based_expr!(3 + 2 sqrt 2));
    }
}