// I would use inari, but because of lack of support for rounding modes, it doesn't work on the web.

use std::{cmp::Ordering, fmt::{Debug, Display}, iter::Sum, ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign}};

use malachite::{base::{num::{arithmetic::traits::{Abs, AbsDiff, CeilingSqrt, DivRound, FloorSqrt, ShrRound, Sign}, basic::traits::Zero, conversion::traits::RoundingFrom}, rounding_modes::RoundingMode}, Integer, Natural};
use num::{Signed, Zero as _};

use crate::{rat::Rat, SqrtExpr};

/// Closed conservative interval between two real numbers.
/// 
/// !!! DO NOT DIVIDE THIS TYPE BY ZERO !!!
/// !!! FIXED-POINT TYPES DON'T HAVE AN INFINITY !!!
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Interval<T: IntervalContent> {
    bounds: Option<[T; 2]>,
    // Whether the interval possibly contains an invalid value.
    // For example, sqrt([-1, 1]) gives [0, 1] with a possible invalid value.
    maybe_invalid: bool,
}

impl<T: IntervalContent> Interval<T> {
    /// The empty interval
    pub const EMPTY: Interval<T> = Interval { bounds: None, maybe_invalid: false };

    pub fn new(bounds: [T; 2]) -> Self {
        Interval {
            bounds: Some(bounds),
            maybe_invalid: false
        }
    }

    fn with_invalid(mut self) -> Self {
        self.maybe_invalid = true;
        self
    }

    fn empty_if_overflowed(bounds: [T; 2]) -> Option<[T; 2]> {
        if T::overflowed(&bounds) { None } else { Some(bounds) }
    }

    fn contains_zero(&self) -> bool {
        match &self.bounds {
            None => false,
            Some(bounds) => bounds[0].true_sign() <= 0 && bounds[1].true_sign() >= 0
        }
    }

    /// Compares this interval to zero.
    /// Returns an Ordering if there's a definite answer and None otherwise.
    pub fn cmp_zero(&self) -> Option<Ordering> {
        match &self.bounds {
            None => None,
            Some(bounds) => if bounds[0].true_sign() > 0 {
                Some(Ordering::Greater)
            } else if bounds[1].true_sign() < 0 {
                Some(Ordering::Less)
            } else if bounds[0].true_sign() == 0 && bounds[1].true_sign() == 0 {
                Some(Ordering::Equal)
            } else {
                None
            }
        }
    }

    // Contains a value that is less than 0
    fn maybe_less_than_zero(&self) -> bool {
        match &self.bounds {
            None => false,
            Some(bounds) => bounds[0].true_sign() < 0
        }
    }

    pub fn sqrt(self) -> Interval<T> {
        let invalid = self.maybe_less_than_zero();
        Self {
            bounds: self.bounds.and_then(|[a, b]| {
                if b.true_sign() < 0 { None }
                else if a.true_sign() < 0 { Some([a.zero(), b.sqrt_up_v()]) }
                else { Some([a.sqrt_down_v(), b.sqrt_up_v()]) }
            }),
            maybe_invalid: self.maybe_invalid || invalid
        }
    }

    pub fn from_natural(nat: Natural, context: &T) -> Self {
        Interval::new([T::from_integer_down(nat.clone().into(), context), T::from_integer_up(nat.into(), context)])
    }

    pub fn from_integer(int: Integer, context: &T) -> Self {
        Interval::new([T::from_integer_down(int.clone(), context), T::from_integer_up(int, context)])
    }

    pub fn from_rational(rat: Rat, context: &T) -> Self {
        let negative = rat.is_negative();
        let (numer, denom) = rat.into_numerator_and_denominator();
        let result = Self::from_natural(numer, context) / Self::from_natural(denom, context);
        if negative { -result } else { result }
    }

    pub fn from_sqrt_expr(sqrt: SqrtExpr, context: &T) -> Self {
        match sqrt {
            SqrtExpr::Int(int) => Self::from_integer(int, context),
            SqrtExpr::Sum(terms) =>
                terms.into_iter().map(|(coeff, sqrt)| Self::from_integer(coeff, context) * Self::from_sqrt_expr(sqrt, context)).sum()
        }.sqrt()
    }
}

impl Interval<f64> {
    /// The entire f64 space (except infinities).
    pub const ENTIRE: Interval<f64> = Interval { bounds: Some([f64::NEG_INFINITY, f64::INFINITY]), maybe_invalid: false };
}

impl<T: IntervalContent> AddAssign<Interval<T>> for Interval<T> {
    fn add_assign(&mut self, rhs: Interval<T>) {
        self.bounds = std::mem::take(&mut self.bounds).zip(rhs.bounds).and_then(|([a0, a1], [b0, b1])|
            Self::empty_if_overflowed([a0.add_down_vv(b0), a1.add_up_vv(b1)]));
        self.maybe_invalid |= rhs.maybe_invalid;
    }
}

impl<T: IntervalContent> AddAssign<&Interval<T>> for Interval<T> {
    fn add_assign(&mut self, rhs: &Interval<T>) {
        self.bounds = std::mem::take(&mut self.bounds).zip(rhs.bounds.as_ref()).and_then(|([a0, a1], [b0, b1])|
            Self::empty_if_overflowed([a0.add_down_vr(b0), a1.add_up_vr(b1)]));
        self.maybe_invalid |= rhs.maybe_invalid;
    }
}

impl<T: IntervalContent> Add<Interval<T>> for Interval<T> {
    type Output = Interval<T>;

    fn add(mut self, rhs: Interval<T>) -> Self::Output {
        self += rhs;
        self
    }
}

impl<T: IntervalContent> Add<&Interval<T>> for Interval<T> {
    type Output = Interval<T>;

    fn add(mut self, rhs: &Interval<T>) -> Self::Output {
        self += rhs;
        self
    }
}

impl<T: IntervalContent> Add<Interval<T>> for &Interval<T> {
    type Output = Interval<T>;

    fn add(self, rhs: Interval<T>) -> Self::Output {
        rhs + self
    }
}

impl<T: IntervalContent> Add<&Interval<T>> for &Interval<T> {
    type Output = Interval<T>;

    fn add(self, rhs: &Interval<T>) -> Self::Output {
        Interval {
            bounds: self.bounds.as_ref().zip(rhs.bounds.as_ref()).and_then(|([a0, a1], [b0, b1])|
                Interval::empty_if_overflowed([a0.add_down_rr(b0), a1.add_up_rr(b1)])),
            maybe_invalid: self.maybe_invalid || rhs.maybe_invalid
        }
    }
}

impl<T: IntervalContent> SubAssign<Interval<T>> for Interval<T> {
    fn sub_assign(&mut self, rhs: Interval<T>) {
        self.bounds = std::mem::take(&mut self.bounds).zip(rhs.bounds).and_then(|([a0, a1], [b0, b1])|
            Self::empty_if_overflowed([a0.sub_down_vv(b1), a1.sub_up_vv(b0)]));
        self.maybe_invalid |= rhs.maybe_invalid;
    }
}

impl<T: IntervalContent> SubAssign<&Interval<T>> for Interval<T> {
    fn sub_assign(&mut self, rhs: &Interval<T>) {
        self.bounds = std::mem::take(&mut self.bounds).zip(rhs.bounds.as_ref()).and_then(|([a0, a1], [b0, b1])|
            Self::empty_if_overflowed([a0.sub_down_vr(b1), a1.sub_up_vr(b0)]));
        self.maybe_invalid |= rhs.maybe_invalid;
    }
}

impl<T: IntervalContent> Sub<Interval<T>> for Interval<T> {
    type Output = Interval<T>;

    fn sub(mut self, rhs: Interval<T>) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<T: IntervalContent> Sub<&Interval<T>> for Interval<T> {
    type Output = Interval<T>;

    fn sub(mut self, rhs: &Interval<T>) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<T: IntervalContent> Sub<Interval<T>> for &Interval<T> {
    type Output = Interval<T>;

    fn sub(self, rhs: Interval<T>) -> Self::Output {
        Interval {
            bounds: self.bounds.as_ref().zip(rhs.bounds).and_then(|([a0, a1], [b0, b1])|
                Interval::empty_if_overflowed([a0.sub_down_rv(b1), a1.sub_up_rv(b0)])),
            maybe_invalid: self.maybe_invalid || rhs.maybe_invalid
        }
    }
}

impl<T: IntervalContent> Sub<&Interval<T>> for &Interval<T> {
    type Output = Interval<T>;

    fn sub(self, rhs: &Interval<T>) -> Self::Output {
        Interval {
            bounds: self.bounds.as_ref().zip(rhs.bounds.as_ref()).and_then(|([a0, a1], [b0, b1])|
                Interval::empty_if_overflowed([a0.sub_down_rr(b1), a1.sub_up_rr(b0)])),
            maybe_invalid: self.maybe_invalid || rhs.maybe_invalid
        }
    }
}

impl<T: IntervalContent> Neg for Interval<T> {
    type Output = Interval<T>;

    fn neg(self) -> Self::Output {
        Interval {
            bounds: self.bounds.map(|[a, b]| [b.neg_v(), a.neg_v()]),
            maybe_invalid: self.maybe_invalid
        }
    }
}

impl<T: IntervalContent> Neg for &Interval<T> {
    type Output = Interval<T>;

    fn neg(self) -> Self::Output {
        Interval {
            bounds: self.bounds.as_ref().map(|[a, b]| [b.neg_r(), a.neg_r()]),
            maybe_invalid: self.maybe_invalid
        }
    }
}

impl<T: IntervalContent> MulAssign<Interval<T>> for Interval<T> {
    fn mul_assign(&mut self, rhs: Interval<T>) {
        self.bounds = std::mem::take(&mut self.bounds).zip(rhs.bounds).and_then(|([a, b], [c, d])| {
            Interval::empty_if_overflowed(if b.true_sign() <= 0 {
                if d.true_sign() <= 0      { [b.mul_down_vv(d), a.mul_up_vv(c)] }
                else if c.true_sign() >= 0 { [a.mul_down_vv(d), b.mul_up_vv(c)] }
                else                     { [a.mul_down_rv(d), a.mul_up_vv(c)] }
            } else if a.true_sign() >= 0 {
                if d.true_sign() <= 0      { [b.mul_down_vv(c), a.mul_up_vv(d)] }
                else if c.true_sign() >= 0 { [a.mul_down_vv(c), b.mul_up_vv(d)] }
                else                     { [b.mul_down_rv(c), b.mul_up_vv(d)] }
            } else {
                if d.true_sign() <= 0      { [b.mul_down_vr(&c), a.mul_up_vv(c)] }
                else if c.true_sign() >= 0 { [a.mul_down_vr(&d), b.mul_up_vv(d)] }
                else { [a.mul_down_rr(&d).min_vv(b.mul_down_rr(&c)), a.mul_up_vv(c).max_vv(b.mul_up_vv(d)) ]}
            })
        });
        self.maybe_invalid |= rhs.maybe_invalid;
    }
}

impl<T: IntervalContent> MulAssign<&Interval<T>> for Interval<T> {
    fn mul_assign(&mut self, rhs: &Interval<T>) {
        self.bounds = std::mem::take(&mut self.bounds).zip(rhs.bounds.as_ref()).and_then(|([a, b], [c, d])| {
            Interval::empty_if_overflowed(if b.true_sign() <= 0 {
                if d.true_sign() <= 0      { [b.mul_down_vr(d), a.mul_up_vr(c)] }
                else if c.true_sign() >= 0 { [a.mul_down_vr(d), b.mul_up_vr(c)] }
                else                     { [a.mul_down_rr(d), a.mul_up_vr(c)] }
            } else if a.true_sign() >= 0 {
                if d.true_sign() <= 0      { [b.mul_down_vr(c), a.mul_up_vr(d)] }
                else if c.true_sign() >= 0 { [a.mul_down_vr(c), b.mul_up_vr(d)] }
                else                     { [b.mul_down_rr(c), b.mul_up_vr(d)] }
            } else {
                if d.true_sign() <= 0      { [b.mul_down_vr(&c), a.mul_up_vr(c)] }
                else if c.true_sign() >= 0 { [a.mul_down_vr(&d), b.mul_up_vr(d)] }
                else { [a.mul_down_rr(&d).min_vv(b.mul_down_rr(&c)), a.mul_up_vr(c).max_vv(b.mul_up_vr(d)) ]}
            })
        });
        self.maybe_invalid |= rhs.maybe_invalid;
    }
}

impl<T: IntervalContent> Mul<Interval<T>> for Interval<T> {
    type Output = Interval<T>;

    fn mul(mut self, rhs: Interval<T>) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<T: IntervalContent> Mul<&Interval<T>> for Interval<T> {
    type Output = Interval<T>;

    fn mul(mut self, rhs: &Interval<T>) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<T: IntervalContent> Mul<Interval<T>> for &Interval<T> {
    type Output = Interval<T>;

    fn mul(self, rhs: Interval<T>) -> Self::Output {
        rhs * self
    }
}

impl<T: IntervalContent> Mul<&Interval<T>> for &Interval<T> {
    type Output = Interval<T>;

    fn mul(self, rhs: &Interval<T>) -> Self::Output {
        Interval {
            bounds: self.bounds.as_ref().zip(rhs.bounds.as_ref()).and_then(|([a, b], [c, d])| {
                Interval::empty_if_overflowed(if b.true_sign() <= 0 {
                    if d.true_sign() <= 0      { [b.mul_down_rr(d), a.mul_up_rr(c)] }
                    else if c.true_sign() >= 0 { [a.mul_down_rr(d), b.mul_up_rr(c)] }
                    else                     { [a.mul_down_rr(d), a.mul_up_rr(c)] }
                } else if a.true_sign() >= 0 {
                    if d.true_sign() <= 0      { [b.mul_down_rr(c), a.mul_up_rr(d)] }
                    else if c.true_sign() >= 0 { [a.mul_down_rr(c), b.mul_up_rr(d)] }
                    else                     { [b.mul_down_rr(c), b.mul_up_rr(d)] }
                } else {
                    if d.true_sign() <= 0      { [b.mul_down_rr(&c), a.mul_up_rr(c)] }
                    else if c.true_sign() >= 0 { [a.mul_down_rr(&d), b.mul_up_rr(d)] }
                    else { [a.mul_down_rr(&d).min_vv(b.mul_down_rr(&c)), a.mul_up_rr(c).max_vv(b.mul_up_rr(d)) ]}
                })
            }),
            maybe_invalid: self.maybe_invalid || rhs.maybe_invalid
        }
    }
}

impl<T: IntervalContent> DivAssign<Interval<T>> for Interval<T> {
    fn div_assign(&mut self, rhs: Interval<T>) {
        let div_by_0 = rhs.contains_zero();
        self.bounds = std::mem::take(&mut self.bounds).zip(rhs.bounds).and_then(|([a, b], [c, d])| {
            if c.true_sign() == 0 && d.true_sign() == 0 { return None; }
            Interval::empty_if_overflowed(
                if a.true_sign() == 0 && b.true_sign() == 0 {
                    [a, b]
                } else if d.true_sign() < 0 {
                    if b.true_sign() <= 0      { [b.div_down_vv(c),  a.div_up_vv(d)] }
                    else if a.true_sign() >= 0 { [b.div_down_vv(d),  a.div_up_vv(c)] }
                    else                       { [b.div_down_vr(&d), a.div_up_vv(d)] }
                } else if c.true_sign() > 0 {
                    if b.true_sign() <= 0      { [a.div_down_vv(c),  b.div_up_vv(d)] }
                    else if a.true_sign() >= 0 { [a.div_down_vv(d),  b.div_up_vv(c)] }
                    else                       { [a.div_down_vr(&c), b.div_up_vv(c)] }
                } else if d.true_sign() == 0 {
                    if b.true_sign() <= 0      { [b.div_down_vv(c), T::pos_inf()] }
                    else if a.true_sign() >= 0 { [T::neg_inf(), a.div_up_vv(c)] }
                    else                       { [T::neg_inf(), T::pos_inf()] }
                } else if c.true_sign() == 0 {
                    if b.true_sign() <= 0      { [T::neg_inf(), b.div_up_vv(d)] }
                    else if a.true_sign() >= 0 { [a.div_down_vv(d), T::pos_inf()] }
                    else                       { [T::neg_inf(), T::pos_inf()] }
                } else {
                    [T::neg_inf(), T::pos_inf()]
                }
            )
        });
        self.maybe_invalid |= rhs.maybe_invalid || div_by_0;
    }
}

impl<T: IntervalContent> DivAssign<&Interval<T>> for Interval<T> {
    fn div_assign(&mut self, rhs: &Interval<T>) {
        let div_by_0 = rhs.contains_zero();
        self.bounds = std::mem::take(&mut self.bounds).zip(rhs.bounds.as_ref()).and_then(|([a, b], [c, d])| {
            if c.true_sign() == 0 && d.true_sign() == 0 { return None; }
            Interval::empty_if_overflowed(
                if a.true_sign() == 0 && b.true_sign() == 0 {
                    [a, b]
                } else if d.true_sign() < 0 {
                    if b.true_sign() <= 0      { [b.div_down_vr(&c), a.div_up_vr(&d)] }
                    else if a.true_sign() >= 0 { [b.div_down_vr(&d), a.div_up_vr(&c)] }
                    else                       { [b.div_down_vr(&d), a.div_up_vr(&d)] }
                } else if c.true_sign() > 0 {
                    if b.true_sign() <= 0      { [a.div_down_vr(&c), b.div_up_vr(&d)] }
                    else if a.true_sign() >= 0 { [a.div_down_vr(&d), b.div_up_vr(&c)] }
                    else                       { [a.div_down_vr(&c), b.div_up_vr(&c)] }
                } else if d.true_sign() == 0 {
                    if b.true_sign() <= 0      { [b.div_down_vr(&c), T::pos_inf()] }
                    else if a.true_sign() >= 0 { [T::neg_inf(), a.div_up_vr(&c)] }
                    else                       { [T::neg_inf(), T::pos_inf()] }
                } else if c.true_sign() == 0 {
                    if b.true_sign() <= 0      { [T::neg_inf(), b.div_up_vr(&d)] }
                    else if a.true_sign() >= 0 { [a.div_down_vr(&d), T::pos_inf()] }
                    else                       { [T::neg_inf(), T::pos_inf()] }
                } else {
                    [T::neg_inf(), T::pos_inf()]
                }
            )
        });
        self.maybe_invalid |= rhs.maybe_invalid || div_by_0;
    }
}

impl<T: IntervalContent> Div<Interval<T>> for Interval<T> {
    type Output = Interval<T>;

    fn div(mut self, rhs: Interval<T>) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<T: IntervalContent> Div<&Interval<T>> for Interval<T> {
    type Output = Interval<T>;

    fn div(mut self, rhs: &Interval<T>) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<T: IntervalContent> Div<Interval<T>> for &Interval<T> {
    type Output = Interval<T>;

    fn div(self, rhs: Interval<T>) -> Self::Output {
        let div_by_0 = rhs.contains_zero();
        Interval {
            bounds: self.bounds.as_ref().zip(rhs.bounds).and_then(|([a, b], [c, d])| {
                if c.true_sign() == 0 && d.true_sign() == 0 { return None; }
                Interval::empty_if_overflowed(
                    if a.true_sign() == 0 && b.true_sign() == 0 {
                        [a.clone(), b.clone()]
                    } else if d.true_sign() < 0 {
                        if b.true_sign() <= 0      { [b.div_down_rv(c),  a.div_up_rv(d)] }
                        else if a.true_sign() >= 0 { [b.div_down_rv(d),  a.div_up_rv(c)] }
                        else                       { [b.div_down_rr(&d), a.div_up_rv(d)] }
                    } else if c.true_sign() > 0 {
                        if b.true_sign() <= 0      { [a.div_down_rv(c),  b.div_up_rv(d)] }
                        else if a.true_sign() >= 0 { [a.div_down_rv(d),  b.div_up_rv(c)] }
                        else                       { [a.div_down_rr(&c), b.div_up_rv(c)] }
                    } else if d.true_sign() == 0 {
                        if b.true_sign() <= 0      { [b.div_down_rv(c), T::pos_inf()] }
                        else if a.true_sign() >= 0 { [T::neg_inf(), a.div_up_rv(c)] }
                        else                       { [T::neg_inf(), T::pos_inf()] }
                    } else if c.true_sign() == 0 {
                        if b.true_sign() <= 0      { [T::neg_inf(), b.div_up_rv(d)] }
                        else if a.true_sign() >= 0 { [a.div_down_rv(d), T::pos_inf()] }
                        else                       { [T::neg_inf(), T::pos_inf()] }
                    } else {
                        [T::neg_inf(), T::pos_inf()]
                    }
                )
            }),
            maybe_invalid: self.maybe_invalid || rhs.maybe_invalid || div_by_0
        }
    }
}

impl<T: IntervalContent> Div<&Interval<T>> for &Interval<T> {
    type Output = Interval<T>;

    fn div(self, rhs: &Interval<T>) -> Self::Output {
        let div_by_0 = rhs.contains_zero();
        Interval {
            bounds: self.bounds.as_ref().zip(rhs.bounds.as_ref()).and_then(|([a, b], [c, d])| {
                if c.true_sign() == 0 && d.true_sign() == 0 { return None; }
                Interval::empty_if_overflowed(
                    if a.true_sign() == 0 && b.true_sign() == 0 {
                        [a.clone(), b.clone()]
                    } else if d.true_sign() < 0 {
                        if b.true_sign() <= 0      { [b.div_down_rr(&c), a.div_up_rr(&d)] }
                        else if a.true_sign() >= 0 { [b.div_down_rr(&d), a.div_up_rr(&c)] }
                        else                       { [b.div_down_rr(&d), a.div_up_rr(&d)] }
                    } else if c.true_sign() > 0 {
                        if b.true_sign() <= 0      { [a.div_down_rr(&c), b.div_up_rr(&d)] }
                        else if a.true_sign() >= 0 { [a.div_down_rr(&d), b.div_up_rr(&c)] }
                        else                       { [a.div_down_rr(&c), b.div_up_rr(&c)] }
                    } else if d.true_sign() == 0 {
                        if b.true_sign() <= 0      { [b.div_down_rr(&c), T::pos_inf()] }
                        else if a.true_sign() >= 0 { [T::neg_inf(), a.div_up_rr(&c)] }
                        else                       { [T::neg_inf(), T::pos_inf()] }
                    } else if c.true_sign() == 0 {
                        if b.true_sign() <= 0      { [T::neg_inf(), b.div_up_rr(&d)] }
                        else if a.true_sign() >= 0 { [a.div_down_rr(&d), T::pos_inf()] }
                        else                       { [T::neg_inf(), T::pos_inf()] }
                    } else {
                        [T::neg_inf(), T::pos_inf()]
                    }
                )
            }),
            maybe_invalid: self.maybe_invalid || rhs.maybe_invalid || div_by_0
        }
    }
}

impl<T: IntervalContent> Sum for Interval<T> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(|a, b| a + b).expect("Interval needs to have elements")
    }
}

impl<T: IntervalContent> From<T> for Interval<T> {
    fn from(value: T) -> Self {
        let clone = value.clone();
        Self::new([value, clone])
    }
}

impl<T: IntervalContent> From<[T; 2]> for Interval<T> {
    fn from(value: [T; 2]) -> Self {
        Self::new(value)
    }
}

/// A helper trait for numbers that can be put in intervals
pub trait IntervalContent: Clone {
    fn from_integer_down(int: Integer, context: &Self) -> Self;
    fn from_integer_up(int: Integer, context: &Self) -> Self;

    fn add_down_vv(self, rhs: Self) -> Self;
    fn add_down_vr(self, rhs: &Self) -> Self;
    fn add_down_rv(&self, rhs: Self) -> Self;
    fn add_down_rr(&self, rhs: &Self) -> Self;
    fn add_up_vv(self, rhs: Self) -> Self;
    fn add_up_vr(self, rhs: &Self) -> Self;
    fn add_up_rv(&self, rhs: Self) -> Self;
    fn add_up_rr(&self, rhs: &Self) -> Self;

    fn sub_down_vv(self, rhs: Self) -> Self;
    fn sub_down_vr(self, rhs: &Self) -> Self;
    fn sub_down_rv(&self, rhs: Self) -> Self;
    fn sub_down_rr(&self, rhs: &Self) -> Self;
    fn sub_up_vv(self, rhs: Self) -> Self;
    fn sub_up_vr(self, rhs: &Self) -> Self;
    fn sub_up_rv(&self, rhs: Self) -> Self;
    fn sub_up_rr(&self, rhs: &Self) -> Self;

    fn neg_v( self) -> Self;
    fn neg_r(&self) -> Self;

    fn mul_down_vv(self, rhs: Self) -> Self;
    fn mul_down_vr(self, rhs: &Self) -> Self;
    fn mul_down_rv(&self, rhs: Self) -> Self;
    fn mul_down_rr(&self, rhs: &Self) -> Self;
    fn mul_up_vv(self, rhs: Self) -> Self;
    fn mul_up_vr(self, rhs: &Self) -> Self;
    fn mul_up_rv(&self, rhs: Self) -> Self;
    fn mul_up_rr(&self, rhs: &Self) -> Self;

    fn div_down_vv(self, rhs: Self) -> Self;
    fn div_down_vr(self, rhs: &Self) -> Self;
    fn div_down_rv(&self, rhs: Self) -> Self;
    fn div_down_rr(&self, rhs: &Self) -> Self;
    fn div_up_vv(self, rhs: Self) -> Self;
    fn div_up_vr(self, rhs: &Self) -> Self;
    fn div_up_rv(&self, rhs: Self) -> Self;
    fn div_up_rr(&self, rhs: &Self) -> Self;
    fn neg_inf() -> Self;
    fn pos_inf() -> Self;
    fn zero(&self) -> Self; // generates a zero with context

    fn sqrt_down_v(self) -> Self;
    fn sqrt_down_r(&self) -> Self;
    fn sqrt_up_v(self) -> Self;
    fn sqrt_up_r(&self) -> Self;

    fn min_vv(self, rhs: Self) -> Self;
    fn min_vr(self, rhs: &Self) -> Self;
    fn min_rv(&self, rhs: Self) -> Self;
    fn min_rr(&self, rhs: &Self) -> Self;
    fn max_vv(self, rhs: Self) -> Self;
    fn max_vr(self, rhs: &Self) -> Self;
    fn max_rv(&self, rhs: Self) -> Self;
    fn max_rr(&self, rhs: &Self) -> Self;

    /// Unlike signum(), should return 0 if the number is 0, even in types with +0 and -0.
    fn true_sign(&self) -> i32;

    fn overflowed(value: &[Self; 2]) -> bool where Self: Sized;
}


// Does not assume a rounding mode, except that a number is rounded to one of the two closest
// f64s to it.
impl IntervalContent for f64 {
    fn from_integer_down(int: Integer, _context: &Self) -> Self { f64::rounding_from(&int, RoundingMode::Floor).0 }
    fn from_integer_up  (int: Integer, _context: &Self) -> Self { f64::rounding_from(&int, RoundingMode::Ceiling).0 }

    // Prioritizing speed over bound tightness here.
    // It's probably rare that the difference will matter.
    fn add_down_vv(self, rhs: Self) -> Self { f64_down(self + rhs) }
    fn add_down_vr(self, rhs: &Self) -> Self { self.add_down_vv(*rhs) }
    fn add_down_rv(&self, rhs: Self) -> Self { self.add_down_vv(rhs) }
    fn add_down_rr(&self, rhs: &Self) -> Self { self.add_down_vv(*rhs) }
    fn add_up_vv(self, rhs: Self) -> Self { f64_up(self + rhs) }
    fn add_up_vr(self, rhs: &Self) -> Self { self.add_up_vv(*rhs) }
    fn add_up_rv(&self, rhs: Self) -> Self { self.add_up_vv(rhs) }
    fn add_up_rr(&self, rhs: &Self) -> Self { self.add_up_vv(*rhs) }

    fn sub_down_vv(self, rhs: Self) -> Self { f64_down(self - rhs) }
    fn sub_down_vr(self, rhs: &Self) -> Self { self.sub_down_vv(*rhs) }
    fn sub_down_rv(&self, rhs: Self) -> Self { self.sub_down_vv(rhs) }
    fn sub_down_rr(&self, rhs: &Self) -> Self { self.sub_down_vv(*rhs) }
    fn sub_up_vv(self, rhs: Self) -> Self { f64_up(self - rhs) }
    fn sub_up_vr(self, rhs: &Self) -> Self { self.sub_up_vv(*rhs) }
    fn sub_up_rv(&self, rhs: Self) -> Self { self.sub_up_vv(rhs) }
    fn sub_up_rr(&self, rhs: &Self) -> Self { self.sub_up_vv(*rhs) }

    fn neg_v( self) -> Self { -self }
    fn neg_r(&self) -> Self { -*self }

    // Zeros stay zero, other numbers get interval bloat regardless of whether the result
    // can be represented exactly
    fn mul_down_vv(self, rhs: Self) -> Self { if self == 0.0 || rhs == 0.0 { 0.0 } else { f64_down(self * rhs) }}
    fn mul_down_vr(self, rhs: &Self) -> Self { self.mul_down_vv(*rhs) }
    fn mul_down_rv(&self, rhs: Self) -> Self { self.mul_down_vv(rhs) }
    fn mul_down_rr(&self, rhs: &Self) -> Self { self.mul_down_vv(*rhs) }
    fn mul_up_vv(self, rhs: Self) -> Self { if self == 0.0 || rhs == 0.0 { 0.0 } else { f64_up(self * rhs) } }
    fn mul_up_vr(self, rhs: &Self) -> Self { self.mul_up_vv(*rhs) }
    fn mul_up_rv(&self, rhs: Self) -> Self { self.mul_up_vv(rhs) }
    fn mul_up_rr(&self, rhs: &Self) -> Self { self.mul_up_vv(*rhs) }

    // Zeros stay zero, other numbers get interval bloat regardless of whether the result
    // can be represented exactly
    fn div_down_vv(self, rhs: Self) -> Self { if self == 0.0 || rhs.is_infinite() { 0.0 } else { f64_down(self / rhs) }}
    fn div_down_vr(self, rhs: &Self) -> Self { self.div_down_vv(*rhs) }
    fn div_down_rv(&self, rhs: Self) -> Self { self.div_down_vv(rhs) }
    fn div_down_rr(&self, rhs: &Self) -> Self { self.div_down_vv(*rhs) }
    fn div_up_vv(self, rhs: Self) -> Self { if self == 0.0 || rhs.is_infinite() { 0.0 } else { f64_up(self / rhs) } }
    fn div_up_vr(self, rhs: &Self) -> Self { self.div_up_vv(*rhs) }
    fn div_up_rv(&self, rhs: Self) -> Self { self.div_up_vv(rhs) }
    fn div_up_rr(&self, rhs: &Self) -> Self { self.div_up_vv(*rhs) }
    fn neg_inf() -> Self { f64::NEG_INFINITY }
    fn pos_inf() -> Self { f64::INFINITY }
    fn zero(&self) -> Self { 0.0 }
    
    fn sqrt_down_v(self) -> Self { if self == 0.0 { 0.0 } else { f64_down(self.sqrt()) } }
    fn sqrt_down_r(&self) -> Self { if self == &0.0 { 0.0 } else { f64_down(self.sqrt()) } }
    fn sqrt_up_v(self) -> Self { if self == 0.0 { 0.0 } else { f64_up(self.sqrt()) } }
    fn sqrt_up_r(&self) -> Self { if self == &0.0 { 0.0 } else { f64_up(self.sqrt()) } }

    fn min_vv(self, rhs: Self) -> Self { self.min(rhs) }
    fn min_vr(self, rhs: &Self) -> Self { self.min_vv(*rhs) }
    fn min_rv(&self, rhs: Self) -> Self { self.min(rhs) }
    fn min_rr(&self, rhs: &Self) -> Self { self.min_vv(*rhs) }
    fn max_vv(self, rhs: Self) -> Self { self.max(rhs) }
    fn max_vr(self, rhs: &Self) -> Self { self.max_vv(*rhs) }
    fn max_rv(&self, rhs: Self) -> Self { self.max(rhs) }
    fn max_rr(&self, rhs: &Self) -> Self { self.max_vv(*rhs) }

    fn true_sign(&self) -> i32 { if self < &0.0 { -1 } else if self == &0.0 { 0 } else { 1 } }

    fn overflowed(value: &[Self; 2]) -> bool where Self: Sized {
        value[0] == f64::INFINITY || value[1] == f64::NEG_INFINITY
    }
}

/// Get the previous float.
fn f64_down(x: f64) -> f64 {
    x.next_down()
}

/// Get the next float.
fn f64_up(x: f64) -> f64 {
    x.next_up()
}

/// Significand and number of bits after the point.
#[derive(Clone, Debug)]
pub struct Fixed(pub Integer, pub u64);

impl PartialEq for Fixed {
    fn eq(&self, other: &Self) -> bool {
        assert_eq!(self.1, other.1);
        self.0 == other.0
    }
}

impl Eq for Fixed {}

impl PartialOrd for Fixed {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        assert_eq!(self.1, other.1);
        self.0.partial_cmp(&other.0)
    }
}

impl Ord for Fixed {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        assert_eq!(self.1, other.1);
        self.0.cmp(&other.0)
    }
}

impl IntervalContent for Fixed {
    fn from_integer_down(int: Integer, context: &Self) -> Self { Fixed(int << context.1, context.1) }
    fn from_integer_up  (int: Integer, context: &Self) -> Self { Self::from_integer_down(int, context) }

    fn add_down_vv( self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed( self.0 +  rhs.0, self.1) }
    fn add_down_vr( self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed( self.0 + &rhs.0, self.1) }
    fn add_down_rv(&self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(&self.0 +  rhs.0, self.1) }
    fn add_down_rr(&self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(&self.0 + &rhs.0, self.1) }
    fn add_up_vv( self, rhs:  Self) -> Self { self.add_down_vv(rhs) }
    fn add_up_vr( self, rhs: &Self) -> Self { self.add_down_vr(rhs) }
    fn add_up_rv(&self, rhs:  Self) -> Self { self.add_down_rv(rhs) }
    fn add_up_rr(&self, rhs: &Self) -> Self { self.add_down_rr(rhs) }

    fn sub_down_vv( self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed( self.0 -  rhs.0, self.1) }
    fn sub_down_vr( self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed( self.0 - &rhs.0, self.1) }
    fn sub_down_rv(&self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(&self.0 -  rhs.0, self.1) }
    fn sub_down_rr(&self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(&self.0 - &rhs.0, self.1) }
    fn sub_up_vv( self, rhs:  Self) -> Self { self.sub_down_vv(rhs) }
    fn sub_up_vr( self, rhs: &Self) -> Self { self.sub_down_vr(rhs) }
    fn sub_up_rv(&self, rhs:  Self) -> Self { self.sub_down_rv(rhs) }
    fn sub_up_rr(&self, rhs: &Self) -> Self { self.sub_down_rr(rhs) }

    fn neg_v( self) -> Self { Fixed(- self.0, self.1) }
    fn neg_r(&self) -> Self { Fixed(-&self.0, self.1) }

    fn mul_down_vv( self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed( self.0 *  rhs.0 >> self.1, self.1) }
    fn mul_down_vr( self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed( self.0 * &rhs.0 >> self.1, self.1) }
    fn mul_down_rv(&self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(&self.0 *  rhs.0 >> self.1, self.1) }
    fn mul_down_rr(&self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(&self.0 * &rhs.0 >> self.1, self.1) }
    fn mul_up_vv( self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(( self.0 *  rhs.0).shr_round(self.1, RoundingMode::Ceiling).0, self.1) }
    fn mul_up_vr( self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(( self.0 * &rhs.0).shr_round(self.1, RoundingMode::Ceiling).0, self.1) }
    fn mul_up_rv(&self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed((&self.0 *  rhs.0).shr_round(self.1, RoundingMode::Ceiling).0, self.1) }
    fn mul_up_rr(&self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed((&self.0 * &rhs.0).shr_round(self.1, RoundingMode::Ceiling).0, self.1) }

    fn div_down_vv( self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(( self.0 << self.1).div_round( rhs.0, RoundingMode::Floor).0, self.1) }
    fn div_down_vr( self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(( self.0 << self.1).div_round(&rhs.0, RoundingMode::Floor).0, self.1) }
    fn div_down_rv(&self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed((&self.0 << self.1).div_round( rhs.0, RoundingMode::Floor).0, self.1) }
    fn div_down_rr(&self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed((&self.0 << self.1).div_round(&rhs.0, RoundingMode::Floor).0, self.1) }
    fn div_up_vv( self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(( self.0 << self.1).div_round( rhs.0, RoundingMode::Ceiling).0, self.1) }
    fn div_up_vr( self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(( self.0 << self.1).div_round(&rhs.0, RoundingMode::Ceiling).0, self.1) }
    fn div_up_rv(&self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed((&self.0 << self.1).div_round( rhs.0, RoundingMode::Ceiling).0, self.1) }
    fn div_up_rr(&self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed((&self.0 << self.1).div_round(&rhs.0, RoundingMode::Ceiling).0, self.1) }
    fn pos_inf() -> Self { panic!("I said not to divide Interval<> by 0!") }
    fn neg_inf() -> Self { panic!("I said not to divide Interval<> by 0!") }
    fn zero(&self) -> Self { Fixed(Integer::ZERO, self.1) }

    fn sqrt_down_v( self) -> Self { Fixed(( self.0 << self.1).floor_sqrt(), self.1) }
    fn sqrt_down_r(&self) -> Self { Fixed((&self.0 << self.1).floor_sqrt(), self.1) }
    fn sqrt_up_v( self) -> Self { Fixed(( self.0 << self.1).ceiling_sqrt(), self.1) }
    fn sqrt_up_r(&self) -> Self { Fixed((&self.0 << self.1).ceiling_sqrt(), self.1) }

    fn min_vv( self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(if self.0 < rhs.0 { self.0         } else { rhs.0         }, self.1) }
    fn min_vr( self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(if self.0 < rhs.0 { self.0         } else { rhs.0.clone() }, self.1) }
    fn min_rv(&self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(if self.0 < rhs.0 { self.0.clone() } else { rhs.0         }, self.1) }
    fn min_rr(&self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(if self.0 < rhs.0 { self.0.clone() } else { rhs.0.clone() }, self.1) }
    fn max_vv( self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(if self.0 > rhs.0 { self.0         } else { rhs.0         }, self.1) }
    fn max_vr( self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(if self.0 > rhs.0 { self.0         } else { rhs.0.clone() }, self.1) }
    fn max_rv(&self, rhs:  Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(if self.0 > rhs.0 { self.0.clone() } else { rhs.0         }, self.1) }
    fn max_rr(&self, rhs: &Self) -> Self { assert_eq!(self.1, rhs.1); Fixed(if self.0 > rhs.0 { self.0.clone() } else { rhs.0.clone() }, self.1) }

    fn overflowed(_value: &[Self; 2]) -> bool where Self: Sized { false } // Remember, no infinities!

    fn true_sign(&self) -> i32 {
        match self.0.sign() {
            std::cmp::Ordering::Less => -1,
            std::cmp::Ordering::Equal => 0,
            std::cmp::Ordering::Greater => 1
        }
    }
}

#[cfg(test)]
mod test {
    use std::cmp::Ordering;

    use crate::interval::{Fixed, Interval, IntervalContent};

    fn add_test<T: IntervalContent + Clone + PartialEq + std::fmt::Debug>(a: Interval<T>, b: Interval<T>, expected: Interval<T>) {
        assert_eq!(a.clone() + b.clone(), expected);
        assert_eq!(a.clone() + &b, expected);
        assert_eq!(&a + b.clone(), expected);
        assert_eq!(&a + &b, expected);
    }

    fn sub_test<T: IntervalContent + Clone + PartialEq + std::fmt::Debug>(a: Interval<T>, b: Interval<T>, expected: Interval<T>) {
        assert_eq!(a.clone() - b.clone(), expected);
        assert_eq!(a.clone() - &b, expected);
        assert_eq!(&a - b.clone(), expected);
        assert_eq!(&a - &b, expected);
    }

    fn neg_test<T: IntervalContent + Clone + PartialEq + std::fmt::Debug>(a: Interval<T>, expected: Interval<T>) {
        assert_eq!(-a.clone(), expected);
        assert_eq!(-&a, expected);
    }

    fn mul_test<T: IntervalContent + Clone + PartialEq + std::fmt::Debug>(a: Interval<T>, b: Interval<T>, expected: Interval<T>) {
        assert_eq!(a.clone() * b.clone(), expected);
        assert_eq!(a.clone() * &b, expected);
        assert_eq!(&a * b.clone(), expected);
        assert_eq!(&a * &b, expected);
    }

    fn div_bounds_test<T: IntervalContent + Clone + PartialEq + std::fmt::Debug>(a: Interval<T>, b: Interval<T>, expected: Interval<T>) {
        assert_eq!((a.clone() / b.clone()).bounds, expected.bounds);
        assert_eq!((a.clone() / &b).bounds, expected.bounds);
        assert_eq!((&a / b.clone()).bounds, expected.bounds);
        assert_eq!((&a / &b).bounds, expected.bounds);
    }

    fn sqrt_bounds_test<T: IntervalContent + Clone + PartialEq + std::fmt::Debug>(a: Interval<T>, expected: Interval<T>) {
        assert_eq!(a.sqrt().bounds, expected.bounds);
    }

    fn cmp_zero_test<T: IntervalContent + Clone + PartialEq + std::fmt::Debug>(a: Interval<T>, expected: Option<Ordering>) {
        assert_eq!(a.cmp_zero(), expected);
    }

    fn n2i(a: f64, b: f64) -> Interval<f64> {
        Interval::new([a, b])
    }

    fn i2i(a: i64, b: i64, prec: u64) -> Interval<Fixed> {
        Interval::new([Fixed(a.into(), prec), Fixed(b.into(), prec)])
    }

    // libieeep1788 tests, but only applicable ones, and modified to account for rounding
    #[test]
    fn test_interval_add_minimal_float() {
        type I = Interval<f64>;
        add_test::<f64>(I::EMPTY, I::EMPTY, I::EMPTY);
        add_test::<f64>(n2i(-1.0, 1.0), I::EMPTY, I::EMPTY);
        add_test::<f64>(I::EMPTY, n2i(-1.0, 1.0), I::EMPTY);
        add_test::<f64>(I::EMPTY, I::ENTIRE, I::EMPTY);
        add_test::<f64>(I::ENTIRE, I::EMPTY, I::EMPTY);
        add_test::<f64>(I::ENTIRE, n2i(f64::NEG_INFINITY, 1.0), I::ENTIRE);
        add_test::<f64>(I::ENTIRE, n2i(-1.0, 1.0), I::ENTIRE);
        add_test::<f64>(I::ENTIRE, n2i(-1.0, f64::INFINITY), I::ENTIRE);
        add_test::<f64>(I::ENTIRE, I::ENTIRE, I::ENTIRE);
        add_test::<f64>(n2i(f64::NEG_INFINITY, 1.0), I::ENTIRE, I::ENTIRE);
        add_test::<f64>(n2i(-1.0, 1.0), I::ENTIRE, I::ENTIRE);
        add_test::<f64>(n2i(-1.0, f64::INFINITY), I::ENTIRE, I::ENTIRE);
        add_test::<f64>(
            n2i(f64::NEG_INFINITY, 2.0), n2i(f64::NEG_INFINITY, 4.0),
            n2i(f64::NEG_INFINITY, 6.0f64.next_up())
        );
        add_test::<f64>(
            n2i(f64::NEG_INFINITY, 2.0), n2i(3.0, 4.0),
            n2i(f64::NEG_INFINITY, 6.0f64.next_up())
        );
        add_test::<f64>(
            n2i(f64::NEG_INFINITY, 2.0), n2i(3.0, f64::INFINITY),
            I::ENTIRE
        );
        add_test::<f64>(
            n2i(1.0, 2.0), n2i(f64::NEG_INFINITY, 4.0),
            n2i(f64::NEG_INFINITY, 6.0f64.next_up())
        );
        add_test::<f64>(n2i(1.0, 2.0), n2i(3.0, 4.0), n2i(4.0f64.next_down(), 6.0f64.next_up()));
        add_test::<f64>(
            n2i(1.0, 2.0), n2i(3.0, f64::INFINITY),
            n2i(4.0f64.next_down(), f64::INFINITY)
        );
        add_test::<f64>(
            n2i(1.0, f64::INFINITY), n2i(f64::NEG_INFINITY, 4.0),
            I::ENTIRE
        );
        add_test::<f64>(
            n2i(1.0, f64::INFINITY), n2i(3.0, 4.0),
            n2i(4.0f64.next_down(), f64::INFINITY)
        );
        add_test::<f64>(
            n2i(1.0, f64::INFINITY), n2i(3.0, f64::INFINITY),
            n2i(4.0f64.next_down(), f64::INFINITY)
        );
        add_test::<f64>(
            n2i(1.0, 1.7976931348623157e+308), n2i(3.0, 4.0),
            n2i(4.0f64.next_down(), f64::INFINITY)
        );
        add_test::<f64>(
            n2i(-1.7976931348623157e+308, 2.0), n2i(-3.0, 4.0),
            n2i(f64::NEG_INFINITY, 6.0f64.next_up())
        );
        add_test::<f64>(
            n2i(-1.7976931348623157e+308, 2.0), n2i(-3.0, 1.7976931348623157e+308),
            I::ENTIRE
        );
        add_test::<f64>(
            n2i(1.0, 1.7976931348623157e+308), n2i(0.0, 0.0),
            n2i(1.0f64.next_down(), 1.7976931348623157e+308f64.next_up())
        );
        add_test::<f64>(
            n2i(1.0, 1.7976931348623157e+308), n2i(-0.0, -0.0),
            n2i(1.0f64.next_down(), 1.7976931348623157e+308f64.next_up())
        );
        add_test::<f64>(n2i(0.0, 0.0), n2i(-3.0, 4.0), n2i((-3.0f64).next_down(), 4.0f64.next_up()));
        add_test::<f64>(
            n2i(-0.0, -0.0), n2i(-3.0, 1.7976931348623157e+308),
            n2i((-3.0f64).next_down(), 1.7976931348623157e+308f64.next_up())
        );
        add_test::<f64>(
            n2i(1.9999999999999964, 1.9999999999999964), n2i(0.1, 0.1),
            n2i(2.099999999999996, 2.0999999999999965f64.next_up())
        );
        add_test::<f64>(
            n2i(1.9999999999999964, 1.9999999999999964), n2i(-0.1, -0.1),
            n2i(1.8999999999999964f64.next_down(), 1.8999999999999966)
        );
        add_test::<f64>(
            n2i(-1.9999999999999964, 1.9999999999999964), n2i(0.1, 0.1),
            n2i(-1.8999999999999966, 2.0999999999999965f64.next_up())
        );
    }

    #[test]
    fn test_interval_add_minimal_fixed() {
        type I = Interval<Fixed>;
        add_test::<Fixed>(I::EMPTY, I::EMPTY, I::EMPTY);
        add_test::<Fixed>(i2i(-1, 1, 0), I::EMPTY, I::EMPTY);
        add_test::<Fixed>(I::EMPTY, i2i(-1, 1, 0), I::EMPTY);
        add_test::<Fixed>(i2i(1, 2, 0), i2i(3, 4, 0), i2i(4, 6, 0));
        add_test::<Fixed>(i2i(0, 0, 0), i2i(-3, 4, 0), i2i(-3, 4, 0));
        add_test::<Fixed>(
            i2i(13, 13, 4), i2i(1, 1, 4),
            i2i(14, 14, 4)
        );
    }

    #[test]
    fn test_interval_subtract_minimal_float() {
        type I = Interval<f64>;
        sub_test::<f64>(I::EMPTY, I::EMPTY, I::EMPTY);
        sub_test::<f64>(n2i(-1.0, 1.0), I::EMPTY, I::EMPTY);
        sub_test::<f64>(I::EMPTY, n2i(-1.0, 1.0), I::EMPTY);
        sub_test::<f64>(I::EMPTY, I::ENTIRE, I::EMPTY);
        sub_test::<f64>(I::ENTIRE, I::EMPTY, I::EMPTY);
        sub_test::<f64>(I::ENTIRE, n2i(f64::NEG_INFINITY, 1.0), I::ENTIRE);
        sub_test::<f64>(I::ENTIRE, n2i(-1.0, 1.0), I::ENTIRE);
        sub_test::<f64>(I::ENTIRE, n2i(-1.0, f64::INFINITY), I::ENTIRE);
        sub_test::<f64>(I::ENTIRE, I::ENTIRE, I::ENTIRE);
        sub_test::<f64>(n2i(f64::NEG_INFINITY, 1.0), I::ENTIRE, I::ENTIRE);
        sub_test::<f64>(n2i(-1.0, 1.0), I::ENTIRE, I::ENTIRE);
        sub_test::<f64>(n2i(-1.0, f64::INFINITY), I::ENTIRE, I::ENTIRE);
        sub_test::<f64>(
            n2i(f64::NEG_INFINITY, 2.0), n2i(f64::NEG_INFINITY, 4.0),
            I::ENTIRE
        );
        sub_test::<f64>(
            n2i(f64::NEG_INFINITY, 2.0), n2i(3.0, 4.0),
            n2i(f64::NEG_INFINITY, (-1.0f64).next_up())
        );
        sub_test::<f64>(
            n2i(f64::NEG_INFINITY, 2.0), n2i(3.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, (-1.0f64).next_up())
        );
        sub_test::<f64>(
            n2i(1.0, 2.0), n2i(f64::NEG_INFINITY, 4.0),
            n2i((-3.0f64).next_down(), f64::INFINITY)
        );
        sub_test::<f64>(n2i(1.0, 2.0), n2i(3.0, 4.0), n2i((-3.0f64).next_down(), (-1.0f64).next_up()));
        sub_test::<f64>(
            n2i(1.0, 2.0), n2i(3.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, (-1.0f64).next_up())
        );
        sub_test::<f64>(
            n2i(1.0, f64::INFINITY), n2i(f64::NEG_INFINITY, 4.0),
            n2i((-3.0f64).next_down(), f64::INFINITY)
        );
        sub_test::<f64>(
            n2i(1.0, f64::INFINITY), n2i(3.0, 4.0),
            n2i((-3.0f64).next_down(), f64::INFINITY)
        );
        sub_test::<f64>(n2i(1.0, f64::INFINITY), n2i(3.0, f64::INFINITY), I::ENTIRE);
        sub_test::<f64>(
            n2i(1.0, 1.7976931348623157e+308), n2i(-3.0, 4.0),
            n2i((-3.0f64).next_down(), f64::INFINITY)
        );
        sub_test::<f64>(
            n2i(-1.7976931348623157e+308, 2.0), n2i(3.0, 4.0),
            n2i(f64::NEG_INFINITY, (-1.0f64).next_up())
        );
        sub_test::<f64>(
            n2i(-1.7976931348623157e+308, 2.0), n2i(-1.7976931348623157e+308, 4.0),
            I::ENTIRE
        );
        sub_test::<f64>(
            n2i(1.0, 1.7976931348623157e+308), n2i(0.0, 0.0),
            n2i(1.0f64.next_down(), 1.7976931348623157e+308f64.next_up())
        );
        sub_test::<f64>(
            n2i(1.0, 1.7976931348623157e+308), n2i(-0.0, -0.0),
            n2i(1.0f64.next_down(), 1.7976931348623157e+308f64.next_up())
        );
        sub_test::<f64>(n2i(0.0, 0.0), n2i(-3.0, 4.0), n2i((-4.0f64).next_down(), 3.0f64.next_up()));
        sub_test::<f64>(
            n2i(-0.0, -0.0), n2i(-3.0, 1.7976931348623157e+308),
            n2i((-1.7976931348623157e+308f64).next_down(), 3.0f64.next_up())
        );
        sub_test::<f64>(
            n2i(1.9999999999999964, 1.9999999999999964), n2i(0.1, 0.1),
            n2i(1.8999999999999964f64.next_down(), 1.8999999999999966)
        );
        sub_test::<f64>(
            n2i(1.9999999999999964, 1.9999999999999964), n2i(-0.1, -0.1),
            n2i(2.099999999999996, 2.0999999999999965f64.next_up())
        );
        sub_test::<f64>(
            n2i(-1.9999999999999964, 1.9999999999999964), n2i(0.1, 0.1),
            n2i((-2.0999999999999965f64).next_down(), 1.8999999999999966)
        );
    }

    #[test]
    fn test_interval_subtract_minimal_fixed() {
        type I = Interval<Fixed>;
        sub_test::<Fixed>(I::EMPTY, I::EMPTY, I::EMPTY);
        sub_test::<Fixed>(i2i(-1, 1, 0), I::EMPTY, I::EMPTY);
        sub_test::<Fixed>(I::EMPTY, i2i(-1, 1, 0), I::EMPTY);
        sub_test::<Fixed>(i2i(1, 2, 0), i2i(3, 4, 0), i2i(-3, -1, 0));
        sub_test::<Fixed>(i2i(0, 0, 0), i2i(-3, 4, 0), i2i(-4, 3, 0));
        sub_test::<Fixed>(
            i2i(13, 13, 4), i2i(1, 1, 4),
            i2i(12, 12, 4)
        );
    }

    #[test]
    fn test_interval_mul_minimal_float() {
        type I = Interval<f64>;
        mul_test::<f64>(I::EMPTY, I::EMPTY, I::EMPTY);
        mul_test::<f64>(n2i(-1.0, 1.0), I::EMPTY, I::EMPTY);
        mul_test::<f64>(I::EMPTY, n2i(-1.0, 1.0), I::EMPTY);
        mul_test::<f64>(I::EMPTY, I::ENTIRE, I::EMPTY);
        mul_test::<f64>(I::ENTIRE, I::EMPTY, I::EMPTY);
        mul_test::<f64>(n2i(0.0, 0.0), I::EMPTY, I::EMPTY);
        mul_test::<f64>(I::EMPTY, n2i(0.0, 0.0), I::EMPTY);
        mul_test::<f64>(n2i(-0.0, -0.0), I::EMPTY, I::EMPTY);
        mul_test::<f64>(I::EMPTY, n2i(-0.0, -0.0), I::EMPTY);
        mul_test::<f64>(I::ENTIRE, n2i(0.0, 0.0), n2i(0.0, 0.0));
        mul_test::<f64>(I::ENTIRE, n2i(-0.0, -0.0), n2i(0.0, 0.0));
        mul_test::<f64>(I::ENTIRE, n2i(-5.0, -1.0), I::ENTIRE);
        mul_test::<f64>(I::ENTIRE, n2i(-5.0, 3.0), I::ENTIRE);
        mul_test::<f64>(I::ENTIRE, n2i(1.0, 3.0), I::ENTIRE);
        mul_test::<f64>(I::ENTIRE, n2i(f64::NEG_INFINITY, -1.0), I::ENTIRE);
        mul_test::<f64>(I::ENTIRE, n2i(f64::NEG_INFINITY, 3.0), I::ENTIRE);
        mul_test::<f64>(I::ENTIRE, n2i(-5.0, f64::INFINITY), I::ENTIRE);
        mul_test::<f64>(I::ENTIRE, n2i(1.0, f64::INFINITY), I::ENTIRE);
        mul_test::<f64>(I::ENTIRE, I::ENTIRE, I::ENTIRE);
        mul_test::<f64>(n2i(1.0, f64::INFINITY), n2i(0.0, 0.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(1.0, f64::INFINITY), n2i(-0.0, -0.0), n2i(0.0, 0.0));
        mul_test::<f64>(
            n2i(1.0, f64::INFINITY), n2i(-5.0, -1.0),
            n2i(f64::NEG_INFINITY, (-1.0f64).next_up())
        );
        mul_test::<f64>(n2i(1.0, f64::INFINITY), n2i(-5.0, 3.0), I::ENTIRE);
        mul_test::<f64>(
            n2i(1.0, f64::INFINITY), n2i(1.0, 3.0),
            n2i(1.0f64.next_down(), f64::INFINITY)
        );
        mul_test::<f64>(
            n2i(1.0, f64::INFINITY), n2i(f64::NEG_INFINITY, -1.0),
            n2i(f64::NEG_INFINITY, (-1.0f64).next_up())
        );
        mul_test::<f64>(
            n2i(1.0, f64::INFINITY), n2i(f64::NEG_INFINITY, 3.0),
            I::ENTIRE
        );
        mul_test::<f64>(
            n2i(1.0, f64::INFINITY), n2i(-5.0, f64::INFINITY),
            I::ENTIRE
        );
        mul_test::<f64>(
            n2i(1.0, f64::INFINITY), n2i(1.0, f64::INFINITY),
            n2i(1.0f64.next_down(), f64::INFINITY)
        );
        mul_test::<f64>(n2i(1.0, f64::INFINITY), I::ENTIRE, I::ENTIRE);
        mul_test::<f64>(n2i(-1.0, f64::INFINITY), n2i(0.0, 0.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(-1.0, f64::INFINITY), n2i(-0.0, -0.0), n2i(0.0, 0.0));
        mul_test::<f64>(
            n2i(-1.0, f64::INFINITY), n2i(-5.0, -1.0),
            n2i(f64::NEG_INFINITY, 5.0f64.next_up())
        );
        mul_test::<f64>(n2i(-1.0, f64::INFINITY), n2i(-5.0, 3.0), I::ENTIRE);
        mul_test::<f64>(
            n2i(-1.0, f64::INFINITY), n2i(1.0, 3.0),
            n2i((-3.0f64).next_down(), f64::INFINITY)
        );
        mul_test::<f64>(
            n2i(-1.0, f64::INFINITY), n2i(f64::NEG_INFINITY, -1.0),
            I::ENTIRE
        );
        mul_test::<f64>(
            n2i(-1.0, f64::INFINITY), n2i(f64::NEG_INFINITY, 3.0),
            I::ENTIRE
        );
        mul_test::<f64>(
            n2i(-1.0, f64::INFINITY), n2i(-5.0, f64::INFINITY),
            I::ENTIRE
        );
        mul_test::<f64>(
            n2i(-1.0, f64::INFINITY), n2i(1.0, f64::INFINITY),
            I::ENTIRE
        );
        mul_test::<f64>(n2i(-1.0, f64::INFINITY), I::ENTIRE, I::ENTIRE);
        mul_test::<f64>(n2i(f64::NEG_INFINITY, 3.0), n2i(0.0, 0.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(f64::NEG_INFINITY, 3.0), n2i(-0.0, -0.0), n2i(0.0, 0.0));
        mul_test::<f64>(
            n2i(f64::NEG_INFINITY, 3.0), n2i(-5.0, -1.0),
            n2i((-15.0f64).next_down(), f64::INFINITY)
        );
        mul_test::<f64>(n2i(f64::NEG_INFINITY, 3.0), n2i(-5.0, 3.0), I::ENTIRE);
        mul_test::<f64>(
            n2i(f64::NEG_INFINITY, 3.0), n2i(1.0, 3.0),
            n2i(f64::NEG_INFINITY, 9.0f64.next_up())
        );
        mul_test::<f64>(
            n2i(f64::NEG_INFINITY, 3.0), n2i(f64::NEG_INFINITY, -1.0),
            I::ENTIRE
        );
        mul_test::<f64>(
            n2i(f64::NEG_INFINITY, 3.0), n2i(f64::NEG_INFINITY, 3.0),
            I::ENTIRE
        );
        mul_test::<f64>(
            n2i(f64::NEG_INFINITY, 3.0), n2i(-5.0, f64::INFINITY),
            I::ENTIRE
        );
        mul_test::<f64>(
            n2i(f64::NEG_INFINITY, 3.0), n2i(1.0, f64::INFINITY),
            I::ENTIRE
        );
        mul_test::<f64>(n2i(f64::NEG_INFINITY, 3.0), I::ENTIRE, I::ENTIRE);
        mul_test::<f64>(n2i(f64::NEG_INFINITY, -3.0), n2i(0.0, 0.0), n2i(0.0, 0.0));
        mul_test::<f64>(
            n2i(f64::NEG_INFINITY, -3.0), n2i(-0.0, -0.0),
            n2i(0.0, 0.0)
        );
        mul_test::<f64>(
            n2i(f64::NEG_INFINITY, -3.0), n2i(-5.0, -1.0),
            n2i(3.0f64.next_down(), f64::INFINITY)
        );
        mul_test::<f64>(n2i(f64::NEG_INFINITY, -3.0), n2i(-5.0, 3.0), I::ENTIRE);
        mul_test::<f64>(
            n2i(f64::NEG_INFINITY, -3.0), n2i(1.0, 3.0),
            n2i(f64::NEG_INFINITY, (-3.0f64).next_up())
        );
        mul_test::<f64>(
            n2i(f64::NEG_INFINITY, -3.0), n2i(f64::NEG_INFINITY, -1.0),
            n2i(3.0f64.next_down(), f64::INFINITY)
        );
        mul_test::<f64>(
            n2i(f64::NEG_INFINITY, -3.0), n2i(f64::NEG_INFINITY, 3.0),
            I::ENTIRE
        );
        mul_test::<f64>(
            n2i(f64::NEG_INFINITY, -3.0), n2i(-5.0, f64::INFINITY),
            I::ENTIRE
        );
        mul_test::<f64>(
            n2i(f64::NEG_INFINITY, -3.0), n2i(1.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, (-3.0f64).next_up())
        );
        mul_test::<f64>(n2i(f64::NEG_INFINITY, -3.0), I::ENTIRE, I::ENTIRE);
        mul_test::<f64>(n2i(0.0, 0.0), n2i(0.0, 0.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(0.0, 0.0), n2i(-0.0, -0.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(0.0, 0.0), n2i(-5.0, -1.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(0.0, 0.0), n2i(-5.0, 3.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(0.0, 0.0), n2i(1.0, 3.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(0.0, 0.0), n2i(f64::NEG_INFINITY, -1.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(0.0, 0.0), n2i(f64::NEG_INFINITY, 3.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(0.0, 0.0), n2i(-5.0, f64::INFINITY), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(0.0, 0.0), n2i(1.0, f64::INFINITY), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(0.0, 0.0), I::ENTIRE, n2i(0.0, 0.0));
        mul_test::<f64>(n2i(-0.0, -0.0), n2i(0.0, 0.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(-0.0, -0.0), n2i(-0.0, -0.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(-0.0, -0.0), n2i(-5.0, -1.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(-0.0, -0.0), n2i(-5.0, 3.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(-0.0, -0.0), n2i(1.0, 3.0), n2i(0.0, 0.0));
        mul_test::<f64>(
            n2i(-0.0, -0.0), n2i(f64::NEG_INFINITY, -1.0),
            n2i(0.0, 0.0)
        );
        mul_test::<f64>(n2i(-0.0, -0.0), n2i(f64::NEG_INFINITY, 3.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(-0.0, -0.0), n2i(-5.0, f64::INFINITY), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(-0.0, -0.0), n2i(1.0, f64::INFINITY), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(-0.0, -0.0), I::ENTIRE, n2i(0.0, 0.0));
        mul_test::<f64>(n2i(1.0, 5.0), n2i(0.0, 0.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(1.0, 5.0), n2i(-0.0, -0.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(1.0, 5.0), n2i(-5.0, -1.0), n2i((-25.0f64).next_down(), (-1.0f64).next_up()));
        mul_test::<f64>(n2i(1.0, 5.0), n2i(-5.0, 3.0), n2i((-25.0f64).next_down(), 15.0f64.next_up()));
        mul_test::<f64>(n2i(1.0, 5.0), n2i(1.0, 3.0), n2i(1.0f64.next_down(), 15.0f64.next_up()));
        mul_test::<f64>(
            n2i(1.0, 5.0), n2i(f64::NEG_INFINITY, -1.0),
            n2i(f64::NEG_INFINITY, (-1.0f64).next_up())
        );
        mul_test::<f64>(
            n2i(1.0, 5.0), n2i(f64::NEG_INFINITY, 3.0),
            n2i(f64::NEG_INFINITY, 15.0f64.next_up())
        );
        mul_test::<f64>(
            n2i(1.0, 5.0), n2i(-5.0, f64::INFINITY),
            n2i((-25.0f64).next_down(), f64::INFINITY)
        );
        mul_test::<f64>(
            n2i(1.0, 5.0), n2i(1.0, f64::INFINITY),
            n2i(1.0f64.next_down(), f64::INFINITY)
        );
        mul_test::<f64>(n2i(1.0, 5.0), I::ENTIRE, I::ENTIRE);
        mul_test::<f64>(n2i(-1.0, 5.0), n2i(0.0, 0.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(-1.0, 5.0), n2i(-0.0, -0.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(-1.0, 5.0), n2i(-5.0, -1.0), n2i((-25.0f64).next_down(), 5.0f64.next_up()));
        //min max
        mul_test::<f64>(n2i(-1.0, 5.0), n2i(-5.0, 3.0), n2i((-25.0f64).next_down(), 15.0f64.next_up()));
        mul_test::<f64>(n2i(-10.0, 2.0), n2i(-5.0, 3.0), n2i((-30.0f64).next_down(), 50.0f64.next_up()));
        mul_test::<f64>(n2i(-1.0, 5.0), n2i(-1.0, 10.0), n2i((-10.0f64).next_down(), 50.0f64.next_up()));
        mul_test::<f64>(n2i(-2.0, 2.0), n2i(-5.0, 3.0), n2i((-10.0f64).next_down(), 10.0f64.next_up()));
        //end min max
        mul_test::<f64>(n2i(-1.0, 5.0), n2i(1.0, 3.0), n2i((-3.0f64).next_down(), 15.0f64.next_up()));
        mul_test::<f64>(n2i(-1.0, 5.0), n2i(f64::NEG_INFINITY, -1.0), I::ENTIRE);
        mul_test::<f64>(n2i(-1.0, 5.0), n2i(f64::NEG_INFINITY, 3.0), I::ENTIRE);
        mul_test::<f64>(n2i(-1.0, 5.0), n2i(-5.0, f64::INFINITY), I::ENTIRE);
        mul_test::<f64>(n2i(-1.0, 5.0), n2i(1.0, f64::INFINITY), I::ENTIRE);
        mul_test::<f64>(n2i(-1.0, 5.0), I::ENTIRE, I::ENTIRE);
        mul_test::<f64>(n2i(-10.0, -5.0), n2i(0.0, 0.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(-10.0, -5.0), n2i(-0.0, -0.0), n2i(0.0, 0.0));
        mul_test::<f64>(n2i(-10.0, -5.0), n2i(-5.0, -1.0), n2i(5.0f64.next_down(), 50.0f64.next_up()));
        mul_test::<f64>(n2i(-10.0, -5.0), n2i(-5.0, 3.0), n2i((-30.0f64).next_down(), 50.0f64.next_up()));
        mul_test::<f64>(n2i(-10.0, -5.0), n2i(1.0, 3.0), n2i((-30.0f64).next_down(), (-5.0f64).next_up()));
        mul_test::<f64>(
            n2i(-10.0, -5.0), n2i(f64::NEG_INFINITY, -1.0),
            n2i(5.0f64.next_down(), f64::INFINITY)
        );
        mul_test::<f64>(
            n2i(-10.0, -5.0), n2i(f64::NEG_INFINITY, 3.0),
            n2i((-30.0f64).next_down(), f64::INFINITY)
        );
        mul_test::<f64>(
            n2i(-10.0, -5.0), n2i(-5.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 50.0f64.next_up())
        );
        mul_test::<f64>(
            n2i(-10.0, -5.0), n2i(1.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, (-5.0f64).next_up())
        );
        mul_test::<f64>(n2i(-10.0, -5.0), I::ENTIRE, I::ENTIRE);
        mul_test::<f64>(
            n2i(0.1, 1.9999999999999964), n2i(-1.9999999999999964, f64::INFINITY),
            n2i(-3.9999999999999862, f64::INFINITY)
        );
        mul_test::<f64>(
            n2i(-0.1, 1.9999999999999964), n2i(-1.9999999999999964, -0.1),
            n2i(-3.9999999999999862, 0.19999999999999968)
        );
        mul_test::<f64>(
            n2i(-0.1, 0.1), n2i(-1.9999999999999964, 0.1),
            n2i(-0.19999999999999968, 0.19999999999999968)
        );
        mul_test::<f64>(
            n2i(-1.9999999999999964, -0.1), n2i(0.1, 1.9999999999999964),
            n2i(-3.9999999999999862, -0.01)
        );
    }
     
    #[test]
    fn test_interval_mul_minimal_fixed() {
        type I = Interval<Fixed>;
        mul_test::<Fixed>(I::EMPTY, I::EMPTY, I::EMPTY);
        mul_test::<Fixed>(i2i(-1, 1, 0), I::EMPTY, I::EMPTY);
        mul_test::<Fixed>(I::EMPTY, i2i(-1, 1, 0), I::EMPTY);
        mul_test::<Fixed>(i2i(0, 0, 0), I::EMPTY, I::EMPTY);
        mul_test::<Fixed>(I::EMPTY, i2i(0, 0, 0), I::EMPTY);
        mul_test::<Fixed>(i2i(-0, -0, 0), I::EMPTY, I::EMPTY);
        mul_test::<Fixed>(I::EMPTY, i2i(-0, -0, 0), I::EMPTY);
        mul_test::<Fixed>(i2i(0, 0, 0), i2i(0, 0, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(0, 0, 0), i2i(-0, -0, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(0, 0, 0), i2i(-5, -1, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(0, 0, 0), i2i(-5, 3, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(0, 0, 0), i2i(1, 3, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(-0, -0, 0), i2i(0, 0, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(-0, -0, 0), i2i(-0, -0, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(-0, -0, 0), i2i(-5, -1, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(-0, -0, 0), i2i(-5, 3, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(-0, -0, 0), i2i(1, 3, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(1, 5, 0), i2i(0, 0, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(1, 5, 0), i2i(-0, -0, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(1, 5, 0), i2i(-5, -1, 0), i2i(-25, -1, 0));
        mul_test::<Fixed>(i2i(1, 5, 0), i2i(-5, 3, 0), i2i(-25, 15, 0));
        mul_test::<Fixed>(i2i(1, 5, 0), i2i(1, 3, 0), i2i(1, 15, 0));
        mul_test::<Fixed>(i2i(-1, 5, 0), i2i(0, 0, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(-1, 5, 0), i2i(-0, -0, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(-1, 5, 0), i2i(-5, -1, 0), i2i(-25, 5, 0));
        //min max
        mul_test::<Fixed>(i2i(-1, 5, 0), i2i(-5, 3, 0), i2i(-25, 15, 0));
        mul_test::<Fixed>(i2i(-10, 2, 0), i2i(-5, 3, 0), i2i(-30, 50, 0));
        mul_test::<Fixed>(i2i(-1, 5, 0), i2i(-1, 10, 0), i2i(-10, 50, 0));
        mul_test::<Fixed>(i2i(-2, 2, 0), i2i(-5, 3, 0), i2i(-10, 10, 0));
        //end min max
        mul_test::<Fixed>(i2i(-1, 5, 0), i2i(1, 3, 0), i2i(-3, 15, 0));
        mul_test::<Fixed>(i2i(-10, -5, 0), i2i(0, 0, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(-10, -5, 0), i2i(-0, -0, 0), i2i(0, 0, 0));
        mul_test::<Fixed>(i2i(-10, -5, 0), i2i(-5, -1, 0), i2i(5, 50, 0));
        mul_test::<Fixed>(i2i(-10, -5, 0), i2i(-5, 3, 0), i2i(-30, 50, 0));
        mul_test::<Fixed>(i2i(-10, -5, 0), i2i(1, 3, 0), i2i(-30, -5, 0));
        mul_test::<Fixed>(
            i2i(5, 5, 4), i2i(6, 6, 4),
            i2i(1, 2, 4)
        );
        mul_test::<Fixed>(
            i2i(-5, -5, 4), i2i(6, 6, 4),
            i2i(-2, -1, 4)
        );
    }
        
    #[test]
    fn test_interval_div_minimal_float() {
        type I = Interval<f64>;
        div_bounds_test::<f64>(I::EMPTY, I::EMPTY, I::EMPTY);
        div_bounds_test::<f64>(n2i(-1.0, 1.0), I::EMPTY, I::EMPTY);
        div_bounds_test::<f64>(I::EMPTY, n2i(-1.0, 1.0), I::EMPTY);
        div_bounds_test::<f64>(I::EMPTY, n2i(0.1, 1.0), I::EMPTY);
        div_bounds_test::<f64>(I::EMPTY, n2i(-1.0, -0.1), I::EMPTY);
        div_bounds_test::<f64>(I::EMPTY, I::ENTIRE, I::EMPTY);
        div_bounds_test::<f64>(I::ENTIRE, I::EMPTY, I::EMPTY);
        div_bounds_test::<f64>(n2i(0.0, 0.0), I::EMPTY, I::EMPTY);
        div_bounds_test::<f64>(I::EMPTY, n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(-0.0, -0.0), I::EMPTY, I::EMPTY);
        div_bounds_test::<f64>(I::EMPTY, n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(I::ENTIRE, n2i(-5.0, -3.0), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, n2i(3.0, 5.0), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, n2i(f64::NEG_INFINITY, -3.0), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, n2i(3.0, f64::INFINITY), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(I::ENTIRE, n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(I::ENTIRE, n2i(-3.0, 0.0), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, n2i(-3.0, -0.0), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, n2i(0.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, n2i(f64::NEG_INFINITY, 0.0), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, n2i(-0.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, n2i(f64::NEG_INFINITY, -0.0), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, n2i(f64::NEG_INFINITY, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, n2i(-3.0, f64::INFINITY), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, n2i(0.0, f64::INFINITY), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, n2i(-0.0, f64::INFINITY), I::ENTIRE);
        div_bounds_test::<f64>(I::ENTIRE, I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, -15.0), n2i(-5.0, -3.0), n2i(3.0f64.next_down(), 10.0f64.next_up()));
        div_bounds_test::<f64>(n2i(-30.0, -15.0), n2i(3.0, 5.0), n2i((-10.0f64).next_down(), (-3.0f64).next_up()));
        div_bounds_test::<f64>(
            n2i(-30.0, -15.0), n2i(f64::NEG_INFINITY, -3.0),
            n2i(0.0, 10.0f64.next_up())
        );
        div_bounds_test::<f64>(n2i(-30.0, -15.0), n2i(3.0, f64::INFINITY), n2i((-10.0f64).next_down(), 0.0));
        div_bounds_test::<f64>(n2i(-30.0, -15.0), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(-30.0, -15.0), n2i(-3.0, 0.0), n2i(5.0f64.next_down(), f64::INFINITY));
        div_bounds_test::<f64>(n2i(-30.0, -15.0), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(-30.0, -15.0), n2i(-3.0, -0.0), n2i(5.0f64.next_down(), f64::INFINITY));
        div_bounds_test::<f64>(n2i(-30.0, -15.0), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(-30.0, -15.0), n2i(0.0, 3.0),
            n2i(f64::NEG_INFINITY, (-5.0f64).next_up())
        );
        div_bounds_test::<f64>(
            n2i(-30.0, -15.0), n2i(f64::NEG_INFINITY, 0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(-30.0, -15.0), n2i(-0.0, 3.0),
            n2i(f64::NEG_INFINITY, (-5.0f64).next_up())
        );
        div_bounds_test::<f64>(
            n2i(-30.0, -15.0), n2i(f64::NEG_INFINITY, -0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(-30.0, -15.0), n2i(f64::NEG_INFINITY, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, -15.0), n2i(-3.0, f64::INFINITY), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(-30.0, -15.0), n2i(0.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(-30.0, -15.0), n2i(-0.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(-30.0, -15.0), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(-5.0, -3.0), n2i((-5.0f64).next_down(), 10.0f64.next_up()));
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(3.0, 5.0), n2i((-10.0f64).next_down(), 5.0f64.next_up()));
        div_bounds_test::<f64>(
            n2i(-30.0, 15.0), n2i(f64::NEG_INFINITY, -3.0),
            n2i((-5.0f64).next_down(), 10.0f64.next_up())
        );
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(3.0, f64::INFINITY), n2i((-10.0f64).next_down(), 5.0f64.next_up()));
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(-3.0, 0.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(-3.0, -0.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(0.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(f64::NEG_INFINITY, 0.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(-0.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(f64::NEG_INFINITY, -0.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(f64::NEG_INFINITY, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(-3.0, f64::INFINITY), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(0.0, f64::INFINITY), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 15.0), n2i(-0.0, f64::INFINITY), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 15.0), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(n2i(15.0, 30.0), n2i(-5.0, -3.0), n2i((-10.0f64).next_down(), (-3.0f64).next_up()));
        div_bounds_test::<f64>(n2i(15.0, 30.0), n2i(3.0, 5.0), n2i(3.0f64.next_down(), 10.0f64.next_up()));
        div_bounds_test::<f64>(
            n2i(15.0, 30.0), n2i(f64::NEG_INFINITY, -3.0),
            n2i((-10.0f64).next_down(), 0.0)
        );
        div_bounds_test::<f64>(n2i(15.0, 30.0), n2i(3.0, f64::INFINITY), n2i(0.0, 10.0f64.next_up()));
        div_bounds_test::<f64>(n2i(15.0, 30.0), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(15.0, 30.0), n2i(-3.0, 0.0),
            n2i(f64::NEG_INFINITY, (-5.0f64).next_up())
        );
        div_bounds_test::<f64>(n2i(15.0, 30.0), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(15.0, 30.0), n2i(-3.0, -0.0),
            n2i(f64::NEG_INFINITY, (-5.0f64).next_up())
        );
        div_bounds_test::<f64>(n2i(15.0, 30.0), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(15.0, 30.0), n2i(0.0, 3.0), n2i(5.0f64.next_down(), f64::INFINITY));
        div_bounds_test::<f64>(
            n2i(15.0, 30.0), n2i(f64::NEG_INFINITY, 0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(15.0, 30.0), n2i(-0.0, 3.0), n2i(5.0f64.next_down(), f64::INFINITY));
        div_bounds_test::<f64>(
            n2i(15.0, 30.0), n2i(f64::NEG_INFINITY, -0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(15.0, 30.0), n2i(f64::NEG_INFINITY, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(15.0, 30.0), n2i(-3.0, f64::INFINITY), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(15.0, 30.0), n2i(0.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(15.0, 30.0), n2i(-0.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(15.0, 30.0), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(-5.0, -3.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(3.0, 5.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(f64::NEG_INFINITY, -3.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(3.0, f64::INFINITY), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(-3.0, 0.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(-3.0, -0.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(-3.0, 3.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(0.0, 3.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(f64::NEG_INFINITY, 0.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(-0.0, 3.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(f64::NEG_INFINITY, -0.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(f64::NEG_INFINITY, 3.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(-3.0, f64::INFINITY), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(0.0, f64::INFINITY), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), n2i(-0.0, f64::INFINITY), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 0.0), I::ENTIRE, n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(-5.0, -3.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(3.0, 5.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(
            n2i(-0.0, -0.0), n2i(f64::NEG_INFINITY, -3.0),
            n2i(0.0, 0.0)
        );
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(3.0, f64::INFINITY), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(-3.0, 0.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(-3.0, -0.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(-3.0, 3.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(0.0, 3.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(f64::NEG_INFINITY, 0.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(-0.0, 3.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(
            n2i(-0.0, -0.0), n2i(f64::NEG_INFINITY, -0.0),
            n2i(0.0, 0.0)
        );
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(f64::NEG_INFINITY, 3.0), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(-3.0, f64::INFINITY), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(0.0, f64::INFINITY), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(-0.0, -0.0), n2i(-0.0, f64::INFINITY), n2i(0.0, 0.0));
        div_bounds_test::<f64>(n2i(-0.0, -0.0), I::ENTIRE, n2i(0.0, 0.0));
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -15.0), n2i(-5.0, -3.0),
            n2i(3.0f64.next_down(), f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -15.0), n2i(3.0, 5.0),
            n2i(f64::NEG_INFINITY, (-3.0f64).next_up())
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -15.0), n2i(f64::NEG_INFINITY, -3.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -15.0), n2i(3.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, -15.0), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -15.0), n2i(-3.0, 0.0),
            n2i(5.0f64.next_down(), f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, -15.0), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -15.0), n2i(-3.0, -0.0),
            n2i(5.0f64.next_down(), f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, -15.0), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -15.0), n2i(0.0, 3.0),
            n2i(f64::NEG_INFINITY, (-5.0f64).next_up())
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -15.0), n2i(f64::NEG_INFINITY, 0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -15.0), n2i(-0.0, 3.0),
            n2i(f64::NEG_INFINITY, (-5.0f64).next_up())
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -15.0), n2i(f64::NEG_INFINITY, -0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -15.0), n2i(f64::NEG_INFINITY, 3.0),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -15.0), n2i(-3.0, f64::INFINITY),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -15.0), n2i(0.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -15.0), n2i(-0.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, -15.0), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 15.0), n2i(-5.0, -3.0),
            n2i((-5.0f64).next_down(), f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 15.0), n2i(3.0, 5.0),
            n2i(f64::NEG_INFINITY, 5.0f64.next_up())
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 15.0), n2i(f64::NEG_INFINITY, -3.0),
            n2i((-5.0f64).next_down(), f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 15.0), n2i(3.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 5.0f64.next_up())
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, 15.0), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, 15.0), n2i(-3.0, 0.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, 15.0), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, 15.0), n2i(-3.0, -0.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, 15.0), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, 15.0), n2i(0.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 15.0), n2i(f64::NEG_INFINITY, 0.0),
            I::ENTIRE
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, 15.0), n2i(-0.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 15.0), n2i(f64::NEG_INFINITY, -0.0),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 15.0), n2i(f64::NEG_INFINITY, 3.0),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 15.0), n2i(-3.0, f64::INFINITY),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 15.0), n2i(0.0, f64::INFINITY),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 15.0), n2i(-0.0, f64::INFINITY),
            I::ENTIRE
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, 15.0), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(-15.0, f64::INFINITY), n2i(-5.0, -3.0),
            n2i(f64::NEG_INFINITY, 5.0f64.next_up())
        );
        div_bounds_test::<f64>(
            n2i(-15.0, f64::INFINITY), n2i(3.0, 5.0),
            n2i((-5.0f64).next_down(), f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(-15.0, f64::INFINITY), n2i(f64::NEG_INFINITY, -3.0),
            n2i(f64::NEG_INFINITY, 5.0f64.next_up())
        );
        div_bounds_test::<f64>(
            n2i(-15.0, f64::INFINITY), n2i(3.0, f64::INFINITY),
            n2i((-5.0f64).next_down(), f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(-15.0, f64::INFINITY), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(-15.0, f64::INFINITY), n2i(-3.0, 0.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-15.0, f64::INFINITY), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(-15.0, f64::INFINITY), n2i(-3.0, -0.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-15.0, f64::INFINITY), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-15.0, f64::INFINITY), n2i(0.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(-15.0, f64::INFINITY), n2i(f64::NEG_INFINITY, 0.0),
            I::ENTIRE
        );
        div_bounds_test::<f64>(n2i(-15.0, f64::INFINITY), n2i(-0.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(-15.0, f64::INFINITY), n2i(f64::NEG_INFINITY, -0.0),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(-15.0, f64::INFINITY), n2i(f64::NEG_INFINITY, 3.0),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(-15.0, f64::INFINITY), n2i(-3.0, f64::INFINITY),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(-15.0, f64::INFINITY), n2i(0.0, f64::INFINITY),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(-15.0, f64::INFINITY), n2i(-0.0, f64::INFINITY),
            I::ENTIRE
        );
        div_bounds_test::<f64>(n2i(-15.0, f64::INFINITY), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(15.0, f64::INFINITY), n2i(-5.0, -3.0),
            n2i(f64::NEG_INFINITY, (-3.0f64).next_up())
        );
        div_bounds_test::<f64>(
            n2i(15.0, f64::INFINITY), n2i(3.0, 5.0),
            n2i(3.0f64.next_down(), f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(15.0, f64::INFINITY), n2i(f64::NEG_INFINITY, -3.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(15.0, f64::INFINITY), n2i(3.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(15.0, f64::INFINITY), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(15.0, f64::INFINITY), n2i(-3.0, 0.0),
            n2i(f64::NEG_INFINITY, (-5.0f64).next_up())
        );
        div_bounds_test::<f64>(n2i(15.0, f64::INFINITY), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(15.0, f64::INFINITY), n2i(-3.0, -0.0),
            n2i(f64::NEG_INFINITY, (-5.0f64).next_up())
        );
        div_bounds_test::<f64>(n2i(15.0, f64::INFINITY), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(15.0, f64::INFINITY), n2i(0.0, 3.0),
            n2i(5.0f64.next_down(), f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(15.0, f64::INFINITY), n2i(f64::NEG_INFINITY, 0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(15.0, f64::INFINITY), n2i(-0.0, 3.0),
            n2i(5.0f64.next_down(), f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(15.0, f64::INFINITY), n2i(f64::NEG_INFINITY, -0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(15.0, f64::INFINITY), n2i(f64::NEG_INFINITY, 3.0),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(15.0, f64::INFINITY), n2i(-3.0, f64::INFINITY),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(15.0, f64::INFINITY), n2i(0.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(15.0, f64::INFINITY), n2i(-0.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(15.0, f64::INFINITY), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 0.0), n2i(-5.0, -3.0), n2i(0.0, 10.0f64.next_up()));
        div_bounds_test::<f64>(n2i(-30.0, 0.0), n2i(3.0, 5.0), n2i((-10.0f64).next_down(), 0.0));
        div_bounds_test::<f64>(
            n2i(-30.0, 0.0), n2i(f64::NEG_INFINITY, -3.0),
            n2i(0.0, 10.0f64.next_up())
        );
        div_bounds_test::<f64>(n2i(-30.0, 0.0), n2i(3.0, f64::INFINITY), n2i((-10.0f64).next_down(), 0.0));
        div_bounds_test::<f64>(n2i(-30.0, 0.0), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(-30.0, 0.0), n2i(-3.0, 0.0), n2i(0.0, f64::INFINITY));
        div_bounds_test::<f64>(n2i(-30.0, 0.0), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(-30.0, 0.0), n2i(-3.0, -0.0), n2i(0.0, f64::INFINITY));
        div_bounds_test::<f64>(n2i(-30.0, 0.0), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 0.0), n2i(0.0, 3.0), n2i(f64::NEG_INFINITY, 0.0));
        div_bounds_test::<f64>(
            n2i(-30.0, 0.0), n2i(f64::NEG_INFINITY, 0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(-30.0, 0.0), n2i(-0.0, 3.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(-30.0, 0.0), n2i(f64::NEG_INFINITY, -0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(-30.0, 0.0), n2i(f64::NEG_INFINITY, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, 0.0), n2i(-3.0, f64::INFINITY), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(-30.0, 0.0), n2i(0.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(-30.0, 0.0), n2i(-0.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(-30.0, 0.0), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, -0.0), n2i(-5.0, -3.0), n2i(0.0, 10.0f64.next_up()));
        div_bounds_test::<f64>(n2i(-30.0, -0.0), n2i(3.0, 5.0), n2i((-10.0f64).next_down(), 0.0));
        div_bounds_test::<f64>(
            n2i(-30.0, -0.0), n2i(f64::NEG_INFINITY, -3.0),
            n2i(0.0, 10.0f64.next_up())
        );
        div_bounds_test::<f64>(n2i(-30.0, -0.0), n2i(3.0, f64::INFINITY), n2i((-10.0f64).next_down(), 0.0));
        div_bounds_test::<f64>(n2i(-30.0, -0.0), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(-30.0, -0.0), n2i(-3.0, 0.0), n2i(0.0, f64::INFINITY));
        div_bounds_test::<f64>(n2i(-30.0, -0.0), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(-30.0, -0.0), n2i(-3.0, -0.0), n2i(0.0, f64::INFINITY));
        div_bounds_test::<f64>(n2i(-30.0, -0.0), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(-30.0, -0.0), n2i(0.0, 3.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(-30.0, -0.0), n2i(f64::NEG_INFINITY, 0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(-30.0, -0.0), n2i(-0.0, 3.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(-30.0, -0.0), n2i(f64::NEG_INFINITY, -0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(-30.0, -0.0), n2i(f64::NEG_INFINITY, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-30.0, -0.0), n2i(-3.0, f64::INFINITY), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(-30.0, -0.0), n2i(0.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(-30.0, -0.0), n2i(-0.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(-30.0, -0.0), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(n2i(0.0, 30.0), n2i(-5.0, -3.0), n2i((-10.0f64).next_down(), 0.0));
        div_bounds_test::<f64>(n2i(0.0, 30.0), n2i(3.0, 5.0), n2i(0.0, 10.0f64.next_up()));
        div_bounds_test::<f64>(
            n2i(0.0, 30.0), n2i(f64::NEG_INFINITY, -3.0),
            n2i((-10.0f64).next_down(), 0.0)
        );
        div_bounds_test::<f64>(n2i(0.0, 30.0), n2i(3.0, f64::INFINITY), n2i(0.0, 10.0f64.next_up()));
        div_bounds_test::<f64>(n2i(0.0, 30.0), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(n2i(0.0, 30.0), n2i(-3.0, 0.0), n2i(f64::NEG_INFINITY, 0.0));
        div_bounds_test::<f64>(n2i(0.0, 30.0), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(0.0, 30.0), n2i(-3.0, -0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(0.0, 30.0), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(0.0, 30.0), n2i(0.0, 3.0), n2i(0.0, f64::INFINITY));
        div_bounds_test::<f64>(
            n2i(0.0, 30.0), n2i(f64::NEG_INFINITY, 0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(0.0, 30.0), n2i(-0.0, 3.0), n2i(0.0, f64::INFINITY));
        div_bounds_test::<f64>(
            n2i(0.0, 30.0), n2i(f64::NEG_INFINITY, -0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(0.0, 30.0), n2i(f64::NEG_INFINITY, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(0.0, 30.0), n2i(-3.0, f64::INFINITY), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(0.0, 30.0), n2i(0.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(0.0, 30.0), n2i(-0.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(0.0, 30.0), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(n2i(-0.0, 30.0), n2i(-5.0, -3.0), n2i((-10.0f64).next_down(), 0.0));
        div_bounds_test::<f64>(n2i(-0.0, 30.0), n2i(3.0, 5.0), n2i(0.0, 10.0f64.next_up()));
        div_bounds_test::<f64>(
            n2i(-0.0, 30.0), n2i(f64::NEG_INFINITY, -3.0),
            n2i((-10.0f64).next_down(), 0.0)
        );
        div_bounds_test::<f64>(n2i(-0.0, 30.0), n2i(3.0, f64::INFINITY), n2i(0.0, 10.0f64.next_up()));
        div_bounds_test::<f64>(n2i(-0.0, 30.0), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(-0.0, 30.0), n2i(-3.0, 0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(-0.0, 30.0), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(-0.0, 30.0), n2i(-3.0, -0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(-0.0, 30.0), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-0.0, 30.0), n2i(0.0, 3.0), n2i(0.0, f64::INFINITY));
        div_bounds_test::<f64>(
            n2i(-0.0, 30.0), n2i(f64::NEG_INFINITY, 0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(-0.0, 30.0), n2i(-0.0, 3.0), n2i(0.0, f64::INFINITY));
        div_bounds_test::<f64>(
            n2i(-0.0, 30.0), n2i(f64::NEG_INFINITY, -0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(-0.0, 30.0), n2i(f64::NEG_INFINITY, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(n2i(-0.0, 30.0), n2i(-3.0, f64::INFINITY), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(-0.0, 30.0), n2i(0.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(-0.0, 30.0), n2i(-0.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(-0.0, 30.0), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 0.0), n2i(-5.0, -3.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 0.0), n2i(3.0, 5.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 0.0), n2i(f64::NEG_INFINITY, -3.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 0.0), n2i(3.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, 0.0), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 0.0), n2i(-3.0, 0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, 0.0), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 0.0), n2i(-3.0, -0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, 0.0), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 0.0), n2i(0.0, 3.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 0.0), n2i(f64::NEG_INFINITY, 0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 0.0), n2i(-0.0, 3.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 0.0), n2i(f64::NEG_INFINITY, -0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 0.0), n2i(f64::NEG_INFINITY, 3.0),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 0.0), n2i(-3.0, f64::INFINITY),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 0.0), n2i(0.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, 0.0), n2i(-0.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, 0.0), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -0.0), n2i(-5.0, -3.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -0.0), n2i(3.0, 5.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -0.0), n2i(f64::NEG_INFINITY, -3.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -0.0), n2i(3.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, -0.0), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -0.0), n2i(-3.0, 0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, -0.0), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -0.0), n2i(-3.0, -0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, -0.0), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -0.0), n2i(0.0, 3.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -0.0), n2i(f64::NEG_INFINITY, 0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -0.0), n2i(-0.0, 3.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -0.0), n2i(f64::NEG_INFINITY, -0.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -0.0), n2i(f64::NEG_INFINITY, 3.0),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -0.0), n2i(-3.0, f64::INFINITY),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -0.0), n2i(0.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -0.0), n2i(-0.0, f64::INFINITY),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(f64::NEG_INFINITY, -0.0), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(0.0, f64::INFINITY), n2i(-5.0, -3.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(0.0, f64::INFINITY), n2i(3.0, 5.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(0.0, f64::INFINITY), n2i(f64::NEG_INFINITY, -3.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(0.0, f64::INFINITY), n2i(3.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(0.0, f64::INFINITY), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(0.0, f64::INFINITY), n2i(-3.0, 0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(0.0, f64::INFINITY), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(0.0, f64::INFINITY), n2i(-3.0, -0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(0.0, f64::INFINITY), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(0.0, f64::INFINITY), n2i(0.0, 3.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(0.0, f64::INFINITY), n2i(f64::NEG_INFINITY, 0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(0.0, f64::INFINITY), n2i(-0.0, 3.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(0.0, f64::INFINITY), n2i(f64::NEG_INFINITY, -0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(0.0, f64::INFINITY), n2i(f64::NEG_INFINITY, 3.0),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(0.0, f64::INFINITY), n2i(-3.0, f64::INFINITY),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(0.0, f64::INFINITY), n2i(0.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(0.0, f64::INFINITY), n2i(-0.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(0.0, f64::INFINITY), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(-0.0, f64::INFINITY), n2i(-5.0, -3.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(-0.0, f64::INFINITY), n2i(3.0, 5.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(-0.0, f64::INFINITY), n2i(f64::NEG_INFINITY, -3.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(-0.0, f64::INFINITY), n2i(3.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(-0.0, f64::INFINITY), n2i(0.0, 0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(-0.0, f64::INFINITY), n2i(-3.0, 0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(-0.0, f64::INFINITY), n2i(-0.0, -0.0), I::EMPTY);
        div_bounds_test::<f64>(
            n2i(-0.0, f64::INFINITY), n2i(-3.0, -0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(n2i(-0.0, f64::INFINITY), n2i(-3.0, 3.0), I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(-0.0, f64::INFINITY), n2i(0.0, 3.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(-0.0, f64::INFINITY), n2i(f64::NEG_INFINITY, 0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(-0.0, f64::INFINITY), n2i(-0.0, 3.0),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(-0.0, f64::INFINITY), n2i(f64::NEG_INFINITY, -0.0),
            n2i(f64::NEG_INFINITY, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(-0.0, f64::INFINITY), n2i(f64::NEG_INFINITY, 3.0),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(-0.0, f64::INFINITY), n2i(-3.0, f64::INFINITY),
            I::ENTIRE
        );
        div_bounds_test::<f64>(
            n2i(-0.0, f64::INFINITY), n2i(0.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(
            n2i(-0.0, f64::INFINITY), n2i(-0.0, f64::INFINITY),
            n2i(0.0, f64::INFINITY)
        );
        div_bounds_test::<f64>(n2i(-0.0, f64::INFINITY), I::ENTIRE, I::ENTIRE);
        div_bounds_test::<f64>(
            n2i(-2.0, -1.0), n2i(-10.0, -3.0),
            n2i(0.09999999999999999, 0.6666666666666667)
        );
        div_bounds_test::<f64>(
            n2i(-2.0, -1.0), n2i(0.0, 10.0),
            n2i(f64::NEG_INFINITY, -0.09999999999999999)
        );
        div_bounds_test::<f64>(
            n2i(-2.0, -1.0), n2i(-0.0, 10.0),
            n2i(f64::NEG_INFINITY, -0.09999999999999999)
        );
        div_bounds_test::<f64>(n2i(-1.0, 2.0), n2i(10.0, f64::INFINITY), n2i((-0.1f64).next_down(), 0.2f64.next_up()));
        div_bounds_test::<f64>(
            n2i(1.0, 3.0), n2i(f64::NEG_INFINITY, -10.0),
            n2i(-0.30000000000000004, 0.0)
        );
        div_bounds_test::<f64>(
            n2i(f64::NEG_INFINITY, -1.0), n2i(1.0, 3.0),
            n2i(f64::NEG_INFINITY, (-0.3333333333333333f64).next_up())
        );
    }

    #[test]
    fn test_interval_sqrt_minimal_float() {
        type I = Interval<f64>;
        sqrt_bounds_test::<f64>(I::EMPTY, I::EMPTY);
        sqrt_bounds_test::<f64>(I::ENTIRE, n2i(0.0, f64::INFINITY));
        sqrt_bounds_test::<f64>(n2i(f64::NEG_INFINITY, -5e-324), I::EMPTY);
        sqrt_bounds_test::<f64>(n2i(-1.0, 1.0), n2i(0.0, 1.0f64.next_up()));
        sqrt_bounds_test::<f64>(n2i(0.0, 1.0), n2i(0.0, 1.0f64.next_up()));
        sqrt_bounds_test::<f64>(n2i(-0.0, 1.0), n2i(0.0, 1.0f64.next_up()));
        sqrt_bounds_test::<f64>(n2i(-5.0, 25.0), n2i(0.0, 5.0f64.next_up()));
        sqrt_bounds_test::<f64>(n2i(0.0, 25.0), n2i(0.0, 5.0f64.next_up()));
        sqrt_bounds_test::<f64>(n2i(-0.0, 25.0), n2i(0.0, 5.0f64.next_up()));
        sqrt_bounds_test::<f64>(n2i(-5.0, f64::INFINITY), n2i(0.0, f64::INFINITY));
        sqrt_bounds_test::<f64>(
            n2i(0.1, 0.1),
            n2i(0.31622776601683794f64.next_down(), 0.316227766016838)
        );
        sqrt_bounds_test::<f64>(
            n2i(-1.9999999999999964, 0.1),
            n2i(0.0, 0.316227766016838)
        );
        sqrt_bounds_test::<f64>(
            n2i(0.1, 1.9999999999999964),
            n2i(0.31622776601683794f64.next_down(), 1.4142135623730938f64.next_up())
        );
    }
 
    #[test]
    fn test_interval_div_minimal_fixed() {
        type I = Interval<Fixed>;
        // I'm lazy, so just test rounding
        div_bounds_test::<Fixed>(
            i2i(16, 16, 4), i2i(48, 48, 4),
            i2i(5, 6, 4)
        );
        div_bounds_test::<Fixed>(
            i2i(-16, -16, 4), i2i(48, 48, 4),
            i2i(-6, -5, 4)
        );
        div_bounds_test::<Fixed>(
            i2i(16, 16, 4), i2i(-48, -48, 4),
            i2i(-6, -5, 4)
        );
        div_bounds_test::<Fixed>(
            i2i(-16, -16, 4), i2i(-48, -48, 4),
            i2i(5, 6, 4)
        );
    }

    #[test]
    fn test_interval_sqrt_minimal_fixed() {
        type I = Interval<Fixed>;
        sqrt_bounds_test::<Fixed>(I::EMPTY, I::EMPTY);
        sqrt_bounds_test::<Fixed>(i2i(-1, 1, 0), i2i(0, 1, 0));
        sqrt_bounds_test::<Fixed>(i2i(0, 1, 0), i2i(0, 1, 0));
        sqrt_bounds_test::<Fixed>(i2i(-0, 1, 0), i2i(0, 1, 0));
        sqrt_bounds_test::<Fixed>(i2i(-5, 25, 0), i2i(0, 5, 0));
        sqrt_bounds_test::<Fixed>(i2i(0, 25, 0), i2i(0, 5, 0));
        sqrt_bounds_test::<Fixed>(i2i(-0, 25, 0), i2i(0, 5, 0));
        sqrt_bounds_test::<Fixed>(
            i2i(32, 32, 4),
            i2i(22, 23, 4)
        );
    }

    #[test]
    fn test_interval_neg_minimal_float() {
        type I = Interval<f64>;
        neg_test::<f64>(n2i(1.0, 2.0), n2i(-2.0, -1.0));
        neg_test::<f64>(I::EMPTY, I::EMPTY);
        neg_test::<f64>(I::ENTIRE, I::ENTIRE);
        neg_test::<f64>(n2i(1.0, f64::INFINITY), n2i(f64::NEG_INFINITY, -1.0));
        neg_test::<f64>(n2i(f64::NEG_INFINITY, 1.0), n2i(-1.0, f64::INFINITY));
        neg_test::<f64>(n2i(0.0, 2.0), n2i(-2.0, 0.0));
        neg_test::<f64>(n2i(-0.0, 2.0), n2i(-2.0, 0.0));
        neg_test::<f64>(n2i(-2.0, 0.0), n2i(0.0, 2.0));
        neg_test::<f64>(n2i(-2.0, -0.0), n2i(0.0, 2.0));
        neg_test::<f64>(n2i(0.0, -0.0), n2i(0.0, 0.0));
        neg_test::<f64>(n2i(-0.0, -0.0), n2i(0.0, 0.0));
    }

    #[test]
    fn test_interval_cmp_minimal_float() {
        type I = Interval<f64>;
        cmp_zero_test::<f64>(I::EMPTY, None);
        cmp_zero_test::<f64>(n2i(-2.0, -1.0), Some(Ordering::Less));
        cmp_zero_test::<f64>(n2i(-2.0, 0.0), None);
        cmp_zero_test::<f64>(n2i(-2.0, 2.0), None);
        cmp_zero_test::<f64>(n2i(0.0, 2.0), None);
        cmp_zero_test::<f64>(n2i(1.0, 2.0), Some(Ordering::Greater));
        cmp_zero_test::<f64>(n2i(0.0, 0.0), Some(Ordering::Equal));
    }
}