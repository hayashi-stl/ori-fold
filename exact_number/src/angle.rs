use std::{mem, ops::{Add, AddAssign, Neg, Sub, SubAssign}};

use malachite::base::num::arithmetic::traits::{CheckedDiv, NegAssign};
use num::Signed;

use crate::BasedExpr;

/// A fold angle represented exactly
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Angle {
    /// Specifically, ⌊angle/π + 1/2⌋ - 1/2
    turn_value: i32,
    tan: Option<BasedExpr>
}

macro_rules! angle {
    ($turn:literal, -inf) => { $crate::angle::Angle::new($turn, None) };
    ($turn:literal, $($tt:tt)*) => { $crate::angle::Angle::new($turn, Some($crate::based_expr!($($tt)*))) };
}

impl Angle {
    pub const DEG_360: Angle = Angle::new(2, Some(BasedExpr::BASELESS_ZERO));
    pub const DEG_180: Angle = Angle::new(1, Some(BasedExpr::BASELESS_ZERO));
    pub const DEG_90: Angle = Angle::new(1, None);
    pub const DEG_45: Angle = Angle::new(0, Some(BasedExpr::BASELESS_ONE));
    pub const ZERO: Angle = Angle::new(0, Some(BasedExpr::BASELESS_ZERO));

    /// The `turn value` is ⌊angle/π + 1/2⌋ - 1/2.
    /// The `tan` is the tangent of the angle, or `None` if the tangent is undefined.
    pub const fn new(turn_value: i32, tan: Option<BasedExpr>) -> Self {
        Self { turn_value, tan }
    }

    pub fn atan2(y: BasedExpr, x: BasedExpr) -> Self {
        let x_sign = x.cmp_zero();
        let y_sign = y.cmp_zero();
        let turn_value = if x_sign.is_lt() && y_sign.is_le() {
            -1
        } else if x_sign.is_gt() || y_sign.is_le() {
            0
        } else {
            1
        };
        let tan = if x_sign.is_eq() && y_sign.is_eq() { Some(x.into_zero()) } else { y.checked_div(x) };
        Self::new(turn_value, tan)
    }

    /// Returns the tangent of the angle. If undefined, returns `None`.
    pub fn tan(&self) -> Option<BasedExpr> {
        self.tan.clone()
    }

    /// Converts this into its turn value $\lfloor angle/π + 1/2\rfloor - 1/2$ and its tangent.
    pub fn into_turn_value_tan(self) -> (i32, Option<BasedExpr>) {
        (self.turn_value, self.tan)
    }
}

impl Default for Angle {
    fn default() -> Self { Self::ZERO }
}

impl NegAssign for Angle {
    fn neg_assign(&mut self) {
        self.turn_value.neg_assign();
        if let Some(tan) = &mut self.tan {
            tan.neg_assign();
        } else {
            self.turn_value += 1;
        }
    }
}

impl Neg for Angle {
    type Output = Angle;

    fn neg(mut self) -> Self::Output {
        self.neg_assign();
        self
    }
}

impl Neg for &Angle {
    type Output = Angle;

    fn neg(self) -> Self::Output {
        -self.clone()
    }
}

impl AddAssign<Angle> for Angle {
    fn add_assign(&mut self, mut rhs: Angle) {
        if self.tan == None { mem::swap(self, &mut rhs); }
        match (self.tan.as_mut(), rhs.tan) {
            (None, None) => {
                self.turn_value += rhs.turn_value - 1;
                self.tan = Some(BasedExpr::BASELESS_ZERO);
            },
            (None, Some(_)) => unreachable!(),
            (Some(tan), None) => {
                self.turn_value += rhs.turn_value;
                if tan.cmp_zero().is_lt() { self.turn_value -= 1; }
                self.tan = (-BasedExpr::BASELESS_ONE).checked_div(mem::take(tan));
            },
            (Some(tan_a), Some(tan_b)) => {
                self.turn_value += rhs.turn_value;
                let cmp_a = tan_a.cmp_zero();
                let cmp_b = tan_b.cmp_zero();
                let maybe_overflow = cmp_a.is_gt() && cmp_b.is_gt();
                let maybe_underflow = cmp_a.is_lt() && cmp_b.is_lt();
                self.tan = (&*tan_a + &tan_b).checked_div(BasedExpr::BASELESS_ONE - mem::take(tan_a) * tan_b);
                if maybe_overflow  && self.tan < Some(BasedExpr::BASELESS_ZERO) { self.turn_value += 1; }
                if maybe_underflow && self.tan > Some(BasedExpr::BASELESS_ZERO) { self.turn_value -= 1; }
            },
        }
    }
}

impl AddAssign<&Angle> for Angle {
    fn add_assign(&mut self, mut rhs: &Angle) {
        match (self.tan.as_mut(), rhs.tan.as_ref()) {
            (None, None) => {
                self.turn_value += rhs.turn_value - 1;
                self.tan = Some(BasedExpr::BASELESS_ZERO);
            },
            (None, Some(tan)) => {
                self.turn_value += rhs.turn_value;
                if tan.cmp_zero().is_lt() { self.turn_value -= 1; }
                self.tan = (-BasedExpr::BASELESS_ONE).checked_div(tan);
            },
            (Some(tan), None) => {
                self.turn_value += rhs.turn_value;
                if tan.cmp_zero().is_lt() { self.turn_value -= 1; }
                self.tan = (-BasedExpr::BASELESS_ONE).checked_div(mem::take(tan));
            },
            (Some(tan_a), Some(tan_b)) => {
                self.turn_value += rhs.turn_value;
                let cmp_a = tan_a.cmp_zero();
                let cmp_b = tan_b.cmp_zero();
                let maybe_overflow = cmp_a.is_gt() && cmp_b.is_gt();
                let maybe_underflow = cmp_a.is_lt() && cmp_b.is_lt();
                self.tan = (&*tan_a + tan_b).checked_div(BasedExpr::BASELESS_ONE - mem::take(tan_a) * tan_b);
                if maybe_overflow  && self.tan < Some(BasedExpr::BASELESS_ZERO) { self.turn_value += 1; }
                if maybe_underflow && self.tan > Some(BasedExpr::BASELESS_ZERO) { self.turn_value -= 1; }
            },
        }
    }
}

impl Add<Angle> for Angle {
    type Output = Angle;

    fn add(mut self, rhs: Angle) -> Self::Output {
        self += rhs;
        self
    }
}

impl Add<&Angle> for Angle {
    type Output = Angle;

    fn add(mut self, rhs: &Angle) -> Self::Output {
        self += rhs;
        self
    }
}

impl Add<Angle> for &Angle {
    type Output = Angle;

    fn add(self, mut rhs: Angle) -> Self::Output {
        rhs += self;
        rhs
    }
}

impl Add<&Angle> for &Angle {
    type Output = Angle;

    fn add(self, rhs: &Angle) -> Self::Output {
        if self.tan.is_none() && rhs.tan.is_some() { return rhs + self; }
        match (self.tan.as_ref(), rhs.tan.as_ref()) {
            (None, None) => {
                let turn_value = self.turn_value + rhs.turn_value - 1;
                let tan = Some(BasedExpr::BASELESS_ZERO);
                Angle::new(turn_value, tan)
            },
            (None, Some(_)) => unreachable!(),
            (Some(tan), None) => {
                let mut turn_value = self.turn_value + rhs.turn_value;
                if tan.cmp_zero().is_lt() { turn_value -= 1; }
                let tan = (-BasedExpr::BASELESS_ONE).checked_div(tan);
                Angle::new(turn_value, tan)
            },
            (Some(tan_a), Some(tan_b)) => {
                let mut turn_value = self.turn_value + rhs.turn_value;
                let cmp_a = tan_a.cmp_zero();
                let cmp_b = tan_b.cmp_zero();
                let maybe_overflow = cmp_a.is_gt() && cmp_b.is_gt();
                let maybe_underflow = cmp_a.is_lt() && cmp_b.is_lt();
                let tan = (tan_a + tan_b).checked_div(BasedExpr::BASELESS_ONE - tan_a * tan_b);
                if maybe_overflow  && tan < Some(BasedExpr::BASELESS_ZERO) { turn_value += 1; }
                if maybe_underflow && tan > Some(BasedExpr::BASELESS_ZERO) { turn_value -= 1; }
                Angle::new(turn_value, tan)
            },
        }
    }
}

impl SubAssign<Angle> for Angle {
    fn sub_assign(&mut self, rhs: Angle) {
        *self += -rhs;
    }
}

impl SubAssign<&Angle> for Angle {
    fn sub_assign(&mut self, rhs: &Angle) {
        *self += -rhs;
    }
}

impl Sub<Angle> for Angle {
    type Output = Angle;
    fn sub(mut self, rhs: Angle) -> Self::Output {
        self -= rhs;
        self
    }
}

impl Sub<&Angle> for Angle {
    type Output = Angle;
    fn sub(mut self, rhs: &Angle) -> Self::Output {
        self -= rhs;
        self
    }
}

impl Sub<Angle> for &Angle {
    type Output = Angle;
    fn sub(self, rhs: Angle) -> Self::Output {
        self + -rhs
    }
}

impl Sub<&Angle> for &Angle {
    type Output = Angle;
    fn sub(self, rhs: &Angle) -> Self::Output {
        self + -rhs
    }
}

#[cfg(test)]
mod test {
    use super::Angle;

    fn angle_neg_test(a: Angle, expected: Angle) {
        assert_eq!(-a.clone(), expected);
        assert_eq!(-&a, expected);
    }

    fn angle_add_test(a: Angle, b: Angle, expected: Angle) {
        assert_eq!(a.clone() + b.clone(), expected);
        assert_eq!(a.clone() + &b, expected);
        assert_eq!(&a + b.clone(), expected);
        assert_eq!(&a + &b, expected);
        assert_eq!(b.clone() + a.clone(), expected);
        assert_eq!(b.clone() + &a, expected);
        assert_eq!(&b + a.clone(), expected);
        assert_eq!(&b + &a, expected);
    }

    #[test]
    fn test_angle_neg() {
        angle_neg_test(angle!(0, 0), angle!(0, 0));
        angle_neg_test(angle!(0, -inf), angle!(1, -inf));
        angle_neg_test(angle!(1, -inf), angle!(0, -inf));
        angle_neg_test(angle!(0, 1), angle!(0, -1));
        angle_neg_test(angle!(2, 1), angle!(-2, -1));
        angle_neg_test(angle!(1, 0 - sqrt 3), angle!(-1, 0 + sqrt 3));
    }

    #[test]
    fn test_angle_add() {
        angle_add_test(angle!(0, 0), angle!(0, 0), angle!(0, 0));
        angle_add_test(angle!(1, 0), angle!(2, 0), angle!(3, 0));
        angle_add_test(angle!(1, -inf), angle!(2, 0), angle!(3, -inf));
        angle_add_test(angle!(1, -inf), angle!(2, -inf), angle!(2, 0));
        angle_add_test(angle!(1, 1/2), angle!(2, 1/2), angle!(3, 4/3));
        angle_add_test(angle!(1, -1/2), angle!(2, -1/2), angle!(3, -4/3));
        angle_add_test(angle!(2, 1), angle!(2, 1), angle!(5, -inf));
        angle_add_test(angle!(2, 1 + 0 sqrt 2), angle!(2, 1 + sqrt 2), angle!(5, -1 - sqrt 2));
        angle_add_test(angle!(2, -1), angle!(2, -1), angle!(4, -inf));
        angle_add_test(angle!(2, -1 + 0 sqrt 2), angle!(2, -1 - sqrt 2), angle!(3, 1 + sqrt 2));
    }
}