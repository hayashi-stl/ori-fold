use std::{mem, ops::{Add, AddAssign, Neg, Sub, SubAssign}};
use std::f32::consts::PI as PI_F32;
use std::f64::consts::PI as PI_F64;

use malachite::base::num::arithmetic::traits::{CheckedDiv, NegAssign};

use crate::BasedExpr;

/// A fold angle represented exactly
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Angle {
    /// Specifically, ⌊angle/π + 1/2⌋ - 1/2
    turn_value: i32,
    tan: Option<BasedExpr>
}

#[allow(unused_macros)]
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

    /// Gets the angle of the vector `[x, y]`. Returns `None` if that is the zero vector.
    pub fn from_xy<X, Y>(x: X, y: Y) -> Option<Self> where X: IntoAngle<Y, Output = Self> {
        x.into_angle(y)
    }

    /// Returns the tangent of the angle. If undefined, returns `None`.
    pub fn tan(&self) -> Option<BasedExpr> {
        self.tan.clone()
    }

    /// Converts this into its turn value $\lfloor angle/π + 1/2\rfloor - 1/2$ and its tangent.
    pub fn into_turn_value_tan(self) -> (i32, Option<BasedExpr>) {
        (self.turn_value, self.tan)
    }

    /// Converts this into an f64.
    /// 
    /// Beware, precision is unspecified and corresponds to the
    /// precision of [`f64::atan2`].
    pub fn to_f64_unspecified(&self) -> f64 {
        let turn = self.turn_value as f64 * PI_F64;
        let atan = self.tan.as_ref()
            .map(|tan| tan.round_to_nearest_f64().atan())
            .unwrap_or(-PI_F64 * 0.5);
        turn + atan
    }
}

pub trait IntoAngle<Y = Self> {
    type Output;
    /// Finds the angle of the vector `[self, y]`,
    /// with the range [-π, π), returning `None` if `self` = `y` = 0.
    /// 
    /// Like `atan2`, but with x first and y second.
    fn into_angle(self, y: Y) -> Option<Self::Output>;
}

impl IntoAngle<f32> for f32 {
    type Output = f32;
    fn into_angle(self, y: Self) -> Option<Self::Output> {
        if self == 0.0 && y == 0.0 { None } else { Some(if self < 0.0 && y == 0.0 { -PI_F32 } else { y.atan2(self) }) }
    }
}
impl IntoAngle<&f32> for f32 { type Output = f32; fn into_angle(self, y: &f32)  -> Option<Self::Output> { self.into_angle(*y) } }
impl IntoAngle<f32> for &f32 { type Output = f32; fn into_angle(self, y: f32)   -> Option<Self::Output> { (*self).into_angle(y) } }
impl IntoAngle<&f32> for &f32 { type Output = f32; fn into_angle(self, y: &f32) -> Option<Self::Output> { (*self).into_angle(*y) } }

impl IntoAngle<f64> for f64 {
    type Output = f64;
    fn into_angle(self, y: Self) -> Option<Self::Output> {
        if self == 0.0 && y == 0.0 { None } else { Some(if self < 0.0 && y == 0.0 { -PI_F64 } else { y.atan2(self) }) }
    }
}
impl IntoAngle<&f64> for f64 { type Output = f64; fn into_angle(self, y: &f64)  -> Option<Self::Output> { self.into_angle(*y) } }
impl IntoAngle<f64> for &f64 { type Output = f64; fn into_angle(self, y: f64)   -> Option<Self::Output> { (*self).into_angle(y) } }
impl IntoAngle<&f64> for &f64 { type Output = f64; fn into_angle(self, y: &f64) -> Option<Self::Output> { (*self).into_angle(*y) } }

fn turn_and_00(x: &BasedExpr, y: &BasedExpr) -> (i32, bool) {
    let x_sign = x.cmp_zero();
    let y_sign = y.cmp_zero();
    let turn_value = if x_sign.is_lt() && y_sign.is_le() {
        -1
    } else if x_sign.is_gt() || y_sign.is_le() {
        0
    } else {
        1
    };
    (turn_value, x_sign.is_eq() && y_sign.is_eq())
}

impl IntoAngle<BasedExpr> for BasedExpr {
    type Output = Angle;
    fn into_angle(self, y: BasedExpr) -> Option<Self::Output> {
        let x = self;
        let (turn_value, zeros) = turn_and_00(&x, &y);
        let tan = if zeros { None } else { Some(y.checked_div(x)) };
        tan.map(|tan| Angle::new(turn_value, tan))
    }
}

impl IntoAngle<&BasedExpr> for BasedExpr {
    type Output = Angle;
    fn into_angle(self, y: &BasedExpr) -> Option<Self::Output> {
        let x = self;
        let (turn_value, zeros) = turn_and_00(&x, &y);
        let tan = if zeros { None } else { Some(y.checked_div(x)) };
        tan.map(|tan| Angle::new(turn_value, tan))
    }
}

impl IntoAngle<BasedExpr> for &BasedExpr {
    type Output = Angle;
    fn into_angle(self, y: BasedExpr) -> Option<Self::Output> {
        let x = self;
        let (turn_value, zeros) = turn_and_00(&x, &y);
        let tan = if zeros { None } else { Some(y.checked_div(x)) };
        tan.map(|tan| Angle::new(turn_value, tan))
    }
}

impl IntoAngle<&BasedExpr> for &BasedExpr {
    type Output = Angle;
    fn into_angle(self, y: &BasedExpr) -> Option<Self::Output> {
        let x = self;
        let (turn_value, zeros) = turn_and_00(&x, &y);
        let tan = if zeros { None } else { Some(y.checked_div(x)) };
        tan.map(|tan| Angle::new(turn_value, tan))
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
    fn add_assign(&mut self, rhs: &Angle) {
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
    use approx::assert_relative_eq;

    use crate::{BasedExpr, angle::IntoAngle, based_expr};

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

    fn into_angle_test(x: BasedExpr, y: BasedExpr, expected: Option<Angle>) {
        assert_eq!(x.clone().into_angle(y.clone()), expected);
        assert_eq!(x.clone().into_angle(&y), expected);
        assert_eq!((&x).into_angle(y.clone()), expected);
        assert_eq!((&x).into_angle(&y), expected);

        // Also check the float implementations
        let x = x.round_to_nearest_f64();
        let y = y.round_to_nearest_f64();
        let expected = expected.map(|e| e.to_f64_unspecified());
        if let Some(expected) = expected {
            assert_relative_eq!(x.into_angle(y).unwrap(), expected);
            assert_relative_eq!(x.into_angle(&y).unwrap(), expected);
            assert_relative_eq!((&x).into_angle(y).unwrap(), expected);
            assert_relative_eq!((&x).into_angle(&y).unwrap(), expected);

            let x = x as f32;
            let y = y as f32;
            let expected = expected as f32;
            assert_relative_eq!(x.into_angle(y).unwrap(), expected);
            assert_relative_eq!(x.into_angle(&y).unwrap(), expected);
            assert_relative_eq!((&x).into_angle(y).unwrap(), expected);
            assert_relative_eq!((&x).into_angle(&y).unwrap(), expected);
        } else {
            assert_eq!(x.into_angle(y), None);
            assert_eq!(x.into_angle(&y), None);
            assert_eq!((&x).into_angle(y), None);
            assert_eq!((&x).into_angle(&y), None);

            let x = x as f32;
            let y = y as f32;
            assert_eq!(x.into_angle(y), None);
            assert_eq!(x.into_angle(&y), None);
            assert_eq!((&x).into_angle(y), None);
            assert_eq!((&x).into_angle(&y), None);
        }
    }

    #[test]
    fn test_into_angle() {
        into_angle_test(based_expr!(0), based_expr!(0), None);
        into_angle_test(based_expr!(1 + 0 sqrt 2), based_expr!(0 + 0 sqrt 2), Some(angle!(0, 0 + 0 sqrt 2))); // basis check
        into_angle_test(based_expr!(1/16), based_expr!(0), Some(angle!(0, 0)));
        into_angle_test(based_expr!(50), based_expr!(0), Some(angle!(0, 0)));

        // Edge cases
        into_angle_test(based_expr!(-1), based_expr!( 0), Some(angle!(-1, 0)));
        into_angle_test(based_expr!(-1), based_expr!(-1), Some(angle!(-1, 1)));
        into_angle_test(based_expr!( 0), based_expr!(-1), Some(angle!(0, -inf)));
        into_angle_test(based_expr!( 1), based_expr!(-1), Some(angle!(0, -1)));
        into_angle_test(based_expr!( 1), based_expr!( 0), Some(angle!(0, 0)));
        into_angle_test(based_expr!( 1), based_expr!( 1), Some(angle!(0, 1)));
        into_angle_test(based_expr!( 0), based_expr!( 1), Some(angle!(1, -inf)));
        into_angle_test(based_expr!(-1), based_expr!( 1), Some(angle!(1, -1)));

        // Other cases (each slice)
        into_angle_test(based_expr!(-5), based_expr!(-1), Some(angle!(-1, 1/5)));
        into_angle_test(based_expr!(-4), based_expr!(-6), Some(angle!(-1, 3/2)));
        into_angle_test(based_expr!( 1), based_expr!(-5), Some(angle!(0, -5)));
        into_angle_test(based_expr!( 6), based_expr!(-4), Some(angle!(0, -2/3)));
        into_angle_test(based_expr!( 5), based_expr!( 1), Some(angle!(0, 1/5)));
        into_angle_test(based_expr!( 4), based_expr!( 6), Some(angle!(0, 3/2)));
        into_angle_test(based_expr!(-1), based_expr!( 5), Some(angle!(1, -5)));
        into_angle_test(based_expr!(-6), based_expr!( 4), Some(angle!(1, -2/3)));
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