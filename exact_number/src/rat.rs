// Workaround to make Rational implement num_traits stuff

use std::{cmp::Ordering, fmt::Display, iter::{Product, Sum}, ops::{Add, AddAssign, Deref, DerefMut, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Shl, ShlAssign, Shr, ShrAssign, Sub, SubAssign}, str::FromStr};

use malachite::{base::{named::Named, num::{arithmetic::traits::{Abs, AbsAssign, AbsDiff, AbsDiffAssign, Ceiling, CeilingAssign, CeilingLogBase, CeilingLogBase2, CeilingLogBasePowerOf2, CheckedDiv, CheckedLogBase, CheckedLogBase2, CheckedLogBasePowerOf2, CheckedRoot, CheckedSqrt, Floor, FloorAssign, FloorLogBase, FloorLogBase2, FloorLogBasePowerOf2, IsPowerOf2, NegAssign, NextPowerOf2, NextPowerOf2Assign, Pow as _, PowAssign, PowerOf2, Reciprocal, ReciprocalAssign, RoundToMultipleAssign, RoundToMultipleOfPowerOf2, RoundToMultipleOfPowerOf2Assign, Sign, Square, SquareAssign}, basic::traits::Two, comparison::traits::{EqAbs, OrdAbs, PartialOrdAbs}, conversion::{string::options::FromSciStringOptions, traits::{ConvertibleFrom, FromSciString, IsInteger, RoundingFrom, SciMantissaAndExponent, ToSci}}, logic::traits::SignificantBits}, rounding_modes::RoundingMode}, natural::conversion::primitive_int_from_natural::SignedFromNaturalError, rational::{arithmetic::traits::{Approximate, ApproximateAssign, DenominatorsInClosedInterval}, conversion::{continued_fraction::{convergents::RationalConvergents, to_continued_fraction::RationalContinuedFraction}, from_primitive_float::RationalFromPrimitiveFloatError, integer_from_rational::IntegerFromRationalError, natural_from_rational::NaturalFromRationalError, primitive_float_from_rational::FloatConversionError, primitive_int_from_rational::{SignedFromRationalError, UnsignedFromRationalError}, traits::{ContinuedFraction, Convergents}}, Rational}, Integer, Natural};
use malachite::base::num::arithmetic::traits::RoundToMultiple;
use malachite::base::num::basic::traits::{Zero as _, One as _};
use num::{pow::Pow, Num, One, Zero};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct Rat(pub Rational);

impl Rat {
    pub const ZERO: Rat = Rat(Rational::ZERO);
    pub const ONE: Rat = Rat(Rational::ONE);
    pub const TWO: Rat = Rat(Rational::TWO);
}

impl Deref for Rat {
    type Target = Rational;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Rat {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Abs for Rat { fn abs(self) -> Self::Output { Rat(self.0.abs()) } type Output = Rat; }
impl Abs for &Rat { fn abs(self) -> Self::Output { Rat((&self.0).abs()) } type Output = Rat; }
impl AbsAssign for Rat { fn abs_assign(&mut self) { self.0.abs_assign(); } }
impl AbsDiff<Rat> for Rat { fn abs_diff(self, other: Rat) -> Self::Output { Rat(self.0.abs_diff(other.0)) } type Output = Rat;}
impl AbsDiff<&Rat> for Rat { fn abs_diff(self, other: &Rat) -> Self::Output { Rat(self.0.abs_diff(&other.0)) } type Output = Rat;}
impl AbsDiff<Rat> for &Rat { fn abs_diff(self, other: Rat) -> Self::Output { Rat((&self.0).abs_diff(other.0)) } type Output = Rat;}
impl AbsDiff<&Rat> for &Rat { fn abs_diff(self, other: &Rat) -> Self::Output { Rat((&self.0).abs_diff(&other.0)) } type Output = Rat;}
impl AbsDiffAssign<Rat> for Rat { fn abs_diff_assign(&mut self, other: Rat) { self.0.abs_diff_assign(other.0); }}
impl AbsDiffAssign<&Rat> for Rat { fn abs_diff_assign(&mut self, other: &Rat) { self.0.abs_diff_assign(&other.0); }}

impl Add<Rat> for Rat { fn add(self, other: Rat) -> Self::Output { Rat(self.0.add(other.0)) } type Output = Rat;}
impl Add<&Rat> for Rat { fn add(self, other: &Rat) -> Self::Output { Rat(self.0.add(&other.0)) } type Output = Rat;}
impl Add<Rat> for &Rat { fn add(self, other: Rat) -> Self::Output { Rat((&self.0).add(other.0)) } type Output = Rat;}
impl Add<&Rat> for &Rat { fn add(self, other: &Rat) -> Self::Output { Rat((&self.0).add(&other.0)) } type Output = Rat;}
impl AddAssign<Rat> for Rat { fn add_assign(&mut self, other: Rat) { self.0.add_assign(other.0); }}
impl AddAssign<&Rat> for Rat { fn add_assign(&mut self, other: &Rat) { self.0.add_assign(&other.0); }}

impl Approximate for Rat { fn approximate(self, max_denominator: &Natural) -> Rational { self.0.approximate(max_denominator) } }
impl Approximate for &Rat { fn approximate(self, max_denominator: &Natural) -> Rational { (&self.0).approximate(max_denominator) } }
impl ApproximateAssign for Rat { fn approximate_assign(&mut self, max_denominator: &Natural) { self.0.approximate_assign(max_denominator) } }

impl Ceiling for Rat { fn ceiling(self) -> Self::Output { self.0.ceiling() } type Output = Integer; }
impl Ceiling for &Rat { fn ceiling(self) -> Self::Output { (&self.0).ceiling() } type Output = Integer; }
impl CeilingAssign for Rat { fn ceiling_assign(&mut self) { self.0.ceiling_assign(); } }
impl CeilingLogBase<&Rat> for &Rat  { fn ceiling_log_base(self, base: &Rat) -> Self::Output { self.0.ceiling_log_base(&base.0) } type Output = i64; }
impl CeilingLogBase2 for &Rat { fn ceiling_log_base_2(self) -> Self::Output { self.0.ceiling_log_base_2() } type Output = i64; }
impl CeilingLogBasePowerOf2<i64> for &Rat { fn ceiling_log_base_power_of_2(self, pow: i64) -> Self::Output { self.0.ceiling_log_base_power_of_2(pow) } type Output = i64;}

impl CheckedDiv<Rat> for Rat { fn checked_div(self, other: Rat) -> Option<Self::Output> { self.0.checked_div(other.0).map(Rat) } type Output = Rat;}
impl CheckedDiv<&Rat> for Rat { fn checked_div(self, other: &Rat) -> Option<Self::Output> { self.0.checked_div(&other.0).map(Rat) } type Output = Rat;}
impl CheckedDiv<Rat> for &Rat { fn checked_div(self, other: Rat) -> Option<Self::Output> { (&self.0).checked_div(other.0).map(Rat) } type Output = Rat;}
impl CheckedDiv<&Rat> for &Rat { fn checked_div(self, other: &Rat) -> Option<Self::Output> { (&self.0).checked_div(&other.0).map(Rat) } type Output = Rat;}
impl CheckedLogBase<&Rat> for &Rat  { fn checked_log_base(self, base: &Rat) -> Option<Self::Output> { self.0.checked_log_base(&base.0) } type Output = i64; }
impl CheckedLogBase2 for &Rat { fn checked_log_base_2(self) -> Option<Self::Output> { self.0.checked_log_base_2() } type Output = i64; }
impl CheckedLogBasePowerOf2<i64> for &Rat { fn checked_log_base_power_of_2(self, pow: i64) -> Option<Self::Output> { self.0.checked_log_base_power_of_2(pow) } type Output = i64;}
impl CheckedRoot<i64> for Rat { fn checked_root(self, pow: i64) -> Option<Self::Output> { self.0.checked_root(pow).map(Rat) } type Output = Rat; }
impl CheckedRoot<i64> for &Rat { fn checked_root(self, pow: i64) -> Option<Self::Output> { (&self.0).checked_root(pow).map(Rat) } type Output = Rat; }
impl CheckedRoot<u64> for Rat { fn checked_root(self, pow: u64) -> Option<Self::Output> { self.0.checked_root(pow).map(Rat) } type Output = Rat; }
impl CheckedRoot<u64> for &Rat { fn checked_root(self, pow: u64) -> Option<Self::Output> { (&self.0).checked_root(pow).map(Rat) } type Output = Rat; }
impl CheckedSqrt for Rat { fn checked_sqrt(self) -> Option<Self::Output> { self.0.checked_sqrt().map(Rat) } type Output = Rat; }
impl CheckedSqrt for &Rat { fn checked_sqrt(self) -> Option<Self::Output> { (&self.0).checked_sqrt().map(Rat) } type Output = Rat; }

impl ContinuedFraction for Rat { fn continued_fraction(self) -> (Integer, Self::CF) { self.0.continued_fraction() } type CF = RationalContinuedFraction; }
impl ContinuedFraction for &Rat { fn continued_fraction(self) -> (Integer, Self::CF) { (&self.0).continued_fraction() } type CF = RationalContinuedFraction; }
impl Convergents for Rat { fn convergents(self) -> Self::C { self.0.convergents() } type C = RationalConvergents; }
impl Convergents for &Rat { fn convergents(self) -> Self::C { (&self.0).convergents() } type C = RationalConvergents; }

impl ConvertibleFrom<&Rat> for i8    { fn convertible_from(value: &Rat) -> bool { i8   ::convertible_from(&value.0) } }
impl ConvertibleFrom<&Rat> for i16   { fn convertible_from(value: &Rat) -> bool { i16  ::convertible_from(&value.0) } }
impl ConvertibleFrom<&Rat> for i32   { fn convertible_from(value: &Rat) -> bool { i32  ::convertible_from(&value.0) } }
impl ConvertibleFrom<&Rat> for i64   { fn convertible_from(value: &Rat) -> bool { i64  ::convertible_from(&value.0) } }
impl ConvertibleFrom<&Rat> for i128  { fn convertible_from(value: &Rat) -> bool { i128 ::convertible_from(&value.0) } }
impl ConvertibleFrom<&Rat> for isize { fn convertible_from(value: &Rat) -> bool { isize::convertible_from(&value.0) } }
impl ConvertibleFrom<&Rat> for u8    { fn convertible_from(value: &Rat) -> bool { u8   ::convertible_from(&value.0) } }
impl ConvertibleFrom<&Rat> for u16   { fn convertible_from(value: &Rat) -> bool { u16  ::convertible_from(&value.0) } }
impl ConvertibleFrom<&Rat> for u32   { fn convertible_from(value: &Rat) -> bool { u32  ::convertible_from(&value.0) } }
impl ConvertibleFrom<&Rat> for u64   { fn convertible_from(value: &Rat) -> bool { u64  ::convertible_from(&value.0) } }
impl ConvertibleFrom<&Rat> for u128  { fn convertible_from(value: &Rat) -> bool { u128 ::convertible_from(&value.0) } }
impl ConvertibleFrom<&Rat> for usize { fn convertible_from(value: &Rat) -> bool { usize::convertible_from(&value.0) } }
impl ConvertibleFrom<Rat> for f32 { fn convertible_from(value: Rat) -> bool { f32::convertible_from(value.0) } }
impl ConvertibleFrom<Rat> for f64 { fn convertible_from(value: Rat) -> bool { f64::convertible_from(value.0) } }
impl ConvertibleFrom<f32> for Rat { fn convertible_from(value: f32) -> bool { Rational::convertible_from(value) } }
impl ConvertibleFrom<f64> for Rat { fn convertible_from(value: f64) -> bool { Rational::convertible_from(value) } }

impl Default for Rat { fn default() -> Self { Rat(Rational::default()) } }
impl Display for Rat { fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { self.0.fmt(f) } }

impl Div<Rat> for Rat { fn div(self, other: Rat) -> Self::Output { Rat(self.0.div(other.0)) } type Output = Rat;}
impl Div<&Rat> for Rat { fn div(self, other: &Rat) -> Self::Output { Rat(self.0.div(&other.0)) } type Output = Rat;}
impl Div<Rat> for &Rat { fn div(self, other: Rat) -> Self::Output { Rat((&self.0).div(other.0)) } type Output = Rat;}
impl Div<&Rat> for &Rat { fn div(self, other: &Rat) -> Self::Output { Rat((&self.0).div(&other.0)) } type Output = Rat;}
impl DivAssign<Rat> for Rat { fn div_assign(&mut self, other: Rat) { self.0.div_assign(other.0); }}
impl DivAssign<&Rat> for Rat { fn div_assign(&mut self, other: &Rat) { self.0.div_assign(&other.0); }}

impl EqAbs for Rat { fn eq_abs(&self, other: &Self) -> bool { self.0.eq_abs(&other.0) } }
impl EqAbs<Integer> for Rat { fn eq_abs(&self, other: &Integer) -> bool { self.0.eq_abs(other) } }
impl EqAbs<Natural> for Rat { fn eq_abs(&self, other: &Natural) -> bool { self.0.eq_abs(other) } }
impl EqAbs<u8   > for Rat { fn eq_abs(&self, other: &u8   ) -> bool { self.0.eq_abs(other) } }
impl EqAbs<u16  > for Rat { fn eq_abs(&self, other: &u16  ) -> bool { self.0.eq_abs(other) } }
impl EqAbs<u32  > for Rat { fn eq_abs(&self, other: &u32  ) -> bool { self.0.eq_abs(other) } }
impl EqAbs<u64  > for Rat { fn eq_abs(&self, other: &u64  ) -> bool { self.0.eq_abs(other) } }
impl EqAbs<u128 > for Rat { fn eq_abs(&self, other: &u128 ) -> bool { self.0.eq_abs(other) } }
impl EqAbs<usize> for Rat { fn eq_abs(&self, other: &usize) -> bool { self.0.eq_abs(other) } }
impl EqAbs<i8   > for Rat { fn eq_abs(&self, other: &i8   ) -> bool { self.0.eq_abs(other) } }
impl EqAbs<i16  > for Rat { fn eq_abs(&self, other: &i16  ) -> bool { self.0.eq_abs(other) } }
impl EqAbs<i32  > for Rat { fn eq_abs(&self, other: &i32  ) -> bool { self.0.eq_abs(other) } }
impl EqAbs<i64  > for Rat { fn eq_abs(&self, other: &i64  ) -> bool { self.0.eq_abs(other) } }
impl EqAbs<i128 > for Rat { fn eq_abs(&self, other: &i128 ) -> bool { self.0.eq_abs(other) } }
impl EqAbs<isize> for Rat { fn eq_abs(&self, other: &isize) -> bool { self.0.eq_abs(other) } }
impl EqAbs<f32  > for Rat { fn eq_abs(&self, other: &f32  ) -> bool { self.0.eq_abs(other) } }
impl EqAbs<f64  > for Rat { fn eq_abs(&self, other: &f64  ) -> bool { self.0.eq_abs(other) } }
impl EqAbs<Rat> for Integer { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for Natural { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for u8    { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for u16   { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for u32   { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for u64   { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for u128  { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for usize { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for i8    { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for i16   { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for i32   { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for i64   { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for i128  { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for isize { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for f32   { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }
impl EqAbs<Rat> for f64   { fn eq_abs(&self, other: &Rat) -> bool { self.eq_abs(&other.0) } }

impl Floor for Rat { fn floor(self) -> Self::Output { self.0.floor() } type Output = Integer; }
impl Floor for &Rat { fn floor(self) -> Self::Output { (&self.0).floor() } type Output = Integer; }
impl FloorAssign for Rat { fn floor_assign(&mut self) { self.0.floor_assign(); } }
impl FloorLogBase<&Rat> for &Rat  { fn floor_log_base(self, base: &Rat) -> Self::Output { self.0.floor_log_base(&base.0) } type Output = i64; }
impl FloorLogBase2 for &Rat { fn floor_log_base_2(self) -> Self::Output { self.0.floor_log_base_2() } type Output = i64; }
impl FloorLogBasePowerOf2<i64> for &Rat { fn floor_log_base_power_of_2(self, pow: i64) -> Self::Output { self.0.floor_log_base_power_of_2(pow) } type Output = i64;}

impl From<Rational> for Rat { fn from(value: Rational) -> Self { Rat(value) } }
impl From<&Integer> for Rat { fn from(value: &Integer) -> Self { Rational::from(value).into() } }
impl From<&Natural> for Rat { fn from(value: &Natural) -> Self { Rational::from(value).into() } }
impl From<Integer> for Rat { fn from(value: Integer) -> Self { Rational::from(value).into() } }
impl From<Natural> for Rat { fn from(value: Natural) -> Self { Rational::from(value).into() } }
impl From<bool> for Rat { fn from(value: bool) -> Self { Rational::from(value).into() } }
impl From<i8   > for Rat { fn from(value: i8   ) -> Self { Rational::from(value).into() } }
impl From<i16  > for Rat { fn from(value: i16  ) -> Self { Rational::from(value).into() } }
impl From<i32  > for Rat { fn from(value: i32  ) -> Self { Rational::from(value).into() } }
impl From<i64  > for Rat { fn from(value: i64  ) -> Self { Rational::from(value).into() } }
impl From<i128 > for Rat { fn from(value: i128 ) -> Self { Rational::from(value).into() } }
impl From<isize> for Rat { fn from(value: isize) -> Self { Rational::from(value).into() } }
impl From<u8   > for Rat { fn from(value: u8   ) -> Self { Rational::from(value).into() } }
impl From<u16  > for Rat { fn from(value: u16  ) -> Self { Rational::from(value).into() } }
impl From<u32  > for Rat { fn from(value: u32  ) -> Self { Rational::from(value).into() } }
impl From<u64  > for Rat { fn from(value: u64  ) -> Self { Rational::from(value).into() } }
impl From<u128 > for Rat { fn from(value: u128 ) -> Self { Rational::from(value).into() } }
impl From<usize> for Rat { fn from(value: usize) -> Self { Rational::from(value).into() } }
impl FromSciString for Rat { fn from_sci_string_with_options(s: &str, options: FromSciStringOptions) -> Option<Self> {
    Rational::from_sci_string_with_options(s, options).map(Rat) } }
impl FromStr for Rat { fn from_str(s: &str) -> Result<Self, Self::Err> { Rational::from_str(s).map(Rat) } type Err = (); }

impl IsInteger for &Rat { fn is_integer(self) -> bool { self.0.is_integer() } }
impl IsPowerOf2 for Rat { fn is_power_of_2(&self) -> bool { self.0.is_power_of_2() } }

impl Mul<Rat> for Rat { fn mul(self, other: Rat) -> Self::Output { Rat(self.0.mul(other.0)) } type Output = Rat;}
impl Mul<&Rat> for Rat { fn mul(self, other: &Rat) -> Self::Output { Rat(self.0.mul(&other.0)) } type Output = Rat;}
impl Mul<Rat> for &Rat { fn mul(self, other: Rat) -> Self::Output { Rat((&self.0).mul(other.0)) } type Output = Rat;}
impl Mul<&Rat> for &Rat { fn mul(self, other: &Rat) -> Self::Output { Rat((&self.0).mul(&other.0)) } type Output = Rat;}
impl MulAssign<Rat> for Rat { fn mul_assign(&mut self, other: Rat) { self.0.mul_assign(other.0); }}
impl MulAssign<&Rat> for Rat { fn mul_assign(&mut self, other: &Rat) { self.0.mul_assign(&other.0); }}

impl Named for Rat { const NAME: &'static str = "Rat"; }

impl Neg for Rat { fn neg(self) -> Self::Output { Rat(self.0.neg()) } type Output = Rat; }
impl Neg for &Rat { fn neg(self) -> Self::Output { Rat((&self.0).neg()) } type Output = Rat; }
impl NegAssign for Rat { fn neg_assign(&mut self) { self.0.neg_assign(); } }

impl NextPowerOf2 for Rat { fn next_power_of_2(self) -> Self::Output { Rat(self.0.next_power_of_2()) } type Output = Rat; }
impl NextPowerOf2 for &Rat { fn next_power_of_2(self) -> Self::Output { Rat((&self.0).next_power_of_2()) } type Output = Rat; }
impl NextPowerOf2Assign for Rat { fn next_power_of_2_assign(&mut self) { self.0.next_power_of_2_assign(); } }

impl OrdAbs for Rat { fn cmp_abs(&self, other: &Self) -> std::cmp::Ordering { self.0.cmp_abs(&other.0) } }

impl PartialEq<Integer> for Rat { fn eq(&self, other: &Integer) -> bool { self.0.eq(other) } }
impl PartialEq<Natural> for Rat { fn eq(&self, other: &Natural) -> bool { self.0.eq(other) } }
impl PartialEq<u8   > for Rat { fn eq(&self, other: &u8   ) -> bool { self.0.eq(other) } }
impl PartialEq<u16  > for Rat { fn eq(&self, other: &u16  ) -> bool { self.0.eq(other) } }
impl PartialEq<u32  > for Rat { fn eq(&self, other: &u32  ) -> bool { self.0.eq(other) } }
impl PartialEq<u64  > for Rat { fn eq(&self, other: &u64  ) -> bool { self.0.eq(other) } }
impl PartialEq<u128 > for Rat { fn eq(&self, other: &u128 ) -> bool { self.0.eq(other) } }
impl PartialEq<usize> for Rat { fn eq(&self, other: &usize) -> bool { self.0.eq(other) } }
impl PartialEq<i8   > for Rat { fn eq(&self, other: &i8   ) -> bool { self.0.eq(other) } }
impl PartialEq<i16  > for Rat { fn eq(&self, other: &i16  ) -> bool { self.0.eq(other) } }
impl PartialEq<i32  > for Rat { fn eq(&self, other: &i32  ) -> bool { self.0.eq(other) } }
impl PartialEq<i64  > for Rat { fn eq(&self, other: &i64  ) -> bool { self.0.eq(other) } }
impl PartialEq<i128 > for Rat { fn eq(&self, other: &i128 ) -> bool { self.0.eq(other) } }
impl PartialEq<isize> for Rat { fn eq(&self, other: &isize) -> bool { self.0.eq(other) } }
impl PartialEq<f32  > for Rat { fn eq(&self, other: &f32  ) -> bool { self.0.eq(other) } }
impl PartialEq<f64  > for Rat { fn eq(&self, other: &f64  ) -> bool { self.0.eq(other) } }
impl PartialEq<Rat> for Integer { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for Natural { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for u8    { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for u16   { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for u32   { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for u64   { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for u128  { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for usize { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for i8    { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for i16   { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for i32   { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for i64   { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for i128  { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for isize { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for f32   { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }
impl PartialEq<Rat> for f64   { fn eq(&self, other: &Rat) -> bool { self.eq(&other.0) } }

impl PartialOrd<Integer> for Rat { fn partial_cmp(&self, other: &Integer) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<Natural> for Rat { fn partial_cmp(&self, other: &Natural) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<u8   > for Rat { fn partial_cmp(&self, other: &u8   ) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<u16  > for Rat { fn partial_cmp(&self, other: &u16  ) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<u32  > for Rat { fn partial_cmp(&self, other: &u32  ) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<u64  > for Rat { fn partial_cmp(&self, other: &u64  ) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<u128 > for Rat { fn partial_cmp(&self, other: &u128 ) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<usize> for Rat { fn partial_cmp(&self, other: &usize) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<i8   > for Rat { fn partial_cmp(&self, other: &i8   ) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<i16  > for Rat { fn partial_cmp(&self, other: &i16  ) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<i32  > for Rat { fn partial_cmp(&self, other: &i32  ) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<i64  > for Rat { fn partial_cmp(&self, other: &i64  ) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<i128 > for Rat { fn partial_cmp(&self, other: &i128 ) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<isize> for Rat { fn partial_cmp(&self, other: &isize) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<f32  > for Rat { fn partial_cmp(&self, other: &f32  ) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<f64  > for Rat { fn partial_cmp(&self, other: &f64  ) -> Option<Ordering> { self.0.partial_cmp(other) } }
impl PartialOrd<Rat> for Integer { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for Natural { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for u8    { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for u16   { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for u32   { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for u64   { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for u128  { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for usize { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for i8    { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for i16   { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for i32   { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for i64   { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for i128  { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for isize { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for f32   { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }
impl PartialOrd<Rat> for f64   { fn partial_cmp(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp(&other.0) } }

impl PartialOrdAbs for Rat { fn partial_cmp_abs(&self, other: &Self) -> Option<std::cmp::Ordering> { self.0.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Integer> for Rat { fn partial_cmp_abs(&self, other: &Integer) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<Natural> for Rat { fn partial_cmp_abs(&self, other: &Natural) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<u8   > for Rat { fn partial_cmp_abs(&self, other: &u8   ) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<u16  > for Rat { fn partial_cmp_abs(&self, other: &u16  ) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<u32  > for Rat { fn partial_cmp_abs(&self, other: &u32  ) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<u64  > for Rat { fn partial_cmp_abs(&self, other: &u64  ) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<u128 > for Rat { fn partial_cmp_abs(&self, other: &u128 ) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<usize> for Rat { fn partial_cmp_abs(&self, other: &usize) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<i8   > for Rat { fn partial_cmp_abs(&self, other: &i8   ) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<i16  > for Rat { fn partial_cmp_abs(&self, other: &i16  ) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<i32  > for Rat { fn partial_cmp_abs(&self, other: &i32  ) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<i64  > for Rat { fn partial_cmp_abs(&self, other: &i64  ) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<i128 > for Rat { fn partial_cmp_abs(&self, other: &i128 ) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<isize> for Rat { fn partial_cmp_abs(&self, other: &isize) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<f32  > for Rat { fn partial_cmp_abs(&self, other: &f32  ) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<f64  > for Rat { fn partial_cmp_abs(&self, other: &f64  ) -> Option<Ordering> { self.0.partial_cmp_abs(other) } }
impl PartialOrdAbs<Rat> for Integer { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for Natural { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for u8    { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for u16   { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for u32   { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for u64   { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for u128  { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for usize { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for i8    { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for i16   { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for i32   { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for i64   { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for i128  { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for isize { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for f32   { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }
impl PartialOrdAbs<Rat> for f64   { fn partial_cmp_abs(&self, other: &Rat) -> Option<Ordering> { self.partial_cmp_abs(&other.0) } }

impl Pow<i64> for Rat { fn pow(self, other: i64) -> Self::Output { Rat(self.0.pow(other)) } type Output = Rat;}
impl Pow<i64> for &Rat { fn pow(self, other: i64) -> Self::Output { Rat((&self.0).pow(other)) } type Output = Rat;}
impl PowAssign<i64> for Rat { fn pow_assign(&mut self, other: i64) { self.0.pow_assign(other); }}
impl PowAssign<u64> for Rat { fn pow_assign(&mut self, other: u64) { self.0.pow_assign(other); }}

impl PowerOf2<i64> for Rat { fn power_of_2(pow: i64) -> Self { Rational::power_of_2(pow).into() } }
impl PowerOf2<u64> for Rat { fn power_of_2(pow: u64) -> Self { Rational::power_of_2(pow).into() } }

impl Product for Rat { fn product<I: Iterator<Item = Self>>(iter: I) -> Self { Rational::product(iter.map(|r| r.0)).into() } }
impl<'a> Product<&'a Rat> for Rat { fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self { Rational::product(iter.map(|r| &r.0)).into() } }

impl Reciprocal for Rat { fn reciprocal(self) -> Self::Output { Rat(self.0.reciprocal()) } type Output = Rat; }
impl Reciprocal for &Rat { fn reciprocal(self) -> Self::Output { Rat((&self.0).reciprocal()) } type Output = Rat; }
impl ReciprocalAssign for Rat { fn reciprocal_assign(&mut self) { self.0.reciprocal_assign(); } }

fn first_to_rat<T>(tup: (Rational, T)) -> (Rat, T) {
    (Rat(tup.0), tup.1)
}
impl RoundToMultiple<Rat> for Rat { fn round_to_multiple(self, other: Rat, mode: RoundingMode) -> (Self::Output, Ordering) { first_to_rat(self.0.round_to_multiple(other.0, mode)) } type Output = Rat;}
impl RoundToMultiple<&Rat> for Rat { fn round_to_multiple(self, other: &Rat, mode: RoundingMode) -> (Self::Output, Ordering) { first_to_rat(self.0.round_to_multiple(&other.0, mode)) } type Output = Rat;}
impl RoundToMultiple<Rat> for &Rat { fn round_to_multiple(self, other: Rat, mode: RoundingMode) -> (Self::Output, Ordering) { first_to_rat((&self.0).round_to_multiple(other.0, mode)) } type Output = Rat;}
impl RoundToMultiple<&Rat> for &Rat { fn round_to_multiple(self, other: &Rat, mode: RoundingMode) -> (Self::Output, Ordering) { first_to_rat((&self.0).round_to_multiple(&other.0, mode)) } type Output = Rat;}
impl RoundToMultipleAssign<Rat> for Rat { fn round_to_multiple_assign(&mut self, other: Rat, mode: RoundingMode) -> Ordering { self.0.round_to_multiple_assign(other.0, mode) }}
impl RoundToMultipleAssign<&Rat> for Rat { fn round_to_multiple_assign(&mut self, other: &Rat, mode: RoundingMode) -> Ordering { self.0.round_to_multiple_assign(&other.0, mode) }}
impl RoundToMultipleOfPowerOf2<i64> for Rat { fn round_to_multiple_of_power_of_2(self, other: i64, mode: RoundingMode) -> (Self::Output, Ordering) { first_to_rat(self.0.round_to_multiple_of_power_of_2(other, mode)) } type Output = Rat;}
impl RoundToMultipleOfPowerOf2<i64> for &Rat { fn round_to_multiple_of_power_of_2(self, other: i64, mode: RoundingMode) -> (Self::Output, Ordering) { first_to_rat((&self.0).round_to_multiple_of_power_of_2(other, mode)) } type Output = Rat;}
impl RoundToMultipleOfPowerOf2Assign<i64> for Rat { fn round_to_multiple_of_power_of_2_assign(&mut self, other: i64, mode: RoundingMode) -> Ordering { self.0.round_to_multiple_of_power_of_2_assign(other, mode) }}

impl RoundingFrom<&Rat> for i8    { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { i8   ::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for i16   { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { i16  ::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for i32   { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { i32  ::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for i64   { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { i64  ::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for i128  { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { i128 ::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for isize { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { isize::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for u8    { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { u8   ::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for u16   { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { u16  ::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for u32   { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { u32  ::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for u64   { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { u64  ::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for u128  { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { u128 ::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for usize { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { usize::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for f32   { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { f32  ::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for f64   { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { f64  ::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for Integer { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { Integer::rounding_from(&value.0, rm) } }
impl RoundingFrom<&Rat> for Natural { fn rounding_from(value: &Rat, rm: RoundingMode) -> (Self, Ordering) { Natural::rounding_from(&value.0, rm) } }
impl RoundingFrom<Rat> for f32 { fn rounding_from(value: Rat, rm: RoundingMode) -> (Self, Ordering) { f32::rounding_from(value.0, rm) } }
impl RoundingFrom<Rat> for f64 { fn rounding_from(value: Rat, rm: RoundingMode) -> (Self, Ordering) { f64::rounding_from(value.0, rm) } }
impl RoundingFrom<Rat> for Integer { fn rounding_from(value: Rat, rm: RoundingMode) -> (Self, Ordering) { Integer::rounding_from(value.0, rm) } }
impl RoundingFrom<Rat> for Natural { fn rounding_from(value: Rat, rm: RoundingMode) -> (Self, Ordering) { Natural::rounding_from(value.0, rm) } }

impl Shl<i8   > for Rat { fn shl(self, rhs: i8   ) -> Self::Output { Rat(self.0.shl(rhs)) } type Output = Rat; }
impl Shl<i16  > for Rat { fn shl(self, rhs: i16  ) -> Self::Output { Rat(self.0.shl(rhs)) } type Output = Rat; }
impl Shl<i32  > for Rat { fn shl(self, rhs: i32  ) -> Self::Output { Rat(self.0.shl(rhs)) } type Output = Rat; }
impl Shl<i64  > for Rat { fn shl(self, rhs: i64  ) -> Self::Output { Rat(self.0.shl(rhs)) } type Output = Rat; }
impl Shl<i128 > for Rat { fn shl(self, rhs: i128 ) -> Self::Output { Rat(self.0.shl(rhs)) } type Output = Rat; }
impl Shl<isize> for Rat { fn shl(self, rhs: isize) -> Self::Output { Rat(self.0.shl(rhs)) } type Output = Rat; }
impl Shl<u8   > for Rat { fn shl(self, rhs: u8   ) -> Self::Output { Rat(self.0.shl(rhs)) } type Output = Rat; }
impl Shl<u16  > for Rat { fn shl(self, rhs: u16  ) -> Self::Output { Rat(self.0.shl(rhs)) } type Output = Rat; }
impl Shl<u32  > for Rat { fn shl(self, rhs: u32  ) -> Self::Output { Rat(self.0.shl(rhs)) } type Output = Rat; }
impl Shl<u64  > for Rat { fn shl(self, rhs: u64  ) -> Self::Output { Rat(self.0.shl(rhs)) } type Output = Rat; }
impl Shl<u128 > for Rat { fn shl(self, rhs: u128 ) -> Self::Output { Rat(self.0.shl(rhs)) } type Output = Rat; }
impl Shl<usize> for Rat { fn shl(self, rhs: usize) -> Self::Output { Rat(self.0.shl(rhs)) } type Output = Rat; }
impl Shl<i8   > for &Rat { fn shl(self, rhs: i8   ) -> Self::Output { Rat((&self.0).shl(rhs)) } type Output = Rat; }
impl Shl<i16  > for &Rat { fn shl(self, rhs: i16  ) -> Self::Output { Rat((&self.0).shl(rhs)) } type Output = Rat; }
impl Shl<i32  > for &Rat { fn shl(self, rhs: i32  ) -> Self::Output { Rat((&self.0).shl(rhs)) } type Output = Rat; }
impl Shl<i64  > for &Rat { fn shl(self, rhs: i64  ) -> Self::Output { Rat((&self.0).shl(rhs)) } type Output = Rat; }
impl Shl<i128 > for &Rat { fn shl(self, rhs: i128 ) -> Self::Output { Rat((&self.0).shl(rhs)) } type Output = Rat; }
impl Shl<isize> for &Rat { fn shl(self, rhs: isize) -> Self::Output { Rat((&self.0).shl(rhs)) } type Output = Rat; }
impl Shl<u8   > for &Rat { fn shl(self, rhs: u8   ) -> Self::Output { Rat((&self.0).shl(rhs)) } type Output = Rat; }
impl Shl<u16  > for &Rat { fn shl(self, rhs: u16  ) -> Self::Output { Rat((&self.0).shl(rhs)) } type Output = Rat; }
impl Shl<u32  > for &Rat { fn shl(self, rhs: u32  ) -> Self::Output { Rat((&self.0).shl(rhs)) } type Output = Rat; }
impl Shl<u64  > for &Rat { fn shl(self, rhs: u64  ) -> Self::Output { Rat((&self.0).shl(rhs)) } type Output = Rat; }
impl Shl<u128 > for &Rat { fn shl(self, rhs: u128 ) -> Self::Output { Rat((&self.0).shl(rhs)) } type Output = Rat; }
impl Shl<usize> for &Rat { fn shl(self, rhs: usize) -> Self::Output { Rat((&self.0).shl(rhs)) } type Output = Rat; }
impl ShlAssign<i8   > for Rat { fn shl_assign(&mut self, rhs: i8   ) { self.0.shl_assign(rhs) } }
impl ShlAssign<i16  > for Rat { fn shl_assign(&mut self, rhs: i16  ) { self.0.shl_assign(rhs) } }
impl ShlAssign<i32  > for Rat { fn shl_assign(&mut self, rhs: i32  ) { self.0.shl_assign(rhs) } }
impl ShlAssign<i64  > for Rat { fn shl_assign(&mut self, rhs: i64  ) { self.0.shl_assign(rhs) } }
impl ShlAssign<i128 > for Rat { fn shl_assign(&mut self, rhs: i128 ) { self.0.shl_assign(rhs) } }
impl ShlAssign<isize> for Rat { fn shl_assign(&mut self, rhs: isize) { self.0.shl_assign(rhs) } }
impl ShlAssign<u8   > for Rat { fn shl_assign(&mut self, rhs: u8   ) { self.0.shl_assign(rhs) } }
impl ShlAssign<u16  > for Rat { fn shl_assign(&mut self, rhs: u16  ) { self.0.shl_assign(rhs) } }
impl ShlAssign<u32  > for Rat { fn shl_assign(&mut self, rhs: u32  ) { self.0.shl_assign(rhs) } }
impl ShlAssign<u64  > for Rat { fn shl_assign(&mut self, rhs: u64  ) { self.0.shl_assign(rhs) } }
impl ShlAssign<u128 > for Rat { fn shl_assign(&mut self, rhs: u128 ) { self.0.shl_assign(rhs) } }
impl ShlAssign<usize> for Rat { fn shl_assign(&mut self, rhs: usize) { self.0.shl_assign(rhs) } }
impl Shr<i8   > for Rat { fn shr(self, rhs: i8   ) -> Self::Output { Rat(self.0.shr(rhs)) } type Output = Rat; }
impl Shr<i16  > for Rat { fn shr(self, rhs: i16  ) -> Self::Output { Rat(self.0.shr(rhs)) } type Output = Rat; }
impl Shr<i32  > for Rat { fn shr(self, rhs: i32  ) -> Self::Output { Rat(self.0.shr(rhs)) } type Output = Rat; }
impl Shr<i64  > for Rat { fn shr(self, rhs: i64  ) -> Self::Output { Rat(self.0.shr(rhs)) } type Output = Rat; }
impl Shr<i128 > for Rat { fn shr(self, rhs: i128 ) -> Self::Output { Rat(self.0.shr(rhs)) } type Output = Rat; }
impl Shr<isize> for Rat { fn shr(self, rhs: isize) -> Self::Output { Rat(self.0.shr(rhs)) } type Output = Rat; }
impl Shr<u8   > for Rat { fn shr(self, rhs: u8   ) -> Self::Output { Rat(self.0.shr(rhs)) } type Output = Rat; }
impl Shr<u16  > for Rat { fn shr(self, rhs: u16  ) -> Self::Output { Rat(self.0.shr(rhs)) } type Output = Rat; }
impl Shr<u32  > for Rat { fn shr(self, rhs: u32  ) -> Self::Output { Rat(self.0.shr(rhs)) } type Output = Rat; }
impl Shr<u64  > for Rat { fn shr(self, rhs: u64  ) -> Self::Output { Rat(self.0.shr(rhs)) } type Output = Rat; }
impl Shr<u128 > for Rat { fn shr(self, rhs: u128 ) -> Self::Output { Rat(self.0.shr(rhs)) } type Output = Rat; }
impl Shr<usize> for Rat { fn shr(self, rhs: usize) -> Self::Output { Rat(self.0.shr(rhs)) } type Output = Rat; }
impl Shr<i8   > for &Rat { fn shr(self, rhs: i8   ) -> Self::Output { Rat((&self.0).shr(rhs)) } type Output = Rat; }
impl Shr<i16  > for &Rat { fn shr(self, rhs: i16  ) -> Self::Output { Rat((&self.0).shr(rhs)) } type Output = Rat; }
impl Shr<i32  > for &Rat { fn shr(self, rhs: i32  ) -> Self::Output { Rat((&self.0).shr(rhs)) } type Output = Rat; }
impl Shr<i64  > for &Rat { fn shr(self, rhs: i64  ) -> Self::Output { Rat((&self.0).shr(rhs)) } type Output = Rat; }
impl Shr<i128 > for &Rat { fn shr(self, rhs: i128 ) -> Self::Output { Rat((&self.0).shr(rhs)) } type Output = Rat; }
impl Shr<isize> for &Rat { fn shr(self, rhs: isize) -> Self::Output { Rat((&self.0).shr(rhs)) } type Output = Rat; }
impl Shr<u8   > for &Rat { fn shr(self, rhs: u8   ) -> Self::Output { Rat((&self.0).shr(rhs)) } type Output = Rat; }
impl Shr<u16  > for &Rat { fn shr(self, rhs: u16  ) -> Self::Output { Rat((&self.0).shr(rhs)) } type Output = Rat; }
impl Shr<u32  > for &Rat { fn shr(self, rhs: u32  ) -> Self::Output { Rat((&self.0).shr(rhs)) } type Output = Rat; }
impl Shr<u64  > for &Rat { fn shr(self, rhs: u64  ) -> Self::Output { Rat((&self.0).shr(rhs)) } type Output = Rat; }
impl Shr<u128 > for &Rat { fn shr(self, rhs: u128 ) -> Self::Output { Rat((&self.0).shr(rhs)) } type Output = Rat; }
impl Shr<usize> for &Rat { fn shr(self, rhs: usize) -> Self::Output { Rat((&self.0).shr(rhs)) } type Output = Rat; }
impl ShrAssign<i8   > for Rat { fn shr_assign(&mut self, rhs: i8   ) { self.0.shr_assign(rhs) } }
impl ShrAssign<i16  > for Rat { fn shr_assign(&mut self, rhs: i16  ) { self.0.shr_assign(rhs) } }
impl ShrAssign<i32  > for Rat { fn shr_assign(&mut self, rhs: i32  ) { self.0.shr_assign(rhs) } }
impl ShrAssign<i64  > for Rat { fn shr_assign(&mut self, rhs: i64  ) { self.0.shr_assign(rhs) } }
impl ShrAssign<i128 > for Rat { fn shr_assign(&mut self, rhs: i128 ) { self.0.shr_assign(rhs) } }
impl ShrAssign<isize> for Rat { fn shr_assign(&mut self, rhs: isize) { self.0.shr_assign(rhs) } }
impl ShrAssign<u8   > for Rat { fn shr_assign(&mut self, rhs: u8   ) { self.0.shr_assign(rhs) } }
impl ShrAssign<u16  > for Rat { fn shr_assign(&mut self, rhs: u16  ) { self.0.shr_assign(rhs) } }
impl ShrAssign<u32  > for Rat { fn shr_assign(&mut self, rhs: u32  ) { self.0.shr_assign(rhs) } }
impl ShrAssign<u64  > for Rat { fn shr_assign(&mut self, rhs: u64  ) { self.0.shr_assign(rhs) } }
impl ShrAssign<u128 > for Rat { fn shr_assign(&mut self, rhs: u128 ) { self.0.shr_assign(rhs) } }
impl ShrAssign<usize> for Rat { fn shr_assign(&mut self, rhs: usize) { self.0.shr_assign(rhs) } }

impl Sign for Rat { fn sign(&self) -> Ordering { self.0.sign() } }

impl SignificantBits for &Rat { fn significant_bits(self) -> u64 { self.0.significant_bits() } }

impl Square for Rat { fn square(self) -> Self::Output { Rat(self.0.square()) } type Output = Rat; }
impl Square for &Rat { fn square(self) -> Self::Output { Rat((&self.0).square()) } type Output = Rat; }
impl SquareAssign for Rat { fn square_assign(&mut self) { self.0.square_assign(); } }

impl Sub<Rat> for Rat { fn sub(self, other: Rat) -> Self::Output { Rat(self.0.sub(other.0)) } type Output = Rat;}
impl Sub<&Rat> for Rat { fn sub(self, other: &Rat) -> Self::Output { Rat(self.0.sub(&other.0)) } type Output = Rat;}
impl Sub<Rat> for &Rat { fn sub(self, other: Rat) -> Self::Output { Rat((&self.0).sub(other.0)) } type Output = Rat;}
impl Sub<&Rat> for &Rat { fn sub(self, other: &Rat) -> Self::Output { Rat((&self.0).sub(&other.0)) } type Output = Rat;}
impl SubAssign<Rat> for Rat { fn sub_assign(&mut self, other: Rat) { self.0.sub_assign(other.0); }}
impl SubAssign<&Rat> for Rat { fn sub_assign(&mut self, other: &Rat) { self.0.sub_assign(&other.0); }}

impl Sum for Rat { fn sum<I: Iterator<Item = Self>>(iter: I) -> Self { Rational::sum(iter.map(|r| r.0)).into() } }
impl<'a> Sum<&'a Rat> for Rat { fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self { Rational::sum(iter.map(|r| &r.0)).into() } }

impl TryFrom<&Rat> for i8    { fn try_from(value: &Rat) -> Result<Self, Self::Error> { i8   ::try_from(&value.0) } type Error = SignedFromRationalError; }
impl TryFrom<&Rat> for i16   { fn try_from(value: &Rat) -> Result<Self, Self::Error> { i16  ::try_from(&value.0) } type Error = SignedFromRationalError; }
impl TryFrom<&Rat> for i32   { fn try_from(value: &Rat) -> Result<Self, Self::Error> { i32  ::try_from(&value.0) } type Error = SignedFromRationalError; }
impl TryFrom<&Rat> for i64   { fn try_from(value: &Rat) -> Result<Self, Self::Error> { i64  ::try_from(&value.0) } type Error = SignedFromRationalError; }
impl TryFrom<&Rat> for i128  { fn try_from(value: &Rat) -> Result<Self, Self::Error> { i128 ::try_from(&value.0) } type Error = SignedFromRationalError; }
impl TryFrom<&Rat> for isize { fn try_from(value: &Rat) -> Result<Self, Self::Error> { isize::try_from(&value.0) } type Error = SignedFromRationalError; }
impl TryFrom<&Rat> for u8    { fn try_from(value: &Rat) -> Result<Self, Self::Error> { u8   ::try_from(&value.0) } type Error = UnsignedFromRationalError; }
impl TryFrom<&Rat> for u16   { fn try_from(value: &Rat) -> Result<Self, Self::Error> { u16  ::try_from(&value.0) } type Error = UnsignedFromRationalError; }
impl TryFrom<&Rat> for u32   { fn try_from(value: &Rat) -> Result<Self, Self::Error> { u32  ::try_from(&value.0) } type Error = UnsignedFromRationalError; }
impl TryFrom<&Rat> for u64   { fn try_from(value: &Rat) -> Result<Self, Self::Error> { u64  ::try_from(&value.0) } type Error = UnsignedFromRationalError; }
impl TryFrom<&Rat> for u128  { fn try_from(value: &Rat) -> Result<Self, Self::Error> { u128 ::try_from(&value.0) } type Error = UnsignedFromRationalError; }
impl TryFrom<&Rat> for usize { fn try_from(value: &Rat) -> Result<Self, Self::Error> { usize::try_from(&value.0) } type Error = UnsignedFromRationalError; }
impl TryFrom<&Rat> for f32 { fn try_from(value: &Rat) -> Result<Self, Self::Error> { f32::try_from(&value.0) } type Error = FloatConversionError; }
impl TryFrom<&Rat> for f64 { fn try_from(value: &Rat) -> Result<Self, Self::Error> { f64::try_from(&value.0) } type Error = FloatConversionError; }
impl TryFrom<&Rat> for Integer { fn try_from(value: &Rat) -> Result<Self, Self::Error> { Integer::try_from(&value.0) } type Error = IntegerFromRationalError; }
impl TryFrom<&Rat> for Natural { fn try_from(value: &Rat) -> Result<Self, Self::Error> { Natural::try_from(&value.0) } type Error = NaturalFromRationalError; }
impl TryFrom<Rat> for f32 { fn try_from(value: Rat) -> Result<Self, Self::Error> { f32::try_from(value.0) } type Error = FloatConversionError; }
impl TryFrom<Rat> for f64 { fn try_from(value: Rat) -> Result<Self, Self::Error> { f64::try_from(value.0) } type Error = FloatConversionError; }
impl TryFrom<Rat> for Integer { fn try_from(value: Rat) -> Result<Self, Self::Error> { Integer::try_from(value.0) } type Error = IntegerFromRationalError; }
impl TryFrom<Rat> for Natural { fn try_from(value: Rat) -> Result<Self, Self::Error> { Natural::try_from(value.0) } type Error = NaturalFromRationalError; }
impl TryFrom<f32> for Rat { fn try_from(value: f32) -> Result<Self, Self::Error> { Rational::try_from(value).map(Rat) } type Error = RationalFromPrimitiveFloatError; }
impl TryFrom<f64> for Rat { fn try_from(value: f64) -> Result<Self, Self::Error> { Rational::try_from(value).map(Rat) } type Error = RationalFromPrimitiveFloatError; }

impl SciMantissaAndExponent<f32, i64, Rat> for &Rat {
    fn from_sci_mantissa_and_exponent(sci_mantissa: f32, sci_exponent: i64) -> Option<Rat> {
        Rational::from_sci_mantissa_and_exponent(sci_mantissa, sci_exponent).map(Rat)
    }
    fn sci_mantissa_and_exponent(self) -> (f32, i64) {
        (&self.0).sci_mantissa_and_exponent()
    }
}

impl SciMantissaAndExponent<f32, i64, Rat> for Rat {
    fn from_sci_mantissa_and_exponent(sci_mantissa: f32, sci_exponent: i64) -> Option<Rat> {
        Rational::from_sci_mantissa_and_exponent(sci_mantissa, sci_exponent).map(Rat)
    }
    fn sci_mantissa_and_exponent(self) -> (f32, i64) {
        self.0.sci_mantissa_and_exponent()
    }
}

impl SciMantissaAndExponent<f64, i64, Rat> for &Rat {
    fn from_sci_mantissa_and_exponent(sci_mantissa: f64, sci_exponent: i64) -> Option<Rat> {
        Rational::from_sci_mantissa_and_exponent(sci_mantissa, sci_exponent).map(Rat)
    }
    fn sci_mantissa_and_exponent(self) -> (f64, i64) {
        (&self.0).sci_mantissa_and_exponent()
    }
}

impl SciMantissaAndExponent<f64, i64, Rat> for Rat {
    fn from_sci_mantissa_and_exponent(sci_mantissa: f64, sci_exponent: i64) -> Option<Rat> {
        Rational::from_sci_mantissa_and_exponent(sci_mantissa, sci_exponent).map(Rat)
    }
    fn sci_mantissa_and_exponent(self) -> (f64, i64) {
        self.0.sci_mantissa_and_exponent()
    }
}

impl ToSci for Rat {
    fn fmt_sci(&self, f: &mut std::fmt::Formatter, options: malachite::base::num::conversion::string::options::ToSciOptions) -> std::fmt::Result {
        self.0.fmt_sci(f, options)
    }

    fn fmt_sci_valid(&self, options: malachite::base::num::conversion::string::options::ToSciOptions) -> bool {
        self.0.fmt_sci_valid(options)
    }
}

impl Zero for Rat {
    fn zero() -> Self {
        Self(Rational::ZERO)
    }

    fn is_zero(&self) -> bool {
        self == &Self::zero()
    }
}

impl One for Rat {
    fn one() -> Self {
        Self(Rational::ONE)
    }

    fn is_one(&self) -> bool {
        self == &Self::one()
    }
}

impl RemAssign<Self> for Rat {
    fn rem_assign(&mut self, rhs: Self) {
        self.0 -= (&self.0).round_to_multiple(rhs.0, RoundingMode::Floor).0
    }
}

impl RemAssign<&Self> for Rat {
    fn rem_assign(&mut self, rhs: &Self) {
        self.0 -= (&self.0).round_to_multiple(&rhs.0, RoundingMode::Floor).0
    }
}

impl Rem<Rat> for Rat {
    type Output = Rat;

    fn rem(self, rhs: Rat) -> Self::Output {
        Rat(&self.0 - (&self.0).round_to_multiple(rhs.0, RoundingMode::Floor).0)
    }
}

impl Rem<&Rat> for Rat {
    type Output = Rat;

    fn rem(self, rhs: &Rat) -> Self::Output {
        Rat(&self.0 - (&self.0).round_to_multiple(&rhs.0, RoundingMode::Floor).0)
    }
}

impl<'a> Rem<Rat> for &'a Rat {
    type Output = Rat;

    fn rem(self, rhs: Rat) -> Self::Output {
        Rat(&self.0 - (&self.0).round_to_multiple(rhs.0, RoundingMode::Floor).0)
    }
}

impl<'a> Rem<&Rat> for &'a Rat {
    type Output = Rat;

    fn rem(self, rhs: &Rat) -> Self::Output {
        Rat(&self.0 - (&self.0).round_to_multiple(&rhs.0, RoundingMode::Floor).0)
    }
}

impl Num for Rat {
    type FromStrRadixErr = ();

    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        let mut options = FromSciStringOptions::default();
        options.set_base(radix as u8);
        Rational::from_sci_string_simplest_with_options(str, options).ok_or(()).map(|r| Rat(r))
    }
}
