
use std::sync::LazyLock;

use malachite::{base::{num::conversion::traits::{IsInteger, RoundingFrom}, rounding_modes::RoundingMode}, Integer, Natural};
use nalgebra::DVector;
use num::Num;

use crate::{interval::Interval, rat::Rat, BasedExpr};

impl BasedExpr {
    /// Rounds this expression to the integer type with the given rounding mode.
    pub fn round_to_integer(&self, mode: RoundingMode) -> Integer {
        match self {
            BasedExpr::Baseless(a) => Integer::rounding_from(a, mode).0,
            BasedExpr::Based(a, basis) => {
                if a.iter().skip(1).all(|a| a == &Rat::ZERO) {
                    Integer::rounding_from(&a[0], mode).0.into()
                } else {
                    // Not an integer
                    Self::narrow_down(&a, &basis, Some(|k: Interval<f64>| k.definite_round(mode)), |k| k.definite_round(mode)).into()
                }
            }
        }
    }

    /// Rounds this expression to the f64 with the given rounding mode.
    pub fn round_to_f64(&self, mode: RoundingMode) -> f64 {
        match self {
            BasedExpr::Baseless(a) => f64::rounding_from(a, mode).0,
            BasedExpr::Based(a, basis) => {
                if a.iter().skip(1).all(|a| a == &Rat::ZERO) {
                    f64::rounding_from(&a[0], mode).0
                } else {
                    // Not an exact f64
                    Self::narrow_down(&a, &basis,
                            None::<fn(Interval<f64>) -> Option<f64>>, |k| k.definite_round_into_f64(mode))
                }
            }
        }
    }

    pub fn is_integer(&self) -> bool {
        match self {
            BasedExpr::Baseless(a) => a.is_integer(),
            BasedExpr::Based(a, _) => {
                a[0].is_integer() && a.iter().skip(1).all(|a| a == &Rat::ZERO)
            }
        }
    }

    pub fn try_to_integer(&self) -> Option<Integer> {
        match self {
            BasedExpr::Baseless(a) => a.try_into().ok(),
            BasedExpr::Based(a, _basis) => {
                if a.iter().skip(1).all(|a| a == &Rat::ZERO) {
                    (&a[0]).try_into().ok()
                } else {
                    None
                }
            }
        }
    }

    pub fn try_into_integer(self) -> Option<Integer> {
        match self {
            BasedExpr::Baseless(a) => a.try_into().ok(),
            BasedExpr::Based(mut a, _basis) => {
                if a.iter().skip(1).all(|a| a == &Rat::ZERO) {
                    std::mem::take(&mut a[0]).try_into().ok()
                } else {
                    None
                }
            }
        }
    }

    /// Rounds this expression to an integer with the given rounding mode,
    /// while preserving the basis if it exists.
    pub fn round(&self, mode: RoundingMode) -> Self {
        let int = self.round_to_integer(mode);
        match self {
            BasedExpr::Baseless(_) => BasedExpr::Baseless(int.into()),
            BasedExpr::Based(a, basis) => {
                let mut coeffs = DVector::repeat(a.len(), Rat::ZERO);
                coeffs[0] = int.into();
                BasedExpr::Based(coeffs, basis.clone())
            }
        }
    }


}

/// Converts a number into a baseless BasedExpr.
/// This trait forces you to think about the fact that you don't have a basis yet.
pub trait BaselessFrom<T> {
    /// Converts `t` to a *baseless* `BasedExpr
    fn baseless_from(t: T) -> Self;
}

pub trait IntoBaseless where
    Self: Sized,
    BasedExpr: BaselessFrom<Self>
{
    /// Converts `self` to a *baseless* `BasedExpr
    fn into_baseless(self) -> BasedExpr;
}

impl<T> IntoBaseless for T where BasedExpr: BaselessFrom<T> {
    fn into_baseless(self) -> BasedExpr {
        BasedExpr::baseless_from(self)
    }
}

/// Tries to convert a number into a baseless BasedExpr.
/// This trait forces you to think about the fact that you don't have a basis yet.
pub trait TryBaselessFrom<T>: Sized {
    /// The type returned in the event of a conversion error.
    type Error;
    /// Performs the conversion.
    fn try_baseless_from(t: T) -> Result<Self, Self::Error>;
}

impl<T> TryBaselessFrom<T> for BasedExpr where BasedExpr: BaselessFrom<T> {
    type Error = ();

    fn try_baseless_from(t: T) -> Result<Self, Self::Error> {
        Ok(t.into_baseless())
    }
}

/// Tries to convert a number into a baseless BasedExpr.
/// This trait forces you to think about the fact that you don't have a basis yet.
pub trait TryIntoBaseless where
    Self: Sized,
    BasedExpr: TryBaselessFrom<Self>
{
    /// The type returned in the event of a conversion error.
    type Error;
    /// Performs the conversion.
    fn try_into_baseless(t: Self) -> Result<BasedExpr, Self::Error>;
}

impl<T> TryIntoBaseless for T where BasedExpr: TryBaselessFrom<T> {
    type Error = <BasedExpr as TryBaselessFrom<T>>::Error;

    fn try_into_baseless(t: Self) -> Result<BasedExpr, Self::Error> {
        BasedExpr::try_baseless_from(t)
    }
}

impl BaselessFrom<bool > for BasedExpr { fn baseless_from(t: bool ) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<i8   > for BasedExpr { fn baseless_from(t: i8   ) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<i16  > for BasedExpr { fn baseless_from(t: i16  ) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<i32  > for BasedExpr { fn baseless_from(t: i32  ) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<i64  > for BasedExpr { fn baseless_from(t: i64  ) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<i128 > for BasedExpr { fn baseless_from(t: i128 ) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<isize> for BasedExpr { fn baseless_from(t: isize) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<u8   > for BasedExpr { fn baseless_from(t: u8   ) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<u16  > for BasedExpr { fn baseless_from(t: u16  ) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<u32  > for BasedExpr { fn baseless_from(t: u32  ) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<u64  > for BasedExpr { fn baseless_from(t: u64  ) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<u128 > for BasedExpr { fn baseless_from(t: u128 ) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<usize> for BasedExpr { fn baseless_from(t: usize) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<Natural> for BasedExpr { fn baseless_from(t: Natural) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<Integer> for BasedExpr { fn baseless_from(t: Integer) -> BasedExpr { Self::Baseless(t.into()) } }
impl BaselessFrom<Rat>     for BasedExpr { fn baseless_from(t: Rat    ) -> BasedExpr { Self::Baseless(t) } }

impl TryBaselessFrom<f32> for BasedExpr {
    type Error = (); // lazy for now
    fn try_baseless_from(t: f32) -> Result<Self, Self::Error> { Rat::try_from(t).map_err(|_e| ()).map(Self::Baseless) }
}
impl TryBaselessFrom<f64> for BasedExpr {
    type Error = (); // lazy for now
    fn try_baseless_from(t: f64) -> Result<Self, Self::Error> { Rat::try_from(t).map_err(|_e| ()).map(Self::Baseless) }
}

impl TryFrom<&BasedExpr> for i8    { fn try_from(t: &BasedExpr) -> Result<i8   , ()>{ t.try_to_integer().and_then(|i| (&i).try_into().ok()).ok_or(()) } type Error = (); }
impl TryFrom<&BasedExpr> for i16   { fn try_from(t: &BasedExpr) -> Result<i16  , ()>{ t.try_to_integer().and_then(|i| (&i).try_into().ok()).ok_or(()) } type Error = (); }
impl TryFrom<&BasedExpr> for i32   { fn try_from(t: &BasedExpr) -> Result<i32  , ()>{ t.try_to_integer().and_then(|i| (&i).try_into().ok()).ok_or(()) } type Error = (); }
impl TryFrom<&BasedExpr> for i64   { fn try_from(t: &BasedExpr) -> Result<i64  , ()>{ t.try_to_integer().and_then(|i| (&i).try_into().ok()).ok_or(()) } type Error = (); }
impl TryFrom<&BasedExpr> for i128  { fn try_from(t: &BasedExpr) -> Result<i128 , ()>{ t.try_to_integer().and_then(|i| (&i).try_into().ok()).ok_or(()) } type Error = (); }
impl TryFrom<&BasedExpr> for isize { fn try_from(t: &BasedExpr) -> Result<isize, ()>{ t.try_to_integer().and_then(|i| (&i).try_into().ok()).ok_or(()) } type Error = (); }
impl TryFrom<&BasedExpr> for u8    { fn try_from(t: &BasedExpr) -> Result<u8   , ()>{ t.try_to_integer().and_then(|i| (&i).try_into().ok()).ok_or(()) } type Error = (); }
impl TryFrom<&BasedExpr> for u16   { fn try_from(t: &BasedExpr) -> Result<u16  , ()>{ t.try_to_integer().and_then(|i| (&i).try_into().ok()).ok_or(()) } type Error = (); }
impl TryFrom<&BasedExpr> for u32   { fn try_from(t: &BasedExpr) -> Result<u32  , ()>{ t.try_to_integer().and_then(|i| (&i).try_into().ok()).ok_or(()) } type Error = (); }
impl TryFrom<&BasedExpr> for u64   { fn try_from(t: &BasedExpr) -> Result<u64  , ()>{ t.try_to_integer().and_then(|i| (&i).try_into().ok()).ok_or(()) } type Error = (); }
impl TryFrom<&BasedExpr> for u128  { fn try_from(t: &BasedExpr) -> Result<u128 , ()>{ t.try_to_integer().and_then(|i| (&i).try_into().ok()).ok_or(()) } type Error = (); }
impl TryFrom<&BasedExpr> for usize { fn try_from(t: &BasedExpr) -> Result<usize, ()>{ t.try_to_integer().and_then(|i| (&i).try_into().ok()).ok_or(()) } type Error = (); }
impl TryFrom<&BasedExpr> for Natural { fn try_from(t: &BasedExpr) -> Result<Natural, ()>{ t.try_to_integer().and_then(|i| (&i).try_into().ok()).ok_or(()) } type Error = (); }
impl TryFrom<&BasedExpr> for Integer { fn try_from(t: &BasedExpr) -> Result<Integer, ()>{ t.try_to_integer().ok_or(()) } type Error = (); }
impl TryFrom<BasedExpr> for Natural { fn try_from(t: BasedExpr) -> Result<Natural, ()>{ t.try_into_integer().and_then(|i| (&i).try_into().ok()).ok_or(()) } type Error = (); }
impl TryFrom<BasedExpr> for Integer { fn try_from(t: BasedExpr) -> Result<Integer, ()>{ t.try_into_integer().ok_or(()) } type Error = (); }

// Oof, not a regular expression
//static REGEX: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"(?x)
//(
//    \s*
//    (?P<sign> [+-])\s*
//    (?P<numerator> [0-9A-Za-z_]+)\s*
//    (/ \s* (?P<denominator> [0-9A-Za-z_]+)\s*)?
//    ((√ | \s sqrt))?
//)
//").unwrap());

#[cfg(test)]
mod test {
    use malachite::Integer;

    use crate::{based_expr, rat::Rat, BasedExpr};

    fn integer_test(a: BasedExpr, expected: Option<Integer>) {
        assert_eq!(a.is_integer(), expected.is_some());
        assert_eq!(a.try_to_integer(), expected);
        assert_eq!(a.try_into_integer(), expected);
    }

    #[test]
    fn test_based_expr_integer() {
        let _lock = crate::TEST_LOCK.read();

        integer_test(BasedExpr::Baseless(3.into()), Some(3.into()));
        integer_test(BasedExpr::Baseless(Rat::from_signeds(-1, 3)), None);

        integer_test(based_expr!(-3), Some((-3i64).into()));
        integer_test(based_expr!(-1/3), None);
        integer_test(based_expr!(100 + 0 sqrt 2), Some(100.into()));
        integer_test(based_expr!(100 + 1 sqrt 2), None);
        integer_test(based_expr!(5/6 + 0 sqrt 2), None);
    }
}