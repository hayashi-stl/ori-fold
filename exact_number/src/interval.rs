// I would use inari, but because of lack of support for rounding modes, it doesn't work on the web.

use std::ops::{Add, AddAssign, Sub, SubAssign};

/// Closed conservative interval between two real numbers.
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

    fn empty_if_overflowed(bounds: [T; 2]) -> Option<[T; 2]> {
        if T::overflowed(&bounds) { None } else { Some(bounds) }
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

impl<T: IntervalContent + Clone> From<T> for Interval<T> {
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
pub trait IntervalContent {
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

    fn overflowed(value: &[Self; 2]) -> bool where Self: Sized;
}


// Does not assume a rounding mode, except that a number is rounded to one of the two closest
// f64s to it.
impl IntervalContent for f64 {
    // Prioritizing speed over bound tightness here.
    // It's probably rare that the difference will matter.
    fn add_down_vv(self, rhs: Self) -> Self {
        f64_down(self + rhs)
    }

    fn add_down_vr(self, rhs: &Self) -> Self {
        self.add_down_vv(*rhs)
    }

    fn add_down_rv(&self, rhs: Self) -> Self {
        self.add_down_vv(rhs)
    }

    fn add_down_rr(&self, rhs: &Self) -> Self {
        self.add_down_vv(*rhs)
    }

    fn add_up_vv(self, rhs: Self) -> Self {
        f64_up(self + rhs)
    }

    fn add_up_vr(self, rhs: &Self) -> Self {
        self.add_up_vv(*rhs)
    }

    fn add_up_rv(&self, rhs: Self) -> Self {
        self.add_up_vv(rhs)
    }

    fn add_up_rr(&self, rhs: &Self) -> Self {
        self.add_up_vv(*rhs)
    }

    fn sub_down_vv(self, rhs: Self) -> Self {
        f64_down(self - rhs)
    }

    fn sub_down_vr(self, rhs: &Self) -> Self {
        self.sub_down_vv(*rhs)
    }

    fn sub_down_rv(&self, rhs: Self) -> Self {
        self.sub_down_vv(rhs)
    }

    fn sub_down_rr(&self, rhs: &Self) -> Self {
        self.sub_down_vv(*rhs)
    }

    fn sub_up_vv(self, rhs: Self) -> Self {
        f64_up(self - rhs)
    }

    fn sub_up_vr(self, rhs: &Self) -> Self {
        self.sub_up_vv(*rhs)
    }

    fn sub_up_rv(&self, rhs: Self) -> Self {
        self.sub_up_vv(rhs)
    }

    fn sub_up_rr(&self, rhs: &Self) -> Self {
        self.sub_up_vv(*rhs)
    }

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

#[cfg(test)]
mod test {
    use crate::interval::{Interval, IntervalContent};

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

    fn n2i(a: f64, b: f64) -> Interval<f64> {
        Interval::new([a, b])
    }

    // libieeep1788 tests, but only applicable ones, and modified to account for rounding
    #[test]
    fn test_interval_add_minimal() {
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
    fn minimal_sub_test() {
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
}