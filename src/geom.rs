use std::{cmp::Ordering, iter::Sum, mem, ops::{Mul, Neg, Sub}};

use duplicate::duplicate_item;
use exact_number::{Angle, BasedExpr, angle::IntoAngle, malachite::base::num::arithmetic::traits::Sign};
use nalgebra::{Affine2, ClosedSubAssign, Const, DefaultAllocator, Dim, DimNameAdd, DimNameSum, Dyn, Matrix2, MatrixView2xX, OVector, RealField, SVector, Scalar, Storage, TAffine, Transform, U1, U2, Vector, Vector2, VectorView, VectorView2, allocator::Allocator, matrix, vector};
use num_traits::{Num, NumAssign, NumAssignRef, NumRef, RefNum, Signed};
use robust_geometry as robust;

pub type VectorViewDyn<'s, T, D> = VectorView<'s, T, D, U1, Dyn>;
pub type VectorView2Dyn<'s, T> = VectorView2<'s, T, U1, Dyn>;
pub type MatrixView2Dyn<'s, T> = MatrixView2xX<'s, T, U1, Dyn>;

pub trait NumEx: Default + PartialOrd + Num + NumRef + NumAssignRef + Neg + Scalar + NumAssign + Signed + Neg<Output = Self> {}
impl<T> NumEx for T where T: Default + PartialOrd + Num + NumRef + NumAssignRef + Neg + Scalar + NumAssign + Signed + Neg<Output = Self> {}

/// A wrapper for floats, that implements total equality and ordering
/// and hashing.
/// 
/// Taken from https://docs.rs/float-ord/0.3.2/src/float_ord/lib.rs.html,
/// and modified to have From and Into implementations and to compare 0.0 equal to -0.0.
#[derive(Clone, Copy, Debug)]
#[repr(transparent)]
pub struct FloatOrd<T>(pub T);

macro_rules! float_ord_impl {
    ($f:ident, $i:ident, $n:expr) => {
        impl FloatOrd<$f> {
            fn convert(self) -> $i {
                let u = self.0.to_bits();
                let bit = 1 << ($n - 1);
                if u & bit == 0 {
                    u | bit
                } else {
                    // Make sure that 0 == -0
                    if u == bit { u } else { !u }
                }
            }
        }
        impl PartialEq for FloatOrd<$f> {
            fn eq(&self, other: &Self) -> bool {
                self.convert() == other.convert()
            }
        }
        impl Eq for FloatOrd<$f> {}
        impl PartialOrd for FloatOrd<$f> {
            fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
                self.convert().partial_cmp(&other.convert())
            }
        }
        impl Ord for FloatOrd<$f> {
            fn cmp(&self, other: &Self) -> Ordering {
                self.convert().cmp(&other.convert())
            }
        }
        impl std::hash::Hash for FloatOrd<$f> {
            fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
                self.convert().hash(state);
            }
        }

        impl From<$f> for FloatOrd<$f> {
            fn from(val: $f) -> Self {
                Self(val)
            }
        }

        impl From<FloatOrd<$f>> for $f {
            fn from(val: FloatOrd<$f>) -> Self {
                val.0
            }
        }
    }
}

float_ord_impl!(f32, u32, 32);
float_ord_impl!(f64, u64, 64);

/// Stores an angle as a directed line segment for robust comparison with other angles.
/// 
/// Angles compare such that counterclockwise order = increasing order.
/// There is a branch cut in the direction [-1, ε] for infinitesimal ε.
#[derive(Clone, Copy, Debug)]
pub struct AngleF64(pub [Vector2<f64>; 2]);

impl AngleF64 {
    fn nonzero_equiv(self) -> Self {
        if self.0[0] == self.0[1] {
            AngleF64([vector![0.0, 0.0], vector![1.0, 0.0]])
        } else {
            self
        }
    }

    fn yx_cmp(self) -> Ordering {
        [self.0[1].y, self.0[1].x].map(FloatOrd).cmp(&[self.0[0].y, self.0[0].x].map(FloatOrd))
    }

    /// Compares two angles that are known to be parallel/antiparallel to each other.
    /// using a specific coordinate.
    fn cmp_parallel(self, other: Self) -> Ordering {
        use Ordering::*;
        match (self.yx_cmp(), other.yx_cmp()) {
            (Less, Less) | (Greater, Greater) => Equal,
            (Less, Greater) => Less,
            (Greater, Less) => Greater,
            _ => panic!("the vectors aren't parallel"),
        }
    }

    /// Returns False if this is to the right of the branch cut
    /// and True if it's to the left.
    fn is_branch_cut_left(self) -> bool {
        self.yx_cmp().is_ge()
    }
}

impl PartialEq for AngleF64 {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}

impl Eq for AngleF64 {}

impl PartialOrd for AngleF64 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for AngleF64 {
    fn cmp(&self, other: &Self) -> Ordering {
        let mut a = self.nonzero_equiv();
        let mut b = other.nonzero_equiv();
        let mut orient = robust::cross_2d(b.0[0], b.0[1], a.0[0], a.0[1]);
        if orient == 0.0 {
            // angles are parallel/antiparallel
            a.cmp_parallel(b)
        } else {
            if orient > 0.0 { mem::swap(&mut a, &mut b); }
            if a.is_branch_cut_left() && !b.is_branch_cut_left() { orient = -orient; }
            FloatOrd(orient).cmp(&FloatOrd(0.0))
        }
    }
}

pub trait IntoOrdAngle<Y = Self>: {
    type LVec<'a>;
    type RVec<'a>;
    type Output: Ord;
    /// Finds the angle of the vector `[self, y]`,
    /// with the range [-π, π), returning `None` if `self` = `y` = 0.
    /// 
    /// Like `atan2`, but with x first and y second.
    fn into_ord_angle(self, y: Y) -> Option<Self::Output>;
}

#[duplicate_item(
    f32_a  f32_b  vec_f32_a                 vec_f32_b                 val_a  val_b;
    [ f32] [ f32] [           Vector2<f32>] [           Vector2<f32>] [ self] [ y];
    [ f32] [&f32] [           Vector2<f32>] [VectorView2Dyn<'a, f32>] [ self] [*y];
    [&f32] [ f32] [VectorView2Dyn<'a, f32>] [           Vector2<f32>] [*self] [ y];
    [&f32] [&f32] [VectorView2Dyn<'a, f32>] [VectorView2Dyn<'a, f32>] [*self] [*y];
)]
impl IntoOrdAngle<f32_b> for f32_a {
    type LVec<'a> = vec_f32_a;
    type RVec<'a> = vec_f32_b;
    type Output = AngleF64;
    fn into_ord_angle(self, y: f32_b) -> Option<Self::Output> {
        (val_a as f64).into_ord_angle(val_b as f64)
    }
}

#[duplicate_item(
    f64_a  f64_b  vec_f64_a                 vec_f64_b                 val_a  val_b;
    [ f64] [ f64] [           Vector2<f64>] [           Vector2<f64>] [ self] [ y];
    [ f64] [&f64] [           Vector2<f64>] [VectorView2Dyn<'a, f64>] [ self] [*y];
    [&f64] [ f64] [VectorView2Dyn<'a, f64>] [           Vector2<f64>] [*self] [ y];
    [&f64] [&f64] [VectorView2Dyn<'a, f64>] [VectorView2Dyn<'a, f64>] [*self] [*y];
)]
impl IntoOrdAngle<f64_b> for f64_a {
    type LVec<'a> = vec_f64_a;
    type RVec<'a> = vec_f64_b;
    type Output = AngleF64;
    fn into_ord_angle(self, y: f64_b) -> Option<Self::Output> {
        if val_a == 0.0 && val_b == 0.0 { None } else { Some(AngleF64([Vector2::zeros(), vector![val_a, val_b]])) }
    }
}

#[duplicate_item(
    ty_a         ty_b         vec_ty_a                        vec_ty_b;
    [ BasedExpr] [ BasedExpr] [           Vector2<BasedExpr>] [           Vector2<BasedExpr>];
    [ BasedExpr] [&BasedExpr] [           Vector2<BasedExpr>] [VectorView2Dyn<'a, BasedExpr>];
    [&BasedExpr] [ BasedExpr] [VectorView2Dyn<'a, BasedExpr>] [           Vector2<BasedExpr>];
    [&BasedExpr] [&BasedExpr] [VectorView2Dyn<'a, BasedExpr>] [VectorView2Dyn<'a, BasedExpr>];
)]
impl IntoOrdAngle<ty_b> for ty_a {
    type LVec<'a> = vec_ty_a;
    type RVec<'a> = vec_ty_b;
    type Output = Angle;
    fn into_ord_angle(self, y: ty_b) -> Option<Self::Output> {
        self.into_angle(y)
    }
}

pub trait IntoOrdAngleOp: Sized
    + for<'b, 'c> IntoOrdAngle<Self, LVec<'b> = Vector2<Self>, RVec<'c> = Vector2<Self>>
    + for<'a, 'b, 'c> IntoOrdAngle<&'a Self, Output = <Self as IntoOrdAngle>::Output, LVec<'b> = Vector2<Self>, RVec<'c> = VectorView2Dyn<'c, Self>> {}
impl<T: Sized
    + for<'b, 'c> IntoOrdAngle<Self, LVec<'b> = Vector2<Self>, RVec<'c> = Vector2<Self>>
    + for<'a, 'b, 'c> IntoOrdAngle<&'a Self, Output = <Self as IntoOrdAngle>::Output, LVec<'b> = Vector2<Self>, RVec<'c> = VectorView2Dyn<'c, Self>>>
    IntoOrdAngleOp for T {}

pub trait RefIntoOrdAngleOp<Base: IntoOrdAngleOp>:
    for<'b, 'c> IntoOrdAngle<Base, Output = <Base as IntoOrdAngle>::Output, LVec<'b> = VectorView2Dyn<'b, Base>, RVec<'c> = Vector2<Base>> +
    for<'a, 'b, 'c> IntoOrdAngle<&'a Base, Output = <Base as IntoOrdAngle>::Output, LVec<'b> = VectorView2Dyn<'b, Base>, RVec<'c> = VectorView2Dyn<'c, Base>> {}
impl<Base: IntoOrdAngleOp, T:
    for<'b, 'c> IntoOrdAngle<Base, Output = <Base as IntoOrdAngle>::Output, LVec<'b> = VectorView2Dyn<'b, Base>, RVec<'c> = Vector2<Base>> +
    for<'a, 'b, 'c> IntoOrdAngle<&'a Base, Output = <Base as IntoOrdAngle>::Output, LVec<'b> = VectorView2Dyn<'b, Base>, RVec<'c> = VectorView2Dyn<'c, Base>>>
    RefIntoOrdAngleOp<Base> for T {}

/// Operations involving angles in some sense, like 2D determinant and conversion into an angle
pub trait AngleOps: NumEx + IntoOrdAngle {
    type Det: Ord;
    type OrdAngle: Ord;
    /// Finds the angle from vector `v0` to vector `v1`.
    /// with the range [-π, π), returning `None` if `v0` = `v1`.
    /// 
    /// This gives an angle that compares equal to `(other - self).x.into_ord_angle(other - self).y`.
    fn ord_angle_to<S0: Storage<Self, U2>, S1: Storage<Self, U2>>(v0: Vector<Self, U2, S0>, v1: Vector<Self, U2, S1>) -> Option<Self::OrdAngle>;

    /// Computes the following determinant:
    /// ```notrust
    /// | b.x - a.x    d.x - c.x |
    /// | b.y - a.y    d.y - c.y |
    /// ```
    fn det_2d<Sa, Sb, Sc, Sd>(a: Vector<Self, U2, Sa>, b: Vector<Self, U2, Sb>, c: Vector<Self, U2, Sc>, d: Vector<Self, U2, Sd>) -> Self::Det where
            Sa: Storage<Self, U2>,
            Sb: Storage<Self, U2>,
            Sc: Storage<Self, U2>,
            Sd: Storage<Self, U2>,;

    /// Computes the following determinant, comparing it to 0:
    /// ```notrust
    /// | b.x - a.x    d.x - c.x |
    /// | b.y - a.y    d.y - c.y |
    /// ```
    fn cmp_det_2d<Sa, Sb, Sc, Sd>(a: Vector<Self, U2, Sa>, b: Vector<Self, U2, Sb>, c: Vector<Self, U2, Sc>, d: Vector<Self, U2, Sd>) -> Ordering where
            Sa: Storage<Self, U2>,
            Sb: Storage<Self, U2>,
            Sc: Storage<Self, U2>,
            Sd: Storage<Self, U2>,;
}

impl AngleOps for f32 {
    type OrdAngle = AngleF64;
    type Det = FloatOrd<f64>;

    fn ord_angle_to<S0: Storage<Self, U2>, S1: Storage<Self, U2>>(v0: Vector<Self, U2, S0>, v1: Vector<Self, U2, S1>) -> Option<Self::OrdAngle> {
        f64::ord_angle_to(v0.map(|c| c as f64), v1.map(|c| c as f64))
    }

    fn det_2d<Sa, Sb, Sc, Sd>(a: Vector<Self, U2, Sa>, b: Vector<Self, U2, Sb>, c: Vector<Self, U2, Sc>, d: Vector<Self, U2, Sd>) -> Self::Det where
                Sa: Storage<Self, U2>,
                Sb: Storage<Self, U2>,
                Sc: Storage<Self, U2>,
                Sd: Storage<Self, U2>, {
        f64::det_2d(a.map(|x| x as f64), b.map(|x| x as f64), c.map(|x| x as f64), d.map(|x| x as f64))
    }

    fn cmp_det_2d<Sa, Sb, Sc, Sd>(a: Vector<Self, U2, Sa>, b: Vector<Self, U2, Sb>, c: Vector<Self, U2, Sc>, d: Vector<Self, U2, Sd>) -> Ordering where
                Sa: Storage<Self, U2>,
                Sb: Storage<Self, U2>,
                Sc: Storage<Self, U2>,
                Sd: Storage<Self, U2>, {
        Self::det_2d(a, b, c, d).cmp(&FloatOrd(0.0))
    }
}

impl AngleOps for f64 {
    type OrdAngle = AngleF64;
    type Det = FloatOrd<f64>;

    fn ord_angle_to<S0: Storage<Self, U2>, S1: Storage<Self, U2>>(v0: Vector<Self, U2, S0>, v1: Vector<Self, U2, S1>) -> Option<Self::OrdAngle> {
        if v0 == v1 { None } else { Some(AngleF64([v0.into_owned(), v1.into_owned()])) }
    }

    fn det_2d<Sa, Sb, Sc, Sd>(a: Vector<Self, U2, Sa>, b: Vector<Self, U2, Sb>, c: Vector<Self, U2, Sc>, d: Vector<Self, U2, Sd>) -> Self::Det where
                Sa: Storage<Self, U2>,
                Sb: Storage<Self, U2>,
                Sc: Storage<Self, U2>,
                Sd: Storage<Self, U2>, {
        FloatOrd(robust::cross_2d(a.into_owned(), b.into_owned(), c.into_owned(), d.into_owned()))
    }

    fn cmp_det_2d<Sa, Sb, Sc, Sd>(a: Vector<Self, U2, Sa>, b: Vector<Self, U2, Sb>, c: Vector<Self, U2, Sc>, d: Vector<Self, U2, Sd>) -> Ordering where
                Sa: Storage<Self, U2>,
                Sb: Storage<Self, U2>,
                Sc: Storage<Self, U2>,
                Sd: Storage<Self, U2>, {
        Self::det_2d(a, b, c, d).cmp(&FloatOrd(0.0))
    }
}

impl AngleOps for BasedExpr {
    type OrdAngle = Angle;
    type Det = BasedExpr;

    fn ord_angle_to<S0: Storage<Self, U2>, S1: Storage<Self, U2>>(v0: Vector<Self, U2, S0>, v1: Vector<Self, U2, S1>) -> Option<Self::OrdAngle> {
        let [[x, y]] = (v1 - v0).data.0;
        y.into_ord_angle(x)
    }

    fn det_2d<Sa, Sb, Sc, Sd>(a: Vector<Self, U2, Sa>, b: Vector<Self, U2, Sb>, c: Vector<Self, U2, Sc>, d: Vector<Self, U2, Sd>) -> Self::Det where
                Sa: Storage<Self, U2>,
                Sb: Storage<Self, U2>,
                Sc: Storage<Self, U2>,
                Sd: Storage<Self, U2>, {
        (b - a).perp(&(d - c))
    }

    fn cmp_det_2d<Sa, Sb, Sc, Sd>(a: Vector<Self, U2, Sa>, b: Vector<Self, U2, Sb>, c: Vector<Self, U2, Sc>, d: Vector<Self, U2, Sd>) -> Ordering where
                Sa: Storage<Self, U2>,
                Sb: Storage<Self, U2>,
                Sc: Storage<Self, U2>,
                Sd: Storage<Self, U2>, {
        Self::det_2d(a, b, c, d).cmp_zero()
    }
}

/// Converts a type into an orderable type.
/// Unfortunately only implemented for a few types due to lack of the understanding
/// that `f64` does not and will never implement `Ord`.
pub trait IntoOrd: Sized {
    type Output: Ord;
    type OutputRef<'a>: Ord where Self: 'a;
    /// Converts this to an ordered type.
    fn into_ord(self) -> Self::Output;
    /// Converts this to an ordered type by reference.
    /// This is used instead of just implementing `IntoOrd` for reference types
    /// because doing so would propagate reference trait bounds.
    fn to_ord_ref(&self) -> Self::OutputRef<'_>;
    /// Converts this slice to an ordered type by reference
    fn array_to_ord<const N: usize>(array: [Self; N]) -> [Self::Output; N];
    /// Converts this slice to an ordered type by reference
    fn slice_to_ord(slice: &[Self]) -> &[Self::Output];
}

#[duplicate_item(
    f_ty; [f32]; [f64];
)]
impl IntoOrd for f_ty {
    type Output = FloatOrd<f_ty>;
    type OutputRef<'a> = FloatOrd<f_ty>;
    fn into_ord(self) -> Self::Output { FloatOrd(self) }
    fn to_ord_ref(&self) -> Self::OutputRef<'_> { FloatOrd(*self) }
    fn array_to_ord<const N: usize>(array: [Self; N]) -> [Self::Output; N] {
        array.map(Self::into_ord)
    }
    fn slice_to_ord(slice: &[Self]) -> &[Self::Output] {
        // SAFETY: FloatOrd<f_ty> is repr(transparent) over f_ty
        unsafe { mem::transmute(slice) }
    }
}

impl IntoOrd for BasedExpr {
    type Output = BasedExpr;
    type OutputRef<'a> = &'a BasedExpr;
    fn into_ord(self) -> Self::Output { self }
    fn to_ord_ref(&self) -> Self::OutputRef<'_> { self }
    fn array_to_ord<const N: usize>(array: [Self; N]) -> [Self::Output; N] { array }
    fn slice_to_ord(slice: &[Self]) -> &[Self::Output] { slice }
}

pub const fn pow2i(exp: i32) -> f64 {
    if exp < (f64::MIN_EXP - 1) - (f64::MANTISSA_DIGITS - 1) as i32 {
        panic!("exponent is too small")
    } else if exp < f64::MIN_EXP - 1 {
        f64::from_bits(1 << (exp - ((f64::MIN_EXP - 1) - (f64::MANTISSA_DIGITS - 1) as i32)))
    } else if exp <= f64::MAX_EXP - 1 {
        f64::from_bits(((exp - (f64::MIN_EXP - 2)) as u64) << (f64::MANTISSA_DIGITS - 1))
    } else {
        panic!("exponent is too big")
    }
}

/// Rounds the float to one that's a multiple of 2^`lsb`, rounding ties to even.
pub fn with_max_lsb(this: f64, lsb: i32) -> f64 {
    let min_mag = pow2i(lsb + f64::MANTISSA_DIGITS as i32 - 1);
    let prec = pow2i(lsb);
    if this.abs() < min_mag {
        (this / prec).round_ties_even() * prec
    } else {
        this
    }
}

/// Gets the angle of a vector, even with exact coordinates. Returns None for the zero vector.
/// Use this instead of [`Vector2::angle`].
pub fn angle<T: IntoAngle>(vector: Vector2<T>) -> Option<T::Output> {
    let [[x, y]] = vector.data.0;
    x.into_angle(y)
}

/// Gets the angle of a vector, even with exact coordinates. Returns None for hte zero vector.
/// Use this instead of [`Vector2::angle`].
pub fn angle_ref<T>(vector: VectorView2Dyn<'_, T>) -> Option<<&T as IntoAngle>::Output> where
    for<'a> &'a T: IntoAngle
{
    let slice = vector.data.into_slice();
    (&slice[0]).into_angle(&slice[1])
}

/// Gets the orderable angle of a vector, even with exact coordinates.
/// Returns None for the zero vector.
pub fn ord_angle<T: IntoOrdAngle>(vector: Vector2<T>) -> Option<T::Output> {
    let [[x, y]] = vector.data.0;
    x.into_ord_angle(y)
}

/// Gets the orderable angle of a vector, even with exact coordinates.
/// Returns None for the zero vector.
pub fn ord_angle_ref<T>(vector: VectorView2Dyn<'_, T>) -> Option<<&T as IntoOrdAngle>::Output> where
    for<'a> &'a T: IntoOrdAngle
{
    let slice = vector.data.into_slice();
    (&slice[0]).into_ord_angle(&slice[1])
}

/// Sorts the coordinates counterclockwise.
pub fn sort_by_angle<T, S, F: FnMut(&T) -> Vector2<S>>(points: &mut [T], origin: VectorView2Dyn<S>, mut mapping: F) where 
    S: Scalar + ClosedSubAssign + IntoOrdAngle,
{
    points.sort_by_key(|p| ord_angle(mapping(p) - &origin));
}

/// Sorts the coordinates counterclockwise.
/// 
/// Unlike `sort_by_angle_field`, the return value of `mapping` *cannot* borrow from the argument,
/// but *can* from elsewhere. See https://github.com/rust-lang/rust/issues/34162
pub fn sort_by_angle_ref<'a, T, S, F: FnMut(&T) -> VectorView2Dyn<'a, S>>(points: &mut [T], origin: VectorView2Dyn<S>, mut mapping: F) where 
    S: Scalar + ClosedSubAssign + IntoOrdAngle,
{
    points.sort_by_key(|p| ord_angle(mapping(p) - &origin));
}

/// Sorts the coordinates counterclockwise.
/// 
/// Unlike `sort_by_angle_ref`, the return value of `mapping` *can* borrow from the argument,
/// but *not* from elsewhere. See https://github.com/rust-lang/rust/issues/34162
pub fn sort_by_angle_field<T, S, F: FnMut(&T) -> VectorView2Dyn<S>>(points: &mut [T], origin: VectorView2Dyn<S>, mut mapping: F) where 
    S: Scalar + ClosedSubAssign + IntoOrdAngle,
{
    points.sort_by_key(|p| ord_angle(mapping(p) - &origin));
}

/// Returns twice the signed area of the polygon in 2D defined by the input points.
//pub fn twice_signed_area_ref<'a, T, S, F: FnMut(&T) -> VectorView2Dyn<'a, S>>(points: &[T], mut mapping: F) -> T where
//    S: 'a
//{
//    todo!()
//    //points.column_iter().zip(points.column_iter().cycle().skip(1))
//    //    .map(|(v0, v1)| &v0.x * &v1.y - &v1.x * &v0.y)
//    //    .sum()
//}

/// Returns the orientation of the 2D polygon defined by the input points.
/// +1 for counterclockwise, -1 for clockwise
/// via computing sum of signed areas of triangles formed with origin
/// Assumes the polygon is simple (no self-touching) and nonempty
pub fn polygon_orientation<T, S, F: FnMut(&T) -> Vector2<S>>(points: &[T], mut mapping: F) -> i32 where
    S: IntoOrd + Scalar + ClosedSubAssign + AngleOps,
{
    let left_bottom = points.iter().enumerate()
        .min_by_key(|(_, point)| S::array_to_ord({let [arr] = mapping(*point).data.0; arr}))
        .expect("polygon must have points").0;
    let p0 = mapping(&points[(left_bottom + points.len() - 1) % points.len()]);
    let p1 = mapping(&points[left_bottom]);
    let p2 = mapping(&points[(left_bottom + 1) % points.len()]);
    match S::cmp_det_2d(p0.as_view::<U2, U1, U1, Dyn>(), p1, p0.as_view::<U2, U1, U1, Dyn>(), p2) {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

/// Returns the orientation of the 2D polygon defined by the input points.
/// +1 for counterclockwise, -1 for clockwise
/// via computing sum of signed areas of triangles formed with origin
/// Assumes the polygon is simple (no self-touching) and nonempty
pub fn polygon_orientation_ref<'a, T, S, F: FnMut(&T) -> VectorView2Dyn<'a, S>>(points: &[T], mut mapping: F) -> i32 where
    S: 'a + IntoOrd + Scalar + ClosedSubAssign + AngleOps,
{
    let left_bottom = points.iter().enumerate()
        .min_by_key(|(_, point)| S::slice_to_ord(mapping(*point).data.into_slice()))
        .expect("polygon must have points").0;
    let p0 = mapping(&points[(left_bottom + points.len() - 1) % points.len()]);
    let p1 = mapping(&points[left_bottom]);
    let p2 = mapping(&points[(left_bottom + 1) % points.len()]);
    match S::cmp_det_2d(p0.clone(), p1, p0.clone(), p2) {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

/// Returns the orientation of the 2D polygon defined by the input points.
/// +1 for counterclockwise, -1 for clockwise
/// via computing sum of signed areas of triangles formed with origin
/// Assumes the polygon is simple (no self-touching) and nonempty
pub fn polygon_orientation_field<T, S, F: FnMut(&T) -> VectorView2Dyn<S>>(points: &[T], mut mapping: F) -> i32 where
    S: IntoOrd + Scalar + ClosedSubAssign + AngleOps,
{
    let left_bottom = points.iter().enumerate()
        .min_by_key(|(_, point)| S::slice_to_ord(mapping(*point).data.into_slice()))
        .expect("polygon must have points").0;
    let p0 = mapping(&points[(left_bottom + points.len() - 1) % points.len()]);
    let p1 = mapping(&points[left_bottom]);
    let p2 = mapping(&points[(left_bottom + 1) % points.len()]);
    match S::cmp_det_2d(p0.clone(), p1, p0.clone(), p2) {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

/// Rotates a vector 90° counterclockwise and gets the result
pub fn perp_ccw<T: Scalar + Neg<Output = T>>(v: Vector2<T>) -> Vector2<T> {
    let [[x, y]] = v.data.0;
    vector![-y, x]
}

/// Rotates a vector 90° clockwise and gets the result
pub fn perp_cw<T: Scalar + Neg<Output = T>>(v: Vector2<T>) -> Vector2<T> {
    let [[x, y]] = v.data.0;
    vector![y, -x]
}

/// Tries to reflect vector `v` across the line through the origin and `a`
/// Returns `None` if `a` == `b`.
pub fn try_reflect<T: NumEx>(v: VectorView2Dyn<T>, a: VectorView2Dyn<T>) -> Option<Vector2<T>> {
    if a == Vector2::zeros() { None? }
    let n = perp_ccw(a.into());
    let two = T::one() + T::one();
    Some(&v - &n * (v.dot(&n) * two) / n.dot(&n))
}

/// Reflects vector `v` across the line through the origin and `a`
/// 
/// # Panics
/// Panics if `a` == `b`.
pub fn reflect<T: NumEx>(v: VectorView2Dyn<T>, a: VectorView2Dyn<T>) -> Vector2<T> {
    try_reflect(v, a).unwrap()
}

/// Tries to reflect point `p` across the line through `a` and `b`
/// Returns `None` if `a` == `b`.
pub fn try_reflect_line<T: NumEx>(p: VectorView2Dyn<T>, a: VectorView2Dyn<T>, b: VectorView2Dyn<T>) -> Option<Vector2<T>> {
    try_reflect((p - &a).as_view(), (b - &a).as_view()).map(|result| result + a)
}

/// Reflects point `p` across the line through `a` and `b`
/// 
/// # Panics
/// Panics if `a` == `b`.
pub fn reflect_line<T: NumEx>(p: VectorView2Dyn<T>, a: VectorView2Dyn<T>, b: VectorView2Dyn<T>) -> Vector2<T> {
    try_reflect_line(p, a, b).unwrap()
}

/// Tries to get the matrix for reflecting across the line `a` and `b`
/// Returns `None` if `a` == `b`.
pub fn try_reflect_line_matrix<T: NumEx + RealField>(a: VectorView2Dyn<T>, b: VectorView2Dyn<T>) -> Option<Affine2<T>> where
    for<'a> &'a T: RefNum<T>
{
    if a == b { None? }
    let v = b - &a;
    let dot = v.dot(&a);
    let len2 = v.dot(&v);
    let two = T::one() + T::one();
    Some(Affine2::from_matrix_unchecked(matrix![
        &v.x * &v.x / &len2 * &two - T::one(), &v.x * &v.y / &len2 * &two,            (&a.x - &v.x * &dot / &len2) * &two;
        &v.y * &v.x / &len2 * &two,            &v.y * &v.y / &len2 * &two - T::one(), (&a.y - &v.y * &dot / &len2) * &two;
        T::zero(),                             T::zero(),                             T::one()
    ]))
}

/// Gets the matrix for reflecting across the line `a` and `b`.
/// 
/// # Panics
/// Panics if `a` == `b`.
pub fn reflect_line_matrix<T: NumEx + RealField>(a: VectorView2Dyn<T>, b: VectorView2Dyn<T>) -> Affine2<T> where
    for<'a> &'a T: RefNum<T>
{
    try_reflect_line_matrix(a, b).unwrap()
}

/// Applies an affine transform to a point `p`.
/// 
/// This is different from `nalgebra::geometry::Transform::transform_point` because
/// it doesn't require creating a point
pub fn transform<T, const D: usize>(transform: &Transform<T, TAffine, D>, p: VectorView<'_, T, Const<D>, U1, Dyn>) -> SVector<T, D> where
    T: NumEx + RealField,
    Const<D>: DimNameAdd<U1>,
    DefaultAllocator: Allocator<DimNameSum<Const<D>, U1>, DimNameSum<Const<D>, U1>>
{
    transform.matrix().fixed_view::<D, D>(0, 0) * p + transform.matrix().fixed_view::<D, 1>(0, D)
}

/// Gets the parameters t1 and t2 for the intersection
/// of two parameterized lines specified by t=0 and t=1 points
/// Returns `None` if the lines are parallel
/// 
/// Warning: this does not use epsilon comparison. Not relevant for exact coordinates,
/// but be careful when using floating-point numbers
pub fn parametric_line_intersect<T: NumEx + RealField>
    (line_a: [VectorView2Dyn<T>; 2], line_b: [VectorView2Dyn<T>; 2]) -> Option<Vector2<T>>
{
    let vec_a = &line_a[1] - &line_a[0];
    let vec_b = &line_b[1] - &line_b[0];
    let start_diff = &line_b[0] - &line_a[0];
    let mtx = Matrix2::from_columns(&[vec_a, -vec_b]);
    mtx.try_inverse().map(|inv| inv * start_diff)
}

/// The result of intersecting two lines
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LineIntersection<T> {
    /// The lines don't intersect (a.k.a. they're parallel)
    None,
    /// An intersection between non-parallel segments
    Intersection(Vector2<T>),
    /// The lines are collinear
    Collinear,
}

/// The result of intersecting two segments
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SegmentIntersection<T> {
    /// The segments don't intersect
    None,
    /// An intersection between non-parallel segments
    Intersection(Vector2<T>),
    /// A shared section between two collinear segments.
    /// The order of the points is the order within the first segment.
    Collinear([Vector2<T>; 2]),
}

/// Gets the intersection between two lines `line_a` and `line_b`
/// 
/// Warning: this does not use epsilon comparison. Not relevant for exact coordinates,
/// but be careful when using floating-point numbers
pub fn line_intersect<T: NumEx + RealField>(line_a: [VectorView2Dyn<T>; 2], line_b: [VectorView2Dyn<T>; 2]) -> LineIntersection<T> {
    if let Some(mut t) = parametric_line_intersect(line_a.clone(), line_b.clone()) {
        LineIntersection::Intersection((&line_a[1] - &line_a[0]) * mem::take(&mut t.x) + &line_a[0])
    } else {
        // The lines are parallel.
        let vec_a = &line_a[1] - &line_a[0];

        if vec_a == Vector2::zeros() {
            let vec_b = &line_b[1] - &line_b[0];
            if !(&line_b[0] - &line_a[0]).perp(&vec_b).is_zero() { return LineIntersection::None } // not collinear
            return LineIntersection::Collinear
        }

        if !(&line_b[0] - &line_a[0]).perp(&vec_a).is_zero() { return LineIntersection::None } // not collinear in the first place
        LineIntersection::Collinear
    }
}

/// Gets the intersection between two segments `seg_a` and `seg_b`. Note that boundary points are included.
/// 
/// Warning: this does not use epsilon comparison. Not relevant for exact coordinates,
/// but be careful when using floating-point numbers
pub fn segment_intersect<T: NumEx + RealField>(seg_a: [VectorView2Dyn<T>; 2], seg_b: [VectorView2Dyn<T>; 2]) -> SegmentIntersection<T> {
    if let Some(mut t) = parametric_line_intersect(seg_a.clone(), seg_b.clone()) {
        if t.x >= T::zero() && t.x <= T::one() && t.y >= T::zero() && t.y <= T::one() {
            SegmentIntersection::Intersection((&seg_a[1] - &seg_a[0]) * mem::take(&mut t.x) + &seg_a[0])
        } else {
            SegmentIntersection::None
        }
    } else {
        // The lines are parallel.
        let vec_a = &seg_a[1] - &seg_a[0];

        if vec_a == Vector2::zeros() {
            let vec_b = &seg_b[1] - &seg_b[0];
            if !(&seg_b[0] - &seg_a[0]).perp(&vec_b).is_zero() { return SegmentIntersection::None } // not collinear

            if vec_b == Vector2::zeros() {
                return if seg_a[0] == seg_b[0] {
                    SegmentIntersection::Collinear(seg_a.map(|p| p.into_owned()))
                } else { SegmentIntersection::None }
            }
            let len2 = vec_b.dot(&vec_b);
            let t = (&seg_a[0] - &seg_b[0]).dot(&vec_b) / len2;
            return if t >= T::zero() && t <= T::one() {
                SegmentIntersection::Collinear(seg_a.map(|p| p.into_owned()))
            } else { SegmentIntersection::None }
        }

        if !(&seg_b[0] - &seg_a[0]).perp(&vec_a).is_zero() { return SegmentIntersection::None } // not collinear in the first place
        let len2 = vec_a.dot(&vec_a);
        let mut t0 = (&seg_b[0] - &seg_a[0]).dot(&vec_a) / &len2;
        let mut t1 = (&seg_b[1] - &seg_a[0]).dot(&vec_a) / &len2;
        if t0 > t1 { mem::swap(&mut t0, &mut t1) }
        t0 = t0.max(T::zero());
        t1 = t1.min(T::one());
        if t0 <= t1 {
            SegmentIntersection::Collinear([&vec_a * t0 + &seg_a[0], &vec_a * t1 + &seg_a[0]])
        } else {
            SegmentIntersection::None
        }
    }
}

#[cfg(test)]
mod test {
    use std::fmt::Debug;

    use exact_number::{based_expr};
    use nalgebra::{Affine2, ClosedSubAssign, Matrix2xX, Scalar, U2, Vector2, matrix, vector};

    use crate::geom::{AngleOps, IntoOrd, IntoOrdAngle, IntoOrdAngleOp, RefIntoOrdAngleOp, SegmentIntersection, VectorView2Dyn, polygon_orientation, polygon_orientation_field, polygon_orientation_ref, pow2i, reflect_line, reflect_line_matrix, segment_intersect, sort_by_angle, sort_by_angle_field, sort_by_angle_ref, try_reflect_line, try_reflect_line_matrix, with_max_lsb};

    /// Gets the permutation required to take `a` to `b`.
    /// `permutation(a, b).0[i] == i`
    /// `a[permutation(a, b).1[i]] == b[i]`
    fn permutation<T: PartialEq>(a: &[T], b: &[T]) -> (Vec<usize>, Vec<usize>) {
        let indexes = Vec::from_iter(0..a.len());
        let perm = b.iter()
            .map(|b| a.iter().position(|x| x == b).unwrap())
            .collect::<Vec<_>>();
        (indexes, perm)
    }

    #[test]
    fn test_pow2i() {
        // There aren't that many; test every valid exponent. Every single one.
        for exp in (-1074)..=1023 {
            assert_eq!(pow2i(exp), if exp < 0 { 0.5f64.powi(-exp) } else { 2f64.powi(exp) }, "failed for exp={}", exp);
        }
    }

    #[test]
    fn test_with_max_lsb() {
        assert_eq!(with_max_lsb(0.0, -100), 0.0);
        assert_eq!(with_max_lsb(1.0, -100), 1.0);
        assert_eq!(with_max_lsb(pow2i(-100), -100), pow2i(-100));
        assert_eq!(with_max_lsb(pow2i(-101), -100), 0.0); // break ties toward even
        assert_eq!(with_max_lsb(pow2i(-101) * 1.5, -100), pow2i(-100));
        assert_eq!(with_max_lsb(pow2i(-48) + pow2i(-100), -100), pow2i(-48) + pow2i(-100));
    }
    
    #[test]
    fn test_cmp_angle_f64() {
        let angles = vec![
            ( 0, (-1.0f64).into_ord_angle( 0.0f64).unwrap()),
            ( 0, (-2.0f64).into_ord_angle( 0.0f64).unwrap()),
            ( 1, (-2.0f64).into_ord_angle(-1.0f64).unwrap()),
            ( 1, (-4.0f64).into_ord_angle(-2.0f64).unwrap()),
            ( 2, (-1.0f64).into_ord_angle(-1.0f64).unwrap()),
            ( 2, (-2.0f64).into_ord_angle(-2.0f64).unwrap()),
            ( 3, (-1.0f64).into_ord_angle(-2.0f64).unwrap()),
            ( 3, (-2.0f64).into_ord_angle(-4.0f64).unwrap()),
            ( 4, ( 0.0f64).into_ord_angle(-1.0f64).unwrap()),
            ( 4, ( 0.0f64).into_ord_angle(-2.0f64).unwrap()),
            ( 5, ( 1.0f64).into_ord_angle(-2.0f64).unwrap()),
            ( 5, ( 2.0f64).into_ord_angle(-4.0f64).unwrap()),
            ( 6, ( 1.0f64).into_ord_angle(-1.0f64).unwrap()),
            ( 6, ( 2.0f64).into_ord_angle(-2.0f64).unwrap()),
            ( 7, ( 2.0f64).into_ord_angle(-1.0f64).unwrap()),
            ( 7, ( 4.0f64).into_ord_angle(-2.0f64).unwrap()),
            ( 8, ( 1.0f64).into_ord_angle( 0.0f64).unwrap()),
            ( 8, ( 2.0f64).into_ord_angle( 0.0f64).unwrap()),
            ( 9, ( 2.0f64).into_ord_angle( 1.0f64).unwrap()),
            ( 9, ( 4.0f64).into_ord_angle( 2.0f64).unwrap()),
            (10, ( 1.0f64).into_ord_angle( 1.0f64).unwrap()),
            (10, ( 2.0f64).into_ord_angle( 2.0f64).unwrap()),
            (11, ( 1.0f64).into_ord_angle( 2.0f64).unwrap()),
            (11, ( 2.0f64).into_ord_angle( 4.0f64).unwrap()),
            (12, ( 0.0f64).into_ord_angle( 1.0f64).unwrap()),
            (12, ( 0.0f64).into_ord_angle( 2.0f64).unwrap()),
            (13, (-1.0f64).into_ord_angle( 2.0f64).unwrap()),
            (13, (-2.0f64).into_ord_angle( 4.0f64).unwrap()),
            (14, (-1.0f64).into_ord_angle( 1.0f64).unwrap()),
            (14, (-2.0f64).into_ord_angle( 2.0f64).unwrap()),
            (15, (-2.0f64).into_ord_angle( 1.0f64).unwrap()),
            (15, (-4.0f64).into_ord_angle( 2.0f64).unwrap()),
        ];

        for &(ord_a, angle_a) in &angles {
            for &(ord_b, angle_b) in &angles {
                let result = angle_a.cmp(&angle_b);
                let expected = ord_a.cmp(&ord_b);
                assert_eq!(result, expected, "[{}, {}] vs [{}, {}] comparison is wrong: got {result:?}, expected {expected:?}",
                    angle_a.0[1].x, angle_a.0[1].y, angle_b.0[1].x, angle_b.0[1].y);
            }
        }
    }

    fn sort_by_angle_test<T, S>(points: Vec<T>, origin: Vector2<S>, mut mapping: impl FnMut(&T) -> VectorView2Dyn<S>, expected: Vec<T>) where
        T: Clone + Debug + PartialEq,
        S: Scalar + ClosedSubAssign + IntoOrdAngle,
    {
        let (mut indexes, expected_indexes) = permutation(&points, &expected);
        let mut points_a = points.clone();
        let mut points_b = points.clone();
        sort_by_angle_ref(&mut indexes, origin.as_view(), |i| mapping(&points[*i]));
        assert_eq!(indexes, expected_indexes);
        sort_by_angle_field(&mut points_a, origin.as_view(), &mut mapping);
        assert_eq!(points_a, expected);
        sort_by_angle(&mut points_b, origin.as_view(), |t| mapping(t).clone_owned());
        assert_eq!(points_b, expected);
    }

    #[test]
    fn test_sort_by_angle() {
        sort_by_angle_test(vec![], Vector2::zeros(), |x: &Vector2<f64>| x.as_view(), vec![]);

        sort_by_angle_test(vec![vector!(0.0, 0.0)], Vector2::zeros(), |x| x.as_view(), vec![vector!(0.0, 0.0)]);
        sort_by_angle_test(vec![vector!(1.0, 0.0)], Vector2::zeros(), |x| x.as_view(), vec![vector!(1.0, 0.0)]);
        sort_by_angle_test(vec![vector!(1.0, 0.0), vector!(-1.0, 0.0), vector!(0.0, 1.0), vector!(0.0, -1.0)],
            Vector2::zeros(), |x| x.as_view(),
            vec![vector!(-1.0, 0.0), vector!(0.0, -1.0), vector!(1.0, 0.0), vector!(0.0, 1.0)]);
        sort_by_angle_test(vec![vector!(1.0, 0.1), vector!(-1.0, 0.0), vector!(0.0, 1.0), vector!(0.0, -1.0)],
            vector!(2.0, 0.0), |x| x.as_view(),
            vec![vector!(-1.0, 0.0), vector!(0.0, -1.0), vector!(0.0, 1.0), vector!(1.0, 0.1)]);
    }

    fn polygon_orientation_test<T, S>(points: Vec<T>, mut mapping: impl FnMut(&T) -> VectorView2Dyn<S>, expected: i32) where
        T: Clone + Debug + PartialEq,
        S: IntoOrd + Scalar + ClosedSubAssign + AngleOps,
    {
        let indexes = (0..points.len()).collect::<Vec<_>>();
        let result = polygon_orientation_ref(&indexes, |i| mapping(&points[*i]));
        assert_eq!(result, expected);
        let result = polygon_orientation_field(&points, &mut mapping);
        assert_eq!(result, expected);
        let result = polygon_orientation(&points, |t| mapping(t).clone_owned());
        assert_eq!(result, expected);
    }

    // TODO: Reactivate
    //#[test]
    //fn test_twice_signed_area() {
    //    assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
    //        based_expr!(0), based_expr!(0),
    //        based_expr!(2), based_expr!(0),
    //    ]).as_view()), based_expr!(0));
    //    assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
    //        based_expr!(0), based_expr!(0),
    //        based_expr!(0), based_expr!(0),
    //        based_expr!(0), based_expr!(0),
    //    ]).as_view()), based_expr!(0));
    //    assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
    //        based_expr!(5), based_expr!(9),
    //        based_expr!(2), based_expr!(7),
    //        based_expr!(1), based_expr!(7),
    //    ]).as_view()), based_expr!(-2));
    //    assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
    //        based_expr!(0), based_expr!(3),
    //        based_expr!(5), based_expr!(8),
    //        based_expr!(3), based_expr!(1),
    //    ]).as_view()), based_expr!(-25));
    //    assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
    //        4.5, 1.0, // Mixing it up
    //        4.5, 8.0,
    //        2.0, 5.0,
    //    ]).as_view()), 17.5);
    //    assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
    //        1.0, 7.0,
    //        0.0, 4.0,
    //        5.0, 9.0,
    //        2.0, 7.0,
    //    ]).as_view()), 8.0);
    //    assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
    //        0.0, 3.0,
    //        -1.0, 5.0,
    //        3.0, 1.0,
    //        5.0, 8.0
    //    ]).as_view()), 21.0);
    //    assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
    //        based_expr!(9/2), based_expr!(1),
    //        based_expr!(9/2), based_expr!(8),
    //        based_expr!(2), based_expr!(5),
    //        based_expr!(7/2), based_expr!(3),
    //    ]).as_view()), based_expr!(33/2));
    //}

    #[test]
    fn test_polygon_orientation() {
        polygon_orientation_test(vec![
            vector![based_expr!(0), based_expr!(0)],
            vector![based_expr!(2), based_expr!(0)],
        ], |v| v.as_view(), 0);
        polygon_orientation_test(vec![
            vector![based_expr!(0), based_expr!(0)],
            vector![based_expr!(0), based_expr!(0)],
            vector![based_expr!(0), based_expr!(0)],
        ], |v| v.as_view(), 0);
        polygon_orientation_test(vec![
            vector![based_expr!(5), based_expr!(9)],
            vector![based_expr!(2), based_expr!(7)],
            vector![based_expr!(1), based_expr!(7)],
        ], |v| v.as_view(), -1);
        polygon_orientation_test(vec![
            vector![based_expr!(0), based_expr!(3)],
            vector![based_expr!(5), based_expr!(8)],
            vector![based_expr!(3), based_expr!(1)],
        ], |v| v.as_view(), -1);
        polygon_orientation_test(vec![
            vector![4.5, 1.0], // Mixing it up
            vector![4.5, 8.0],
            vector![2.0, 5.0],
        ], |v| v.as_view(), 1);
        polygon_orientation_test(vec![
            vector![1.0, 7.0],
            vector![0.0, 4.0],
            vector![5.0, 9.0],
            vector![2.0, 7.0],
        ], |v| v.as_view(), 1);
        polygon_orientation_test(vec![
            vector![based_expr!(9/2), based_expr!(1)],
            vector![based_expr!(9/2), based_expr!(8)],
            vector![based_expr!(2), based_expr!(5)],
            vector![based_expr!(7/2), based_expr!(3)],
        ], |v| v.as_view(), 1);
    }

    #[test]
    fn test_reflect_line() {
        assert_eq!(try_reflect_line(vector![1.0, 2.0].as_view(), vector![-0.5, 1.5].as_view(), vector![-0.5, 1.5].as_view()), None);

        // Reflections through the origin
        assert_eq!(reflect_line(vector![1.0, 2.0].as_view(), vector![0.0, 0.0].as_view(), vector![1.0, 0.0].as_view()), vector![1.0, -2.0]);
        assert_eq!(reflect_line(vector![1.0, 2.0].as_view(), vector![0.0, 0.0].as_view(), vector![-1.0, 0.0].as_view()), vector![1.0, -2.0]);
        assert_eq!(reflect_line(vector![1.0, 2.0].as_view(), vector![0.0, 0.0].as_view(), vector![0.0, 1.0].as_view()), vector![-1.0, 2.0]);
        assert_eq!(reflect_line(vector![1.0, 2.0].as_view(), vector![0.0, 0.0].as_view(), vector![0.0, -1.0].as_view()), vector![-1.0, 2.0]);

        // Reflections away from the origin
        assert_eq!(reflect_line(vector![1.0, 2.0].as_view(), vector![-3.0, 4.0].as_view(), vector![-2.0, 4.0].as_view()), vector![1.0, 6.0]);
        assert_eq!(reflect_line(vector![1.0, 2.0].as_view(), vector![-3.0, 4.0].as_view(), vector![-4.0, 4.0].as_view()), vector![1.0, 6.0]);
        assert_eq!(reflect_line(vector![1.0, 2.0].as_view(), vector![-3.0, 4.0].as_view(), vector![-3.0, 5.0].as_view()), vector![-7.0, 2.0]);
        assert_eq!(reflect_line(vector![1.0, 2.0].as_view(), vector![-3.0, 4.0].as_view(), vector![-3.0, 3.0].as_view()), vector![-7.0, 2.0]);

        // Reflection that doesn't do anything
        assert_eq!(reflect_line(vector![0.0, 1.0].as_view(), vector![0.0, 0.0].as_view(), vector![0.0, 2.0].as_view()), vector![0.0, 1.0]);
        
        // Reflections across slanted lines
        assert_eq!(reflect_line(vector![1.0, 0.0].as_view(), vector![0.0, 0.0].as_view(), vector![1.0, 1.0].as_view()), vector![0.0, 1.0]);
        assert_eq!(reflect_line(
            vector![based_expr!(1), based_expr!(0)].as_view(),
            vector![based_expr!(0), based_expr!(0)].as_view(),
            vector![based_expr!(2), based_expr!(1)].as_view()
        ), vector![based_expr!(3/5), based_expr!(4/5)]);
    }

    #[test]
    fn test_reflect_line_matrix() {
        assert_eq!(try_reflect_line_matrix(vector![-0.5, 1.5].as_view(), vector![-0.5, 1.5].as_view()), None);

        // Reflections through the origin
        assert_eq!(reflect_line_matrix(vector![0.0, 0.0].as_view(), vector![1.0, 0.0].as_view()), Affine2::from_matrix_unchecked(matrix![
            1.0, 0.0, 0.0;
            0.0, -1.0, 0.0;
            0.0, 0.0, 1.0;
        ]));
        assert_eq!(reflect_line_matrix(vector![0.0, 0.0].as_view(), vector![-1.0, 0.0].as_view()), Affine2::from_matrix_unchecked(matrix![
            1.0, 0.0, 0.0;
            0.0, -1.0, 0.0;
            0.0, 0.0, 1.0;
        ]));
        assert_eq!(reflect_line_matrix(vector![0.0, 0.0].as_view(), vector![0.0, 1.0].as_view()), Affine2::from_matrix_unchecked(matrix![
            -1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 1.0;
        ]));
        assert_eq!(reflect_line_matrix(vector![0.0, 0.0].as_view(), vector![0.0, -1.0].as_view()), Affine2::from_matrix_unchecked(matrix![
            -1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 1.0;
        ]));

        // Reflections away from the origin
        assert_eq!(reflect_line_matrix(vector![-3.0, 4.0].as_view(), vector![-2.0, 4.0].as_view()), Affine2::from_matrix_unchecked(matrix![
            1.0, 0.0, 0.0;
            0.0, -1.0, 8.0;
            0.0, 0.0, 1.0;
        ]));
        assert_eq!(reflect_line_matrix(vector![-3.0, 4.0].as_view(), vector![-3.0, 5.0].as_view()), Affine2::from_matrix_unchecked(matrix![
            -1.0, 0.0, -6.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 1.0;
        ]));

        // Reflections across slanted lines
        assert_eq!(reflect_line_matrix(vector![0.0, 0.0].as_view(), vector![1.0, 1.0].as_view()), Affine2::from_matrix_unchecked(matrix![
            0.0, 1.0, 0.0;
            1.0, 0.0, 0.0;
            0.0, 0.0, 1.0;
        ]));
        assert_eq!(reflect_line_matrix(
            vector![based_expr!(0), based_expr!(0)].as_view(),
            vector![based_expr!(2), based_expr!(1)].as_view()
        ), Affine2::from_matrix_unchecked(matrix![
            based_expr!(3/5), based_expr!(4/5),  based_expr!(0);
            based_expr!(4/5), based_expr!(-3/5), based_expr!(0);
            based_expr!(0),   based_expr!(0),    based_expr!(1);
        ]));
    }

    macro_rules! exact_vec_2 {
        (($($a:tt)*), ($($b:tt)*)$(,)?) => {
            nalgebra::vector![exact_number::based_expr!($($a)*), exact_number::based_expr!($($b)*)]
        };
    }

    macro_rules! exact_view_2 {
        (($($a:tt)*), ($($b:tt)*)$(,)?) => {
            nalgebra::vector![exact_number::based_expr!($($a)*), exact_number::based_expr!($($b)*)].as_view()
        };
    }

    #[test]
    fn test_segment_intersect() {
        // t extremes
        assert_eq!(segment_intersect(
            [exact_view_2![(0), (0)], exact_view_2![(1), (0)]], [exact_view_2![(0), (0)], exact_view_2![(0), (1)]]
        ), SegmentIntersection::Intersection(exact_vec_2![(0), (0)]));
        assert_eq!(segment_intersect(
            [exact_view_2![(0), (0)], exact_view_2![(1), (0)]], [exact_view_2![(1), (0)], exact_view_2![(1), (1)]]
        ), SegmentIntersection::Intersection(exact_vec_2![(1), (0)]));
        assert_eq!(segment_intersect(
            [exact_view_2![(0), (1)], exact_view_2![(1), (1)]], [exact_view_2![(0), (0)], exact_view_2![(0), (1)]]
        ), SegmentIntersection::Intersection(exact_vec_2![(0), (1)]));
        assert_eq!(segment_intersect(
            [exact_view_2![(0), (1)], exact_view_2![(1), (1)]], [exact_view_2![(1), (0)], exact_view_2![(1), (1)]]
        ), SegmentIntersection::Intersection(exact_vec_2![(1), (1)]));

        // Some "normal" intersections
        assert_eq!(segment_intersect(
            [exact_view_2![(1), (1)], exact_view_2![(2), (2)]], [exact_view_2![(1), (2)], exact_view_2![(2), (1)]]
        ), SegmentIntersection::Intersection(exact_vec_2![(3/2), (3/2)]));
        assert_eq!(segment_intersect(
            [exact_view_2![(0), (0)], exact_view_2![(1), (1)]], [exact_view_2![(0), (1/2)], exact_view_2![(1), (0)]]
        ), SegmentIntersection::Intersection(exact_vec_2![(1/3), (1/3)]));

        // out of bounds, so no intersection
        assert_eq!(segment_intersect(
            [exact_view_2![(1), (0)], exact_view_2![(2), (0)]], [exact_view_2![(0), (-1)], exact_view_2![(0), (1)]]
        ), SegmentIntersection::None);
        assert_eq!(segment_intersect(
            [exact_view_2![(-2), (0)], exact_view_2![(-1), (0)]], [exact_view_2![(0), (-1)], exact_view_2![(0), (1)]]
        ), SegmentIntersection::None);
        assert_eq!(segment_intersect(
            [exact_view_2![(-1), (0)], exact_view_2![(1), (0)]], [exact_view_2![(0), (1)], exact_view_2![(0), (2)]]
        ), SegmentIntersection::None);
        assert_eq!(segment_intersect(
            [exact_view_2![(-1), (0)], exact_view_2![(1), (0)]], [exact_view_2![(0), (-2)], exact_view_2![(0), (-1)]]
        ), SegmentIntersection::None);
        
        // parallel lines
        assert_eq!(segment_intersect(
            [exact_view_2![(0), (1)], exact_view_2![(1), (2)]], [exact_view_2![(2), (4)], exact_view_2![(4), (6)]]
        ), SegmentIntersection::None);
        assert_eq!(segment_intersect(
            [exact_view_2![(0), (1)], exact_view_2![(1), (2)]], [exact_view_2![(0), (0)], exact_view_2![(1), (1)]]
        ), SegmentIntersection::None);

        // collinear lines that don't intersect
        assert_eq!(segment_intersect(
            [exact_view_2![(0), (1)], exact_view_2![(1), (2)]], [exact_view_2![(2), (3)], exact_view_2![(3), (4)]]
        ), SegmentIntersection::None);
        assert_eq!(segment_intersect(
            [exact_view_2![(0), (1)], exact_view_2![(1), (2)]], [exact_view_2![(3), (4)], exact_view_2![(2), (3)]]
        ), SegmentIntersection::None);
        assert_eq!(segment_intersect(
            [exact_view_2![(1), (2)], exact_view_2![(0), (1)]], [exact_view_2![(2), (3)], exact_view_2![(3), (4)]]
        ), SegmentIntersection::None);
        assert_eq!(segment_intersect(
            [exact_view_2![(1), (2)], exact_view_2![(0), (1)]], [exact_view_2![(3), (4)], exact_view_2![(2), (3)]]
        ), SegmentIntersection::None);

        // collinear lines that intersect
        // one segment is properly inside the other
        assert_eq!(segment_intersect(
            [exact_view_2![(0), (1)], exact_view_2![(3), (4)]], [exact_view_2![(1), (2)], exact_view_2![(2), (3)]]
        ), SegmentIntersection::Collinear([exact_vec_2![(1), (2)], exact_vec_2![(2), (3)]]));
        assert_eq!(segment_intersect(
            [exact_view_2![(0), (1)], exact_view_2![(3), (4)]], [exact_view_2![(2), (3)], exact_view_2![(1), (2)]]
        ), SegmentIntersection::Collinear([exact_vec_2![(1), (2)], exact_vec_2![(2), (3)]]));
        assert_eq!(segment_intersect(
            [exact_view_2![(3), (4)], exact_view_2![(0), (1)]], [exact_view_2![(1), (2)], exact_view_2![(2), (3)]]
        ), SegmentIntersection::Collinear([exact_vec_2![(2), (3)], exact_vec_2![(1), (2)]]));
        assert_eq!(segment_intersect(
            [exact_view_2![(3), (4)], exact_view_2![(0), (1)]], [exact_view_2![(2), (3)], exact_view_2![(1), (2)]]
        ), SegmentIntersection::Collinear([exact_vec_2![(2), (3)], exact_vec_2![(1), (2)]]));
        // no segment contains the other
        assert_eq!(segment_intersect(
            [exact_view_2![(0), (1)], exact_view_2![(2), (3)]], [exact_view_2![(1), (2)], exact_view_2![(3), (4)]]
        ), SegmentIntersection::Collinear([exact_vec_2![(1), (2)], exact_vec_2![(2), (3)]]));
        // the segments intersect at a point
        assert_eq!(segment_intersect(
            [exact_view_2![(0), (1)], exact_view_2![(2), (3)]], [exact_view_2![(2), (3)], exact_view_2![(3), (4)]]
        ), SegmentIntersection::Collinear([exact_vec_2![(2), (3)], exact_vec_2![(2), (3)]]));

        // point intersections
        assert_eq!(segment_intersect(
            [exact_view_2![(0), (0)], exact_view_2![(0), (2)]], [exact_view_2![(1), (1)], exact_view_2![(1), (1)]]
        ), SegmentIntersection::None);
        assert_eq!(segment_intersect(
            [exact_view_2![(1), (1)], exact_view_2![(1), (1)]], [exact_view_2![(2), (0)], exact_view_2![(2), (2)]]
        ), SegmentIntersection::None);
        assert_eq!(segment_intersect(
            [exact_view_2![(1), (1)], exact_view_2![(1), (1)]], [exact_view_2![(2), (1)], exact_view_2![(2), (1)]]
        ), SegmentIntersection::None);

        assert_eq!(segment_intersect(
            [exact_view_2![(0), (0)], exact_view_2![(0), (2)]], [exact_view_2![(0), (1)], exact_view_2![(0), (1)]]
        ), SegmentIntersection::Collinear([exact_vec_2![(0), (1)], exact_vec_2![(0), (1)]]));
        assert_eq!(segment_intersect(
            [exact_view_2![(1), (1)], exact_view_2![(1), (1)]], [exact_view_2![(1), (0)], exact_view_2![(1), (2)]]
        ), SegmentIntersection::Collinear([exact_vec_2![(1), (1)], exact_vec_2![(1), (1)]]));
        assert_eq!(segment_intersect(
            [exact_view_2![(1), (1)], exact_view_2![(1), (1)]], [exact_view_2![(1), (1)], exact_view_2![(1), (1)]]
        ), SegmentIntersection::Collinear([exact_vec_2![(1), (1)], exact_vec_2![(1), (1)]]));
    }
}