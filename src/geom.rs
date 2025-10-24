use std::{cmp::Ordering, iter::Sum, ops::{Mul, Sub, Neg}};

use exact_number::{malachite::base::num::arithmetic::traits::{Sign, CheckedDiv}, rat::Rat, BasedExpr, Angle};
use float_ord::FloatOrd;
use nalgebra::{allocator::Allocator, matrix, vector, Affine2, ArrayStorage, ClosedSubAssign, Const, DefaultAllocator, DimNameAdd, DimNameSum, Dyn, MatrixView2xX, SVector, RealField, Scalar, TAffine, ToTypenum, Transform, Vector, Vector2, VectorView, VectorView2, U1};
use num_traits::{Num, NumAssign, NumAssignRef, NumRef, RefNum, Signed, Zero};

pub type VectorView2Dyn<'s, T> = VectorView2<'s, T, U1, Dyn>;
pub type MatrixView2Dyn<'s, T> = MatrixView2xX<'s, T, U1, Dyn>;

pub trait NumEx: Default + PartialOrd + Num + NumRef + NumAssignRef + Scalar + NumAssign + Signed + Neg<Output = Self> {}
impl<T> NumEx for T where T: Default + PartialOrd + Num + NumRef + NumAssignRef + Scalar + NumAssign + Signed + Neg<Output = Self> {}

///// An angle, either exact or approximate
//#[derive(Clone, Debug, PartialEq, PartialOrd)]
//pub enum AngleRef {
//    Exact(Angle),
//    Approx(f64)
//}

/// A trait that allows getting the angle of a 2D vector of this type.
/// This is often used for sorting by angle.
pub trait AngleRep: Sized {
    type Output: Ord;
    /// The angle representative, representing the angle of this vector.
    /// 
    /// # Requirements
    /// * The `[0, 0]` vector must return 0.
    /// * The `[1, 0]` direction must return 0.
    /// * The `[-1, 0]` direction must return either the minimum or maximum value this function returns,
    ///     and it must be deterministic.
    /// * If vectors `[-1, 0]`, `a`, and `b` are ordered counterclockwise
    ///     and neither `a` nor `b` is `[0, 0]` or pointing in the -x direction, then
    ///     `angle_rep(a)` ≤ `angle_rep(b)`, with strict inequality
    ///     where vectors aren't pointing in the same direction.
    fn angle_rep(self) -> Self::Output;
}

impl AngleRep for Vector2<f32> {
    type Output = FloatOrd<f32>;

    /// The `[-1, 0]` direction returns the maximum value because that's how `f32::atan2` is defined
    fn angle_rep(self) -> Self::Output {
        FloatOrd(self.y.atan2(self.x))
    }
}

impl AngleRep for Vector2<f64> {
    type Output = FloatOrd<f64>;
    
    /// The `[-1, 0]` direction returns the maximum value because that's how `f64::atan2` is defined
    fn angle_rep(self) -> Self::Output {
        FloatOrd(self.y.atan2(self.x))
    }
}

impl AngleRep for Vector2<BasedExpr> {
    type Output = Angle;

    /// The `[-1, 0]` direction returns the minimum value for representation convience.
    fn angle_rep(self) -> Self::Output {
        let [[x, y]] = self.data.0;
        Angle::atan2(y, x)
        //let x = &self.x;
        //let y = &self.y;
        //if x.is_zero() && y.is_zero() { return x.clone() }

        //let (y, flip_y) = if y.is_negative() { (&-y, true) } else { (y, false) };
        //let (x, flip_x) = if x.is_negative() { (&-x, true) } else { (x, false) };
        //let (x, y, transpose) = if y > x { (y, x, true) } else { (x, y, false) };
        //let mut result = y / x;
        //if transpose { result = BasedExpr::Baseless(2.into()) - result }
        //if flip_x { result = BasedExpr::Baseless(4.into()) - result }
        //if flip_y { result = -result }
        //result
    }
}

/// Sorts the coordinates by angle increasing, according to `AngleRep::angle_rep`.
pub fn sort_by_angle<T, S, F: FnMut(&T) -> Vector2<S>>(points: &mut [T], origin: &T, mut mapping: F) where 
    S: Scalar + ClosedSubAssign,
    Vector2<S>: AngleRep
{
    let origin = mapping(origin);
    points.sort_by_key(|p| (mapping(p) - &origin).angle_rep());
}

/// Sorts the coordinates by angle increasing, according to `AngleRep::angle_rep`.
/// 
/// Unlike `sort_by_angle_field`, the return value of `mapping` *cannot* borrow from the argument,
/// but *can* from elsewhere. See https://github.com/rust-lang/rust/issues/34162
pub fn sort_by_angle_ref<'a, T, S, F: FnMut(&T) -> VectorView2Dyn<'a, S>>(points: &mut [T], origin: &T, mut mapping: F) where 
    S: Scalar + ClosedSubAssign,
    Vector2<S>: AngleRep
{
    let origin = mapping(origin);
    points.sort_by_key(|p| (mapping(p) - &origin).angle_rep());
}

/// Sorts the coordinates by angle increasing, according to `AngleRep::angle_rep`.
/// 
/// Unlike `sort_by_angle_ref`, the return value of `mapping` *can* borrow from the argument,
/// but *not* from elsewhere. See https://github.com/rust-lang/rust/issues/34162
pub fn sort_by_angle_field<T, S, F: FnMut(&T) -> VectorView2Dyn<S>>(points: &mut [T], origin: &T, mut mapping: F) where 
    S: Scalar + ClosedSubAssign,
    Vector2<S>: AngleRep
{
    let origin = mapping(origin);
    points.sort_by_key(|p| (mapping(p) - &origin).angle_rep());
}

/// Returns twice the signed area of the polygon in 2D defined by the input points.
pub fn twice_signed_area<T>(points: MatrixView2Dyn<T>) -> T where
    T: Scalar + Sub<Output = T> + Sum,
    for<'a> &'a T: Mul<&'a T, Output = T>,
{
    points.column_iter().zip(points.column_iter().cycle().skip(1))
        .map(|(v0, v1)| &v0.x * &v1.y - &v1.x * &v0.y)
        .sum()
}

/// Returns the orientation of the 2D polygon defined by the input points.
/// +1 for counterclockwise, -1 for clockwise
/// via computing sum of signed areas of triangles formed with origin
pub fn polygon_orientation<T>(points: MatrixView2Dyn<T>) -> i32 where
    T: Scalar + Sub<Output = T> + Sum + Sign,
    for<'a> &'a T: Mul<&'a T, Output = T>,
{
    match twice_signed_area(points).sign() {
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

#[cfg(test)]
mod test {
    use std::fmt::Debug;

    use exact_number::based_expr;
    use nalgebra::{matrix, vector, Affine2, ClosedSubAssign, Matrix2xX, Scalar, Vector2};

    use crate::geom::{polygon_orientation, reflect_line, reflect_line_matrix, sort_by_angle, sort_by_angle_field, sort_by_angle_ref, try_reflect_line, try_reflect_line_matrix, twice_signed_area, AngleRep, Angle, VectorView2Dyn};

    macro_rules! assert_lt {
        ($left:expr, $right:expr) => {
            match (&$left, &$right) {
                (left_val, right_val) => {
                    if !(*left_val < *right_val) {
                        panic!("assert_lt failed: {:?} is not less than {:?}", &*left_val, &*right_val);
                    }
                }
            }
        };
    }

    #[test]
    fn test_angle_rep() {
        assert_eq!(vector![based_expr!(0), based_expr!(0)].angle_rep(), Angle::new(0, Some(based_expr!(0))));
        assert_eq!(vector![based_expr!(0 + 0 sqrt 2), based_expr!(0 + 0 sqrt 2)].angle_rep(), Angle::new(0, Some(based_expr!(0 + 0 sqrt 2)))); // basis check
        assert_eq!(vector![based_expr!(1), based_expr!(0)].angle_rep(), Angle::new(0, Some(based_expr!(0))));
        assert_eq!(vector![based_expr!(1/16), based_expr!(0)].angle_rep(), Angle::new(0, Some(based_expr!(0))));
        assert_eq!(vector![based_expr!(50), based_expr!(0)].angle_rep(), Angle::new(0, Some(based_expr!(0))));

        assert_eq!(vector![based_expr!(13), based_expr!(-25)].angle_rep(), vector!(based_expr!(26), based_expr!(-50)).angle_rep());

        // Edge cases
        assert_lt!(vector![based_expr!(-1), based_expr!( 0)].angle_rep(), vector!(based_expr!(-1), based_expr!(-1)).angle_rep());
        assert_lt!(vector![based_expr!(-1), based_expr!(-1)].angle_rep(), vector!(based_expr!( 0), based_expr!(-1)).angle_rep());
        assert_lt!(vector![based_expr!( 0), based_expr!(-1)].angle_rep(), vector!(based_expr!( 1), based_expr!(-1)).angle_rep());
        assert_lt!(vector![based_expr!( 1), based_expr!(-1)].angle_rep(), vector!(based_expr!( 1), based_expr!( 0)).angle_rep());
        assert_lt!(vector![based_expr!( 1), based_expr!( 0)].angle_rep(), vector!(based_expr!( 1), based_expr!( 1)).angle_rep());
        assert_lt!(vector![based_expr!( 1), based_expr!( 1)].angle_rep(), vector!(based_expr!( 0), based_expr!( 1)).angle_rep());
        assert_lt!(vector![based_expr!( 0), based_expr!( 1)].angle_rep(), vector!(based_expr!(-1), based_expr!( 1)).angle_rep());
        //assert_lt!(vector![based_expr!(-1), based_expr!( 1)].angle_rep(), vector!(based_expr!(-1), based_expr!( 0)).angle_rep());

        // Other cases (each slice)
        assert_lt!(vector![based_expr!(-5), based_expr!(-1)].angle_rep(), vector!(based_expr!(-4), based_expr!(-2)).angle_rep());
        assert_lt!(vector![based_expr!(-4), based_expr!(-6)].angle_rep(), vector!(based_expr!(-2), based_expr!(-6)).angle_rep());
        assert_lt!(vector![based_expr!( 1), based_expr!(-5)].angle_rep(), vector!(based_expr!( 2), based_expr!(-4)).angle_rep());
        assert_lt!(vector![based_expr!( 6), based_expr!(-4)].angle_rep(), vector!(based_expr!( 6), based_expr!(-2)).angle_rep());
        assert_lt!(vector![based_expr!( 5), based_expr!( 1)].angle_rep(), vector!(based_expr!( 4), based_expr!( 2)).angle_rep());
        assert_lt!(vector![based_expr!( 4), based_expr!( 6)].angle_rep(), vector!(based_expr!( 2), based_expr!( 6)).angle_rep());
        assert_lt!(vector![based_expr!(-1), based_expr!( 5)].angle_rep(), vector!(based_expr!(-2), based_expr!( 4)).angle_rep());
        assert_lt!(vector![based_expr!(-6), based_expr!( 4)].angle_rep(), vector!(based_expr!(-6), based_expr!( 2)).angle_rep());
    }

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

    fn sort_by_angle_test<T, S>(mut points: Vec<T>, origin: T, mut mapping: impl FnMut(&T) -> VectorView2Dyn<S>, expected: Vec<T>) where
        T: Clone + Debug + PartialEq,
        S: Scalar + ClosedSubAssign,
        Vector2<S>: AngleRep
    {
        let (mut indexes, expected_indexes) = permutation(&points, &expected);
        let mut points_a = points.clone();
        let mut points_b = points.clone();
        points.push(origin.clone());
        sort_by_angle_ref(&mut indexes, &(points.len() - 1), |i| mapping(&points[*i]));
        assert_eq!(indexes, expected_indexes);
        sort_by_angle_field(&mut points_a, &origin, &mut mapping);
        assert_eq!(points_a, expected);
        sort_by_angle(&mut points_b, &origin, |t| mapping(t).clone_owned());
        assert_eq!(points_b, expected);
    }

    #[test]
    fn test_sort_by_angle() {
        sort_by_angle_test(vec![], Vector2::zeros(), |x: &Vector2<f64>| x.as_view(), vec![]);

        sort_by_angle_test(vec![vector!(0.0, 0.0)], Vector2::zeros(), |x| x.as_view(), vec![vector!(0.0, 0.0)]);
        sort_by_angle_test(vec![vector!(1.0, 0.0)], Vector2::zeros(), |x| x.as_view(), vec![vector!(1.0, 0.0)]);
        sort_by_angle_test(vec![vector!(1.0, 0.0), vector!(-1.0, 0.0), vector!(0.0, 1.0), vector!(0.0, -1.0)],
            Vector2::zeros(), |x| x.as_view(),
            vec![vector!(0.0, -1.0), vector!(1.0, 0.0), vector!(0.0, 1.0), vector!(-1.0, 0.0)]);
        sort_by_angle_test(vec![vector!(1.0, 0.1), vector!(-1.0, 0.0), vector!(0.0, 1.0), vector!(0.0, -1.0)],
            vector!(2.0, 0.0), |x| x.as_view(),
            vec![vector!(0.0, -1.0), vector!(0.0, 1.0), vector!(1.0, 0.1), vector!(-1.0, 0.0)]);
    }

    #[test]
    fn test_twice_signed_area() {
        assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
            based_expr!(0), based_expr!(0),
            based_expr!(2), based_expr!(0),
        ]).as_view()), based_expr!(0));
        assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
            based_expr!(0), based_expr!(0),
            based_expr!(0), based_expr!(0),
            based_expr!(0), based_expr!(0),
        ]).as_view()), based_expr!(0));
        assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
            based_expr!(5), based_expr!(9),
            based_expr!(2), based_expr!(7),
            based_expr!(1), based_expr!(7),
        ]).as_view()), based_expr!(-2));
        assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
            based_expr!(0), based_expr!(3),
            based_expr!(5), based_expr!(8),
            based_expr!(3), based_expr!(1),
        ]).as_view()), based_expr!(-25));
        assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
            4.5, 1.0, // Mixing it up
            4.5, 8.0,
            2.0, 5.0,
        ]).as_view()), 17.5);
        assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
            1.0, 7.0,
            0.0, 4.0,
            5.0, 9.0,
            2.0, 7.0,
        ]).as_view()), 8.0);
        assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
            0.0, 3.0,
            -1.0, 5.0,
            3.0, 1.0,
            5.0, 8.0
        ]).as_view()), 21.0);
        assert_eq!(twice_signed_area(Matrix2xX::from_vec(vec![
            based_expr!(9/2), based_expr!(1),
            based_expr!(9/2), based_expr!(8),
            based_expr!(2), based_expr!(5),
            based_expr!(7/2), based_expr!(3),
        ]).as_view()), based_expr!(33/2));
    }

    #[test]
    fn test_polygon_orientation() {
        assert_eq!(polygon_orientation(Matrix2xX::from_vec(vec![
            based_expr!(0), based_expr!(0),
            based_expr!(2), based_expr!(0),
        ]).as_view()), 0);
        assert_eq!(polygon_orientation(Matrix2xX::from_vec(vec![
            based_expr!(0), based_expr!(0),
            based_expr!(0), based_expr!(0),
            based_expr!(0), based_expr!(0),
        ]).as_view()), 0);
        assert_eq!(polygon_orientation(Matrix2xX::from_vec(vec![
            based_expr!(5), based_expr!(9),
            based_expr!(2), based_expr!(7),
            based_expr!(1), based_expr!(7),
        ]).as_view()), -1);
        assert_eq!(polygon_orientation(Matrix2xX::from_vec(vec![
            based_expr!(0), based_expr!(3),
            based_expr!(5), based_expr!(8),
            based_expr!(3), based_expr!(1),
        ]).as_view()), -1);
        assert_eq!(polygon_orientation(Matrix2xX::from_vec(vec![
            4.5, 1.0, // Mixing it up
            4.5, 8.0,
            2.0, 5.0,
        ]).as_view()), 1);
        assert_eq!(polygon_orientation(Matrix2xX::from_vec(vec![
            1.0, 7.0,
            0.0, 4.0,
            5.0, 9.0,
            2.0, 7.0,
        ]).as_view()), 1);
        assert_eq!(polygon_orientation(Matrix2xX::from_vec(vec![
            0.0, 3.0,
            -1.0, 5.0,
            3.0, 1.0,
            5.0, 8.0
        ]).as_view()), 1);
        assert_eq!(polygon_orientation(Matrix2xX::from_vec(vec![
            based_expr!(9/2), based_expr!(1),
            based_expr!(9/2), based_expr!(8),
            based_expr!(2), based_expr!(5),
            based_expr!(7/2), based_expr!(3),
        ]).as_view()), 1);
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
}