use std::{cmp::Ordering, iter::Sum, ops::{Mul, Sub}};

use exact_number::{malachite::base::num::arithmetic::traits::Sign, rat::Rat, BasedExpr};
use float_ord::FloatOrd;
use nalgebra::{ClosedSubAssign, Dyn, MatrixView2xX, Scalar, Vector2, VectorView2, U1};
use num_traits::{Num, Signed, Zero};

type VectorView2Dyn<'s, T> = VectorView2<'s, T, U1, Dyn>;
type MatrixView2Dyn<'s, T> = MatrixView2xX<'s, T, U1, Dyn>;

/// Convenience function for constructing `Vector2`
pub fn vec2<T: Scalar>(x: T, y: T) -> Vector2<T> {
    Vector2::from_vec(vec![x, y])
}

/// A trait that allows getting the angle of a 2D vector of this type.
/// This is often used for sorting by angle.
pub trait AngleRep: Sized {
    type Elem: Ord;
    /// The angle representative, representing the angle of this vector.
    /// 
    /// # Requirements
    /// * The `[0, 0]` vector must return 0.
    /// * The `[1, 0]` direction must return 0.
    /// * If vectors `[-1, 0]`, `a`, and `b` are ordered counterclockwise and none of them are `[0, 0]`, then
    ///     `angle_rep(a)` ≤ `angle_rep(b)` ≤ `angle_rep([-1, 0])`, with strict inequality
    ///     where vectors aren't pointing in the same direction.
    fn angle_rep(&self) -> Self::Elem;
}

impl AngleRep for Vector2<f32> {
    type Elem = FloatOrd<f32>;

    fn angle_rep(&self) -> Self::Elem {
        FloatOrd(self.y.atan2(self.x))
    }
}

impl AngleRep for Vector2<f64> {
    type Elem = FloatOrd<f64>;
    
    fn angle_rep(&self) -> Self::Elem {
        FloatOrd(self.y.atan2(self.x))
    }
}

impl AngleRep for Vector2<BasedExpr> {
    type Elem = BasedExpr;

    fn angle_rep(&self) -> Self::Elem {
        let x = &self.x;
        let y = &self.y;
        if x.is_zero() && y.is_zero() { return x.clone() }

        let (y, flip_y) = if y.is_negative() { (&-y, true) } else { (y, false) };
        let (x, flip_x) = if x.is_negative() { (&-x, true) } else { (x, false) };
        let (x, y, transpose) = if y > x { (y, x, true) } else { (x, y, false) };
        let mut result = y / x;
        if transpose { result = BasedExpr::Baseless(2.into()) - result }
        if flip_x { result = BasedExpr::Baseless(4.into()) - result }
        if flip_y { result = -result }
        result
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

#[cfg(test)]
mod test {
    use std::fmt::Debug;

    use exact_number::based_expr;
    use nalgebra::{ClosedSubAssign, Matrix2xX, Scalar, Vector2};

    use crate::geom::{polygon_orientation, sort_by_angle, sort_by_angle_field, sort_by_angle_ref, twice_signed_area, vec2, AngleRep, VectorView2Dyn};

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
        assert_eq!(vec2(based_expr!(0), based_expr!(0)).angle_rep(), based_expr!(0));
        assert_eq!(vec2(based_expr!(0 + 0 sqrt 2), based_expr!(0 + 0 sqrt 2)).angle_rep(), based_expr!(0 + 0 sqrt 2)); // basis check
        assert_eq!(vec2(based_expr!(1), based_expr!(0)).angle_rep(), based_expr!(0));
        assert_eq!(vec2(based_expr!(1/16), based_expr!(0)).angle_rep(), based_expr!(0));
        assert_eq!(vec2(based_expr!(50), based_expr!(0)).angle_rep(), based_expr!(0));

        assert_eq!(vec2(based_expr!(13), based_expr!(-25)).angle_rep(), vec2(based_expr!(26), based_expr!(-50)).angle_rep());

        // Edge cases
        assert_lt!(vec2(based_expr!(-1), based_expr!(-1)).angle_rep(), vec2(based_expr!( 0), based_expr!(-1)).angle_rep());
        assert_lt!(vec2(based_expr!( 0), based_expr!(-1)).angle_rep(), vec2(based_expr!( 1), based_expr!(-1)).angle_rep());
        assert_lt!(vec2(based_expr!( 1), based_expr!(-1)).angle_rep(), vec2(based_expr!( 1), based_expr!( 0)).angle_rep());
        assert_lt!(vec2(based_expr!( 1), based_expr!( 0)).angle_rep(), vec2(based_expr!( 1), based_expr!( 1)).angle_rep());
        assert_lt!(vec2(based_expr!( 1), based_expr!( 1)).angle_rep(), vec2(based_expr!( 0), based_expr!( 1)).angle_rep());
        assert_lt!(vec2(based_expr!( 0), based_expr!( 1)).angle_rep(), vec2(based_expr!(-1), based_expr!( 1)).angle_rep());
        assert_lt!(vec2(based_expr!(-1), based_expr!( 1)).angle_rep(), vec2(based_expr!(-1), based_expr!( 0)).angle_rep());

        // Other cases (each slice)
        assert_lt!(vec2(based_expr!(-5), based_expr!(-1)).angle_rep(), vec2(based_expr!(-4), based_expr!(-2)).angle_rep());
        assert_lt!(vec2(based_expr!(-4), based_expr!(-6)).angle_rep(), vec2(based_expr!(-2), based_expr!(-6)).angle_rep());
        assert_lt!(vec2(based_expr!( 1), based_expr!(-5)).angle_rep(), vec2(based_expr!( 2), based_expr!(-4)).angle_rep());
        assert_lt!(vec2(based_expr!( 6), based_expr!(-4)).angle_rep(), vec2(based_expr!( 6), based_expr!(-2)).angle_rep());
        assert_lt!(vec2(based_expr!( 5), based_expr!( 1)).angle_rep(), vec2(based_expr!( 4), based_expr!( 2)).angle_rep());
        assert_lt!(vec2(based_expr!( 4), based_expr!( 6)).angle_rep(), vec2(based_expr!( 2), based_expr!( 6)).angle_rep());
        assert_lt!(vec2(based_expr!(-1), based_expr!( 5)).angle_rep(), vec2(based_expr!(-2), based_expr!( 4)).angle_rep());
        assert_lt!(vec2(based_expr!(-6), based_expr!( 4)).angle_rep(), vec2(based_expr!(-6), based_expr!( 2)).angle_rep());
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

        sort_by_angle_test(vec![vec2(0.0, 0.0)], Vector2::zeros(), |x| x.as_view(), vec![vec2(0.0, 0.0)]);
        sort_by_angle_test(vec![vec2(1.0, 0.0)], Vector2::zeros(), |x| x.as_view(), vec![vec2(1.0, 0.0)]);
        sort_by_angle_test(vec![vec2(1.0, 0.0), vec2(-1.0, 0.0), vec2(0.0, 1.0), vec2(0.0, -1.0)],
            Vector2::zeros(), |x| x.as_view(),
            vec![vec2(0.0, -1.0), vec2(1.0, 0.0), vec2(0.0, 1.0), vec2(-1.0, 0.0)]);
        sort_by_angle_test(vec![vec2(1.0, 0.1), vec2(-1.0, 0.0), vec2(0.0, 1.0), vec2(0.0, -1.0)],
            vec2(2.0, 0.0), |x| x.as_view(),
            vec![vec2(0.0, -1.0), vec2(0.0, 1.0), vec2(1.0, 0.1), vec2(-1.0, 0.0)]);
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
}