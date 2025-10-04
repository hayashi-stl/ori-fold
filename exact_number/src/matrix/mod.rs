// Some functions in Matrix claim to require ComplexField when they really don't.
// This fixes that.

mod lu;

use std::ops::Neg;

use nalgebra::{allocator::Allocator, constraint::{DimEq, ShapeConstraint}, DefaultAllocator, Dim, Matrix, OMatrix, Scalar, SquareMatrix, StorageMut};
use num::{traits::RefNum, Num};

use crate::rat::NumEx;

pub trait SquareMatrixEx {
    type T;
    type D: Dim;
    type S;

    /// Attempts to invert this square matrix in-place.
    ///
    /// Returns `true` if the inversion succeeded, `false` otherwise.
    ///
    /// # Behavior
    ///
    /// - For small dimensions (`n < 5`), the matrix is left unchanged if the inversion fails.
    /// - For dimensions `n >= 5`, the matrix may be **partially modified** even if the inversion fails,
    ///   because LU decomposition is used and it modifies the matrix in-place.
    ///
    /// If you need to preserve the original matrix regardless of success or failure,
    /// consider using [`Self::try_inverse`] instead.
    ///
    /// # Panics
    ///
    /// Panics if `self` isn’t a square matrix.
    fn try_inverse_mut(&mut self) -> bool where
        Self::T: NumEx,
        for<'a> &'a Self::T: RefNum<Self::T> + Neg<Output = Self::T>,
        Self::S: StorageMut<Self::T, Self::D, Self::D>,
        DefaultAllocator: Allocator<Self::D, Self::D>;

        /// Computes the solution of the linear system `self . x = b` where `x` is the unknown and only
        /// the lower-triangular part of `self` (including the diagonal) is considered not-zero.
        #[must_use = "Did you mean to use solve_lower_triangular_mut()?"]
        #[inline]
        pub fn solve_lower_triangular<R2: Dim, C2: Dim, S2>(
            &self,
            b: &Matrix<T, R2, C2, S2>,
        ) -> Option<OMatrix<T, R2, C2>>
        where
            S2: Storage<T, R2, C2>,
            DefaultAllocator: Allocator<R2, C2>,
            ShapeConstraint: SameNumberOfRows<R2, D>;

    
    /// Computes the solution of the linear system `self . x = b` where `x` is the unknown and only
    /// the upper-triangular part of `self` (including the diagonal) is considered not-zero.
    #[must_use = "Did you mean to use solve_upper_triangular_mut()?"]
    #[inline]
    pub fn solve_upper_triangular<R2: Dim, C2: Dim, S2>(
        &self,
        b: &Matrix<T, R2, C2, S2>,
    ) -> Option<OMatrix<T, R2, C2>>
    where
        S2: Storage<T, R2, C2>,
        DefaultAllocator: Allocator<R2, C2>,
        ShapeConstraint: SameNumberOfRows<R2, D>;

    /// Solves the linear system `self . x = b` where `x` is the unknown and only the
    /// lower-triangular part of `self` (including the diagonal) is considered not-zero.
    pub fn solve_lower_triangular_mut<R2: Dim, C2: Dim, S2>(
        &self,
        b: &mut Matrix<T, R2, C2, S2>,
    ) -> bool
    where
        S2: StorageMut<T, R2, C2>,
        ShapeConstraint: SameNumberOfRows<R2, D>;
}

impl<T, D: Dim, S> SquareMatrixEx for Matrix<T, D, D, S> {
    type T = T;
    type D = D;
    type S = S;

    fn try_inverse_mut(&mut self) -> bool where
            Self::T: NumEx,
            for<'a> &'a Self::T: RefNum<Self::T> + Neg<Output = T>,
            Self::S: StorageMut<Self::T, Self::D, Self::D>,
            DefaultAllocator: Allocator<Self::D, Self::D>
    {
        assert!(self.is_square(), "Unable to invert a non-square matrix.");

        let dim = self.shape().0;

        unsafe {
            match dim {
                0 => true,
                1 => {
                    let determinant = self.get_unchecked((0, 0));
                    if determinant.is_zero() {
                        false
                    } else {
                        *self.get_unchecked_mut((0, 0)) = T::one() / determinant;
                        true
                    }
                }
                2 => {
                    let m11 = std::mem::take(self.get_unchecked_mut((0, 0)));
                    let m12 = std::mem::take(self.get_unchecked_mut((0, 1)));
                    let m21 = std::mem::take(self.get_unchecked_mut((1, 0)));
                    let m22 = std::mem::take(self.get_unchecked_mut((1, 1)));

                    let determinant = &m11 * &m22 - &m21 * &m12;

                    if determinant.is_zero() {
                        false
                    } else {
                        *self.get_unchecked_mut((0, 0)) = m22 / &determinant;
                        *self.get_unchecked_mut((0, 1)) = -m12 / &determinant;

                        *self.get_unchecked_mut((1, 0)) = -m21 / &determinant;
                        *self.get_unchecked_mut((1, 1)) = m11 / determinant;

                        true
                    }
                }
                3 => {
                    let m11 = std::mem::take(self.get_unchecked_mut((0, 0)));
                    let m12 = std::mem::take(self.get_unchecked_mut((0, 1)));
                    let m13 = std::mem::take(self.get_unchecked_mut((0, 2)));

                    let m21 = std::mem::take(self.get_unchecked_mut((1, 0)));
                    let m22 = std::mem::take(self.get_unchecked_mut((1, 1)));
                    let m23 = std::mem::take(self.get_unchecked_mut((1, 2)));

                    let m31 = std::mem::take(self.get_unchecked_mut((2, 0)));
                    let m32 = std::mem::take(self.get_unchecked_mut((2, 1)));
                    let m33 = std::mem::take(self.get_unchecked_mut((2, 2)));

                    let minor_m12_m23 = &m22 * &m33 - &m32 * &m23;
                    let minor_m11_m23 = &m21 * &m33 - &m31 * &m23;
                    let minor_m11_m22 = &m21 * &m32 - &m31 * &m22;

                    let determinant = &m11 * &minor_m12_m23
                        - &m12 * &minor_m11_m23
                        + &m13 * &minor_m11_m22;

                    if determinant.is_zero() {
                        false
                    } else {
                        *self.get_unchecked_mut((0, 0)) = minor_m12_m23 / &determinant;
                        *self.get_unchecked_mut((0, 1)) = (&m13 * &m32
                            - &m33 * &m12)
                            / &determinant;
                        *self.get_unchecked_mut((0, 2)) = (&m12 * &m23
                            - &m22 * &m13)
                            / &determinant;

                        *self.get_unchecked_mut((1, 0)) = -minor_m11_m23 / &determinant;
                        *self.get_unchecked_mut((1, 1)) =
                            (&m11 * m33 - &m31 * &m13) / &determinant;
                        *self.get_unchecked_mut((1, 2)) =
                            (m13 * &m21 - m23 * &m11) / &determinant;

                        *self.get_unchecked_mut((2, 0)) = minor_m11_m22 / &determinant;
                        *self.get_unchecked_mut((2, 1)) =
                            (&m12 * m31 - m32 * &m11) / &determinant;
                        *self.get_unchecked_mut((2, 2)) = (m11 * m22 - m21 * m12) / determinant;

                        true
                    }
                }
                4 => {
                    let oself = self.clone_owned();
                    do_inverse4(&oself, self)
                }
                _ => {
                    unimplemented!()
                    //let oself = self.clone_owned();
                    //lu::try_invert_to(oself, self)
                }
            }
        } 
    }
}

// NOTE: this is an extremely efficient, loop-unrolled matrix inverse from MESA (MIT licensed).
fn do_inverse4<T: NumEx, D: Dim, S: StorageMut<T, D, D>>(
    m: &OMatrix<T, D, D>,
    out: &mut SquareMatrix<T, D, S>,
) -> bool
where
    for<'a> &'a T: RefNum<T> + Neg<Output = T>,
    DefaultAllocator: Allocator<D, D>,
{
    let m = m.as_slice();

    let cofactor00 = &m[5] * &m[10] * &m[15]
        - &m[5] * &m[11] * &m[14]
        - &m[9] * &m[6] * &m[15]
        + &m[9] * &m[7] * &m[14]
        + &m[13] * &m[6] * &m[11]
        - &m[13] * &m[7] * &m[10];

    let cofactor01 = -&m[4] * &m[10] * &m[15]
        + &m[4] * &m[11] * &m[14]
        + &m[8] * &m[6] * &m[15]
        - &m[8] * &m[7] * &m[14]
        - &m[12] * &m[6] * &m[11]
        + &m[12] * &m[7] * &m[10];

    let cofactor02 = &m[4] * &m[9] * &m[15]
        - &m[4] * &m[11] * &m[13]
        - &m[8] * &m[5] * &m[15]
        + &m[8] * &m[7] * &m[13]
        + &m[12] * &m[5] * &m[11]
        - &m[12] * &m[7] * &m[9];

    let cofactor03 = -&m[4] * &m[9] * &m[14]
        + &m[4] * &m[10] * &m[13]
        + &m[8] * &m[5] * &m[14]
        - &m[8] * &m[6] * &m[13]
        - &m[12] * &m[5] * &m[10]
        + &m[12] * &m[6] * &m[9];

    let det = &m[0] * &cofactor00
        + &m[1] * &cofactor01
        + &m[2] * &cofactor02
        + &m[3] * &cofactor03;

    if det.is_zero() {
        return false;
    }
    out[(0, 0)] = cofactor00;

    out[(1, 0)] = -&m[1] * &m[10] * &m[15]
        + &m[1] * &m[11] * &m[14]
        + &m[9] * &m[2] * &m[15]
        - &m[9] * &m[3] * &m[14]
        - &m[13] * &m[2] * &m[11]
        + &m[13] * &m[3] * &m[10];

    out[(2, 0)] = &m[1] * &m[6] * &m[15]
        - &m[1] * &m[7] * &m[14]
        - &m[5] * &m[2] * &m[15]
        + &m[5] * &m[3] * &m[14]
        + &m[13] * &m[2] * &m[7]
        - &m[13] * &m[3] * &m[6];

    out[(3, 0)] = -&m[1] * &m[6] * &m[11]
        + &m[1] * &m[7] * &m[10]
        + &m[5] * &m[2] * &m[11]
        - &m[5] * &m[3] * &m[10]
        - &m[9] * &m[2] * &m[7]
        + &m[9] * &m[3] * &m[6];

    out[(0, 1)] = cofactor01;

    out[(1, 1)] = &m[0] * &m[10] * &m[15]
        - &m[0] * &m[11] * &m[14]
        - &m[8] * &m[2] * &m[15]
        + &m[8] * &m[3] * &m[14]
        + &m[12] * &m[2] * &m[11]
        - &m[12] * &m[3] * &m[10];

    out[(2, 1)] = -&m[0] * &m[6] * &m[15]
        + &m[0] * &m[7] * &m[14]
        + &m[4] * &m[2] * &m[15]
        - &m[4] * &m[3] * &m[14]
        - &m[12] * &m[2] * &m[7]
        + &m[12] * &m[3] * &m[6];

    out[(3, 1)] = &m[0] * &m[6] * &m[11]
        - &m[0] * &m[7] * &m[10]
        - &m[4] * &m[2] * &m[11]
        + &m[4] * &m[3] * &m[10]
        + &m[8] * &m[2] * &m[7]
        - &m[8] * &m[3] * &m[6];

    out[(0, 2)] = cofactor02;

    out[(1, 2)] = -&m[0] * &m[9] * &m[15]
        + &m[0] * &m[11] * &m[13]
        + &m[8] * &m[1] * &m[15]
        - &m[8] * &m[3] * &m[13]
        - &m[12] * &m[1] * &m[11]
        + &m[12] * &m[3] * &m[9];

    out[(2, 2)] = &m[0] * &m[5] * &m[15]
        - &m[0] * &m[7] * &m[13]
        - &m[4] * &m[1] * &m[15]
        + &m[4] * &m[3] * &m[13]
        + &m[12] * &m[1] * &m[7]
        - &m[12] * &m[3] * &m[5];

    out[(0, 3)] = cofactor03;

    out[(3, 2)] = -&m[0] * &m[5] * &m[11]
        + &m[0] * &m[7] * &m[9]
        + &m[4] * &m[1] * &m[11]
        - &m[4] * &m[3] * &m[9]
        - &m[8] * &m[1] * &m[7]
        + &m[8] * &m[3] * &m[5];

    out[(1, 3)] = &m[0] * &m[9] * &m[14]
        - &m[0] * &m[10] * &m[13]
        - &m[8] * &m[1] * &m[14]
        + &m[8] * &m[2] * &m[13]
        + &m[12] * &m[1] * &m[10]
        - &m[12] * &m[2] * &m[9];

    out[(2, 3)] = -&m[0] * &m[5] * &m[14]
        + &m[0] * &m[6] * &m[13]
        + &m[4] * &m[1] * &m[14]
        - &m[4] * &m[2] * &m[13]
        - &m[12] * &m[1] * &m[6]
        + &m[12] * &m[2] * &m[5];

    out[(3, 3)] = &m[0] * &m[5] * &m[10]
        - &m[0] * &m[6] * &m[9]
        - &m[4] * &m[1] * &m[10]
        + &m[4] * &m[2] * &m[9]
        + &m[8] * &m[1] * &m[6]
        - &m[8] * &m[2] * &m[5];

    let inv_det = T::one() / det;

    for j in 0..4 {
        for i in 0..4 {
            out[(i, j)] *= &inv_det;
        }
    }
    true
}