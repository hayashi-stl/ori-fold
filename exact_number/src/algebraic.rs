// Algebraics. Used for checking solutions from the PSLQ algorithm for correctness.
// Unfortunately, now we're using 3 bigint libraries because they're all incompatible with each other
// and used by different necessary libraries
// * malachite::Integer (because malachite_q::Rational exists. And it's faster than num_rational. And it was first.)
// * dashu_int::IBig (because dashu_float is used for arbitrary precision floats for pslq )
// * num_bigint::BigInt (because algebraics, used for checking pslq solutions, uses it, and there's no alternative)

use std::cmp::Ordering;
use std::iter::once;

use algebraics::RealAlgebraicNumber;
use malachite::Integer;
use malachite::base::num::arithmetic::traits::{Sign as _};
use num::Zero;
use num::{bigint::Sign, BigInt, BigUint};
use num::traits::Pow;

use crate::{SqrtExpr, SqrtExprSum};

fn integer_to_algebraic(int: &Integer) -> RealAlgebraicNumber {
    let words = int.unsigned_abs_ref().limbs()
        .flat_map(|limb| [limb as u32, (limb >> 32) as u32])
        .collect::<Vec<_>>();
    let bigint = BigInt::from_slice(if int.sign() == Ordering::Less {Sign::Minus} else {Sign::Plus}, &words);
    bigint.into()
}

fn terms_to_algebraic(terms: &SqrtExprSum) -> RealAlgebraicNumber {
    dbg!(terms.iter().map(|(coeff, sqrt)| sqrt.to_algebraic() * integer_to_algebraic(coeff))
        .fold(RealAlgebraicNumber::zero(), |acc, curr| acc + curr))
}

impl SqrtExpr {
    fn to_algebraic(&self) -> RealAlgebraicNumber {
        match self {
            Self::Int(int) => integer_to_algebraic(int).pow((1, 2)),
            Self::Sum(terms) => terms_to_algebraic(terms).pow((1, 2)),
        }
    }
}

pub(crate) fn check_combination(coeffs: &[Integer], basis: &[SqrtExpr]) -> bool {
    coeffs.iter().zip(basis.iter())
        .map(|(coeff, b)| {
            let coeff = integer_to_algebraic(coeff);
            let b = b.to_algebraic();
            coeff * b
        })
        .fold(RealAlgebraicNumber::zero(), |acc, curr| acc + curr)
        .is_zero()
}

pub(crate) fn check_combination_product(coeffs: &[Integer], basis: &[SqrtExpr], factor_a: &SqrtExpr, factor_b: &SqrtExpr) -> bool {
    coeffs.iter().zip(basis.iter())
        .map(|(coeff, b)| {
            let coeff = integer_to_algebraic(coeff);
            let b = b.to_algebraic();
            coeff * b
        })
        .chain(once({
            let coeff = integer_to_algebraic(coeffs.last().unwrap());
            let b = factor_a.to_algebraic() * factor_b.to_algebraic();
            coeff * b
        }))
        .fold(RealAlgebraicNumber::zero(), |acc, curr| {
            let acc = dbg!(acc);
            let curr = dbg!(curr);
            dbg!(acc + curr)
        })
        .is_zero()
}

#[cfg(test)]
mod test {
    use algebraics::RealAlgebraicNumber;
    use num::traits::{Pow, Zero};

    use crate::{algebraic::{check_combination, check_combination_product}, sqrt_expr};

    #[test]
    fn test_to_algebraic_perfect_square() {
        let input = sqrt_expr!(4);
        let expected = RealAlgebraicNumber::from(2);
        let result = input.to_algebraic();
        assert_eq!(result, expected);
    }

    #[test]
    fn test_to_algebraic_not_perfect_square() {
        let input = sqrt_expr!(17);
        let expected = RealAlgebraicNumber::from(17).pow((1, 2));
        let result = input.to_algebraic();
        assert_eq!(result, expected);
    }

    #[test]
    fn test_to_algebraic_complicated() {
        let input = sqrt_expr!(5 - 2 sqrt 3);
        let expected = (RealAlgebraicNumber::from(5) - RealAlgebraicNumber::from(12).pow((1, 2))).pow((1, 2));
        let result = input.to_algebraic();
        assert_eq!(result, expected);
    }

    #[test]
    fn test_check_basis_trivial() {
        let basis = [sqrt_expr!(1), sqrt_expr!(1)];
        let coeffs = [1.into(), (-1).into()];
        assert!(check_combination(&coeffs, &basis));
    }

    #[test]
    fn test_check_basis_trivial_fail() {
        let basis = [sqrt_expr!(1), sqrt_expr!(1)];
        let coeffs = [1.into(), 1.into()];
        assert!(!check_combination(&coeffs, &basis));
    }

    #[test]
    fn test_check_basis_simple_perpendicular() {
        let basis = [sqrt_expr!(36), sqrt_expr!(49)];
        let coeffs = [7.into(), (-6).into()];
        assert!(check_combination(&coeffs, &basis));
    }

    #[test]
    fn test_check_basis_simple_perpendicular_fail() {
        let basis = [sqrt_expr!(36), sqrt_expr!(49)];
        let coeffs = [6.into(), (-7).into()];
        assert!(!check_combination(&coeffs, &basis));
    }

    #[test]
    fn test_check_basis_sqrt() {
        let basis = [sqrt_expr!(1), sqrt_expr!(2), sqrt_expr!(18 - 8 sqrt 2)];
        let coeffs = [4.into(), (-1).into(), (-1).into()];
        assert!(check_combination(&coeffs, &basis));
    }

    #[test]
    fn test_check_basis_product_past_fail() {
        let basis = vec![sqrt_expr!(1), sqrt_expr!(2), sqrt_expr!(2 + sqrt 2), sqrt_expr!(2 - sqrt 2)];
        let a = sqrt_expr!(2);
        let b = sqrt_expr!(2 - sqrt 2);
        let coeffs = [0.into(), 0.into(), 1.into(), (-1).into(), (-1).into()];
        assert!(check_combination_product(&coeffs, &basis, &a, &b));
    }

    #[test]
    fn test_algebraics_0_3_0_fail() {
        // This fails on algebraics version 0.3.0
        let two = RealAlgebraicNumber::from(2);
        let a = (&two + (&two).pow((1, 2))).pow((1, 2));
        let b = (&two - (&two).pow((1, 2))).pow((1, 2));
        let result = a - b;
        assert_ne!(result, RealAlgebraicNumber::zero());
    }
}