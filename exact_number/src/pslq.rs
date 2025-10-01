// Implementation of PSLQ from https://www.davidhbailey.com/dhbtalks/dhb-carma-20100824.pdf,

use dashu_float::FBig;
use dashu_float::ops::{SquareRoot, Abs};
use dashu_int::{IBig, Sign, UBig};
use dashu_int::ops::UnsignedAbs;
use malachite::base::num::basic::traits::Zero;
use malachite::Natural;
use malachite::{rational::Rational, Integer};
use nalgebra::{DMatrix, DVector};

use crate::algebraic::{check_combination, check_combination_target};
use crate::{SqrtExpr, SqrtExprSum};

#[derive(Clone, Debug)]
pub enum IntegerRelationError {
    ZeroInBasis,
    SizeOneBasis,
    DivisionByZero,
    OutOfPrecision,
    WasWrongGaveUp,
}

fn checked_div_v(a: FBig, b: FBig) -> Result<FBig, IntegerRelationError> {
    if b == FBig::<dashu_float::round::mode::Zero, 2>::ZERO { Err(IntegerRelationError::DivisionByZero)?; }
    Ok(a / b)
}

fn checked_div_r(a: &FBig, b: &FBig) -> Result<FBig, IntegerRelationError> {
    if b == &FBig::<dashu_float::round::mode::Zero, 2>::ZERO { Err(IntegerRelationError::DivisionByZero)?; }
    Ok(a / b)
}

fn checked_round(a: FBig) -> Result<FBig, IntegerRelationError> {
    // Skip the panicking
    #[cfg(debug_assertions)]
    {
        // Taken from internal code
        if a.repr().exponent() + (a.repr().digits_ub() as isize) >= -2 {
            let (fract, exponent) = if a.repr().exponent() + (a.repr().digits_ub() as isize) < -1 {
                (a.repr().significand().clone(), -(a.context().precision() as isize))
            } else {
                a.fract().into_repr().into_parts()
            };
            println!("int: {}, fract: {}, exp: {}, true exp: {}", a.trunc(), fract, exponent, a.repr().exponent());
            if exponent <= 0 && fract.clone().unsigned_abs() >= UBig::from_word(2).pow(-exponent as usize) { Err(IntegerRelationError::OutOfPrecision)?; }
        }
    }
    Ok(a.round())
}

/// Tries to find integer coefficients q_i such that 0 = sum_{i=0}^{|basis|} q_i*basis, where not all q_i are 0.
/// 0 is not allowed to appear in `basis`.
/// Fails if it can't find one with the precision it's given.
/// The result cannot be trusted as is and should be checked.
/// 
/// See https://www.davidhbailey.com/dhbtalks/dhb-carma-20100824.pdf
fn maybe_integer_relation(basis: &[FBig], gamma: FBig) -> Result<Vec<Integer>, IntegerRelationError>  {
    if basis.len() == 0 { return Ok(vec![]); }
    if basis.len() == 1 { return Err(IntegerRelationError::SizeOneBasis); }

    let precision = basis.iter().map(|b| b.precision()).max().unwrap();
    let zero = FBig::ZERO.with_precision(precision).value();
    let one = FBig::ONE.with_precision(precision).value();
    let two = &one + &one;
    let max_magnitude = two.powi((precision - 2).into());
    let detection = two.powi(32.into());

    // Initialization
    // Step 1
    let x = DVector::from_iterator(basis.len(), basis.iter().cloned());
    let mut A =
        DMatrix::<FBig>::repeat(basis.len(), basis.len(), zero.clone());
    for i in 0..basis.len() {
        A[(i, i)] = one.clone();
    }
    let mut B = A.clone();

    // Step 2
    let mut s = DVector::repeat(basis.len(), zero.clone());
    s[basis.len() - 1] = basis.last().unwrap().sqr();
    for k in (0..(basis.len() - 1)).rev() {
        s[k] = &s[k + 1] + basis[k].sqr();
    }
    for value in s.data.as_mut_slice().iter_mut() {
        *value = value.sqrt();
    }
    //println!("Step 2, s: {}", s);
    let t = checked_div_r(&one, &s[0])?;
    let mut y = x * t.clone();
    //println!("Step 2, y: {}", y);
    s *= t.clone();

    // Step 3: Initial H
    let mut H =
        DMatrix::<FBig>::repeat(basis.len(), basis.len() - 1, zero.clone());
    for j in 0..H.ncols() {
        H[(j, j)] = checked_div_r(&s[j + 1], &s[j])?;
        for i in (j + 1)..H.nrows() {
            H[(i, j)] = checked_div_v(-&y[i] * &y[j], &s[j] * &s[j + 1])?;
        }
    }
    //println!("Step 3, H: {}", H);

    // Step 4: Reduce H
    for i in 1..H.nrows() {
        for j in (0..=(i - 1)).rev() {
            let t = checked_round(checked_div_r(&H[(i, j)], &H[(j, j)])?)?;
            y[j] = &y[j] + &t * &y[i];
            for k in 0..=j {
                H[(i, k)] = &H[(i, k)] - &t * &H[(j, k)];
            }
            A.set_row(i, &(A.row(i) - A.row(j) * t.clone()));
            B.set_column(j, &(B.column(j) + B.column(i) * t.clone()));
        }
    }
    //println!("Step 4, H: {} A: {} B: {} y: {}", H, A, B, y);

    let mut prev_min_y = y.iter().cloned().enumerate().min_by_key(|(_, y)| y.clone().abs()).unwrap();
    let mut num_iters = 0;
    // Iteration
    loop {
        // Step 6: Termination test
        if A.iter().map(|a| a.clone().abs()).max().unwrap() > max_magnitude {
            return Err(IntegerRelationError::OutOfPrecision); // Precision exhausted
        }
        let min_y = y.iter().cloned().map(|a| a.abs()).enumerate().min_by_key(|(_, y)| y.clone()).unwrap();
        if prev_min_y.1 == zero || prev_min_y.1 > &min_y.1 * &detection {
            return Ok(B.column(min_y.0).iter().map(|b| {
                let int = b.round().to_int().unwrap();
                let (sign, words) = int.as_sign_words();
                let result = Integer::from_sign_and_abs(sign == Sign::Negative,
                    Natural::from_limbs_asc(words));
                result
            }).collect::<Vec<_>>());
        }
        prev_min_y = min_y;

        // Step 1: Find max
        let m = (0..(basis.len() - 1)).max_by_key(|i| H[(*i, *i)].clone().abs() * gamma.powi((*i).into())).unwrap();
        //println!("***Iteration {}***", num_iters);
        //println!("Step 1 (Find max), m: {}", m);

        // Step 2: Exchange
        y.swap_rows(m, m + 1);
        A.swap_rows(m, m + 1);
        H.swap_rows(m, m + 1);
        B.swap_columns(m, m + 1);
        //println!("Step 2 (Swap), H: {} A: {} B: {} y: {}", H, A, B, y);

        // Step 3: Remove corner on H diagonal
        if m < basis.len() - 2 {
            let t0 = (H[(m, m)].sqr() + H[(m, m + 1)].sqr()).sqrt();
            if t0 == zero { Err(IntegerRelationError::DivisionByZero)?; }
            let t1 = &H[(m, m)] / &t0;
            let t2 = &H[(m, m + 1)] / &t0;
            for i in m..basis.len() {
                let t3 = H[(i, m)].clone();
                let t4 = H[(i, m + 1)].clone();
                H[(i, m)] = &t1 * &t3 + &t2 * &t4;
                H[(i, m + 1)] = -&t2 * &t3 + &t1 * &t4;
            }
        }
        //println!("Step 3 (Remove corner), H: {}", H);

        // Step 4: Reduce H
        for i in (m + 1)..basis.len() {
            for j in (0..=(i - 1).min(m + 1)).rev() {
                if H[(j, j)] == zero { Err(IntegerRelationError::DivisionByZero)?; }
                let t = checked_round(checked_div_r(&H[(i, j)], &H[(j, j)])?)?;
                y[j] = &y[j] + &t * &y[i];
                for k in 0..=j {
                    H[(i, k)] = &H[(i, k)] - &t * &H[(j, k)];
                }
                A.set_row(i, &(A.row(i) - A.row(j) * t.clone()));
                B.set_column(j, &(B.column(j) + B.column(i) * t.clone()));
            }
        }
        //println!("Step 4 (Reduce), H: {} A: {} B: {} y: {}", H, A, B, y);

        // Step 5: Norm bound
        //let norm_bound = (0..(basis.len() - 1)).map(|j| &H[(j, j)]).max().unwrap();
        num_iters += 1;
    }
}

pub(crate) fn checked_integer_relation(basis: &[SqrtExpr]) -> Result<Vec<Integer>, IntegerRelationError> {
    let precision = 75;
    let fbig_basis = basis.iter().map(|b| b.to_fbig(precision)).collect::<Vec<_>>();

    let coeffs = maybe_integer_relation(&fbig_basis, FBig::<dashu_float::round::mode::Zero, 2>::from(2))?;
    if !check_combination(&coeffs, basis) {
        Err(IntegerRelationError::WasWrongGaveUp)?;
    }

    Ok(coeffs)
}

pub(crate) fn checked_integer_relation_target(basis: &[SqrtExpr], target: &SqrtExprSum) -> Result<Vec<Integer>, IntegerRelationError> {
    let precision = 75;
    let mut fbig_basis = basis.iter().map(|b| b.to_fbig(precision)).collect::<Vec<_>>();
    fbig_basis.push(terms_to_fbig(target, precision));

    let coeffs = maybe_integer_relation(&fbig_basis, FBig::<dashu_float::round::mode::Zero, 2>::from(2))?;
    if !check_combination_target(&coeffs, basis, target) {
        Err(IntegerRelationError::WasWrongGaveUp)?;
    }

    Ok(coeffs)
}

fn integer_to_fbig(int: &Integer, precision: usize) -> FBig {
    let words = int.unsigned_abs_ref().to_limbs_asc();
    let ibig = IBig::from_parts(if int < &Integer::ZERO {Sign::Negative} else {Sign::Positive}, UBig::from_words(&words));
    FBig::from(ibig).with_precision(precision).value()
}

fn terms_to_fbig(terms: &SqrtExprSum, precision: usize) -> FBig {
    terms.iter().map(|(coeff, sqrt)| integer_to_fbig(coeff, precision) * sqrt.to_fbig(precision)).sum::<FBig>().sqrt()
}

impl SqrtExpr {
    fn to_fbig(&self, precision: usize) -> FBig {
        match self {
            Self::Int(i) => integer_to_fbig(i, precision).sqrt(),
            Self::Sum(terms) => terms_to_fbig(terms, precision).sqrt()
        }
    }
}

#[cfg(test)]
mod test {
    use dashu_float::{round::Rounding, FBig};
    use dashu_float::round::mode;
    use malachite::Integer;

    use crate::pslq::maybe_integer_relation;
    use crate::{sqrt_expr, SqrtExpr};

    fn assert_integer_relation(input: Vec<SqrtExpr>, expected: Option<Vec<i64>>) {
        let basis = input.iter().map(|b| b.to_fbig(50)).collect::<Vec<_>>();
        assert_integer_relation_fbig(basis, expected)
    }

    fn assert_integer_relation_fbig(basis: Vec<FBig>, expected: Option<Vec<i64>>) {
        let expected = expected.map(|v| v.into_iter().map(|i| Integer::from(i)).collect::<Vec<_>>());
        let mut result = maybe_integer_relation(&basis, FBig::<mode::Zero, 2>::from(2));

        match expected {
            None => assert_eq!(result.ok(), None),

            Some(mut expected) => {
                let mut result = result.unwrap();
                if result.last().unwrap() > &Integer::from(0u64) {
                    result.iter_mut().for_each(|r| *r = -&*r );
                }
                assert_eq!(result, expected)
            }
        }
    }

    #[test]
    fn test_sqrt_to_fbig_int_perfect_square() {
        let input = sqrt_expr!(4);
        let result = input.to_fbig(256);
        let expected = FBig::<mode::Zero, 2>::from(2);
        assert_eq!(result, expected)
    }

    #[test]
    fn test_sqrt_to_fbig_int_not_perfect_square() {
        let input = sqrt_expr!(2);
        let result = input.to_fbig(256);
        let expected_start = "1.414213562373095048801688724209698078569671875376948073176679737990732";
        let start = &result.to_decimal().value().to_string()[0..expected_start.len()];
        assert_eq!(start, expected_start);
    }

    #[test]
    fn test_sqrt_to_fbig_complicated() {
        let input = sqrt_expr!(10 - 2 sqrt 5);
        let result = input.to_fbig(256);
        let expected_start = "2.35114100916989251667482381855629107439060975057258396428908992302911";
        let start = &result.to_decimal().value().to_string()[0..expected_start.len()];
        assert_eq!(start, expected_start);
    }

    #[test]
    fn test_relation_trivial() {
        assert_integer_relation(
            vec![sqrt_expr!(1), sqrt_expr!(1)],
            Some(vec![1, -1]),
        );
    }

    #[test]
    fn test_relation_simple_perpendicular() {
        assert_integer_relation(
            vec![sqrt_expr!(4), sqrt_expr!(25)],
            Some(vec![5, -2]),
        );
    }

    #[test]
    fn test_relation_false_sqrt() {
        assert_integer_relation(
            vec![sqrt_expr!(1), sqrt_expr!(2)],
            None
        );
    }

    #[test]
    fn test_relation_true_sqrt() {
        assert_integer_relation(
            vec![sqrt_expr!(1), sqrt_expr!(2), sqrt_expr!(18 - 8 sqrt 2)],
            Some(vec![4, -1, -1])
        );
    }

    #[test]
    fn test_relation_true_typoed_fractional_sqrt() {
        // Make sure this doesn't panic on a round(), 'cause it used to
        assert_integer_relation(
            vec![sqrt_expr!(4), sqrt_expr!(20), sqrt_expr!(15 - 8 sqrt 2)],
            None
        );
    }

    #[test]
    fn test_relation_true_fractional_sqrt() {
        assert_integer_relation(
            vec![sqrt_expr!(4), sqrt_expr!(20), sqrt_expr!(30 + 10 sqrt 5)],
            Some(vec![5, 1, -2])
        );
    }

    #[test]
    fn test_relation_unneeded_component() {
        assert_integer_relation(
            vec![sqrt_expr!(4), sqrt_expr!(20), sqrt_expr!(3), sqrt_expr!(30 + 10 sqrt 5)],
            Some(vec![5, 1, 0, -2])
        );
    }

    #[test]
    fn test_relation_tan_11_25_basis() {
        assert_integer_relation(
            vec![sqrt_expr!(1), sqrt_expr!(2), sqrt_expr!(2 - sqrt 2), sqrt_expr!(2 + sqrt 2)],
            None
        );
    }

    #[test]
    fn test_relation_tan_11_25_basis_tricky() {
        assert_integer_relation(
            vec![sqrt_expr!(1), sqrt_expr!(2), sqrt_expr!(2 - sqrt 2), sqrt_expr!(2 + sqrt 2), sqrt_expr!(4 + 2 sqrt 2)],
            Some(vec![0, 0, 1, 1, -1])
        );
    }

    #[test]
    fn test_relation_tan_18_basis_tricky() {
        assert_integer_relation(
            vec![sqrt_expr!(1), sqrt_expr!(5), sqrt_expr!(10 - 2 sqrt 5), sqrt_expr!(10 + 2 sqrt 5), sqrt_expr!(50 - 10 sqrt 5)],
            Some(vec![0, 0, -1, 2, -1])
        );
    }

    #[test]
    #[ignore]
    fn test_relation_frrt_3_minus_frrt_2() {
        // A classic.
        let fourth = FBig::<mode::Zero, 2>::from_parts(1.into(), -2);
        let precision = 350;
        let a = FBig::from(3).with_precision(precision).unwrap().powf(&fourth)
            + FBig::from(2).with_precision(precision).unwrap().powf(&fourth);
        let mut curr = FBig::from(1).with_precision(precision).unwrap();
        let powers = std::iter::from_fn(|| {
            let result = curr.clone();
            curr *= &a;
            Some(result)
        }).take(17).collect::<Vec<_>>();

        assert_integer_relation_fbig(
            powers,
            Some(vec![-1, 0, 0, 0, 3860, 0, 0, 0, 666, 0, 0, 0, 20, 0, 0, 0, -1]),
        );
    }
}