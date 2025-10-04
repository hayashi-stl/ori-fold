// Implementation of PSLQ from https://www.davidhbailey.com/dhbtalks/dhb-carma-20100824.pdf,

use dashu_float::{DBig};
use dashu_float::ops::SquareRoot;
use dashu_int::{IBig, Sign, UBig};
use dashu_int::ops::UnsignedAbs;
use malachite::base::num::basic::traits::{NegativeOne, One, Zero};
use malachite::base::num::arithmetic::traits::{Abs as _, FloorLogBase2, FloorSqrt, Pow, PowerOf2, Square, CheckedDiv};
use malachite::Natural;
use malachite::Integer;
use nalgebra::allocator::Allocator;
use nalgebra::{DMatrix, DVector, DefaultAllocator, Dim, Matrix, Storage, VecStorage};

use crate::algebraic::{check_combination, check_combination_product};
use crate::{SqrtExpr, SqrtExprSum};
use crate::rat::Rat;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IntegerRelationError {
    TinyNumberInBasis(usize), // also returns the index
    SizeOneBasis,
    DivisionByZero,
    OutOfPrecision,
    ErrorNormTooBig,
    OutOfIterations,
    WasWrongGaveUp,
}

type FBig = dashu_float::FBig<dashu_float::round::mode::Zero, 2>;

fn checked_div_v(a: FBig, b: FBig) -> Result<FBig, IntegerRelationError> {
    if b == FBig::ZERO { Err(IntegerRelationError::DivisionByZero)?; }
    Ok(a / b)
}

fn checked_div_r(a: &FBig, b: &FBig) -> Result<FBig, IntegerRelationError> {
    if b == &FBig::ZERO { Err(IntegerRelationError::DivisionByZero)?; }
    Ok(a / b)
}

fn checked_round(a: FBig) -> Result<FBig, IntegerRelationError> {
    // Skip the panicking
    //print!("Rounding {}, ", a);
    let (mut sig, exp) = (a + FBig::from_parts(1.into(), -1)).into_repr().into_parts();
    if exp > 0 {
        sig <<= exp as usize;
    } else {
        sig >>= -exp as usize;
    }
    let result = FBig::from_parts(sig, 0);
    //println!("Round result: {}", result);
    Ok(result)
    //#[cfg(debug_assertions)]
    //{
    //    // Taken from internal code
    //    if a.repr().exponent() + (a.repr().digits_ub() as isize) >= -2 {
    //        let (fract, exponent) = if a.repr().exponent() + (a.repr().digits_ub() as isize) < -1 {
    //            (a.repr().significand().clone(), -(a.context().precision() as isize))
    //        } else {
    //            a.fract().into_repr().into_parts()
    //        };
    //        println!("int: {}, fract: {}, exp: {}, true exp: {}", a.trunc(), fract, exponent, a.repr().exponent());
    //        if exponent <= 0 && fract.clone().unsigned_abs() >= UBig::from_word(2).pow(-exponent as usize) { Err(IntegerRelationError::FloatRoundGlitched)?; }
    //    }
    //}
    //Ok(a.round())
}

fn dec<R: Dim, C: Dim>(mtx: &Matrix<Integer, R, C, VecStorage<Integer, R, C>>, precision: usize) -> Matrix<DBig, R, C, VecStorage<DBig, R, C>> where
    DefaultAllocator: Allocator<R, C, Buffer<DBig> = VecStorage<DBig, R, C>>,
    VecStorage<Integer, R, C>: Storage<Integer, R, C>,
{
    mtx.map(|x| fixed_to_dbig(&x, precision))
}

fn to_fixed(fbig: &FBig, precision: usize) -> Integer {
    let repr = fbig.repr();
    let (sign, significand) = repr.significand().as_sign_words();
    let significand = Integer::from_sign_and_abs(sign == Sign::Positive, Natural::from_limbs_asc(significand));
    let shift = precision as isize + repr.exponent();
    if shift < 0 {
        significand >> -shift
    } else {
        significand << shift
    }
}

fn fixed_to_dbig(x: &Integer, precision: usize) -> DBig {
    let significand = x.unsigned_abs_ref();
    let significand_size = if significand == &Natural::ZERO { 4 } else { significand.floor_log_base_2() + 4 };
    let ubig = UBig::from_words(&significand.to_limbs_asc());
    let ibig  = IBig::from_parts(if x < &Integer::ZERO { Sign::Negative } else { Sign::Positive }, ubig);
    let float = FBig::from_parts(ibig, -(precision as isize));
    // Unlimited precision happens if float is 0
    let float = float.with_precision(significand_size as usize).value();
    float.to_decimal().value()
}

fn sqrt_fixed(x: Integer, precision: usize) -> Integer {
    (x << precision).floor_sqrt()
}

fn sqrt_fixed_r(x: &Integer, precision: usize) -> Integer {
    (x << precision).floor_sqrt()
}

fn round_fixed(x: Integer, precision: usize) -> Integer {
    return ((x + (Integer::power_of_2(precision as u64 - 1))) >> precision) << precision
}

fn round_fixed_r(x: &Integer, precision: usize) -> Integer {
    return ((x + (Integer::power_of_2(precision as u64 - 1))) >> precision) << precision
}

/// Tries to find integer coefficients q_i such that 0 = sum_{i=0}^{|basis|} q_i*basis, where not all q_i are 0.
/// 0 is not allowed to appear in `basis`.
/// Fails if it can't find one with the precision it's given.
/// The result cannot be trusted as is and should be checked, if the basis came from exact numbers.
/// 
/// Taken from https://github.com/mpmath/mpmath/blob/33eace1005bc0a6d5330104cd1e8fa5fbab0aee3/mpmath/identification.py#L19
/// See https://www.davidhbailey.com/dhbtalks/dhb-carma-20100824.pdf
fn maybe_integer_relation(basis: &[FBig], max_coeff: &Integer, max_steps: usize) -> Result<Vec<Integer>, IntegerRelationError>  {
    const VERBOSE: bool = false;
    if basis.len() == 0 { return Ok(vec![]); }
    if basis.len() == 1 {
        return if basis[0] == FBig::ZERO { Ok(vec![Integer::ONE]) } else { Err(IntegerRelationError::SizeOneBasis) };
    }

    let precision = basis.iter().map(|b| b.precision()).max().unwrap();
    // At too low precision, the algorithm becomes meaningless
    if precision < 53 { return Err(IntegerRelationError::OutOfPrecision); }
    // Also, precision should probably be at least 5 * max(2, basis.len()).

    let target = precision * 3 / 4;
    let tolerance = FBig::from_parts(IBig::ONE, -(target as isize));
    let extra = 60;
    let precision = precision + extra;
    let tolerance = to_fixed(&tolerance, precision);

    if VERBOSE {
        println!("PSLQ using prec {} and tol {}", precision, fixed_to_dbig(&tolerance, precision));
    }

    // Convert to fixed-point numbers. The dummy None is added so we can
    // use 1-based indexing. (This just allows us to be consistent with
    // Bailey's indexing. The algorithm is 100 lines long, so debugging
    // a single wrong index can be painful.)
    let x = DVector::from_iterator(basis.len() + 1, std::iter::once(Integer::ZERO).chain(
        basis.iter().map(|xk| to_fixed(xk, precision))
    ));

    // Sanity check the magnitudes
    let (min_x_i, min_x) = x.iter().skip(1).enumerate().min_by_key(|(_, xx)| xx.abs()).unwrap();
    let min_x = min_x.abs();
    if min_x == Integer::ZERO {
        let mut result = vec![Integer::ZERO; basis.len()];
        result[min_x_i] = Integer::ONE;
        return Ok(result); // that's what happens when you have a zero in your basis!
    }
    if min_x < &tolerance / Integer::from(100) { return Err(IntegerRelationError::TinyNumberInBasis(min_x_i)); }

    let gamma = sqrt_fixed(Integer::power_of_2(precision as u64) * Integer::from(4) / Integer::from(3), precision);
    let mut A = DMatrix::repeat(basis.len() + 1, basis.len() + 1, Integer::ZERO);
    let mut B = DMatrix::repeat(basis.len() + 1, basis.len() + 1, Integer::ZERO);
    let mut H = DMatrix::repeat(basis.len() + 1, basis.len() + 1, Integer::ZERO);
    
    // Initialization
    // step 1
    for i in 1..(basis.len() + 1) {
        for j in 1..(basis.len() + 1) {
            A[(i, j)] = Integer::from(i == j) << precision;
            B[(i, j)] = Integer::from(i == j) << precision;
            H[(i, j)] = Integer::ZERO;
        }
    }

    // step 2
    let mut s = DVector::repeat(basis.len() + 1, Integer::ZERO);
    for k in 1..(basis.len() + 1) {
        let mut t = Integer::ZERO;
        for j in k..(basis.len() + 1) {
            t += (&x[j]).square() >> precision;
        }
        s[k] = sqrt_fixed(t, precision);
    }

    let mut y = x.clone();
    let t = s[1].clone();
    for k in 1..(basis.len() + 1) {
        y[k] = (&x[k] << precision) / &t;
        s[k] = (&s[k] << precision) / &t;
    }
    if VERBOSE {
        println!("Init. step 2: s: {}", dec(&s, precision));
    }

    // step 3
    for i in 1..(basis.len() + 1) {
        for j in (i + 1)..basis.len() {
            H[(i, j)] = Integer::ZERO;
        }
        if i <= basis.len() - 1 {
            if s[i] != Integer::ZERO {
                H[(i, i)] = (&s[i + 1] << precision) / &s[i];
            } else {
                H[(i, i)] = Integer::ZERO;
            }
        }
        for j in 1..i {
            let sjj1 = &s[j] * &s[j + 1];
            if sjj1 != Integer::ZERO {
                H[(i, j)] = ((-&y[i] * &y[j]) << precision) / sjj1;
            } else {
                H[(i, j)] = Integer::ZERO;
            }
        }
    }
    if VERBOSE {
        println!("Init. step 3: H: {}, A: {}, B: {}, y: {}", dec(&H, precision), dec(&A, precision), dec(&B, precision), dec(&y, precision));
    }

    // step 4
    for i in 2..(basis.len() + 1) {
        for j in (1..=(i - 1)).rev() {
            let t = if H[(j, j)] != Integer::ZERO {
                round_fixed((&H[(i, j)] << precision) / &H[(j, j)], precision)
            } else {
                continue;
            };
            y[j] = &y[j] + (&t * &y[i] >> precision);
            for k in 1..(j + 1) {
                H[(i, k)] = &H[(i, k)] - (&t * &H[(j, k)] >> precision);
            }
            for k in 1..(basis.len() + 1) {
                A[(i, k)] = &A[(i, k)] - (&t * &A[(j, k)] >> precision);
                B[(k, j)] = &B[(k, j)] + (&t * &B[(k, i)] >> precision);
            }
        }
    }

    if VERBOSE {
        println!("Init. step 4: H: {}, A: {}, B: {}, y: {}", dec(&H, precision), dec(&A, precision), dec(&B, precision), dec(&y, precision));
    }

    // Main algorithm
    let mut curr_iter = 0;
    for rep in 0..max_steps {
        curr_iter = rep;
        // Until a relation is found, the error typically decreases
        // slowly (e.g. a factor 1-10) with each step TODO: we could
        // compare err from two successive iterations. If there is a
        // large drop (several orders of magnitude), that indicates a
        // "high quality" relation was detected. Reporting this to
        // the user somehow might be useful.
        let mut best_error = max_coeff << precision;
        for i in 1..(basis.len() + 1) {
            let err = (&y[i]).abs();
            // Maybe we are done?
            if err < tolerance {
                // We are done if the coefficients are acceptable
                let vector = (1..(basis.len() + 1)).map(|j| round_fixed_r(&B[(j, i)], precision) >> precision).collect::<Vec<_>>();
                if &vector.iter().map(|v| v.abs()).max().unwrap() < max_coeff {
                    return Ok(vector);
                }
            }
            best_error = err.min(best_error);
        }

        // Calculate a lower bound for the norm. We could do this
        // more exactly (using the Euclidean norm) but there is probably
        // no practical benefit.
        let recnorm = H.iter().map(|h| h.abs()).max().unwrap();
        let norm = if recnorm != Integer::ZERO {
            ((Integer::power_of_2(2 * precision as u64) / recnorm) >> precision) / Integer::from(100)
        } else {
            if VERBOSE {
                println!("{}/{}:  Error: {} Norm: infinity",
                    rep, max_steps, fixed_to_dbig(&best_error, precision));
            }
            return Err(IntegerRelationError::ErrorNormTooBig);
        };
        if VERBOSE {
            println!("{}/{}:  Error: {}   Norm: {}",
                rep, max_steps, fixed_to_dbig(&best_error, precision), norm);
        }
        if &norm >= max_coeff {
            return Err(IntegerRelationError::ErrorNormTooBig);
        }

        // Step 1
        let mut m = 0;
        let mut sz_max = Integer::NEGATIVE_ONE;
        for i in 1..basis.len() {
            let sz = ((&gamma).pow(i as u64) * (&H[(i, i)]).abs()) >> (precision * (i - 1));
            if sz > sz_max {
                m = i;
                sz_max = sz;
            }
        }
        if VERBOSE {
            println!("***ITERATION {}***", rep);
            println!("Main step 1 (max): m: {}", m);
        }

        // Step 2
        y.swap_rows(m, m + 1);
        H.swap_rows(m, m + 1);
        A.swap_rows(m, m + 1);
        B.swap_columns(m, m + 1);
        if VERBOSE {
            println!("Main. step 2 (swap): H: {}, A: {}, B: {}, y: {}", dec(&H, precision), dec(&A, precision), dec(&B, precision), dec(&y, precision));
        }
        
        // Step 3
        if m <= basis.len() - 2 {
            let t0 = sqrt_fixed(((&H[(m, m)]).square() + (&H[(m, m + 1)]).square()) >> precision, precision);
            // A zero element probably indicates that the precision has
            // been exhausted. XXX: this could be spurious, due to
            // using fixed-point arithmetic
            if t0 == Integer::ZERO { return Err(IntegerRelationError::DivisionByZero); }
            let t1 = (&H[(m, m)] << precision) / &t0;
            let t2 = (&H[(m, m + 1)] << precision) / &t0;
            for i in m..(basis.len() + 1) {
                let t3 = H[(i, m)].clone();
                let t4 = H[(i, m + 1)].clone();
                H[(i, m)] = (&t1 * &t3 + &t2 * &t4) >> precision;
                H[(i, m + 1)] = (-&t2 * &t3 + &t1 * &t4) >> precision;
            }
        }
        if VERBOSE {
            println!("Main. step 3 (corner): H: {}", dec(&H, precision));
        }

        // Step 4
        for i in (m + 1)..(basis.len() + 1) {
            for j in (1..=(i - 1).min(m + 1)).rev() {
                let t = (&H[(i, j)] << precision).checked_div(&H[(j, j)])
                    .map(|t| round_fixed(t, precision));
                let t = if let Some(t) = t { t } else { return Err(IntegerRelationError::DivisionByZero) };
                y[j] = &y[j] + ((&t * &y[i]) >> precision);
                for k in 1..(j + 1) {
                    H[(i, k)] = &H[(i, k)] - (&t * &H[(j, k)] >> precision);
                }
                for k in 1..(basis.len() + 1) {
                    A[(i, k)] = &A[(i, k)] - (&t * &A[(j, k)] >> precision);
                    B[(k, j)] = &B[(k, j)] + (&t * &B[(k, i)] >> precision);
                }
            }
        }
        if VERBOSE {
            println!("Main. step 4 (reduce): H: {}, A: {}, B: {}, y: {}", dec(&H, precision), dec(&A, precision), dec(&B, precision), dec(&y, precision));
        }
    }

    if VERBOSE {
        println!("CANCELLING after step {}/{}.", curr_iter, max_steps);
    }
    Err(IntegerRelationError::OutOfIterations)
}
//fn maybe_integer_relation(basis: &[FBig], gamma: FBig) -> Result<Vec<Integer>, IntegerRelationError>  {
//    if basis.len() == 0 { return Ok(vec![]); }
//    if basis.len() == 1 { return Err(IntegerRelationError::SizeOneBasis); }
//
//    let precision = basis.iter().map(|b| b.precision()).max().unwrap();
//    let zero = FBig::ZERO.with_precision(precision).value();
//    let one = FBig::ONE.with_precision(precision).value();
//    let two = &one + &one;
//    let max_magnitude = two.powi((precision - 2).into());
//    let detection = two.powi(32.into());
//    // Used for detecting coefficients that get discovered before an iteration is even performed.
//    let epsilon = two.powi((-(precision as isize) + 2).into());
//
//    // Initialization
//    // Step 1
//    let x = DVector::from_iterator(basis.len(), basis.iter().cloned());
//    let mut A =
//        DMatrix::<FBig>::repeat(basis.len(), basis.len(), zero.clone());
//    for i in 0..basis.len() {
//        A[(i, i)] = one.clone();
//    }
//    let mut B = A.clone();
//
//    // Step 2
//    let mut s = DVector::repeat(basis.len(), zero.clone());
//    s[basis.len() - 1] = basis.last().unwrap().sqr();
//    for k in (0..(basis.len() - 1)).rev() {
//        s[k] = &s[k + 1] + basis[k].sqr();
//    }
//    for value in s.data.as_mut_slice().iter_mut() {
//        *value = value.sqrt();
//    }
//    println!("Step 2, s: {}", dec(&s));
//    let t = checked_div_r(&one, &s[0])?;
//    let mut y = x * t.clone();
//    println!("Step 2, y: {}", dec(&y));
//    s *= t.clone();
//
//    // Step 3: Initial H
//    let mut H =
//        DMatrix::<FBig>::repeat(basis.len(), basis.len() - 1, zero.clone());
//    for j in 0..H.ncols() {
//        H[(j, j)] = checked_div_r(&s[j + 1], &s[j])?;
//        for i in (j + 1)..H.nrows() {
//            H[(i, j)] = checked_div_v(-&y[i] * &y[j], &s[j] * &s[j + 1])?;
//        }
//    }
//    println!("Step 3, H: {}", dec(&H));
//
//    // Step 4: Reduce H
//    for i in 1..H.nrows() {
//        for j in (0..=(i - 1)).rev() {
//            let t = checked_round(checked_div_r(&H[(i, j)], &H[(j, j)])?)?;
//            y[j] = &y[j] + &t * &y[i];
//            for k in 0..=j {
//                H[(i, k)] = &H[(i, k)] - &t * &H[(j, k)];
//            }
//            A.set_row(i, &(A.row(i) - A.row(j) * t.clone()));
//            B.set_column(j, &(B.column(j) + B.column(i) * t.clone()));
//        }
//    }
//    println!("Step 4, H: {} A: {} B: {} y: {}, eps: {}", dec(&H), dec(&A), dec(&B), dec(&y), epsilon.to_decimal().value());
//
//    let mut prev_min_y = y.iter().cloned().enumerate().min_by_key(|(_, y)| y.clone().abs()).unwrap();
//    let mut num_iters = 0;
//    // Iteration
//    loop {
//        //TODO: Remove
//        if num_iters > 100 {
//            panic!("Too many iterations!");
//        }
//        // Step 6: Termination test
//        if A.iter().map(|a| a.clone().abs()).max().unwrap() > max_magnitude {
//            return Err(IntegerRelationError::OutOfPrecision); // Precision exhausted
//        }
//        let min_y = y.iter().cloned().map(|a| a.abs()).enumerate().min_by_key(|(_, y)| y.clone()).unwrap();
//        println!("Step 6, min_y i: {} min_y: {} eps: {}, min_y < eps: {}", min_y.0, min_y.1.to_decimal().value(), epsilon.to_decimal().value(), min_y.1 < epsilon);
//        if (num_iters == 0 && min_y.1 < epsilon) || prev_min_y.1 > &min_y.1 * &detection {
//            return Ok(B.column(min_y.0).iter().map(|b| {
//                let int = b.round().to_int().unwrap();
//                let (sign, words) = int.as_sign_words();
//                let result = Integer::from_sign_and_abs(sign == Sign::Negative,
//                    Natural::from_limbs_asc(words));
//                result
//            }).collect::<Vec<_>>());
//        }
//        prev_min_y = min_y;
//
//        // Step 1: Find max
//        let m = (0..(basis.len() - 1)).max_by_key(|i| H[(*i, *i)].clone().abs() * gamma.powi((*i).into())).unwrap();
//        println!("***Iteration {}***", num_iters);
//        println!("Step 1 (Find max), m: {}", m);
//
//        // Step 2: Exchange
//        y.swap_rows(m, m + 1);
//        A.swap_rows(m, m + 1);
//        H.swap_rows(m, m + 1);
//        B.swap_columns(m, m + 1);
//        println!("Step 2 (Swap), H: {} A: {} B: {} y: {}", dec(&H), dec(&A), dec(&B), dec(&y));
//
//        // Step 3: Remove corner on H diagonal
//        if m < basis.len() - 2 {
//            let t0 = (H[(m, m)].sqr() + H[(m, m + 1)].sqr()).sqrt();
//            if t0 == zero { Err(IntegerRelationError::DivisionByZero)?; }
//            let t1 = &H[(m, m)] / &t0;
//            let t2 = &H[(m, m + 1)] / &t0;
//            for i in m..basis.len() {
//                let t3 = H[(i, m)].clone();
//                let t4 = H[(i, m + 1)].clone();
//                H[(i, m)] = &t1 * &t3 + &t2 * &t4;
//                H[(i, m + 1)] = -&t2 * &t3 + &t1 * &t4;
//            }
//        }
//        println!("Step 3 (Remove corner), H: {}", dec(&H));
//
//        // Step 4: Reduce H
//        for i in (m + 1)..basis.len() {
//            for j in (0..=(i - 1).min(m + 1)).rev() {
//                if H[(j, j)] == zero { Err(IntegerRelationError::DivisionByZero)?; }
//                let t = checked_round(checked_div_r(&H[(i, j)], &H[(j, j)])?)?;
//                y[j] = &y[j] + &t * &y[i];
//                for k in 0..=j {
//                    H[(i, k)] = &H[(i, k)] - &t * &H[(j, k)];
//                }
//                A.set_row(i, &(A.row(i) - A.row(j) * t.clone()));
//                B.set_column(j, &(B.column(j) + B.column(i) * t.clone()));
//            }
//        }
//        println!("Step 4 (Reduce), H: {} A: {} B: {} y: {}", dec(&H), dec(&A), dec(&B), dec(&y));
//
//        // Step 5: Norm bound
//        let norm_bound = FBig::ONE / (0..(basis.len() - 1)).map(|j| &H[(j, j)]).max().unwrap();
//        println!("Step 5 (Norm Bound), norm_bound: {}", norm_bound.to_decimal().value());
//        num_iters += 1;
//    }
//}

pub(crate) fn checked_integer_relation(basis: &[SqrtExpr]) -> Result<Vec<Integer>, IntegerRelationError> {
    let precision = 75;// + basis.len().saturating_sub(5) * 20; // A guess
    let fbig_basis = basis.iter().map(|b| b.to_fbig(precision)).collect::<Vec<_>>();

    let coeffs = maybe_integer_relation(&fbig_basis, &1000.into(), 100)?;
    if !check_combination(&coeffs, basis) {
        Err(IntegerRelationError::WasWrongGaveUp)?;
    }

    Ok(coeffs)
}

/// Returns `factor_a` * `factor_b` as a rational linear combination of `basis`
pub(crate) fn checked_integer_relation_product(basis: &[SqrtExpr], factor_a: &SqrtExpr, factor_b: &SqrtExpr) -> Result<Vec<Rat>, IntegerRelationError> {
    let precision = 75;// + basis.len().saturating_sub(4) * 20;
    let mut fbig_basis = basis.iter().map(|b| b.to_fbig(precision)).collect::<Vec<_>>();
    fbig_basis.push(factor_a.to_fbig(precision) * factor_b.to_fbig(precision));
    //println!("fbig_basis: [{}]", fbig_basis.iter().map(|f| format!("{}, ", f.to_decimal().value())).collect::<String>());

    let mut coeffs = maybe_integer_relation(&fbig_basis, &1000.into(), 100)?;
    if !check_combination_product(&coeffs, basis, factor_a, factor_b) {
        Err(IntegerRelationError::WasWrongGaveUp)?;
    }

    let divisor = Rat::from(-coeffs.pop().unwrap());
    let coeffs = coeffs.into_iter().map(|c| Rat::from(c) / &divisor).collect::<Vec<_>>();
    Ok(coeffs)
}

fn integer_to_fbig(int: &Integer, precision: usize) -> FBig {
    let words = int.unsigned_abs_ref().to_limbs_asc();
    let ibig = IBig::from_parts(if int < &Integer::ZERO {Sign::Negative} else {Sign::Positive}, UBig::from_words(&words));
    FBig::from(ibig).with_precision(precision).value()
}

fn terms_to_fbig(terms: &SqrtExprSum, precision: usize) -> FBig {
    terms.iter().map(|(coeff, sqrt)| integer_to_fbig(coeff, precision) * sqrt.to_fbig(precision)).sum::<FBig>()
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
    use std::ops::BitAndAssign;

    use dashu_float::{round::Rounding, FBig};
    use dashu_float::round::mode;
    use malachite::Integer;

    use crate::pslq::{checked_integer_relation_product, maybe_integer_relation};
    use crate::{sqrt_expr, SqrtExpr};
    use crate::rat::Rat;

    fn assert_integer_relation(input: Vec<SqrtExpr>, expected: Option<Vec<i64>>) {
        let basis = input.iter().map(|b| b.to_fbig(53)).collect::<Vec<_>>();
        assert_integer_relation_fbig(basis, expected)
    }

    fn assert_integer_relation_product(input: Vec<SqrtExpr>, a: SqrtExpr, b: SqrtExpr, expected: Option<Vec<[i64; 2]>>) {
        let expected = expected.map(|v| v.into_iter().map(|[n, d]|
            Rat::from_integers(n.into(), d.into())).collect::<Vec<_>>());
        let result = checked_integer_relation_product(&input, &a, &b);
        match expected {
            None => assert_eq!(result.ok(), None),
            Some(mut expected) => assert_eq!(result, Ok(expected))
        }
    }

    fn assert_integer_relation_fbig(basis: Vec<FBig>, expected: Option<Vec<i64>>) {
        let expected = expected.map(|v| v.into_iter().map(|i| Integer::from(i)).collect::<Vec<_>>());
        let mut result = maybe_integer_relation(&basis, &1000.into(), 100);

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

    #[test]
    fn test_relation_product_past_fail() {
        assert_integer_relation_product(
            vec![sqrt_expr!(1), sqrt_expr!(2), sqrt_expr!(2 + sqrt 2), sqrt_expr!(2 - sqrt 2)],
            sqrt_expr!(2),
            sqrt_expr!(2 - sqrt 2),
            Some(vec![[0, 1], [0, 1], [1, 1], [-1, 1]])
        );
    }

    //#[test]
    //#[ignore]
    //fn test_run_pslq() {
    //    let basis = [
    //        sqrt_expr!(1),
    //        sqrt_expr!(2),
    //        sqrt_expr!(2 + sqrt 2),
    //        sqrt_expr!(2 - sqrt 2),
    //        sqrt_expr!(2 + sqrt(2 + sqrt 2)),
    //        sqrt_expr!(2 + sqrt(2 - sqrt 2)),
    //        sqrt_expr!(2 - sqrt(2 + sqrt 2)),
    //        sqrt_expr!(2 - sqrt(2 - sqrt 2)),
    //    ];
    //    let mut basis = basis.into_iter().map(|e| e.to_fbig(800)).collect::<Vec<_>>();
    //    let numer = basis[0].clone() * 4;
    //    let denom = (0..8).map(|i| basis[i].clone() * i).reduce(|a, b| a + b).unwrap();
    //    basis.push(numer / denom);
    //    println!("[");
    //    for b in &basis {
    //        println!("    {}", b.to_decimal().value());
    //    }
    //    println!("]");
    //    let result = maybe_integer_relation(&basis, &100000000.into(), 10000);
    //    panic!("Aaaah! A scary result! {:?}", result);
    //}
}