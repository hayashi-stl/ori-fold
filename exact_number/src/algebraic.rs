// Algebraics. Used for checking solutions from the PSLQ algorithm for correctness.
// Unfortunately, now we're using 3 bigint libraries because they're all incompatible with each other
// and used by different necessary libraries
// * malachite::Integer (because malachite_q::Rational exists)
// * dashu_int::IBig (because dashu_float is used for arbitrary precision floats for pslq )
// * num_bigint::BigInt (because algebraics, used for checking pslq solutions, uses it, and there's no alternative)

//use crate::SqrtExpr;
//
//fn integer_to_algebraic(int: &Integer) -> FBig {
//    let words = int.unsigned_abs_ref().to_limbs_asc();
//    let ibig = IBig::from_parts(if int < &Integer::ZERO {Sign::Negative} else {Sign::Positive}, UBig::from_words(&words));
//    FBig::from(ibig).with_precision(precision).value()
//}
//
//impl SqrtExpr {
//    fn to_algebraic(&self) {
//        match self {
//            Self::Int(int) => 
//        }
//    }
//}