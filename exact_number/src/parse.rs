use std::str::FromStr;

use malachite::base::num::basic::traits::One;
use malachite::{base::num::conversion::traits::FromStringBase, Integer, Natural};
use num::Num;
use peg::error::ParseError;
use peg::str::LineCol;
use peg::RuleResult;

use crate::basis::{BasisError, SqrtExpr};
use crate::BasedExpr;
use crate::rat::Rat;

peg::parser! {
    grammar expr_parser(radix: u8) for str {
        rule _() = quiet!{[' ' | '\t']*}
        rule sign(a: usize) -> bool = s:['+' | '-'] { s == '+' } / s:(#{|input, pos|
            if input[a..pos].chars().all(|c| [' ', '\t'].contains(&c)) {
                RuleResult::Matched(pos, ())
            } else {
                RuleResult::Failed
            }}) { true }
        rule digit() -> char = ['0'..='9' | 'A'..='Z' | 'a'..='z' | '_']
        rule natural(a: usize) -> Natural = n:$(digit()+) {? Natural::from_string_base(radix, &n.replace("_", "")).ok_or("natural") }
        rule integer(a: usize) -> Integer = s:sign(a) _ abs:natural(a) { Integer::from_sign_and_abs(s, abs) }
        rule rational(a: usize) -> Rat = n:integer(a) _ d:("/" _ d:natural(a) { d })? { Rat::from_integers(n, d.map(Into::into).unwrap_or(Integer::ONE)) }

        rule int_sum() -> Vec<(Integer, SqrtExpr)> = _ p:position!() v:(
            sign:sign(p) _ s:sqrt_expr() { (Integer::from_sign_and_abs(sign, Natural::ONE), s) } /
            i:integer(p) _ s:sqrt_expr()? { (i, s.unwrap_or(SqrtExpr::ONE)) }
        )* _ { v }

        pub rule sqrt_expr() -> SqrtExpr = _ ("sqrt" / "√") _
            s:(
                p:position!() i:integer(p) { SqrtExpr::Int(i) } /
                s:sqrt_expr() { SqrtExpr::Sum(vec![(Integer::ONE, s)]) } /
                "(" _ n:((p:position!() i:integer(p) _ ")" { SqrtExpr::Int(i) }) / (s:int_sum() _ ")" { SqrtExpr::Sum(s) })) { n }
            ) _ { s }

        pub rule based_expr() -> Vec<(Rat, SqrtExpr)> = _ p:position!() v:(
            sign:sign(p) _ s:sqrt_expr() { (Rat::from(Integer::from_sign_and_abs(sign, Natural::ONE)), s) } /
            r:rational(p) _  s:sqrt_expr()? { (r, s.unwrap_or(SqrtExpr::ONE)) }
        )* _ { v }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FromStrError {
    Parse(ParseError<LineCol>),
    Basis(BasisError)
}

impl FromStr for SqrtExpr {
    type Err = ParseError<LineCol>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        expr_parser::sqrt_expr(s, 10)
    }
}

impl FromStr for BasedExpr {
    type Err = FromStrError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_str_radix(s, 10)
    }
}

impl Num for BasedExpr {
    type FromStrRadixErr = FromStrError;

    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        match expr_parser::based_expr(str, radix as u8) {
            Ok(vec) => Self::try_from_terms(vec).map_err(|e| FromStrError::Basis(e)),
            Err(err) => Err(FromStrError::Parse(err))
        }
    }
}

#[cfg(test)]
mod test {
    use std::str::FromStr;

    use malachite::Integer;

    use crate::{basis::SqrtExpr, sqrt_expr, BasedExpr, based_expr};

    #[test]
    fn test_sqrt_expr_parse() {
        assert_eq!("√1".parse::<SqrtExpr>(), Ok(sqrt_expr!(1)));
        assert_eq!("√ 1".parse::<SqrtExpr>(), Ok(sqrt_expr!(1)));
        assert_eq!("  √1   ".parse::<SqrtExpr>(), Ok(sqrt_expr!(1)));
        assert_eq!("√(1)".parse::<SqrtExpr>(), Ok(sqrt_expr!(1)));
        assert_eq!("√(1) ".parse::<SqrtExpr>(), Ok(sqrt_expr!(1)));
        assert_eq!("√( 1)".parse::<SqrtExpr>(), Ok(sqrt_expr!(1)));
        assert_eq!("√ ( 1 )".parse::<SqrtExpr>(), Ok(sqrt_expr!(1)));
        assert_eq!("sqrt 2".parse::<SqrtExpr>(), Ok(sqrt_expr!(2)));
        assert_eq!("sqrt (2)".parse::<SqrtExpr>(), Ok(sqrt_expr!(2)));
        assert_eq!("sqrt(2)".parse::<SqrtExpr>(), Ok(sqrt_expr!(2)));
        assert_eq!("√√2".parse::<SqrtExpr>(), Ok(sqrt_expr!(sqrt 2)));
        assert_eq!("√(1 +  sqrt 2)".parse::<SqrtExpr>(), Ok(sqrt_expr!(1 + sqrt 2)));
        assert_eq!("√(1 + 1 sqrt 2)".parse::<SqrtExpr>(), Ok(sqrt_expr!(1 + sqrt 2)));
        assert_eq!("sqrt(-12 + 12 sqrt 2)".parse::<SqrtExpr>(), Ok(sqrt_expr!(-12 + 12 sqrt 2))); // needs parenthesis
        assert_eq!("sqrt(123_456_789_012_345_678_901_234_567_890_123_456_789_012_345_677)".parse::<SqrtExpr>(),
            Ok(SqrtExpr::Int("123456789012345678901234567890123456789012345677".parse().unwrap()))); // too big for an i128

        assert_eq!("√(3".parse::<SqrtExpr>().ok(), None); // Unbalanced
        assert_eq!("√3)".parse::<SqrtExpr>().ok(), None); // Unbalanced
        assert_eq!("sqrt(1/2)".parse::<SqrtExpr>().ok(), None); // No rationals allowed inside square roots
        assert_eq!("sqrt -12 + 12 sqrt 2".parse::<SqrtExpr>().ok(), None); // needs parenthesis
        assert_eq!("".parse::<SqrtExpr>().ok(), None); // must start with a square root
    }

    #[test]
    fn test_based_expr_parse() {
        assert_eq!("1".parse::<BasedExpr>(), Ok(based_expr!(1)));
        assert_eq!("+1".parse::<BasedExpr>(), Ok(based_expr!(1)));
        assert_eq!("-1".parse::<BasedExpr>(), Ok(based_expr!(-1)));
        assert_eq!(" 1".parse::<BasedExpr>(), Ok(based_expr!(1)));
        assert_eq!(" 1 ".parse::<BasedExpr>(), Ok(based_expr!(1)));
        assert_eq!("256".parse::<BasedExpr>(), Ok(based_expr!(256)));
        assert_eq!("1+1 sqrt 2".parse::<BasedExpr>(), Ok(based_expr!(1 + sqrt 2)));
        assert_eq!("1+ sqrt 2".parse::<BasedExpr>(), Ok(based_expr!(1 + sqrt 2)));
        assert_eq!("1 - 2 sqrt 2".parse::<BasedExpr>(), Ok(based_expr!(1 - 2 sqrt 2)));
        assert_eq!("1 - 2 sqrt 2 + 3 sqrt 5 + 40 sqrt 10".parse::<BasedExpr>(), Ok(based_expr!(1 - 2 sqrt 2 + 3 sqrt 5 + 40 sqrt 10)));
        assert_eq!("1/2 + 3/4√5".parse::<BasedExpr>(), Ok(based_expr!(1/2 + 3/4 sqrt 5)));
        assert_eq!("-1/2 - 3/4√5".parse::<BasedExpr>(), Ok(based_expr!(-1/2 - 3/4 sqrt 5)));

        assert_eq!("1 sqrt 2 + 1 sqrt 3".parse::<BasedExpr>().ok(), None); // not a basis
        assert_eq!("-/2".parse::<BasedExpr>().ok(), None); // Needs numerator
        assert_eq!("1/-2".parse::<BasedExpr>().ok(), None); // Denominator can't have sign
    }
}