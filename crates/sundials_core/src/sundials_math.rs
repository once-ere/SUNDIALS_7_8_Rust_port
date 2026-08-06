//! Port of `src/sundials/sundials_math.c` + `include/sundials/sundials_math.h`
//! (double-precision branch).

use crate::sundials_types::*;

/// C macro `SUNMIN(A,B)`: `((A) < (B) ? (A) : (B))`.
pub fn SUNMIN<T: PartialOrd>(a: T, b: T) -> T {
    if a < b {
        a
    } else {
        b
    }
}

/// C macro `SUNMAX(A,B)`: `((A) > (B) ? (A) : (B))`.
pub fn SUNMAX<T: PartialOrd>(a: T, b: T) -> T {
    if a > b {
        a
    } else {
        b
    }
}

/// C macro `SUNSQR(A)`: `((A) * (A))`.
pub fn SUNSQR(a: sunrealtype) -> sunrealtype {
    a * a
}

/// C macro `SUNRsqrt(x)`: 0 for `x <= 0`, else `sqrt(x)`.
pub fn SUNRsqrt(x: sunrealtype) -> sunrealtype {
    if x <= 0.0 {
        0.0
    } else {
        x.sqrt()
    }
}

/// C macro `SUNRabs(x)`: `fabs(x)`.
pub fn SUNRabs(x: sunrealtype) -> sunrealtype {
    x.abs()
}

/// C macro `SUNRexp(x)`: `exp(x)`.
pub fn SUNRexp(x: sunrealtype) -> sunrealtype {
    x.exp()
}

/// C macro `SUNRceil(x)`: `ceil(x)`.
pub fn SUNRceil(x: sunrealtype) -> sunrealtype {
    x.ceil()
}

/// C macro `SUNRround(x)`: `round(x)` (halfway cases away from zero).
pub fn SUNRround(x: sunrealtype) -> sunrealtype {
    x.round()
}

/// C macro `SUNRcopysign(x, y)`: `copysign(x, y)`.
pub fn SUNRcopysign(x: sunrealtype, y: sunrealtype) -> sunrealtype {
    x.copysign(y)
}

/// C macro `SUNRsamesign(x, y)`: `signbit(x) == signbit(y)`.
pub fn SUNRsamesign(x: sunrealtype, y: sunrealtype) -> sunbooleantype {
    x.is_sign_negative() == y.is_sign_negative()
}

/// C macro `SUNRdifferentsign(x, y)`: `!SUNRsamesign(x, y)`.
pub fn SUNRdifferentsign(x: sunrealtype, y: sunrealtype) -> sunbooleantype {
    !SUNRsamesign(x, y)
}

/// C macro `SUNRpowerR(base, exponent)`: `pow(base, exponent)`.
pub fn SUNRpowerR(base: sunrealtype, exponent: sunrealtype) -> sunrealtype {
    base.powf(exponent)
}

pub fn SUNIpowerI(base: i32, exponent: i32) -> i32 {
    let mut prod: i32 = 1;
    let mut i = 1;
    while i <= exponent {
        prod *= base;
        i += 1;
    }
    prod
}

pub fn SUNRpowerI(base: sunrealtype, exponent: i32) -> sunrealtype {
    let mut prod: sunrealtype = 1.0;
    let expt = exponent.abs();
    let mut i = 1;
    while i <= expt {
        prod *= base;
        i += 1;
    }
    if exponent < 0 {
        prod = 1.0 / prod;
    }
    prod
}

pub fn SUNRCompare(a: sunrealtype, b: sunrealtype) -> sunbooleantype {
    SUNRCompareTol(a, b, 10.0 * SUN_UNIT_ROUNDOFF)
}

pub fn SUNRCompareTol(a: sunrealtype, b: sunrealtype, tol: sunrealtype) -> sunbooleantype {
    /* If a and b are exactly equal.
     * This also covers the case where a and b are both inf under IEEE 754. */
    if a == b {
        return SUNFALSE;
    }

    let diff = SUNRabs(a - b);
    let norm = SUNMIN(SUNRabs(a + b), SUN_BIG_REAL);

    /* C uses !isless(diff, max) so NaNs compare "not equal" (true);
     * Rust `!(diff < max)` has identical semantics for NaN. */
    !(diff < SUNMAX(10.0 * SUN_UNIT_ROUNDOFF, tol * norm))
}

/// C `SUNStrToReal`: `strtod` semantics — parse the longest valid leading
/// float, ignoring leading whitespace and trailing junk; 0.0 if nothing
/// parses.
pub fn SUNStrToReal(str_: &str) -> sunrealtype {
    let s = str_.trim_start();
    let b = s.as_bytes();
    let mut i = 0usize;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let lower = s[i..].to_ascii_lowercase();
    if lower.starts_with("infinity") {
        return s[..i + 8].parse::<f64>().unwrap_or(0.0);
    }
    if lower.starts_with("inf") {
        return s[..i + 3].parse::<f64>().unwrap_or(0.0);
    }
    if lower.starts_with("nan") {
        return f64::NAN;
    }
    let start_digits = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
    }
    // at least one digit must appear in the mantissa
    if i == start_digits || !s[start_digits..i].bytes().any(|c| c.is_ascii_digit()) {
        return 0.0;
    }
    let mantissa_end = i;
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        let mut j = i + 1;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            j += 1;
        }
        let exp_digits = j;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        if j > exp_digits {
            i = j;
        } else {
            i = mantissa_end;
        }
    }
    s[..i].parse::<f64>().unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn powers() {
        assert_eq!(SUNIpowerI(2, 10), 1024);
        assert_eq!(SUNIpowerI(3, 0), 1);
        assert_eq!(SUNRpowerI(2.0, -2), 0.25);
        assert_eq!(SUNRpowerI(10.0, 3), 1000.0);
    }

    #[test]
    fn compare() {
        assert!(!SUNRCompare(1.0, 1.0));
        assert!(SUNRCompare(1.0, 1.001));
        assert!(SUNRCompare(f64::NAN, 1.0));
        assert!(!SUNRCompare(f64::INFINITY, f64::INFINITY));
    }

    #[test]
    fn str_to_real() {
        assert_eq!(SUNStrToReal("1e-3"), 1e-3);
        assert_eq!(SUNStrToReal("  -2.5rest"), -2.5);
        assert_eq!(SUNStrToReal("junk"), 0.0);
        assert_eq!(SUNStrToReal(".5"), 0.5);
        assert_eq!(SUNStrToReal("3."), 3.0);
        assert_eq!(SUNStrToReal("1e"), 1.0);
        assert_eq!(SUNStrToReal("inf"), f64::INFINITY);
        assert_eq!(SUNStrToReal("-Infinity"), f64::NEG_INFINITY);
    }

    #[test]
    fn sqrt_guard() {
        assert_eq!(SUNRsqrt(-4.0), 0.0);
        assert_eq!(SUNRsqrt(4.0), 2.0);
    }
}
