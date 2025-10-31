//! Contains the geometric tests

use paste::paste;

use crate::Expansion;
use crate::{square, two_product, two_sum};

use nalgebra::{Vector2, Vector3, vector};
type Vec2 = Vector2<f64>;
type Vec3 = Vector3<f64>;

const EPSILON: f64 = f64::EPSILON / 2.0;

const ORIENT_2D_BOUND_A: f64 = (3.0 + 16.0 * EPSILON) * EPSILON;
const ORIENT_2D_BOUND_B: f64 = (2.0 + 12.0 * EPSILON) * EPSILON;
const ORIENT_2D_BOUND_C1: f64 = (3.0 + 8.0 * EPSILON) * EPSILON;
const ORIENT_2D_BOUND_C2: f64 = (9.0 + 64.0 * EPSILON) * EPSILON * EPSILON;

const ORIENT_3D_BOUND_A: f64 = (7.0 + 56.0 * EPSILON) * EPSILON;
const ORIENT_3D_BOUND_B: f64 = (3.0 + 28.0 * EPSILON) * EPSILON;
const ORIENT_3D_BOUND_C1: f64 = (3.0 + 8.0 * EPSILON) * EPSILON;
const ORIENT_3D_BOUND_C2: f64 = (26.0 + 288.0 * EPSILON) * EPSILON * EPSILON;

const IN_CIRCLE_BOUND_A: f64 = (10.0 + 96.0 * EPSILON) * EPSILON;
const IN_CIRCLE_BOUND_B: f64 = (4.0 + 48.0 * EPSILON) * EPSILON;
const IN_CIRCLE_BOUND_C1: f64 = (3.0 + 8.0 * EPSILON) * EPSILON;
const IN_CIRCLE_BOUND_C2: f64 = (44.0 + 576.0 * EPSILON) * EPSILON * EPSILON;

const IN_SPHERE_BOUND_A: f64 = (16.0 + 224.0 * EPSILON) * EPSILON;
const IN_SPHERE_BOUND_B: f64 = (5.0 + 72.0 * EPSILON) * EPSILON;
const IN_SPHERE_BOUND_C1: f64 = (3.0 + 8.0 * EPSILON) * EPSILON;
const IN_SPHERE_BOUND_C2: f64 = (71.0 + 1408.0 * EPSILON) * EPSILON * EPSILON;

const MAGNITUDE_CMP_2D_BOUND_A: f64 = (2.0 + 8.0 * EPSILON) * EPSILON;

const MAGNITUDE_CMP_3D_BOUND_A: f64 = (3.0 + 12.0 * EPSILON) * EPSILON;

const DISTANCE_CMP_3D_BOUND_A: f64 = (6.0 + 32.0 * EPSILON) * EPSILON;

const SIGN_DET_X_X2Y2_BOUND_A: f64 = (5.0 + 32.0 * EPSILON) * EPSILON;

const SIGN_DET_X_X2Y2Z2_BOUND_A: f64 = (6.0 + 32.0 * EPSILON) * EPSILON;

const SIGN_DET_X_Y_X2Y2Z2_BOUND_A: f64 = (9.0 + 64.0 * EPSILON) * EPSILON;

const SIGN_DET_X_Y_X2Y2Z2_PLUS_2X_DET_X_Y_Z_BOUND_A: f64 = (10.0 + 80.0 * EPSILON) * EPSILON;
const SIGN_DET_X_X2Y2Z2_PLUS_2X_DET_X_Y_BOUND_A: f64 = (7.0 + 48.0 * EPSILON) * EPSILON;
const SIGN_DET_X2Y2Z2_PLUS_2X_DET_X_BOUND_A: f64 = (4.0 + 24.0 * EPSILON) * EPSILON;

const SIGN_DET_X_X2Y2_PLUS_2X_DET_X_Y_BOUND_A: f64 = (6.0 + 40.0 * EPSILON) * EPSILON;
const SIGN_DET_X2Y2_PLUS_2X_DET_X_BOUND_A: f64 = (3.0 + 16.0 * EPSILON) * EPSILON;

const COMPLEX_INTERSECT_COMPARE_BOUND: f64 = (13.0 + 128.0 * EPSILON) * EPSILON;

/// Calculates the orientation of points `a`, `b`, `c` in a plane.
/// Returns twice the area of the triangle formed by `a`, `b`, and `c`,
/// which is a positive number if they define a left turn,
/// a negative number if they define a right turn,
/// and 0 if they are collinear.
pub fn orient_2d(a: Vec2, b: Vec2, c: Vec2) -> f64 {
    let diag1 = (a.x - c.x) * (b.y - c.y);
    let diag2 = (a.y - c.y) * (b.x - c.x);
    let det = diag1 - diag2;
    let det_sum = diag1.abs() + diag2.abs();

    if det.abs() >= det_sum * ORIENT_2D_BOUND_A {
        det
    } else {
        orient_2d_adapt(a, b, c, det_sum)
    }
}

fn orient_2d_adapt(a: Vec2, b: Vec2, c: Vec2, det_sum: f64) -> f64 {
    let diag1 = two_product(a.x - c.x, b.y - c.y);
    let diag2 = two_product(a.y - c.y, b.x - c.x);
    let det_hi = diag1 - diag2;
    let det_approx = det_hi.approximate();

    if det_approx.abs() >= det_sum * ORIENT_2D_BOUND_B {
        return det_approx;
    }

    // Calculate correction factor, but don't go for the full exact value yet
    let ax = two_sum(a.x, -c.x);
    let ay = two_sum(a.y, -c.y);
    let bx = two_sum(b.x, -c.x);
    let by = two_sum(b.y, -c.y);

    let det_mid1 = ax[0] * by[1] - ay[0] * bx[1];
    let det_mid2 = ax[1] * by[0] - ay[1] * bx[0];
    let det = det_mid1 + det_mid2 + det_approx;

    if det.abs() >= ORIENT_2D_BOUND_C1 * det_approx.abs() + ORIENT_2D_BOUND_C2 * det_sum {
        return det;
    }

    // Determinant is likely 0; go for exact value
    let det_hi = det_hi.dynamic();
    let det_m1 = (two_product(ax[0], by[1]) - two_product(ay[0], bx[1])).dynamic();
    let det_m2 = (two_product(ax[1], by[0]) - two_product(ay[1], bx[0])).dynamic();
    let det_lo = (two_product(ax[0], by[0]) - two_product(ay[0], bx[0])).dynamic();
    let det = (det_hi + det_m1) + (det_m2 + det_lo);
    det.approximate()
}

/// Computes the following determinant
/// ```notrust
/// | b.x - a.x   b.y - a.y |
/// | d.x - c.x   d.y - c.y |
/// ```
pub fn cross_2d(a: Vec2, b: Vec2, c: Vec2, d: Vec2) -> f64 {
    let diag1 = (b.x - a.x) * (d.y - c.y);
    let diag2 = (b.y - a.y) * (d.x - c.x);
    let det = diag1 - diag2;
    let det_sum = diag1.abs() + diag2.abs();

    if det.abs() >= det_sum * ORIENT_2D_BOUND_A {
        det
    } else {
        cross_2d_adapt(a, b, c, d, det_sum)
    }
}

fn cross_2d_adapt(a: Vec2, b: Vec2, c: Vec2, d: Vec2, det_sum: f64) -> f64 {
    let diag1 = two_product(b.x - a.x, d.y - c.y);
    let diag2 = two_product(b.y - a.y, d.x - c.x);
    let det_hi = diag1 - diag2;
    let det_approx = det_hi.approximate();

    if det_approx.abs() >= det_sum * ORIENT_2D_BOUND_B {
        return det_approx;
    }

    // Calculate correction factor, but don't go for the full exact value yet
    let abx = two_sum(b.x, -a.x);
    let aby = two_sum(b.y, -a.y);
    let cdx = two_sum(d.x, -c.x);
    let cdy = two_sum(d.y, -c.y);

    let det_mid1 = abx[0] * cdy[1] - aby[0] * cdx[1];
    let det_mid2 = abx[1] * cdy[0] - aby[1] * cdx[0];
    let det = det_mid1 + det_mid2 + det_approx;

    if det.abs() >= ORIENT_2D_BOUND_C1 * det_approx.abs() + ORIENT_2D_BOUND_C2 * det_sum {
        return det;
    }

    // Determinant is likely 0; go for exact value
    let det_hi = det_hi.dynamic();
    let det_m1 = (two_product(abx[0], cdy[1]) - two_product(aby[0], cdx[1])).dynamic();
    let det_m2 = (two_product(abx[1], cdy[0]) - two_product(aby[1], cdx[0])).dynamic();
    let det_lo = (two_product(abx[0], cdy[0]) - two_product(aby[0], cdx[0])).dynamic();
    let det = (det_hi + det_m1) + (det_m2 + det_lo);
    det.approximate()
}

macro_rules! sep_let {
    (($($d:tt $var:ident),*) => let $($name:ident : [$($vals:tt)*]),* $(,)? = $($expr:tt)*) => {
        macro_rules! __mac {
            ($($d $var:ident),*) => { $($expr)* };
            ($($d $var:expr),*) => { $($expr)* };
        }

        $(let $name = __mac!($($vals)*);)*
    };
}

// xy cyclic permutation shortcut.
// Simplifies cofactor repetition.
macro_rules! sep_xy {
    (($dx:tt $x:ident, $dy:tt $y:ident, $(, $d:tt $var:ident)*)
        => let $n1:ident : [$($v1:tt)*], $n2:ident : [$($v2:tt)*] = $($expr:tt)*) =>
    {
        sep_let!(($dx $x, $dy $y, $(, $d $var)*) =>
            let $n1: [$x, $y, $($v1)*], $n2: [$y, $x, $($v2)*] = $($expr)*);
    };

    (($dx:tt $x:ident, $dy:tt $y:ident) => let $n1:ident, $n2:ident = $($expr:tt)*) => {
        sep_let!(($dx $x, $dy $y) =>
            let $n1: [$x, $y], $n2: [$y, $x] = $($expr)*);
    };
}

// xyz cyclic permutation shortcut.
// Simplifies cofactor repetition.
macro_rules! sep_xyz {
    (($dx:tt $x:ident, $dy:tt $y:ident, $dz:tt $z:ident $(, $d:tt $var:ident)*)
        => let $n1:ident : [$($v1:tt)*], $n2:ident : [$($v2:tt)*], $n3:ident : [$($v3:tt)*] = $($expr:tt)*) =>
    {
        sep_let!(($dx $x, $dy $y, $dz $z $(, $d $var)*) =>
            let $n1: [$x, $y, $z, $($v1)*], $n2: [$y, $z, $x, $($v2)*], $n3: [$z, $x, $y, $($v3)*] = $($expr)*);
    };

    (($dx:tt $x:ident, $dy:tt $y:ident, $dz:tt $z:ident) => let $n1:ident, $n2:ident, $n3:ident = $($expr:tt)*) => {
        sep_let!(($dx $x, $dy $y, $dz $z) =>
            let $n1: [$x, $y, $z], $n2: [$y, $z, $x], $n3: [$z, $x, $y] = $($expr)*);
    };
}

// xyzw cyclic permutation shortcut.
// Simplifies cofactor repetition.
macro_rules! sep_xyzw {
    (($dx:tt $x:ident, $dy:tt $y:ident, $dz:tt $z:ident, $dw:tt $w:ident $(, $d:tt $var:ident)*)
        => let $n1:ident : [$($v1:tt)*], $n2:ident : [$($v2:tt)*], $n3:ident : [$($v3:tt)*], $n4:ident : [$($v4:tt)*] = $($expr:tt)*) =>
    {
        sep_let!(($dx $x, $dy $y, $dz $z, $dw $w $(, $d $var)*) =>
            let $n1: [$x, $y, $z, $w, $($v1)*], $n2: [$y, $z, $w, $x, $($v2)*], $n3: [$z, $w, $x, $y, $($v3)*], $n4: [$w, $x, $y, $z, $($v4)*]
                = $($expr)*);
    };

    (($dx:tt $x:ident, $dy:tt $y:ident, $dz:tt $z:ident, $dw:tt $w:ident) => let $n1:ident, $n2:ident, $n3:ident, $n4:ident = $($expr:tt)*) => {
        sep_let!(($dx $x, $dy $y, $dz $z, $dw $w) =>
            let $n1: [$x, $y, $z, $w], $n2: [$y, $z, $w, $x], $n3: [$z, $w, $x, $y], $n4: [$w, $x, $y, $z] = $($expr)*);
    };
}

// xyzw half even permutation shortcut.
// Simplifies cofactor repitition.
macro_rules! sep_xyzw6 {
//    (($dx:tt $x:ident, $dy:tt $y:ident, $dz:tt $z:ident, $dw:tt $w:ident $(, $d:tt $var:ident)*)
//        => let $nxy:ident : [$($vxy:tt)*], $nxz:ident : [$($vxz:tt)*], $nxw:ident : [$($vxw:tt)*]
//               $nyz:ident : [$($vyz:tt)*], $nyw:ident : [$($vyw:tt)*], $nyx:ident : [$($vyx:tt)*]
//               $nzw:ident : [$($vzw:tt)*], $nzx:ident : [$($vzx:tt)*], $nzy:ident : [$($vzy:tt)*]
//               $nwx:ident : [$($vwz:tt)*], $nwy:ident : [$($vwy:tt)*], $nwz:ident : [$($vwz:tt)*]
//        = $($expr:tt)*) =>
//    {
//        sep_let!(($dx $x, $dy $y, $dz $z $(, $d $var)*) =>
//            let $nxy: [$x, $y, $z, $w, $($vxy)*], $nxz: [$x, $z, $w, $y, $($vxz)*], $nxw: [ $($v3)*] = $($expr)*);
//    };

    (($dx:tt $x:ident, $dy:tt $y:ident, $dz:tt $z:ident, $dw:tt $w:ident)
        => let $nxy:ident, $nxz:ident, $nxw:ident, $nyz:ident, $nyw:ident, $nzw:ident = $($expr:tt)*) => {
        sep_let!(($dx $x, $dy $y, $dz $z, $dw $w) =>
            let $nxy: [$x, $y, $z, $w], $nxz: [$x, $z, $w, $y], $nxw: [$x, $w, $y, $z],
                $nyz: [$y, $z, $x, $w], $nyw: [$y, $w, $z, $x], $nzw: [$z, $w, $x, $y],
             = $($expr)*);
    };

    //// Because writing out all 12 terms is annoying
    //(($dx:tt $x:ident, $dy:tt $y:ident, $dz:tt $z:ident, $dw:tt $w:ident)
    //    => let $n:ident = $($expr:tt)*) => {
    //    sep_let!(($dx $x, $dy $y, $dz $z) =>
    //        let paste!([<$n$x$y>]): [$x, $y, $z, $w], paste!([<$n$x$z>]): [$x, $z, $w, $y], paste!([<$n$x$w>]): [$x, $w, $y, $z],
    //            paste!([<$n$y$z>]): [$y, $z, $w, $x], paste!([<$n$y$w>]): [$y, $w, $x, $z], paste!([<$n$y$x>]): [$y, $x, $z, $w],
    //            paste!([<$n$z$w>]): [$z, $w, $x, $y], paste!([<$n$z$x>]): [$z, $x, $y, $w], paste!([<$n$z$y>]): [$z, $y, $w, $x],
    //            paste!([<$n$w$x>]): [$w, $x, $y, $z], paste!([<$n$w$y>]): [$w, $y, $z, $x], paste!([<$n$w$z>]): [$w, $z, $x, $y],
    //         = $($expr)*);
    //};
}

// 5 cyclic permutation shortcut.
// Simplifies cofactor repetition.
macro_rules! sep_5 {
    (($dx:tt $x:ident, $dy:tt $y:ident, $dz:tt $z:ident, $dw:tt $w:ident, $dv:tt $u:ident $(, $d:tt $var:ident)*)
        => let $n1:ident : [$($v1:tt)*],
               $n2:ident : [$($v2:tt)*],
               $n3:ident : [$($v3:tt)*],
               $n4:ident : [$($v4:tt)*]
               $n5:ident : [$($v5:tt)*] = $($expr:tt)*) =>
    {
        sep_let!(($dx $x, $dy $y, $dz $z, $dw $w, $du $u $(, $d $var)*) =>
            let $n1: [$x, $y, $z, $w, $u, $($v1)*],
                $n2: [$y, $z, $w, $u, $x, $($v2)*],
                $n3: [$z, $w, $u, $x, $y, $($v3)*],
                $n4: [$w, $u, $x, $y, $z, $($v4)*],
                $n5: [$u, $x, $y, $z, $w, $($v5)*] = $($expr)*);
    };

    (($dx:tt $x:ident, $dy:tt $y:ident, $dz:tt $z:ident, $dw:tt $w:ident, $du:tt $u:ident)
        => let $n1:ident, $n2:ident, $n3:ident, $n4:ident, $n5:ident = $($expr:tt)*) => {
        sep_let!(($dx $x, $dy $y, $dz $z, $dw $w, $du $u) =>
            let $n1: [$x, $y, $z, $w, $u],
                $n2: [$y, $z, $w, $u, $x],
                $n3: [$z, $w, $u, $x, $y],
                $n4: [$w, $u, $x, $y, $z],
                $n5: [$u, $x, $y, $z, $w] = $($expr)*);
    };
}

// These regress performance, unfortunately
//macro_rules! arr_let {
//    (($($d:tt $var:ident),*) => let $name:ident : $([$($vals:tt)*]),* = $($expr:tt)*) => {
//        macro_rules! __mac {
//            ($($d $var:ident),*) => { $($expr)* };
//            ($($d $var:expr),*) => { $($expr)* };
//        }
//
//        let $name = [
//            $(__mac!($($vals)*)),*
//        ];
//    };
//}
//
//macro_rules! arr_xyz {
//    (($dx:tt $x:ident, $dy:tt $y:ident, $dz:tt $z:ident) => let $name:ident = $($expr:tt)*) => {
//        arr_let!(($dx $x, $dy $y, $dz $z) =>
//            let $name: [[0, 1, 2], [1, 2, 0], [2, 0, 1]] = $($expr)*);
//    };
//}

/// Calculates the orientation of points `a`, `b`, `c`, `d` in a space.
/// Returns 6 times the volume of the tetrahedron formed by `a`, `b`, `c`, `d`,
/// which is a positive number if `b`→`c`→`d` defines a left turn when looked at from `a`,
/// a negative number if they define a right turn,
/// and 0 if `a`, `b`, `c`, `d` are coplanar.
pub fn orient_3d(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> f64 {
    // vec![[0, 1, 2], [1, 2, 0], [2, 0, 1]].into_iter().map(|[x, y, z] stuff).sum::<f64>() regressed performance a lot
    sep_xyz!(($x, $y, $z) => let cof1, cof2, cof3 = 
        (a.$z - d.$z) * ((b.$x - d.$x) * (c.$y - d.$y) - (b.$y - d.$y) * (c.$x - d.$x)));
    let det = cof1 + cof2 + cof3;

    sep_xyz!(($x, $y, $z) => let cof1_sum, cof2_sum, cof3_sum = 
        (a.$z - d.$z).abs() * (((b.$x - d.$x) * (c.$y - d.$y)).abs() + ((b.$y - d.$y) * (c.$x - d.$x)).abs()));
    let det_sum = cof1_sum + cof2_sum + cof3_sum;

    if det.abs() >= det_sum * ORIENT_3D_BOUND_A {
        det
    } else {
        orient_3d_adapt(a, b, c, d, det_sum)
    }
}

fn orient_3d_adapt(a: Vec3, b: Vec3, c: Vec3, d: Vec3, det_sum: f64) -> f64 {
    sep_xyz!(($x, $y, $z) => let sub1, sub2, sub3 = (two_product(b.$x - d.$x, c.$y - d.$y) - two_product(b.$y - d.$y, c.$x - d.$x)).dynamic());
    //sep_xyz!(($x, $y, $z, $s) => let cof1: [sub1], cof2: [sub2], cof3: [sub3] = $s.scale_expansion(a.$z - d.$z));
    // Expanding macro until https://github.com/rust-lang/rust-analyzer/issues/20936 gets fixed
    let cof1 = sub1.scale_expansion(a.z - d.z);
    let cof2 = sub2.scale_expansion(a.x - d.x);
    let cof3 = sub3.scale_expansion(a.y - d.y);
    let det_hi = cof1 + cof2 + cof3;
    let det_approx = det_hi.approximate();

    if det_approx.abs() >= det_sum * ORIENT_3D_BOUND_B {
        return det_approx;
    }

    // Correction factor for order ε² error bound
    sep_xyz!(($x, $y, $z) => let ax, ay, az = two_sum(a.$x, -d.$x));
    sep_xyz!(($x, $y, $z) => let bx, by, bz = two_sum(b.$x, -d.$x));
    sep_xyz!(($x, $y, $z) => let cx, cy, cz = two_sum(c.$x, -d.$x));
    sep_xyz!(($x, $y, $z) => let cof1_m1, cof2_m1, cof3_m1 =
        paste!([<a$z>][0] * ([<b$x>][1] * [<c$y>][1] - [<b$y>][1] * [<c$x>][1])));
    sep_xyz!(($x, $y, $z) => let cof1_m2, cof2_m2, cof3_m2 =
        paste!([<a$z>][1] * (([<b$x>][0] * [<c$y>][1] + [<b$x>][1] * [<c$y>][0]) - ([<b$y>][0] * [<c$x>][1] + [<b$y>][1] * [<c$x>][0]))));
    let det = cof1_m1 + cof2_m1 + cof3_m1 + cof1_m2 + cof2_m2 + cof3_m2 + det_approx;

    if det.abs() >= ORIENT_3D_BOUND_C1 * det_approx.abs() + ORIENT_3D_BOUND_C2 * det_sum {
        return det;
    }

    // Exact result time!
    let det_m1 =
        sub1.scale_expansion(az[0]) + sub2.scale_expansion(ax[0]) + sub3.scale_expansion(ay[0]);
    sep_xyz!(($x, $y, $z) => let cof1x, cof1y, cof1z =
        paste!((two_product([<b$x>][0], [<c$y>][1]) + two_product([<b$x>][1], [<c$y>][0]) + two_product([<b$x>][0], [<c$y>][0])).dynamic()));
    sep_xyz!(($x, $y, $z) => let cof2x, cof2y, cof2z =
        paste!((two_product([<b$y>][0], [<c$x>][1]) + two_product([<b$y>][1], [<c$x>][0]) + two_product([<b$y>][0], [<c$x>][0])).dynamic()));
    sep_xyz!(($x, $y, $z) => let cof1, cof2, cof3 = paste!([<cof1$x>] - [<cof2$x>]));
    let det_m2 =
        cof1.scale_expansion(az[1]) + cof2.scale_expansion(ax[1]) + cof3.scale_expansion(ay[1]);
    let det_lo =
        cof1.scale_expansion(az[0]) + cof2.scale_expansion(ax[0]) + cof3.scale_expansion(ay[0]);
    let det = (det_hi + det_m1) + (det_m2 + det_lo);
    det.approximate()
}

/// Computes the following determinant
/// ```notrust
/// | b.x - a.x   b.y - a.y   b.z - a.z |
/// | d.x - c.x   d.y - c.y   d.z - c.z |
/// | f.x - e.x   f.y - e.y   f.z - e.z |
/// ```
pub fn cross_3d(a: Vec3, b: Vec3, c: Vec3, d: Vec3, e: Vec3, f: Vec3) -> f64 {
    // vec![[0, 1, 2], [1, 2, 0], [2, 0, 1]].into_iter().map(|[x, y, z] stuff).sum::<f64>() regressed performance a lot
    sep_xyz!(($x, $y, $z) => let cof1, cof2, cof3 = 
        (b.$z - a.$z) * ((d.$x - c.$x) * (f.$y - e.$y) - (d.$y - c.$y) * (f.$x - e.$x)));
    let det = cof1 + cof2 + cof3;

    sep_xyz!(($x, $y, $z) => let cof1_sum, cof2_sum, cof3_sum = 
        (b.$z - a.$z).abs() * (((d.$x - c.$x) * (f.$y - e.$y)).abs() + ((d.$y - c.$y) * (f.$x - e.$x)).abs()));
    let det_sum = cof1_sum + cof2_sum + cof3_sum;

    if det.abs() >= det_sum * ORIENT_3D_BOUND_A {
        det
    } else {
        cross_3d_adapt(a, b, c, d, e, f, det_sum)
    }
}

fn cross_3d_adapt(a: Vec3, b: Vec3, c: Vec3, d: Vec3, e: Vec3, f: Vec3, det_sum: f64) -> f64 {
    sep_xyz!(($x, $y, $z) => let sub1, sub2, sub3 = (two_product(d.$x - c.$x, f.$y - e.$y) - two_product(d.$y - c.$y, f.$x - e.$x)).dynamic());
    //sep_xyz!(($x, $y, $z, $s) => let cof1: [sub1], cof2: [sub2], cof3: [sub3] = $s.scale_expansion(b.$z - a.$z));
    // Expanding macro until https://github.com/rust-lang/rust-analyzer/issues/20936 gets fixed
    let cof1 = sub1.scale_expansion(b.z - a.z);
    let cof2 = sub2.scale_expansion(b.x - a.x);
    let cof3 = sub3.scale_expansion(b.y - a.y);
    let det_hi = cof1 + cof2 + cof3;
    let det_approx = det_hi.approximate();

    if det_approx.abs() >= det_sum * ORIENT_3D_BOUND_B {
        return det_approx;
    }

    // Correction factor for order ε² error bound
    sep_xyz!(($x, $y, $z) => let ax, ay, az = two_sum(b.$x, -a.$x));
    sep_xyz!(($x, $y, $z) => let bx, by, bz = two_sum(d.$x, -c.$x));
    sep_xyz!(($x, $y, $z) => let cx, cy, cz = two_sum(f.$x, -e.$x));
    sep_xyz!(($x, $y, $z) => let cof1_m1, cof2_m1, cof3_m1 =
        paste!([<a$z>][0] * ([<b$x>][1] * [<c$y>][1] - [<b$y>][1] * [<c$x>][1])));
    sep_xyz!(($x, $y, $z) => let cof1_m2, cof2_m2, cof3_m2 =
        paste!([<a$z>][1] * (([<b$x>][0] * [<c$y>][1] + [<b$x>][1] * [<c$y>][0]) - ([<b$y>][0] * [<c$x>][1] + [<b$y>][1] * [<c$x>][0]))));
    let det = cof1_m1 + cof2_m1 + cof3_m1 + cof1_m2 + cof2_m2 + cof3_m2 + det_approx;

    if det.abs() >= ORIENT_3D_BOUND_C1 * det_approx.abs() + ORIENT_3D_BOUND_C2 * det_sum {
        return det;
    }

    // Exact result time!
    let det_m1 =
        sub1.scale_expansion(az[0]) + sub2.scale_expansion(ax[0]) + sub3.scale_expansion(ay[0]);
    sep_xyz!(($x, $y, $z) => let cof1x, cof1y, cof1z =
        paste!((two_product([<b$x>][0], [<c$y>][1]) + two_product([<b$x>][1], [<c$y>][0]) + two_product([<b$x>][0], [<c$y>][0])).dynamic()));
    sep_xyz!(($x, $y, $z) => let cof2x, cof2y, cof2z =
        paste!((two_product([<b$y>][0], [<c$x>][1]) + two_product([<b$y>][1], [<c$x>][0]) + two_product([<b$y>][0], [<c$x>][0])).dynamic()));
    sep_xyz!(($x, $y, $z) => let cof1, cof2, cof3 = paste!([<cof1$x>] - [<cof2$x>]));
    let det_m2 =
        cof1.scale_expansion(az[1]) + cof2.scale_expansion(ax[1]) + cof3.scale_expansion(ay[1]);
    let det_lo =
        cof1.scale_expansion(az[0]) + cof2.scale_expansion(ax[0]) + cof3.scale_expansion(ay[0]);
    let det = (det_hi + det_m1) + (det_m2 + det_lo);
    det.approximate()
}

/// Returns a positive number if `d` is inside the oriented circle that goes through `a`, `b`, `c`,
/// a negative number if it lies outside,
/// and 0 if `a`, `b`, `c`, `d` are cocircular or `a`, `b`, `c` are collinear.
/// If `a`, `b`, `c` are in counterclockwise order, "inside the circle" is the inside.
/// If `a`, `b`, `c` are in clockwise order, "inside the circle" is the outside.
pub fn in_circle(a: Vec2, b: Vec2, c: Vec2, d: Vec2) -> f64 {
    sep_xyz!(($a, $b, $c) => let sqa, sqb, sqc = ($a.x - d.x) * ($a.x - d.x) + ($a.y - d.y) * ($a.y - d.y));
    sep_xyz!(($a, $b, $c) => let cof1, cof2, cof3 =
        paste!([<sq$a>]) * (($b.x - d.x) * ($c.y - d.y) - ($b.y - d.y) * ($c.x - d.x)));
    let det = cof1 + cof2 + cof3;

    sep_xyz!(($a, $b, $c) => let cof1_sum, cof2_sum, cof3_sum =
        paste!([<sq$a>]) * ((($b.x - d.x) * ($c.y - d.y)).abs() + (($b.y - d.y) * ($c.x - d.x)).abs()));
    let det_sum = cof1_sum + cof2_sum + cof3_sum;

    if det.abs() >= det_sum * IN_CIRCLE_BOUND_A {
        det
    } else {
        in_circle_adapt(a, b, c, d, det_sum)
    }
}

fn in_circle_adapt(a: Vec2, b: Vec2, c: Vec2, d: Vec2, det_sum: f64) -> f64 {
    sep_xyz!(($a, $b, $c) => let sqa, sqb, sqc = (square($a.x - d.x) + square($a.y - d.y)).dynamic());
    sep_xyz!(($a, $b, $c) => let suba, subb, subc = (two_product($b.x - d.x, $c.y - d.y) - two_product($b.y - d.y, $c.x - d.x)).dynamic());
    sep_xyz!(($a, $b, $c) => let cof1, cof2, cof3 = paste!([<sq$a>] * [<sub$a>]));
    let det = cof1 + cof2 + cof3;
    let det_approx = det.approximate();

    if det_approx.abs() >= det_sum * IN_CIRCLE_BOUND_B {
        return det_approx;
    }

    // Correction factor for order ε² error bound
    sep_xyz!(($a, $b, $c) => let xa, xb, xc = two_sum($a.x, -d.x));
    sep_xyz!(($a, $b, $c) => let ya, yb, yc = two_sum($a.y, -d.y));
    sep_xyz!(($a, $b, $c) => let cof1_m1, cof2_m1, cof3_m1 =
        paste!(([<x$a>][0] * [<x$a>][1] + [<y$a>][0] * [<y$a>][1])
            * ([<x$b>][1] * [<y$c>][1] - [<y$b>][1] * [<x$c>][1])));
    sep_xyz!(($a, $b, $c) => let cof1_m2, cof2_m2, cof3_m2 =
        paste!(([<x$a>][1] * [<x$a>][1] + [<y$a>][1] * [<y$a>][1])
            * (([<x$b>][0] * [<y$c>][1] + [<x$b>][1] * [<y$c>][0]) - ([<y$b>][0] * [<x$c>][1] + [<y$b>][1] * [<x$c>][0]))));
    let det = (cof1_m1 + cof2_m1 + cof3_m1) * 2.0 + cof1_m2 + cof2_m2 + cof3_m2 + det_approx;

    if det.abs() >= IN_CIRCLE_BOUND_C1 * det_approx.abs() + IN_CIRCLE_BOUND_C2 * det_sum {
        return det;
    }

    // Exact result time!
    // Let's be lazy
    sep_xyz!(($a, $b, $c) => let xa, xb, xc = paste!([<x$a>].dynamic()));
    sep_xyz!(($a, $b, $c) => let ya, yb, yc = paste!([<y$a>].dynamic()));
    sep_xyz!(($a, $b, $c) => let sqa, sqb, sqc = paste!(&[<x$a>] * &[<x$a>] + &[<y$a>] * &[<y$a>]));
    sep_xyz!(($a, $b, $c) => let cof1, cof2, cof3 = paste!(&[<sq$a>] * (&[<x$b>] * &[<y$c>] - &[<y$b>] * &[<x$c>])));
    let det = cof1 + cof2 + cof3;
    det.highest_magnitude()
}

/// Returns a positive number if `e` is inside the oriented sphere that goes through `a`, `b`, `c`, `d`,
/// a negative number if it lies outside,
/// and 0 if `a`, `b`, `c`, `d`, `e` are cospherical or `a`, `b`, `c`, `d` are coplanar.
/// If `a`, `b`, `c`, `d` are oriented positive, "inside the sphere" is the inside.
/// If `a`, `b`, `c`, `d` are oriented negative, "inside the sphere" is the outside.
pub fn in_sphere(a: Vec3, b: Vec3, c: Vec3, d: Vec3, e: Vec3) -> f64 {
    // Oh boy. 4x4 determinant, *lots* of cofactors.
    sep_xyzw!(($a, $b, $c, $d) => let sqa, sqb, sqc, sqd =
        ($a.x - e.x) * ($a.x - e.x) + ($a.y - e.y) * ($a.y - e.y) + ($a.z - e.z) * ($a.z - e.z));
    sep_xyzw6!(($a, $b, $c, $d) => let cof1, cof2, cof3, cof4, cof5, cof6 =
        paste!((($a.x - e.x) * ($b.y - e.y) - ($a.y - e.y) * ($b.x - e.x)) *
               (($c.z - e.z) * [<sq$d>] - [<sq$c>] * ($d.z - e.z))));
    let det = (cof1 + cof2 + cof3) + (cof4 + cof5 + cof6);

    sep_xyzw6!(($a, $b, $c, $d) => let cof1_sum, cof2_sum, cof3_sum, cof4_sum, cof5_sum, cof6_sum =
        paste!(((($a.x - e.x) * ($b.y - e.y)).abs() + (($a.y - e.y) * ($b.x - e.x)).abs()) *
               ((($c.z - e.z) * [<sq$d>]).abs() + ([<sq$c>] * ($d.z - e.z)).abs())));
    let det_sum = (cof1_sum + cof2_sum + cof3_sum) + (cof4_sum + cof5_sum + cof6_sum);

    if det.abs() >= det_sum * IN_SPHERE_BOUND_A {
        det
    } else {
        in_sphere_adapt(a, b, c, d, e, det_sum)
    }
}

fn in_sphere_adapt(a: Vec3, b: Vec3, c: Vec3, d: Vec3, e: Vec3, det_sum: f64) -> f64 {
    sep_xyzw!(($a, $b, $c, $d) => let sqa, sqb, sqc, sqd = (square($a.x - e.x) + square($a.y - e.y)).dynamic() + square($a.z - e.z).dynamic());
    sep_xyzw6!(($a, $b, $c, $d) => let cof1, cof2, cof3, cof4, cof5, cof6 =
        paste!((two_product(($a.x - e.x), ($b.y - e.y)) - two_product(($a.y - e.y), ($b.x - e.x))).dynamic() *
               ([<sq$d>].scale_expansion($c.z - e.z) - [<sq$c>].scale_expansion($d.z - e.z))));
    let det = (cof1 + cof2 + cof3) + (cof4 + cof5 + cof6); // Balanced summation is faster
    let det_approx = det.approximate();

    if det_approx.abs() >= det_sum * IN_SPHERE_BOUND_B {
        return det_approx;
    }

    sep_xyzw!(($a, $b, $c, $d) => let xa, xb, xc, xd = two_sum($a.x, -e.x));
    sep_xyzw!(($a, $b, $c, $d) => let ya, yb, yc, yd = two_sum($a.y, -e.y));
    sep_xyzw!(($a, $b, $c, $d) => let za, zb, zc, zd = two_sum($a.z, -e.z));

    sep_xyzw!(($a, $b, $c, $d) => let sqa, sqb, sqc, sqd =
        paste!([<x$a>][1] * [<x$a>][1] + [<y$a>][1] * [<y$a>][1] + [<z$a>][1] * [<z$a>][1]));
    sep_xyzw!(($a, $b, $c, $d) => let sra, srb, src, srd =
        paste!(([<x$a>][0] * [<x$a>][1] + [<y$a>][0] * [<y$a>][1] + [<z$a>][0] * [<z$a>][1]) * 2.0));

    sep_xyzw6!(($a, $b, $c, $d) => let c11, c21, c31, c41, c51, c61 =
        paste!((([<x$a>][0] * [<y$b>][1] + [<x$a>][1] * [<y$b>][0]) - ([<y$a>][0] * [<x$b>][1] + [<y$a>][1] * [<x$b>][0])) *
               ([<z$c>][1] * [<sq$d>] - [<sq$c>] * [<z$d>][1])));
    sep_xyzw6!(($a, $b, $c, $d) => let c12, c22, c32, c42, c52, c62 =
        paste!(([<x$a>][1] * [<y$b>][1] - [<y$a>][1] * [<x$b>][1]) *
               (([<z$c>][0] * [<sq$d>] + [<z$c>][1] * [<sr$d>]) - ([<sr$c>] * [<z$d>][1] + [<sq$c>] * [<z$d>][0]))));
    let det = ((c11 + c21 + c31) + (c41 + c51 + c61))
        + ((c12 + c22 + c32) + (c42 + c52 + c62))
        + det_approx;

    if det.abs() >= IN_SPHERE_BOUND_C1 * det_approx.abs() + IN_SPHERE_BOUND_C2 * det_sum {
        det
    } else {
        in_sphere_adapt2(a, b, c, d, e)
    }
}

fn in_sphere_adapt2(a: Vec3, b: Vec3, c: Vec3, d: Vec3, e: Vec3) -> f64 {
    // Exact result time!
    // Calculating the 5x5 determinant because it's less straining on the stack
    sep_5!(($a, $b, $c, $d, $e) => let sqa, sqb, sqc, sqd, sqe =
        (square($a.x) + square($a.y)).dynamic() + square($a.z).dynamic());
    sep_5!(($a, $b, $c, $d, $e) => let ab, bc, cd, de, ea = (two_product($a.x, $b.y) - two_product($b.x, $a.y)).dynamic());
    sep_5!(($a, $b, $c, $d, $e) => let ac, bd, ce, da, eb = (two_product($a.x, $c.y) - two_product($c.x, $a.y)).dynamic());
    sep_5!(($a, $b, $c, $d, $e) => let nab, nbc, ncd, nde, nea = paste!(
        [<$c$d>].scale_expansion($e.z) + [<$d$e>].scale_expansion($c.z) + [<$c$e>].scale_expansion(-$d.z)));
    sep_5!(($a, $b, $c, $d, $e) => let nad, nbe, nca, ndb, nec = paste!(
        [<$b$c>].scale_expansion($e.z) + [<$c$e>].scale_expansion($b.z) + [<$e$b>].scale_expansion($c.z)));
    sep_5!(($a, $b, $c, $d, $e) => let cof1, cof2, cof3, cof4, cof5 = paste!(
        (&[<n$e$a>] * &[<sq$e>] - &[<n$a$d>] * &[<sq$d>]) + (&[<n$c$a>] * &[<sq$c>] - &[<n$a$b>] * &[<sq$b>])));
    let det = (cof1 + cof2) + (cof3 + cof4) + cof5;
    det.highest_magnitude()
}

/// Compares the magnitude of `a` and `b`
/// and returns a positive number if `a`'s magnitude is greater,
/// a negative number if `b`'s magnitude is greater,
/// and 0 if their magnitudes equal.
pub fn magnitude_cmp_2d(a: Vec2, b: Vec2) -> f64 {
    let sqa = a.x * a.x + a.y * a.y;
    let sqb = b.x * b.x + b.y * b.y;
    let det = sqa - sqb;

    if det.abs() >= (sqa + sqb) * MAGNITUDE_CMP_2D_BOUND_A {
        det
    } else {
        magnitude_cmp_2d_adapt(a, b)
    }
}

fn magnitude_cmp_2d_adapt(a: Vec2, b: Vec2) -> f64 {
    let sqa = (square(a.x) + square(a.y)).dynamic();
    let sqb = (square(b.x) + square(b.y)).dynamic();
    (sqa - sqb).highest_magnitude()
}

/// Compares the magnitude of `a` and `b`
/// and returns a positive number if `a`'s magnitude is greater,
/// a negative number if `b`'s magnitude is greater,
/// and 0 if their magnitudes equal.
pub fn magnitude_cmp_3d(a: Vec3, b: Vec3) -> f64 {
    let sqa = a.x * a.x + a.y * a.y + a.z * a.z;
    let sqb = b.x * b.x + b.y * b.y + b.z * b.z;
    let det = sqa - sqb;

    if det.abs() >= (sqa + sqb) * MAGNITUDE_CMP_3D_BOUND_A {
        det
    } else {
        magnitude_cmp_3d_adapt(a, b)
    }
}

fn magnitude_cmp_3d_adapt(a: Vec3, b: Vec3) -> f64 {
    let sqa = (square(a.x) + square(a.y)).dynamic() + square(a.z).dynamic();
    let sqb = (square(b.x) + square(b.y)).dynamic() + square(b.z).dynamic();
    (sqa - sqb).highest_magnitude()
}

/// Compares the distance of `a` and `b` to `c`
/// and returns a positive number if `a`'s distance to `c` is greater,
/// a negative number if `b`'s distance to `c` is greater,
/// and 0 if their distances equal.
pub fn distance_cmp_3d(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    let ac = a - c;
    let bc = b - c;
    let sqa = ac.x * ac.x + ac.y * ac.y + ac.z * ac.z;
    let sqb = bc.x * bc.x + bc.y * bc.y + bc.z * bc.z;
    let det = sqa - sqb;

    if det.abs() >= (sqa + sqb) * DISTANCE_CMP_3D_BOUND_A {
        det
    } else {
        distance_cmp_3d_adapt(a, b, c)
    }
}

fn distance_cmp_3d_adapt(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    sep_xyz!(($x, $y, $z) => let ax, ay, az = two_sum(a.$x, -c.$x).dynamic());
    sep_xyz!(($x, $y, $z) => let bx, by, bz = two_sum(b.$x, -c.$x).dynamic());
    let sqa = &ax * &ax + &ay * &ay + &az * &az;
    let sqb = &bx * &bx + &by * &by + &bz * &bz;
    (sqa - sqb).highest_magnitude()
}

/// Computes the determinant of the following matrix
/// ```notrust
/// ┌─                       ─┐
/// │ a.x   a.x^2 + a.y^2   1 │
/// │ b.x   b.x^2 + b.y^2   1 │
/// │ c.x   c.x^2 + c.y^2   1 │
/// └─                       ─┘
/// ```
/// and returns a number with the same sign as the determinant,
/// or 0 if the determinant is 0
pub fn sign_det_x_x2y2(a: Vec2, b: Vec2, c: Vec2) -> f64 {
    sep_xyz!(($a, $b, $c) => let sqa, sqb, sqc = $a.x * $a.x + $a.y * $a.y);
    sep_xyz!(($a, $b, $c) => let cof1, cof2, cof3 = paste!([<sq$a>] * ($c.x - $b.x)));
    let det = cof1 + cof2 + cof3;

    sep_xyz!(($a, $b, $c) => let cof1_sum, cof2_sum, cof3_sum = paste!([<sq$a>] * ($c.x - $b.x).abs()));
    let det_sum = cof1_sum + cof2_sum + cof3_sum;

    if det.abs() >= det_sum * SIGN_DET_X_X2Y2_BOUND_A {
        det
    } else {
        sign_det_x_x2y2_adapt(a, b, c)
    }
}

fn sign_det_x_x2y2_adapt(a: Vec2, b: Vec2, c: Vec2) -> f64 {
    // Compute exact value at least for now
    sep_xyz!(($a, $b, $c) => let sqa, sqb, sqc = (square($a.x) + square($a.y)).dynamic());
    sep_xyz!(($a, $b, $c) => let cof1, cof2, cof3 = paste!([<sq$a>] * two_sum($c.x, -$b.x).dynamic()));
    let det = cof1 + cof2 + cof3;
    det.highest_magnitude()
}

/// Computes the determinant of the following matrix
/// ```notrust
/// ┌─                               ─┐
/// │ a.x   a.x^2 + a.y^2 + a.z^2   1 │
/// │ b.x   b.x^2 + b.y^2 + b.z^2   1 │
/// │ c.x   c.x^2 + c.y^2 + c.z^2   1 │
/// └─                               ─┘
/// ```
/// and returns a number with the same sign as the determinant,
/// or 0 if the determinant is 0
pub fn sign_det_x_x2y2z2(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    sep_xyz!(($a, $b, $c) => let sqa, sqb, sqc = $a.x * $a.x + $a.y * $a.y + $a.z * $a.z);
    sep_xyz!(($a, $b, $c) => let cof1, cof2, cof3 = paste!([<sq$a>] * ($c.x - $b.x)));
    let det = cof1 + cof2 + cof3;

    sep_xyz!(($a, $b, $c) => let cof1_sum, cof2_sum, cof3_sum = paste!([<sq$a>] * ($c.x - $b.x).abs()));
    let det_sum = cof1_sum + cof2_sum + cof3_sum;

    if det.abs() >= det_sum * SIGN_DET_X_X2Y2Z2_BOUND_A {
        det
    } else {
        sign_det_x_x2y2z2_adapt(a, b, c)
    }
}

fn sign_det_x_x2y2z2_adapt(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    // Compute exact value at least for now
    sep_xyz!(($a, $b, $c) => let sqa, sqb, sqc = (square($a.x) + square($a.y)).dynamic() + square($a.z).dynamic());
    sep_xyz!(($a, $b, $c) => let cof1, cof2, cof3 = paste!([<sq$a>] * two_sum($c.x, -$b.x).dynamic()));
    let det = cof1 + cof2 + cof3;
    det.highest_magnitude()
}

/// Computes the determinant of the following matrix
/// ```notrust
/// ┌─                                     ─┐
/// │ a.x   a.y   a.x^2 + a.y^2 + a.z^2   1 │
/// │ b.x   b.y   b.x^2 + b.y^2 + b.z^2   1 │
/// │ c.x   c.y   c.x^2 + c.y^2 + c.z^2   1 │
/// │ d.x   d.y   d.x^2 + d.y^2 + d.z^2   1 │
/// └─                                     ─┘
/// ```
/// and returns a number with the same sign as the determinant,
/// or 0 if the determinant is 0
pub fn sign_det_x_y_x2y2z2(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> f64 {
    sep_xyzw!(($a, $b, $c, $d) => let sqa, sqb, sqc, sqd = $a.x * $a.x + $a.y * $a.y + $a.z * $a.z);
    sep_xyzw6!(($a, $b, $c, $d) => let cof1, cof2, cof3, cof4, cof5, cof6 =
        paste!(($a.x * $b.y - $a.y * $b.x) * ([<sq$c>] - [<sq$d>])));
    let det = (cof1 + cof2 + cof3) + (cof4 + cof5 + cof6);

    sep_xyzw6!(($a, $b, $c, $d) => let cof1_sum, cof2_sum, cof3_sum, cof4_sum, cof5_sum, cof6_sum =
        paste!((($a.x * $b.y).abs() + ($a.y * $b.x).abs()) * ([<sq$c>] + [<sq$d>])));
    let det_sum = (cof1_sum + cof2_sum + cof3_sum) + (cof4_sum + cof5_sum + cof6_sum);

    if det.abs() >= det_sum * SIGN_DET_X_Y_X2Y2Z2_BOUND_A {
        det
    } else {
        sign_det_x_y_x2y2z2_adapt(a, b, c, d)
    }
}

fn sign_det_x_y_x2y2z2_adapt(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> f64 {
    // Compute exact value at least for now
    sep_xyzw!(($a, $b, $c, $d) => let sqa, sqb, sqc, sqd = (square($a.x) + square($a.y)).dynamic() + square($a.z).dynamic());
    sep_xyzw6!(($a, $b, $c, $d) => let cof1, cof2, cof3, cof4, cof5, cof6 =
        paste!((two_product($a.x, $b.y) - two_product($a.y, $b.x)).dynamic() * (&[<sq$c>] - &[<sq$d>])));
    let det = (cof1 + cof2 + cof3) + (cof4 + cof5 + cof6);
    det.highest_magnitude()
}

/// Computes the following sum of determinants
/// ```notrust
/// │ a.x   a.y   a.x^2 + a.y^2 + a.z^2   1 │     │ i.x   i.y   i.z   1 │
/// │ b.x   b.y   b.x^2 + b.y^2 + b.z^2   1 │ + 2u│ j.x   j.y   j.z   1 │
/// │ c.x   c.y   c.x^2 + c.y^2 + c.z^2   1 │     │ k.x   k.y   k.z   1 │
/// │ d.x   d.y   d.x^2 + d.y^2 + d.z^2   1 │     │ l.x   l.y   l.z   1 │
/// ```
/// and returns a number with the same sign as the sum,
/// or 0 if the sum is 0
pub fn sign_det_x_y_x2y2z2_plus_2x_det_x_y_z(a: Vec3, b: Vec3, c: Vec3, d: Vec3, u: f64, i: Vec3, j: Vec3, k: Vec3, l: Vec3) -> f64 {
    sep_xyzw!(($a, $b, $c, $d) => let sqa, sqb, sqc, sqd = $a.x * $a.x + $a.y * $a.y + $a.z * $a.z);
    sep_xyzw6!(($a, $b, $c, $d) => let cof1, cof2, cof3, cof4, cof5, cof6 =
        paste!(($a.x * $b.y - $a.y * $b.x) * ([<sq$c>] - [<sq$d>])));
    let det1 = (cof1 + cof2 + cof3) + (cof4 + cof5 + cof6);

    sep_xyzw6!(($a, $b, $c, $d) => let cof1_sum, cof2_sum, cof3_sum, cof4_sum, cof5_sum, cof6_sum =
        paste!((($a.x * $b.y).abs() + ($a.y * $b.x).abs()) * ([<sq$c>] + [<sq$d>])));
    let det_sum1 = (cof1_sum + cof2_sum + cof3_sum) + (cof4_sum + cof5_sum + cof6_sum);

    sep_xyzw6!(($i, $j, $k, $l) => let cof1, cof2, cof3, cof4, cof5, cof6 =
        paste!(($i.x * $j.y - $i.y * $j.x) * ($k.z - $l.z)));
    let det2 = ((cof1 + cof2 + cof3) + (cof4 + cof5 + cof6)) * u * 2.0;

    sep_xyzw6!(($i, $j, $k, $l) => let cof1_sum, cof2_sum, cof3_sum, cof4_sum, cof5_sum, cof6_sum =
        paste!((($i.x * $j.y).abs() + ($i.y * $j.x).abs()) * ($k.z - $l.z).abs()));
    let det_sum2 = ((cof1_sum + cof2_sum + cof3_sum) + (cof4_sum + cof5_sum + cof6_sum)) * u.abs() * 2.0;

    if (det1 + det2).abs() >= (det_sum1 + det_sum2) * SIGN_DET_X_Y_X2Y2Z2_PLUS_2X_DET_X_Y_Z_BOUND_A {
        det1 + det2
    } else {
        sign_det_x_y_x2y2z2_plus_2x_det_x_y_z_adapt(a, b, c, d, u, i, j, k, l)
    }
}

fn sign_det_x_y_x2y2z2_plus_2x_det_x_y_z_adapt(a: Vec3, b: Vec3, c: Vec3, d: Vec3, u: f64, i: Vec3, j: Vec3, k: Vec3, l: Vec3) -> f64 {
    // Compute exact value at least for now
    sep_xyzw!(($a, $b, $c, $d) => let sqa, sqb, sqc, sqd = (square($a.x) + square($a.y)).dynamic() + square($a.z).dynamic());
    sep_xyzw6!(($a, $b, $c, $d) => let cof1, cof2, cof3, cof4, cof5, cof6 =
        paste!((two_product($a.x, $b.y) - two_product($a.y, $b.x)).dynamic() * (&[<sq$c>] - &[<sq$d>])));
    let det1 = (cof1 + cof2 + cof3) + (cof4 + cof5 + cof6);

    sep_xyzw6!(($i, $j, $k, $l) => let cof1, cof2, cof3, cof4, cof5, cof6 =
        paste!((two_product($i.x, $j.y) - two_product($i.y, $j.x)).dynamic() * two_sum($k.z, -$l.z).dynamic()));
    let det2 = ((cof1 + cof2 + cof3) + (cof4 + cof5 + cof6)).scale_expansion(u * 2.0);

    (det1 + det2).highest_magnitude()
}

/// Computes the following sum of determinants. 
/// ```notrust
/// │ a.x   a.x^2 + a.y^2 + a.z^2   1 │     │ i.x   i.y   1 │
/// │ b.x   b.x^2 + b.y^2 + b.z^2   1 │ + 2u│ j.x   j.y   1 │
/// │ c.x   c.x^2 + c.y^2 + c.z^2   1 │     │ k.x   k.y   1 │
/// ```
/// and returns a number with the same sign as the sum,
/// or 0 if the sum is 0
pub fn sign_det_x_x2y2z2_plus_2x_det_x_y(a: Vec3, b: Vec3, c: Vec3, u: f64, i: Vec2, j: Vec2, k: Vec2) -> f64 {
    sep_xyz!(($a, $b, $c) => let sqa, sqb, sqc = $a.x * $a.x + $a.y * $a.y + $a.z * $a.z);
    sep_xyz!(($a, $b, $c) => let cof1, cof2, cof3 = paste!([<sq$a>] * ($c.x - $b.x)));
    let det1 = cof1 + cof2 + cof3;

    sep_xyz!(($a, $b, $c) => let cof1_sum, cof2_sum, cof3_sum = paste!([<sq$a>] * ($c.x - $b.x).abs()));
    let det_sum1 = cof1_sum + cof2_sum + cof3_sum;

    sep_xyz!(($i, $j, $k) => let cof1, cof2, cof3 = paste!($i.x * ($j.y - $k.y)));
    let det2 = (cof1 + cof2 + cof3) * u * 2.0;

    sep_xyz!(($i, $j, $k) => let cof1_sum, cof2_sum, cof3_sum = paste!($i.x.abs() * ($j.y - $k.y).abs()));
    let det_sum2 = (cof1_sum + cof2_sum + cof3_sum) * u.abs() * 2.0;

    if (det1 + det2).abs() >= (det_sum1 + det_sum2) * SIGN_DET_X_X2Y2Z2_PLUS_2X_DET_X_Y_BOUND_A {
        det1 + det2
    } else {
        sign_det_x_x2y2z2_plus_2x_det_x_y_adapt(a, b, c, u, i, j, k)
    }
}

fn sign_det_x_x2y2z2_plus_2x_det_x_y_adapt(a: Vec3, b: Vec3, c: Vec3, u: f64, i: Vec2, j: Vec2, k: Vec2) -> f64 {
    // Compute exact value at least for now
    sep_xyz!(($a, $b, $c) => let sqa, sqb, sqc = (square($a.x) + square($a.y)).dynamic() + square($a.z).dynamic());
    sep_xyz!(($a, $b, $c) => let cof1, cof2, cof3 = paste!([<sq$a>] * two_sum($c.x, -$b.x).dynamic()));
    let det1 = cof1 + cof2 + cof3;

    sep_xyz!(($i, $j, $k) => let cof1, cof2, cof3 = paste!(two_sum($j.y, -$k.y).dynamic().scale_expansion($i.x)));
    let det2 = (cof1 + cof2 + cof3).scale_expansion(2.0 * u);

    (det1 + det2).highest_magnitude()
}

/// Computes the following sum of determinants. 
/// ```notrust
/// │ a.x^2 + a.y^2 + a.z^2   1 │ + 2u│ i   1 │
/// │ b.x^2 + b.y^2 + b.z^2   1 │     │ j   1 │
/// ```
/// and returns a number with the same sign as the sum,
/// or 0 if the sum is 0
pub fn sign_det_x2y2z2_plus_2x_det_x(a: Vec3, b: Vec3, u: f64, i: f64, j: f64) -> f64 {
    let sqa = a.x * a.x + a.y * a.y + a.z * a.z;
    let sqb = b.x * b.x + b.y * b.y + b.z * b.z;
    let det1 = sqa - sqb;
    let det2 = (i - j) * u * 2.0;
    let det_sum = sqa + sqb + det2.abs();

    if (det1 + det2).abs() >= det_sum * SIGN_DET_X2Y2Z2_PLUS_2X_DET_X_BOUND_A {
        det1 + det2
    } else {
        sign_det_x2y2z2_plus_2x_det_x_adapt(a, b, u, i, j)
    }
}

fn sign_det_x2y2z2_plus_2x_det_x_adapt(a: Vec3, b: Vec3, u: f64, i: f64, j: f64) -> f64 {
    let sqa = (square(a.x) + square(a.y)).dynamic() + square(a.z).dynamic();
    let sqb = (square(b.x) + square(b.y)).dynamic() + square(b.z).dynamic();
    let diff = two_sum(i, -j).scale_expansion(u * 2.0).dynamic();
    (sqa - sqb + diff).highest_magnitude()
}


/// Computes the following sum of determinants. 
/// ```notrust
/// │ a.x   a.x^2 + a.y^2   1 │     │ i.x   i.y   1 │
/// │ b.x   b.x^2 + b.y^2   1 │ + 2u│ j.x   j.y   1 │
/// │ c.x   c.x^2 + c.y^2   1 │     │ k.x   k.y   1 │
/// ```
/// and returns a number with the same sign as the sum,
/// or 0 if the sum is 0
pub fn sign_det_x_x2y2_plus_2x_det_x_y(a: Vec2, b: Vec2, c: Vec2, u: f64, i: Vec2, j: Vec2, k: Vec2) -> f64 {
    sep_xyz!(($a, $b, $c) => let sqa, sqb, sqc = $a.x * $a.x + $a.y * $a.y);
    sep_xyz!(($a, $b, $c) => let cof1, cof2, cof3 = paste!([<sq$a>] * ($c.x - $b.x)));
    let det1 = cof1 + cof2 + cof3;

    sep_xyz!(($a, $b, $c) => let cof1_sum, cof2_sum, cof3_sum = paste!([<sq$a>] * ($c.x - $b.x).abs()));
    let det_sum1 = cof1_sum + cof2_sum + cof3_sum;

    sep_xyz!(($i, $j, $k) => let cof1, cof2, cof3 = paste!($i.x * ($j.y - $k.y)));
    let det2 = (cof1 + cof2 + cof3) * u * 2.0;

    sep_xyz!(($i, $j, $k) => let cof1_sum, cof2_sum, cof3_sum = paste!($i.x.abs() * ($j.y - $k.y).abs()));
    let det_sum2 = (cof1_sum + cof2_sum + cof3_sum) * u.abs() * 2.0;

    if (det1 + det2).abs() >= (det_sum1 + det_sum2) * SIGN_DET_X_X2Y2_PLUS_2X_DET_X_Y_BOUND_A {
        det1 + det2
    } else {
        sign_det_x_x2y2_plus_2x_det_x_y_adapt(a, b, c, u, i, j, k)
    }
}

fn sign_det_x_x2y2_plus_2x_det_x_y_adapt(a: Vec2, b: Vec2, c: Vec2, u: f64, i: Vec2, j: Vec2, k: Vec2) -> f64 {
    // Compute exact value at least for now
    sep_xyz!(($a, $b, $c) => let sqa, sqb, sqc = (square($a.x) + square($a.y)).dynamic());
    sep_xyz!(($a, $b, $c) => let cof1, cof2, cof3 = paste!([<sq$a>] * two_sum($c.x, -$b.x).dynamic()));
    let det1 = cof1 + cof2 + cof3;

    sep_xyz!(($i, $j, $k) => let cof1, cof2, cof3 = paste!(two_sum($j.y, -$k.y).dynamic().scale_expansion($i.x)));
    let det2 = (cof1 + cof2 + cof3).scale_expansion(2.0 * u);

    (det1 + det2).highest_magnitude()
}

/// Computes the following sum of determinants. 
/// ```notrust
/// │ a.x^2 + a.y^2   1 │ + 2u│ i   1 │
/// │ b.x^2 + b.y^2   1 │     │ j   1 │
/// ```
/// and returns a number with the same sign as the sum,
/// or 0 if the sum is 0
pub fn sign_det_x2y2_plus_2x_det_x(a: Vec2, b: Vec2, u: f64, i: f64, j: f64) -> f64 {
    let sqa = a.x * a.x + a.y * a.y;
    let sqb = b.x * b.x + b.y * b.y;
    let det1 = sqa - sqb;
    let det2 = (i - j) * u * 2.0;
    let det_sum = sqa + sqb + det2.abs();

    if (det1 + det2).abs() >= det_sum * SIGN_DET_X2Y2_PLUS_2X_DET_X_BOUND_A {
        det1 + det2
    } else {
        sign_det_x2y2_plus_2x_det_x_adapt(a, b, u, i, j)
    }
}

fn sign_det_x2y2_plus_2x_det_x_adapt(a: Vec2, b: Vec2, u: f64, i: f64, j: f64) -> f64 {
    let sqa = (square(a.x) + square(a.y)).dynamic();
    let sqb = (square(b.x) + square(b.y)).dynamic();
    let diff = two_sum(i, -j).scale_expansion(u * 2.0).dynamic();
    (sqa - sqb + diff).highest_magnitude()
}

/// Computes the *simple intersection comparison* on coordinate `i`:
/// Given segments *a* = `[a0, a1]` and *b* = `[b0, b1]` and point `p`,
/// where each segment is sorted lexicographically and *b* points clockwise from *a*,
/// returns
/// * if `i == 0`
///     * a negative number if `p` lies left of the intersection of *a* and *b*
///     * a positive number if `p` lies right of the intersection of *a* and *b*
///     * zero if `p` is on the same vertical line as the intersection of *a* and *b*
/// * if `i == 1` (assumes that the function returns zero if `i == 0`)
///     * a negative number if `p` lies below of the intersection of *a* and *b*
///     * a positive number if `p` lies above of the intersection of *a* and *b*
///     * zero if `p` is on the same horizontal line (and thus is the same point) as the intersection of *a* and *b*
/// 
/// This is equivalent to the negative of the following sum:
/// ```notrust
/// (a1[i] - a0[i]) | b1.x - b0.x    b0.x - a0.x | + (p[i] - a0[i]) | a1.x - a0.x    b1.x - b0.x |
///                 | b1.y - b0.y    b0.y - a0.y |                  | a1.y - a0.y    b1.y - b0.y |
/// ```
pub fn simple_intersect_compare_i(p: Vec2, a0: Vec2, a1: Vec2, b0: Vec2, b1: Vec2, i: usize) -> f64 {
    // Implementing as 3D determinant. Optimization is premature
    -cross_3d(
        vector![a0[i], a0.x, a0.y],
        vector![a1[i], a1.x, a1.y],
        vector![0.0, b0.x, b0.y],
        vector![0.0, b1.x, b1.y],
        vector![a0[i], a0.x, a0.y],
        vector![p[i], b0.x, b0.y]
    )
}

/// Computes the *complex intersection comparison* on coordinate `i`:
/// Given segments *a* = `[a0, a1]`, *b* = `[b0, b1]`, *c* = `[c0, c1]`, and *d* = `[d0, d1]`,
/// where each segment is sorted lexicographically, *b* points clockwise from *a*, and *d* points clockwise from *c*,
/// returns (where *x*×*y* is the intersection of segment *x* and segment *y*)
/// * if `i == 0`
///     * a negative number if *a*×*b* lies left of *c*×*d*
///     * a positive number if *a*×*b* lies right of *c*×*d*
///     * zero if *a*×*b* and *c*×*d* are on the same vertical line
/// * if `i == 1` (assumes that the function returns zero if `i == 0`)
///     * a negative number if *a*×*b* lies below of *c*×*d*
///     * a positive number if *a*×*b* lies above of *c*×*d*
///     * zero if *a*×*b* and *c*×*d* are on the same horizontal line (and thus are the same point)
pub fn complex_intersect_compare_i(a0: Vec2, a1: Vec2, b0: Vec2, b1: Vec2, c0: Vec2, c1: Vec2, d0: Vec2, d1: Vec2, i: usize) -> f64 {
    let d_a = a1 - a0;
    let d_b = b1 - b0;
    let d_c = c1 - c0;
    let d_d = d1 - d0;
    let d_ab = b0 - a0;
    let d_cd = d0 - c0;
    let d_ca = a0 - c0;
    let det_ab = d_a.x * d_b.y - d_a.y * d_b.x;
    let det_cd = d_c.x * d_d.y - d_c.y * d_d.x;
    let det_bab = d_b.x * d_ab.y - d_b.y * d_ab.x;
    let det_dcd = d_d.x * d_cd.y - d_d.y * d_cd.x;
    let det_ab_sum = (d_a.x * d_b.y).abs() + (d_a.y * d_b.x).abs();
    let det_cd_sum = (d_c.x * d_d.y).abs() + (d_c.y * d_d.x).abs();
    let det_bab_sum = (d_b.x * d_ab.y).abs() + (d_b.y * d_ab.x).abs();
    let det_dcd_sum = (d_d.x * d_cd.y).abs() + (d_d.y * d_cd.x).abs();
    let result = -d_a[i] * det_bab * det_cd + d_c[i] * det_dcd * det_ab + d_ca[i] * det_ab * det_cd;
    let result_sum = (d_a[i] * det_bab_sum * det_cd_sum).abs() + (d_c[i] * det_dcd_sum * det_ab_sum).abs() + (d_ca[i] * det_ab_sum * det_cd_sum).abs();

    if result.abs() >= result_sum * COMPLEX_INTERSECT_COMPARE_BOUND {
        result
    } else {
        complex_intersect_compare_i_adapt(a0, a1, b0, b1, c0, c1, d0, d1, i)
    }
}
pub fn complex_intersect_compare_i_adapt(a0: Vec2, a1: Vec2, b0: Vec2, b1: Vec2, c0: Vec2, c1: Vec2, d0: Vec2, d1: Vec2, i: usize) -> f64 {
    // Go for exact value for now
    sep_xy!(($x, $y) => let d_ax, d_ay = two_sum(a1.$x, -a0.$x));
    sep_xy!(($x, $y) => let d_bx, d_by = two_sum(b1.$x, -b0.$x));
    sep_xy!(($x, $y) => let d_cx, d_cy = two_sum(c1.$x, -c0.$x));
    sep_xy!(($x, $y) => let d_dx, d_dy = two_sum(d1.$x, -d0.$x));
    sep_xy!(($x, $y) => let d_abx, d_aby = two_sum(b0.$x, -a0.$x));
    sep_xy!(($x, $y) => let d_cdx, d_cdy = two_sum(d0.$x, -c0.$x));
    sep_xy!(($x, $y) => let d_cax, d_cay = two_sum(a0.$x, -c0.$x));

    let det_ab_hi = (two_product(d_ax[1], d_by[1]) - two_product(d_ay[1], d_bx[1])).dynamic();
    let det_cd_hi = (two_product(d_cx[1], d_dy[1]) - two_product(d_cy[1], d_dx[1])).dynamic();
    let det_bab_hi = (two_product(d_bx[1], d_aby[1]) - two_product(d_by[1], d_abx[1])).dynamic();
    let det_dcd_hi = (two_product(d_dx[1], d_cdy[1]) - two_product(d_dy[1], d_cdx[1])).dynamic();
    if [d_ax[0], d_ay[0], d_bx[0], d_by[0], d_cx[0], d_cy[0], d_dx[0], d_dy[0], d_abx[0], d_aby[0], d_cdx[0], d_cdy[0], d_cax[0], d_cay[0]] == [0.0; 14] {
        let det_ab_cd  = &det_ab_hi * &det_cd_hi;
        let det_ab_dcd = &det_ab_hi * &det_dcd_hi;
        let det_cd_bab = &det_cd_hi * &det_bab_hi;
        let d_ai = &[d_ax[1], d_ay[1]][i];
        let d_ci = &[d_cx[1], d_cy[1]][i];
        let d_cai = &[d_cax[1], d_cay[1]][i];
        let result_ab_cd = det_ab_cd.scale_expansion(*d_cai);
        let result_ab_dcd = det_ab_dcd.scale_expansion(*d_ci);
        let result_cd_bab = det_cd_bab.scale_expansion(*d_ai);
        let result = result_ab_cd + result_ab_dcd - result_cd_bab;
        return result.highest_magnitude();
    }

    let det_ab_hi = det_ab_hi;
    let det_ab_m1 = (two_product(d_ax[0], d_by[1]) - two_product(d_ay[0], d_bx[1])).dynamic();
    let det_ab_m2 = (two_product(d_ax[1], d_by[0]) - two_product(d_ay[1], d_bx[0])).dynamic();
    let det_ab_lo = (two_product(d_ax[0], d_by[0]) - two_product(d_ay[0], d_bx[0])).dynamic();
    let det_ab = (det_ab_hi + det_ab_m1) + (det_ab_m2 + det_ab_lo);
    let det_cd_hi = det_cd_hi;
    let det_cd_m1 = (two_product(d_cx[0], d_dy[1]) - two_product(d_cy[0], d_dx[1])).dynamic();
    let det_cd_m2 = (two_product(d_cx[1], d_dy[0]) - two_product(d_cy[1], d_dx[0])).dynamic();
    let det_cd_lo = (two_product(d_cx[0], d_dy[0]) - two_product(d_cy[0], d_dx[0])).dynamic();
    let det_cd = (det_cd_hi + det_cd_m1) + (det_cd_m2 + det_cd_lo);
    let det_bab_hi = det_bab_hi;
    let det_bab_m1 = (two_product(d_bx[0], d_aby[1]) - two_product(d_by[0], d_abx[1])).dynamic();
    let det_bab_m2 = (two_product(d_bx[1], d_aby[0]) - two_product(d_by[1], d_abx[0])).dynamic();
    let det_bab_lo = (two_product(d_bx[0], d_aby[0]) - two_product(d_by[0], d_abx[0])).dynamic();
    let det_bab = (det_bab_hi + det_bab_m1) + (det_bab_m2 + det_bab_lo);
    let det_dcd_hi = det_dcd_hi;
    let det_dcd_m1 = (two_product(d_dx[0], d_cdy[1]) - two_product(d_dy[0], d_cdx[1])).dynamic();
    let det_dcd_m2 = (two_product(d_dx[1], d_cdy[0]) - two_product(d_dy[1], d_cdx[0])).dynamic();
    let det_dcd_lo = (two_product(d_dx[0], d_cdy[0]) - two_product(d_dy[0], d_cdx[0])).dynamic();
    let det_dcd = (det_dcd_hi + det_dcd_m1) + (det_dcd_m2 + det_dcd_lo);

    let det_ab_cd  = &det_ab * &det_cd;
    let det_ab_dcd = &det_ab * &det_dcd;
    let det_cd_bab = &det_cd * &det_bab;
    let d_ai = &[d_ax, d_ay][i];
    let d_ci = &[d_cx, d_cy][i];
    let d_cai = &[d_cax, d_cay][i];
    let result_ab_cd = det_ab_cd.scale_expansion(d_cai[0]) + det_ab_cd.scale_expansion(d_cai[1]);
    let result_ab_dcd = det_ab_dcd.scale_expansion(d_ci[0]) + det_ab_dcd.scale_expansion(d_ci[1]);
    let result_cd_bab = det_cd_bab.scale_expansion(d_ai[0]) + det_cd_bab.scale_expansion(d_ai[1]);
    let result = result_ab_cd + result_ab_dcd - result_cd_bab;
    result.highest_magnitude()
}

#[cfg(test)]
mod test {
    use std::cmp::Ordering;

    use super::*;
    use nalgebra::{Matrix2, Matrix3};
    use rand::distr::{Distribution, Uniform};
    use rand::seq::SliceRandom;
    use rand::Rng;
    use rand::prelude::IndexedRandom;
    use rand_distr::{UnitCircle, UnitSphere};
    use rand_pcg::Pcg64;
    use rug::{Float, Rational};

    type Mtx2 = Matrix2<f64>;
    type Mtx3 = Matrix3<f64>;

    const PCG_STATE: u128 = 0xcafef00dd15ea5e5;
    const PCG_STREAM: u128 = 0xa02bdbf7bb3c0a7ac28fa16a64abf96;

    fn orient_2d_exact(a: Vec2, b: Vec2, c: Vec2) -> Float {
        const PREC: u32 = (f64::MANTISSA_DIGITS + 1) * 2 + 1;
        let ax = Float::with_val(
            PREC,
            &Float::with_val(PREC, a.x) - &Float::with_val(PREC, c.x),
        );
        let ay = Float::with_val(
            PREC,
            &Float::with_val(PREC, a.y) - &Float::with_val(PREC, c.y),
        );
        let bx = Float::with_val(
            PREC,
            &Float::with_val(PREC, b.x) - &Float::with_val(PREC, c.x),
        );
        let by = Float::with_val(
            PREC,
            &Float::with_val(PREC, b.y) - &Float::with_val(PREC, c.y),
        );
        Float::with_val(PREC, &ax * &by - &ay * &bx)
    }

    fn check_orient_2d(a: Vec2, b: Vec2, c: Vec2) {
        let adapt = orient_2d(a, b, c);
        let exact = orient_2d_exact(a, b, c);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}) gave wrong result: {} vs {}",
            a,
            b,
            c,
            adapt,
            exact
        );
    }

    #[test]
    fn test_orient_2d_uniform_random() {
        // Deterministic, portable RNG
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..10000 {
            let vals = dist.sample_iter(&mut rng).take(6).collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            let c = Vec2::new(vals[4], vals[5]);
            check_orient_2d(a, b, c);
        }
    }

    #[test]
    fn test_orient_2d_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..10000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(6)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            let c = Vec2::new(vals[4], vals[5]);
            check_orient_2d(a, b, c);
        }
    }

    #[test]
    fn test_orient_2d_near_collinear() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-1.0, 1.0).unwrap();
        let fac_dist = Uniform::new(-1.0, 2.0).unwrap();

        for _ in 0..10000 {
            let vals = dist.sample_iter(&mut rng).take(4).collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);

            let fac = fac_dist.sample(&mut rng);
            let c = a + (b - a) * fac;

            check_orient_2d(a, b, c);
        }
    }

    #[test]
    fn test_orient_2d_collinear() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new_inclusive(-4096, 4096).unwrap();
        let fac_dist = Uniform::new_inclusive(-4095, 4096).unwrap();

        for _ in 0..10000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(4)
                .map(|x| (x as f64) / 4096.0)
                .collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);

            let fac = fac_dist.sample(&mut rng);
            let c = a + (b - a) * (fac as f64) / 16.0;

            check_orient_2d(a, b, c);
        }
    }

    fn cross_2d_exact(a: Vec2, b: Vec2, c: Vec2, d: Vec2) -> Float {
        const PREC: u32 = (f64::MANTISSA_DIGITS + 1) * 2 + 1;
        let abx = Float::with_val(
            PREC,
            &Float::with_val(PREC, b.x) - &Float::with_val(PREC, a.x),
        );
        let aby = Float::with_val(
            PREC,
            &Float::with_val(PREC, b.y) - &Float::with_val(PREC, a.y),
        );
        let cdx = Float::with_val(
            PREC,
            &Float::with_val(PREC, d.x) - &Float::with_val(PREC, c.x),
        );
        let cdy = Float::with_val(
            PREC,
            &Float::with_val(PREC, d.y) - &Float::with_val(PREC, c.y),
        );
        Float::with_val(PREC, &abx * &cdy - &aby * &cdx)
    }

    fn check_cross_2d(a: Vec2, b: Vec2, c: Vec2, d: Vec2,) {
        let adapt = cross_2d(a, b, c, d);
        let exact = cross_2d_exact(a, b, c, d);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}) gave wrong result: {} vs {}",
            a,
            b,
            c,
            d,
            adapt,
            exact
        );
    }

    #[test]
    fn test_cross_2d_uniform_random() {
        // Deterministic, portable RNG
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..10000 {
            let vals = dist.sample_iter(&mut rng).take(8).collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            let c = Vec2::new(vals[4], vals[5]);
            let d = Vec2::new(vals[6], vals[7]);
            check_cross_2d(a, b, c, d);
        }
    }

    #[test]
    fn test_cross_2d_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..10000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(8)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            let c = Vec2::new(vals[4], vals[5]);
            let d = Vec2::new(vals[6], vals[7]);
            check_cross_2d(a, b, c, d);
        }
    }

    #[test]
    fn test_cross_2d_near_parallel() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-1.0, 1.0).unwrap();
        let fac_dist = Uniform::new(-1.0, 2.0).unwrap();

        for _ in 0..10000 {
            let vals = dist.sample_iter(&mut rng).take(6).collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            let c = Vec2::new(vals[4], vals[5]);

            let fac = fac_dist.sample(&mut rng);
            let d = c + (b - a) * fac;

            check_cross_2d(a, b, c, d);
        }
    }

    #[test]
    fn test_cross_2d_collinear() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new_inclusive(-4096, 4096).unwrap();
        let fac_dist = Uniform::new_inclusive(-4095, 4096).unwrap();

        for _ in 0..10000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(6)
                .map(|x| (x as f64) / 4096.0)
                .collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            let c = Vec2::new(vals[4], vals[5]);

            let fac = fac_dist.sample(&mut rng);
            let d = c + (b - a) * (fac as f64) / 16.0;

            check_cross_2d(a, b, c, d);
        }
    }

    fn orient_3d_exact(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> Float {
        const PREC: u32 = (f64::MANTISSA_DIGITS + 1) * 3 + 3;
        let ax = Float::with_val(
            PREC,
            &Float::with_val(PREC, a.x) - &Float::with_val(PREC, d.x),
        );
        let ay = Float::with_val(
            PREC,
            &Float::with_val(PREC, a.y) - &Float::with_val(PREC, d.y),
        );
        let az = Float::with_val(
            PREC,
            &Float::with_val(PREC, a.z) - &Float::with_val(PREC, d.z),
        );
        let bx = Float::with_val(
            PREC,
            &Float::with_val(PREC, b.x) - &Float::with_val(PREC, d.x),
        );
        let by = Float::with_val(
            PREC,
            &Float::with_val(PREC, b.y) - &Float::with_val(PREC, d.y),
        );
        let bz = Float::with_val(
            PREC,
            &Float::with_val(PREC, b.z) - &Float::with_val(PREC, d.z),
        );
        let cx = Float::with_val(
            PREC,
            &Float::with_val(PREC, c.x) - &Float::with_val(PREC, d.x),
        );
        let cy = Float::with_val(
            PREC,
            &Float::with_val(PREC, c.y) - &Float::with_val(PREC, d.y),
        );
        let cz = Float::with_val(
            PREC,
            &Float::with_val(PREC, c.z) - &Float::with_val(PREC, d.z),
        );
        let xy = Float::with_val(PREC, &bx * &cy - &by * &cx);
        let yz = Float::with_val(PREC, &by * &cz - &bz * &cy);
        let zx = Float::with_val(PREC, &bz * &cx - &bx * &cz);
        let ab = Float::with_val(PREC, &az * &xy + &ax * &yz);
        Float::with_val(PREC, &ab + &ay * &zx)
    }

    fn check_orient_3d(a: Vec3, b: Vec3, c: Vec3, d: Vec3) {
        let adapt = orient_3d(a, b, c, d);
        let exact = orient_3d_exact(a, b, c, d);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}) gave wrong result: {} vs {}",
            a,
            b,
            c,
            d,
            adapt,
            exact
        );
    }

    #[test]
    fn test_orient_3d_uniform_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(12).collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            let d = Vec3::new(vals[9], vals[10], vals[11]);
            check_orient_3d(a, b, c, d);
        }
    }

    #[test]
    fn test_orient_3d_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(12)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            let d = Vec3::new(vals[9], vals[10], vals[11]);
            check_orient_3d(a, b, c, d);
        }
    }

    #[test]
    fn test_orient_3d_near_coplanar() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-1.0, 1.0).unwrap();
        let fac_dist = Uniform::new(-1.0, 2.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(9).collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);

            let fac = fac_dist.sample(&mut rng);
            let fac2 = fac_dist.sample(&mut rng);
            let d = c + (a + (b - a) * fac - c) * fac2;

            check_orient_3d(a, b, c, d);
        }
    }

    #[test]
    fn test_orient_3d_coplanar() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new_inclusive(-4096, 4096).unwrap();
        let fac_dist = Uniform::new_inclusive(-4095, 4096).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(9)
                .map(|x| (x as f64) / 4096.0)
                .collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);

            let fac = fac_dist.sample(&mut rng) as f64 / 16.0;
            let fac2 = fac_dist.sample(&mut rng) as f64 / 16.0;
            let d = c + (a + (b - a) * fac - c) * fac2;

            check_orient_3d(a, b, c, d);
        }
    }

    fn cross_3d_exact(a: Vec3, b: Vec3, c: Vec3, d: Vec3, e: Vec3, f: Vec3) -> Float {
        const PREC: u32 = (f64::MANTISSA_DIGITS + 1) * 3 + 3;
        let ax = Float::with_val(
            PREC,
            &Float::with_val(PREC, b.x) - &Float::with_val(PREC, a.x),
        );
        let ay = Float::with_val(
            PREC,
            &Float::with_val(PREC, b.y) - &Float::with_val(PREC, a.y),
        );
        let az = Float::with_val(
            PREC,
            &Float::with_val(PREC, b.z) - &Float::with_val(PREC, a.z),
        );
        let bx = Float::with_val(
            PREC,
            &Float::with_val(PREC, d.x) - &Float::with_val(PREC, c.x),
        );
        let by = Float::with_val(
            PREC,
            &Float::with_val(PREC, d.y) - &Float::with_val(PREC, c.y),
        );
        let bz = Float::with_val(
            PREC,
            &Float::with_val(PREC, d.z) - &Float::with_val(PREC, c.z),
        );
        let cx = Float::with_val(
            PREC,
            &Float::with_val(PREC, f.x) - &Float::with_val(PREC, e.x),
        );
        let cy = Float::with_val(
            PREC,
            &Float::with_val(PREC, f.y) - &Float::with_val(PREC, e.y),
        );
        let cz = Float::with_val(
            PREC,
            &Float::with_val(PREC, f.z) - &Float::with_val(PREC, e.z),
        );
        let xy = Float::with_val(PREC, &bx * &cy - &by * &cx);
        let yz = Float::with_val(PREC, &by * &cz - &bz * &cy);
        let zx = Float::with_val(PREC, &bz * &cx - &bx * &cz);
        let ab = Float::with_val(PREC, &az * &xy + &ax * &yz);
        Float::with_val(PREC, &ab + &ay * &zx)
    }

    fn check_cross_3d(a: Vec3, b: Vec3, c: Vec3, d: Vec3, e: Vec3, f: Vec3) {
        let adapt = cross_3d(a, b, c, d, e, f);
        let exact = cross_3d_exact(a, b, c, d, e, f);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}, {}, {}) gave wrong result: {} vs {}",
            a,
            b,
            c,
            d,
            e,
            f,
            adapt,
            exact
        );
    }

    #[test]
    fn test_cross_3d_uniform_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(18).collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            let d = Vec3::new(vals[9], vals[10], vals[11]);
            let e = Vec3::new(vals[12], vals[13], vals[14]);
            let f = Vec3::new(vals[15], vals[16], vals[17]);
            check_cross_3d(a, b, c, d, e, f);
        }
    }

    #[test]
    fn test_cross_3d_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(18)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            let d = Vec3::new(vals[9], vals[10], vals[11]);
            let e = Vec3::new(vals[12], vals[13], vals[14]);
            let f = Vec3::new(vals[15], vals[16], vals[17]);
            check_cross_3d(a, b, c, d, e, f);
        }
    }

    #[test]
    fn test_cross_3d_near_coplanar() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-1.0, 1.0).unwrap();
        let fac_dist = Uniform::new(-1.0, 2.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(15).collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            let d = Vec3::new(vals[9], vals[10], vals[11]);
            let e = Vec3::new(vals[12], vals[13], vals[14]);

            let fac = fac_dist.sample(&mut rng);
            let fac2 = fac_dist.sample(&mut rng);
            let f = e + ((b - a) * fac + (d - c) * fac2);

            check_cross_3d(a, b, c, d, e, f);
        }
    }

    #[test]
    fn test_cross_3d_coplanar() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new_inclusive(-4096, 4096).unwrap();
        let fac_dist = Uniform::new_inclusive(-4095, 4096).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(15)
                .map(|x| (x as f64) / 4096.0)
                .collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            let d = Vec3::new(vals[9], vals[10], vals[11]);
            let e = Vec3::new(vals[12], vals[13], vals[14]);

            let fac = fac_dist.sample(&mut rng) as f64 / 16.0;
            let fac2 = fac_dist.sample(&mut rng) as f64 / 16.0;
            let f = e + ((b - a) * fac + (d - c) * fac2);

            check_cross_3d(a, b, c, d, e, f);
        }
    }

    #[test]
    fn test_simple_intersect_compare_i() {
        // This for now just calls cross_3d, so just check to make sure
        // it does the intended thing
        assert_eq!(simple_intersect_compare_i(vector![1.0, 1.0], vector![0.0, 0.0], vector![2.0, 2.0], vector![0.0, 2.0], vector![2.0, 0.0], 0), 0.0);
        assert_eq!(simple_intersect_compare_i(vector![1.0, 1.0], vector![0.0, 0.0], vector![2.0, 2.0], vector![0.0, 2.0], vector![2.0, 0.0], 1), 0.0);
        assert_eq!(simple_intersect_compare_i(vector![1.0, 2.0], vector![0.0, 0.0], vector![2.0, 2.0], vector![0.0, 2.0], vector![2.0, 0.0], 0), 0.0);
        assert!   (simple_intersect_compare_i(vector![1.0, 2.0], vector![0.0, 0.0], vector![2.0, 2.0], vector![0.0, 2.0], vector![2.0, 0.0], 1) > 0.0);
        assert!   (simple_intersect_compare_i(vector![1.0, 0.0], vector![0.0, 0.0], vector![2.0, 2.0], vector![0.0, 2.0], vector![2.0, 0.0], 1) < 0.0);
        assert!   (simple_intersect_compare_i(vector![2.0, 1.0], vector![0.0, 0.0], vector![2.0, 2.0], vector![0.0, 2.0], vector![2.0, 0.0], 0) > 0.0);
        assert!   (simple_intersect_compare_i(vector![0.0, 1.0], vector![0.0, 0.0], vector![2.0, 2.0], vector![0.0, 2.0], vector![2.0, 0.0], 0) < 0.0);
    }

    fn intersect_exact_i(a0: Vec2, a1: Vec2, b0: Vec2, b1: Vec2, i: usize) -> Option<Rational> {
        let dax  = Rational::from_f64(a1.x).unwrap() - Rational::from_f64(a0.x).unwrap();
        let day  = Rational::from_f64(a1.y).unwrap() - Rational::from_f64(a0.y).unwrap();
        let dbx  = Rational::from_f64(b1.x).unwrap() - Rational::from_f64(b0.x).unwrap();
        let dby  = Rational::from_f64(b1.y).unwrap() - Rational::from_f64(b0.y).unwrap();
        let dabx = Rational::from_f64(b0.x).unwrap() - Rational::from_f64(a0.x).unwrap();
        let daby = Rational::from_f64(b0.y).unwrap() - Rational::from_f64(a0.y).unwrap();

        let det = Rational::from(&day * &dbx) - Rational::from(&dby * &dax);
        if &det == Rational::ZERO { return None }

        let ta = (-dby * dabx + dbx * daby) / det;
        Some(if i == 0 {dax} else {day} * ta + Rational::from_f64(a0[i]).unwrap())
    }

    fn check_complex_intersect_compare_i(a0: Vec2, a1: Vec2, b0: Vec2, b1: Vec2, c0: Vec2, c1: Vec2, d0: Vec2, d1: Vec2, i: usize) -> bool {
        let result = complex_intersect_compare_i(a0, a1, b0, b1, c0, c1, d0, d1, i);
        let intersection = intersect_exact_i(a0, a1, b0, b1, i)
            .zip(intersect_exact_i(c0, c1, d0, d1, i))
            .map(|(ia, ic)| ia.cmp(&ic))
            .unwrap_or(Ordering::Equal);
        assert_eq!(result.partial_cmp(&0.0).unwrap(), intersection,
            "{}-{}, {}-{}, {}-{}, {}-{}, {} gives wrong result: {} does not have the same sign as {:?} - {:?}",
            a0, a1, b0, b1, c0, c1, d0, d1, i,
            result, intersect_exact_i(a0, a1, b0, b1, i), intersect_exact_i(c0, c1, d0, d1, i));
        result == 0.0
    }

    #[test]
    fn test_complex_intersect_compare_i_uniform_random() {
        // Deterministic, portable RNG
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..10000 {
            let vals = dist.sample_iter(&mut rng).take(16).collect::<Vec<_>>();
            let mut a0 = Vec2::new(vals[0], vals[1]);
            let mut a1 = Vec2::new(vals[2], vals[3]);
            let mut b0 = Vec2::new(vals[4], vals[5]);
            let mut b1 = Vec2::new(vals[6], vals[7]);
            let mut c0 = Vec2::new(vals[8], vals[9]);
            let mut c1 = Vec2::new(vals[10], vals[11]);
            let mut d0 = Vec2::new(vals[12], vals[13]);
            let mut d1 = Vec2::new(vals[14], vals[15]);
            if cross_2d(a0, a1, b0, b1) > 0.0 { std::mem::swap(&mut a0, &mut b0); std::mem::swap(&mut a1, &mut b1); }
            if cross_2d(c0, c1, d0, d1) > 0.0 { std::mem::swap(&mut c0, &mut d0); std::mem::swap(&mut c1, &mut d1); }
            if check_complex_intersect_compare_i(a0, a1, b0, b1, c0, c1, d0, d1, 0) {
                check_complex_intersect_compare_i(a0, a1, b0, b1, c0, c1, d0, d1, 1);
            }
        }
    }

    #[test]
    fn test_complex_intersect_compare_i_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..10000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(16)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let mut a0 = Vec2::new(vals[0], vals[1]);
            let mut a1 = Vec2::new(vals[2], vals[3]);
            let mut b0 = Vec2::new(vals[4], vals[5]);
            let mut b1 = Vec2::new(vals[6], vals[7]);
            let mut c0 = Vec2::new(vals[8], vals[9]);
            let mut c1 = Vec2::new(vals[10], vals[11]);
            let mut d0 = Vec2::new(vals[12], vals[13]);
            let mut d1 = Vec2::new(vals[14], vals[15]);
            if cross_2d(a0, a1, b0, b1) > 0.0 { std::mem::swap(&mut a0, &mut b0); std::mem::swap(&mut a1, &mut b1); }
            if cross_2d(c0, c1, d0, d1) > 0.0 { std::mem::swap(&mut c0, &mut d0); std::mem::swap(&mut c1, &mut d1); }
            if check_complex_intersect_compare_i(a0, a1, b0, b1, c0, c1, d0, d1, 0) {
                check_complex_intersect_compare_i(a0, a1, b0, b1, c0, c1, d0, d1, 1);
            }
        }
    }

    #[test]
    fn test_complex_intersect_compare_i_near_same_point() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-1.0, 1.0).unwrap();
        let fac_dist = Uniform::new(-1.0, 2.0).unwrap();

        for _ in 0..10000 {
            let vals = dist.sample_iter(&mut rng).take(10).collect::<Vec<_>>();
            let mut a0 = Vec2::new(vals[0], vals[1]);
            let mut b0 = Vec2::new(vals[2], vals[3]);
            let mut c0 = Vec2::new(vals[4], vals[5]);
            let mut d0 = Vec2::new(vals[6], vals[7]);
            let point = Vec2::new(vals[8], vals[9]);

            let fa = fac_dist.sample(&mut rng);
            let fb = fac_dist.sample(&mut rng);
            let fc = fac_dist.sample(&mut rng);
            let fd = fac_dist.sample(&mut rng);
            let mut a1 = a0 + (point - a0) * fa;
            let mut b1 = b0 + (point - b0) * fb;
            let mut c1 = c0 + (point - c0) * fc;
            let mut d1 = d0 + (point - d0) * fd;
            if cross_2d(a0, a1, b0, b1) > 0.0 { std::mem::swap(&mut a0, &mut b0); std::mem::swap(&mut a1, &mut b1); }
            if cross_2d(c0, c1, d0, d1) > 0.0 { std::mem::swap(&mut c0, &mut d0); std::mem::swap(&mut c1, &mut d1); }
            
            if check_complex_intersect_compare_i(a0, a1, b0, b1, c0, c1, d0, d1, 0) {
                check_complex_intersect_compare_i(a0, a1, b0, b1, c0, c1, d0, d1, 1);
            }
        }
    }

    #[test]
    fn test_complex_intersect_compare_i_same_point() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new_inclusive(-4096, 4096).unwrap();
        let fac_dist = Uniform::new_inclusive(-4096, 4096).unwrap();

        for _ in 0..10000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(10)
                .map(|x| (x as f64) / 4096.0)
                .collect::<Vec<_>>();
            let a0 = Vec2::new(vals[0], vals[1]);
            let b0 = Vec2::new(vals[2], vals[3]);
            let c0 = Vec2::new(vals[4], vals[5]);
            let d0 = Vec2::new(vals[6], vals[7]);
            let point = Vec2::new(vals[8], vals[9]);

            let fa = fac_dist.sample(&mut rng);
            let fb = fac_dist.sample(&mut rng);
            let fc = fac_dist.sample(&mut rng);
            let fd = fac_dist.sample(&mut rng);
            let a1 = a0 + (point - a0) * (fa as f64) / 16.0;
            let b1 = b0 + (point - b0) * (fb as f64) / 16.0;
            let c1 = c0 + (point - c0) * (fc as f64) / 16.0;
            let d1 = d0 + (point - d0) * (fd as f64) / 16.0;
            
            if check_complex_intersect_compare_i(a0, a1, b0, b1, c0, c1, d0, d1, 0) {
                check_complex_intersect_compare_i(a0, a1, b0, b1, c0, c1, d0, d1, 1);
            }
        }
    }

    #[test]
    fn test_complex_intersect_compare_i() {
        // Make sure it does the intended thing
        assert_eq!(complex_intersect_compare_i(
            vector![1.0, 1.0], vector![3.0, 3.0], vector![1.0, 3.0], vector![3.0, 1.0],
            vector![2.0, 0.0], vector![2.0, 3.0], vector![0.0, 2.0], vector![5.0, 2.0], 0), 0.0);
        assert_eq!(complex_intersect_compare_i(
            vector![1.0, 1.0], vector![3.0, 3.0], vector![1.0, 3.0], vector![3.0, 1.0],
            vector![2.0, 0.0], vector![2.0, 3.0], vector![0.0, 2.0], vector![5.0, 2.0], 1), 0.0);
        assert_eq!(complex_intersect_compare_i(
            vector![1.0, 1.0], vector![3.0, 4.0], vector![1.0, 4.0], vector![3.0, 1.0],
            vector![2.0, 0.0], vector![2.0, 3.0], vector![0.0, 2.0], vector![5.0, 2.0], 0), 0.0);
        assert!(complex_intersect_compare_i(
            vector![1.0, 1.0], vector![3.0, 4.0], vector![1.0, 4.0], vector![3.0, 1.0],
            vector![2.0, 0.0], vector![2.0, 3.0], vector![0.0, 2.0], vector![5.0, 2.0], 1) > 0.0);
        assert!(complex_intersect_compare_i(
            vector![1.0, 1.0], vector![3.0, 2.0], vector![1.0, 2.0], vector![3.0, 1.0],
            vector![2.0, 0.0], vector![2.0, 3.0], vector![0.0, 2.0], vector![5.0, 2.0], 1) < 0.0);
        assert!(complex_intersect_compare_i(
            vector![1.0, 1.0], vector![3.0, 3.0], vector![1.0, 3.0], vector![4.0, 1.0],
            vector![2.0, 0.0], vector![2.0, 3.0], vector![0.0, 2.0], vector![5.0, 2.0], 0) > 0.0);
        assert!(complex_intersect_compare_i(
            vector![1.0, 1.0], vector![3.0, 3.0], vector![1.0, 3.0], vector![2.0, 1.0],
            vector![2.0, 0.0], vector![2.0, 3.0], vector![0.0, 2.0], vector![5.0, 2.0], 0) < 0.0);
    }

    fn in_circle_exact(a: Vec2, b: Vec2, c: Vec2, d: Vec2) -> Float {
        const PREC: u32 = ((f64::MANTISSA_DIGITS + 1) * 2 + 1) * 2 + 2;
        macro_rules! f {
            ($expr:expr) => {
                Float::with_val(PREC, $expr)
            };
        }

        let ax = f!(&f!(a.x) - &f!(d.x));
        let ay = f!(&f!(a.y) - &f!(d.y));
        let bx = f!(&f!(b.x) - &f!(d.x));
        let by = f!(&f!(b.y) - &f!(d.y));
        let cx = f!(&f!(c.x) - &f!(d.x));
        let cy = f!(&f!(c.y) - &f!(d.y));
        let cof1 = f!(&f!(&ax * &ax + &ay * &ay) * &f!(&bx * &cy - &by * &cx));
        let cof2 = f!(&f!(&bx * &bx + &by * &by) * &f!(&cx * &ay - &cy * &ax));
        let cof3 = f!(&f!(&cx * &cx + &cy * &cy) * &f!(&ax * &by - &ay * &bx));
        f!(&f!(&cof1 + &cof2) + &cof3)
    }

    fn check_in_circle(a: Vec2, b: Vec2, c: Vec2, d: Vec2) {
        let adapt = in_circle(a, b, c, d);
        let exact = in_circle_exact(a, b, c, d);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}) gave wrong result: {} vs {}",
            a,
            b,
            c,
            d,
            adapt,
            exact
        );
    }

    #[test]
    fn test_in_circle_uniform_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(8).collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            let c = Vec2::new(vals[4], vals[5]);
            let d = Vec2::new(vals[6], vals[7]);
            check_in_circle(a, b, c, d);
        }
    }

    #[test]
    fn test_in_circle_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(12)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            let c = Vec2::new(vals[4], vals[5]);
            let d = Vec2::new(vals[6], vals[7]);
            check_in_circle(a, b, c, d);
        }
    }

    #[test]
    fn test_in_circle_near_cocircular() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-1.0, 1.0).unwrap();

        for _ in 0..1000 {
            let vals = UnitCircle.sample_iter(&mut rng).take(4).collect::<Vec<_>>();
            let radius = dist.sample(&mut rng);
            let x = dist.sample(&mut rng);
            let y = dist.sample(&mut rng);
            let offset = Vec2::new(x, y);

            let a = Vec2::new(vals[0][0], vals[0][1]) * radius + offset;
            let b = Vec2::new(vals[1][0], vals[1][1]) * radius + offset;
            let c = Vec2::new(vals[2][0], vals[2][1]) * radius + offset;
            let d = Vec2::new(vals[3][0], vals[3][1]) * radius + offset;

            check_in_circle(a, b, c, d);
        }
    }

    #[test]
    fn test_in_circle_cocircular() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new_inclusive(-4096, 4096).unwrap();
        let abs = Uniform::new_inclusive(0, 4096).unwrap();
        let mtxs = [
            Mtx2::new(1.0, 0.0, 0.0, 1.0),
            Mtx2::new(-1.0, 0.0, 0.0, 1.0),
            Mtx2::new(0.0, -1.0, 1.0, 0.0),
            Mtx2::new(0.0, 1.0, 1.0, 0.0),
            Mtx2::new(-1.0, 0.0, 0.0, -1.0),
            Mtx2::new(1.0, 0.0, 0.0, -1.0),
            Mtx2::new(0.0, 1.0, -1.0, 0.0),
            Mtx2::new(0.0, -1.0, -1.0, 0.0),
        ];

        for _ in 0..1000 {
            let x = abs.sample(&mut rng) as f64 / 4096.0;
            let y = abs.sample(&mut rng) as f64 / 4096.0;
            let v = Vec2::new(x, y);
            let cx = dist.sample(&mut rng) as f64 / 4096.0;
            let cy = dist.sample(&mut rng) as f64 / 4096.0;
            let offset = Vec2::new(cx, cy);

            // Reminder to use choose_multiple for the benchmark to not mix
            // cases that have identical points and are thus more likely to
            // get stopped at the initial floating-point calculation
            let a = mtxs.choose(&mut rng).unwrap() * v + offset;
            let b = mtxs.choose(&mut rng).unwrap() * v + offset;
            let c = mtxs.choose(&mut rng).unwrap() * v + offset;
            let d = mtxs.choose(&mut rng).unwrap() * v + offset;

            check_in_circle(a, b, c, d);
        }
    }

    fn in_sphere_exact(a: Vec3, b: Vec3, c: Vec3, d: Vec3, e: Vec3) -> Float {
        const PREC: u32 = ((f64::MANTISSA_DIGITS + 1) * 2
            + 1
            + f64::MANTISSA_DIGITS
            + 1
            + (f64::MANTISSA_DIGITS + 1) * 2
            + 2)
            + 5;
        macro_rules! f {
            ($expr:expr) => {
                Float::with_val(PREC, $expr)
            };
        }

        let ax = f!(&f!(a.x) - &f!(e.x));
        let ay = f!(&f!(a.y) - &f!(e.y));
        let az = f!(&f!(a.z) - &f!(e.z));
        let bx = f!(&f!(b.x) - &f!(e.x));
        let by = f!(&f!(b.y) - &f!(e.y));
        let bz = f!(&f!(b.z) - &f!(e.z));
        let cx = f!(&f!(c.x) - &f!(e.x));
        let cy = f!(&f!(c.y) - &f!(e.y));
        let cz = f!(&f!(c.z) - &f!(e.z));
        let dx = f!(&f!(d.x) - &f!(e.x));
        let dy = f!(&f!(d.y) - &f!(e.y));
        let dz = f!(&f!(d.z) - &f!(e.z));
        let aq = f!(&f!(&ax * &ax + &ay * &ay) + &az * &az);
        let bq = f!(&f!(&bx * &bx + &by * &by) + &bz * &bz);
        let cq = f!(&f!(&cx * &cx + &cy * &cy) + &cz * &cz);
        let dq = f!(&f!(&dx * &dx + &dy * &dy) + &dz * &dz);
        let cof1 = f!(&f!(&ax * &by - &bx * &ay) * &f!(&cz * &dq - &dz * &cq));
        let cof2 = f!(&f!(&ax * &cy - &cx * &ay) * &f!(&dz * &bq - &bz * &dq));
        let cof3 = f!(&f!(&ax * &dy - &dx * &ay) * &f!(&bz * &cq - &cz * &bq));
        let cof4 = f!(&f!(&bx * &cy - &cx * &by) * &f!(&az * &dq - &dz * &aq));
        let cof5 = f!(&f!(&bx * &dy - &dx * &by) * &f!(&cz * &aq - &az * &cq));
        let cof6 = f!(&f!(&cx * &dy - &dx * &cy) * &f!(&az * &bq - &bz * &aq));
        f!(f!(&f!(&cof1 + &cof2) + &cof3) + f!(&f!(&cof4 + &cof5) + &cof6))
    }

    fn check_in_sphere(a: Vec3, b: Vec3, c: Vec3, d: Vec3, e: Vec3) {
        let adapt = in_sphere(a, b, c, d, e);
        let exact = in_sphere_exact(a, b, c, d, e);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}, {}) gave wrong result: {} vs {}",
            a,
            b,
            c,
            d,
            e,
            adapt,
            exact
        );
    }

    #[test]
    fn test_in_sphere_uniform_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(15).collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            let d = Vec3::new(vals[9], vals[10], vals[11]);
            let e = Vec3::new(vals[12], vals[13], vals[14]);
            check_in_sphere(a, b, c, d, e);
        }
    }

    #[test]
    fn test_in_sphere_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(15)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            let d = Vec3::new(vals[9], vals[10], vals[11]);
            let e = Vec3::new(vals[12], vals[13], vals[14]);
            check_in_sphere(a, b, c, d, e);
        }
    }

    #[test]
    fn test_in_sphere_near_cospherical() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-1.0, 1.0).unwrap();

        for _ in 0..1000 {
            let vals = UnitSphere.sample_iter(&mut rng).take(5).collect::<Vec<_>>();
            let radius = dist.sample(&mut rng);
            let x = dist.sample(&mut rng);
            let y = dist.sample(&mut rng);
            let z = dist.sample(&mut rng);
            let offset = Vec3::new(x, y, z);

            let a = Vec3::new(vals[0][0], vals[0][1], vals[0][2]) * radius + offset;
            let b = Vec3::new(vals[1][0], vals[1][1], vals[1][2]) * radius + offset;
            let c = Vec3::new(vals[2][0], vals[2][1], vals[2][2]) * radius + offset;
            let d = Vec3::new(vals[3][0], vals[3][1], vals[3][2]) * radius + offset;
            let e = Vec3::new(vals[4][0], vals[4][1], vals[4][2]) * radius + offset;

            check_in_sphere(a, b, c, d, e);
        }
    }

    #[test]
    fn test_in_sphere_cospherical() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new_inclusive(-4096, 4096).unwrap();
        let abs = Uniform::new_inclusive(0, 4096).unwrap();
        let flip = Mtx3::new(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        let rot4 = Mtx3::new(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let rot3 = Mtx3::new(0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0);
        let mut mtxs = vec![];

        for i in 0..2 {
            for j in 0..4 {
                for k in 0..3 {
                    mtxs.push(
                        std::iter::repeat(rot3).take(k).product::<Mtx3>()
                            * std::iter::repeat(rot4).take(j).product::<Mtx3>()
                            * std::iter::repeat(flip).take(i).product::<Mtx3>(),
                    );
                }
            }
        }

        for _ in 0..1000 {
            let x = abs.sample(&mut rng) as f64 / 4096.0;
            let y = abs.sample(&mut rng) as f64 / 4096.0;
            let z = abs.sample(&mut rng) as f64 / 4096.0;
            let v = Vec3::new(x, y, z);
            let cx = dist.sample(&mut rng) as f64 / 4096.0;
            let cy = dist.sample(&mut rng) as f64 / 4096.0;
            let cz = dist.sample(&mut rng) as f64 / 4096.0;
            let offset = Vec3::new(cx, cy, cz);

            // Reminder to use choose_multiple for the benchmark to not mix
            // cases that have identical points and are thus more likely to
            // get stopped at the initial floating-point calculation
            let a = mtxs.choose(&mut rng).unwrap() * v + offset;
            let b = mtxs.choose(&mut rng).unwrap() * v + offset;
            let c = mtxs.choose(&mut rng).unwrap() * v + offset;
            let d = mtxs.choose(&mut rng).unwrap() * v + offset;
            let e = mtxs.choose(&mut rng).unwrap() * v + offset;

            check_in_sphere(a, b, c, d, e);
        }
    }

    fn magnitude_cmp_2d_exact(a: Vec2, b: Vec2) -> Float {
        const PREC: u32 = f64::MANTISSA_DIGITS * 2 + 2;
        macro_rules! f {
            ($expr:expr) => {
                Float::with_val(PREC, $expr)
            };
        }

        let ax = f!(a.x);
        let ay = f!(a.y);
        let bx = f!(b.x);
        let by = f!(b.y);
        f!(&f!(&ax * &ax + &ay * &ay) - &f!(&bx * &bx + &by * &by))
    }

    fn check_magnitude_cmp_2d(a: Vec2, b: Vec2) {
        let adapt = magnitude_cmp_2d(a, b);
        let exact = magnitude_cmp_2d_exact(a, b);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}) gave wrong result: {} vs {}",
            a,
            b,
            adapt,
            exact
        );
    }

    #[test]
    fn test_magnitude_cmp_2d_uniform_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(4).collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            check_magnitude_cmp_2d(a, b);
        }
    }

    #[test]
    fn test_magnitude_cmp_2d_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(4)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            check_magnitude_cmp_2d(a, b);
        }
    }

    #[test]
    fn test_magnitude_cmp_2d_near_zero() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-1.0, 1.0).unwrap();

        for _ in 0..1000 {
            let vals = UnitCircle.sample_iter(&mut rng).take(4).collect::<Vec<_>>();
            let radius = dist.sample(&mut rng);

            let a = Vec2::new(vals[0][0], vals[0][1]) * radius;
            let b = Vec2::new(vals[1][0], vals[1][1]) * radius;

            check_magnitude_cmp_2d(a, b);
        }
    }

    #[test]
    fn test_magnitude_cmp_2d_zero() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let abs = Uniform::new_inclusive(0, 4096).unwrap();
        let mtxs = [
            Mtx2::new(1.0, 0.0, 0.0, 1.0),
            Mtx2::new(-1.0, 0.0, 0.0, 1.0),
            Mtx2::new(0.0, -1.0, 1.0, 0.0),
            Mtx2::new(0.0, 1.0, 1.0, 0.0),
            Mtx2::new(-1.0, 0.0, 0.0, -1.0),
            Mtx2::new(1.0, 0.0, 0.0, -1.0),
            Mtx2::new(0.0, 1.0, -1.0, 0.0),
            Mtx2::new(0.0, -1.0, -1.0, 0.0),
        ];

        for _ in 0..1000 {
            let x = abs.sample(&mut rng) as f64 / 4096.0;
            let y = abs.sample(&mut rng) as f64 / 4096.0;
            let v = Vec2::new(x, y);

            let a = mtxs.choose(&mut rng).unwrap() * v;
            let b = mtxs.choose(&mut rng).unwrap() * v;

            check_magnitude_cmp_2d(a, b);
        }
    }

    fn magnitude_cmp_3d_exact(a: Vec3, b: Vec3) -> Float {
        const PREC: u32 = f64::MANTISSA_DIGITS * 2 + 2;
        macro_rules! f {
            ($expr:expr) => {
                Float::with_val(PREC, $expr)
            };
        }

        let ax = f!(a.x);
        let ay = f!(a.y);
        let az = f!(a.z);
        let bx = f!(b.x);
        let by = f!(b.y);
        let bz = f!(b.z);
        f!(&f!(&f!(&ax * &ax + &ay * &ay) + &az * &az)
            - &f!(&f!(&bx * &bx + &by * &by) + &bz * &bz))
    }

    fn check_magnitude_cmp_3d(a: Vec3, b: Vec3) {
        let adapt = magnitude_cmp_3d(a, b);
        let exact = magnitude_cmp_3d_exact(a, b);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}) gave wrong result: {} vs {}",
            a,
            b,
            adapt,
            exact
        );
    }

    #[test]
    fn test_magnitude_cmp_3d_uniform_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(6).collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            check_magnitude_cmp_3d(a, b);
        }
    }

    #[test]
    fn test_magnitude_cmp_3d_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(6)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            check_magnitude_cmp_3d(a, b);
        }
    }

    #[test]
    fn test_magnitude_cmp_3d_near_zero() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-1.0, 1.0).unwrap();

        for _ in 0..1000 {
            let vals = UnitSphere.sample_iter(&mut rng).take(4).collect::<Vec<_>>();
            let radius = dist.sample(&mut rng);

            let a = Vec3::new(vals[0][0], vals[0][1], vals[0][2]) * radius;
            let b = Vec3::new(vals[1][0], vals[1][1], vals[1][2]) * radius;

            check_magnitude_cmp_3d(a, b);
        }
    }

    #[test]
    fn test_magnitude_cmp_3d_zero() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let abs = Uniform::new_inclusive(0, 4096).unwrap();
        let flip = Mtx3::new(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        let rot4 = Mtx3::new(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let rot3 = Mtx3::new(0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0);
        let mut mtxs = vec![];

        for i in 0..2 {
            for j in 0..4 {
                for k in 0..3 {
                    mtxs.push(
                        std::iter::repeat(rot3).take(k).product::<Mtx3>()
                            * std::iter::repeat(rot4).take(j).product::<Mtx3>()
                            * std::iter::repeat(flip).take(i).product::<Mtx3>(),
                    );
                }
            }
        }

        for _ in 0..1000 {
            let x = abs.sample(&mut rng) as f64 / 4096.0;
            let y = abs.sample(&mut rng) as f64 / 4096.0;
            let z = abs.sample(&mut rng) as f64 / 4096.0;
            let v = Vec3::new(x, y, z);

            let a = mtxs.choose(&mut rng).unwrap() * v;
            let b = mtxs.choose(&mut rng).unwrap() * v;

            check_magnitude_cmp_3d(a, b);
        }
    }

    fn sign_det_x_x2y2_exact(a: Vec2, b: Vec2, c: Vec2) -> Float {
        const PREC: u32 = f64::MANTISSA_DIGITS * 2 + 1 + f64::MANTISSA_DIGITS + 3;
        macro_rules! f {
            ($expr:expr) => {
                Float::with_val(PREC, $expr)
            };
        }

        let ax = f!(a.x);
        let ay = f!(a.y);
        let bx = f!(b.x);
        let by = f!(b.y);
        let cx = f!(c.x);
        let cy = f!(c.y);
        let sqa = f!(&ax * &ax + &ay * &ay);
        let sqb = f!(&bx * &bx + &by * &by);
        let sqc = f!(&cx * &cx + &cy * &cy);
        let cof1 = f!(&sqa * f!(&cx - &bx));
        let cof2 = f!(&sqb * f!(&ax - &cx));
        let cof3 = f!(&sqc * f!(&bx - &ax));
        f!(f!(&cof1 + &cof2) + &cof3)
    }

    fn check_sign_det_x_x2y2(a: Vec2, b: Vec2, c: Vec2) {
        let adapt = sign_det_x_x2y2(a, b, c);
        let exact = sign_det_x_x2y2_exact(a, b, c);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}) gave wrong result: {} vs {}",
            a,
            b,
            c,
            adapt,
            exact
        );

        let adapt = sign_det_x_x2y2_adapt(a, b, c);
        let exact = sign_det_x_x2y2_exact(a, b, c);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}) adapted gave wrong result: {} vs {}",
            a,
            b,
            c,
            adapt,
            exact
        );
    }

    #[test]
    fn test_sign_det_x_x2y2_uniform_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(6).collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            let c = Vec2::new(vals[4], vals[5]);
            check_sign_det_x_x2y2(a, b, c);
        }
    }

    #[test]
    fn test_sign_det_x_x2y2_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(6)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            let c = Vec2::new(vals[4], vals[5]);
            check_sign_det_x_x2y2(a, b, c);
        }
    }

    #[test]
    fn test_sign_det_x_x2y2_zero() {
        let a = Vec2::new(1.0, 0.0);
        let b = Vec2::new(2.0, -1.0);
        let c = Vec2::new(3.0, 0.0);
        assert_eq!(sign_det_x_x2y2(a, b, c), 0.0);
    }

    fn sign_det_x_x2y2z2_exact(a: Vec3, b: Vec3, c: Vec3) -> Float {
        const PREC: u32 = f64::MANTISSA_DIGITS * 2 + 2 + f64::MANTISSA_DIGITS + 3;
        macro_rules! f {
            ($expr:expr) => {
                Float::with_val(PREC, $expr)
            };
        }

        let ax = f!(a.x);
        let ay = f!(a.y);
        let az = f!(a.z);
        let bx = f!(b.x);
        let by = f!(b.y);
        let bz = f!(b.z);
        let cx = f!(c.x);
        let cy = f!(c.y);
        let cz = f!(c.z);
        let sqa = f!(&f!(&ax * &ax + &ay * &ay) + &f!(&az * &az));
        let sqb = f!(&f!(&bx * &bx + &by * &by) + &f!(&bz * &bz));
        let sqc = f!(&f!(&cx * &cx + &cy * &cy) + &f!(&cz * &cz));
        let cof1 = f!(&sqa * f!(&cx - &bx));
        let cof2 = f!(&sqb * f!(&ax - &cx));
        let cof3 = f!(&sqc * f!(&bx - &ax));
        f!(f!(&cof1 + &cof2) + &cof3)
    }

    fn check_sign_det_x_x2y2z2(a: Vec3, b: Vec3, c: Vec3) {
        let adapt = sign_det_x_x2y2z2(a, b, c);
        let exact = sign_det_x_x2y2z2_exact(a, b, c);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}) gave wrong result: {} vs {}",
            a,
            b,
            c,
            adapt,
            exact
        );

        let adapt = sign_det_x_x2y2z2_adapt(a, b, c);
        let exact = sign_det_x_x2y2z2_exact(a, b, c);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}) adapted gave wrong result: {} vs {}",
            a,
            b,
            c,
            adapt,
            exact
        );
    }

    #[test]
    fn test_sign_det_x_x2y2z2_uniform_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(9).collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            check_sign_det_x_x2y2z2(a, b, c);
        }
    }

    #[test]
    fn test_sign_det_x_x2y2z2_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(9)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            check_sign_det_x_x2y2z2(a, b, c);
        }
    }

    #[test]
    fn test_sign_det_x_x2y2z2_zero() {
        let a = Vec3::new(1.0, 0.0, 4.0);
        let b = Vec3::new(2.0, -3.0, 3.0);
        let c = Vec3::new(4.0, -4.0, 0.0);
        assert_eq!(sign_det_x_x2y2z2(a, b, c), 0.0);
    }

    fn sign_det_x_y_x2y2z2_exact(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> Float {
        const PREC: u32 = f64::MANTISSA_DIGITS * 2 + 2 + f64::MANTISSA_DIGITS * 2 + 4;
        macro_rules! f {
            ($expr:expr) => {
                Float::with_val(PREC, $expr)
            };
        }

        let ax = f!(a.x);
        let ay = f!(a.y);
        let az = f!(a.z);
        let bx = f!(b.x);
        let by = f!(b.y);
        let bz = f!(b.z);
        let cx = f!(c.x);
        let cy = f!(c.y);
        let cz = f!(c.z);
        let dx = f!(d.x);
        let dy = f!(d.y);
        let dz = f!(d.z);
        let sqa = f!(&f!(&ax * &ax + &ay * &ay) + &f!(&az * &az));
        let sqb = f!(&f!(&bx * &bx + &by * &by) + &f!(&bz * &bz));
        let sqc = f!(&f!(&cx * &cx + &cy * &cy) + &f!(&cz * &cz));
        let sqd = f!(&f!(&dx * &dx + &dy * &dy) + &f!(&dz * &dz));
        let cof1 = f!(&f!(&ax * &by - &bx * &ay) * &f!(&sqc - &sqd));
        let cof2 = f!(&f!(&ax * &cy - &cx * &ay) * &f!(&sqd - &sqb));
        let cof3 = f!(&f!(&ax * &dy - &dx * &ay) * &f!(&sqb - &sqc));
        let cof4 = f!(&f!(&bx * &cy - &cx * &by) * &f!(&sqa - &sqd));
        let cof5 = f!(&f!(&bx * &dy - &dx * &by) * &f!(&sqc - &sqa));
        let cof6 = f!(&f!(&cx * &dy - &dx * &cy) * &f!(&sqa - &sqb));
        f!(f!(&f!(&cof1 + &cof2) + &cof3) + f!(&f!(&cof4 + &cof5) + &cof6))
    }

    fn check_sign_det_x_y_x2y2z2(a: Vec3, b: Vec3, c: Vec3, d: Vec3) {
        let adapt = sign_det_x_y_x2y2z2(a, b, c, d);
        let exact = sign_det_x_y_x2y2z2_exact(a, b, c, d);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}) gave wrong result: {} vs {}",
            a,
            b,
            c,
            d,
            adapt,
            exact
        );

        let adapt = sign_det_x_y_x2y2z2_adapt(a, b, c, d);
        let exact = sign_det_x_y_x2y2z2_exact(a, b, c, d);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}) adapted gave wrong result: {} vs {}",
            a,
            b,
            c,
            d,
            adapt,
            exact
        );
    }

    #[test]
    fn test_sign_det_x_y_x2y2z2_uniform_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(12).collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            let d = Vec3::new(vals[9], vals[10], vals[11]);
            check_sign_det_x_y_x2y2z2(a, b, c, d);
        }
    }

    #[test]
    fn test_sign_det_x_y_x2y2z2_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(12)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            let d = Vec3::new(vals[9], vals[10], vals[11]);
            check_sign_det_x_y_x2y2z2(a, b, c, d);
        }
    }

    #[test]
    fn test_sign_det_x_y_x2y2z2_zero() {
        let a = Vec3::new(1.0, 0.0, 4.0);
        let b = Vec3::new(2.0, -3.0, 3.0);
        let c = Vec3::new(4.0, -4.0, 0.0);
        let d = Vec3::new(2.0, -3.0, -3.0);
        assert_eq!(sign_det_x_y_x2y2z2(a, b, c, d), 0.0);
    }

    fn sign_det_x_y_x2y2z2_plus_2x_det_x_y_z_exact(a: Vec3, b: Vec3, c: Vec3, d: Vec3, u: f64, i: Vec3, j: Vec3, k: Vec3, l: Vec3) -> Float {
        const PREC: u32 = 1000; // lazy; this should be enough
        macro_rules! f {
            ($expr:expr) => {
                Float::with_val(PREC, $expr)
            };
        }
        f!(&sign_det_x_y_x2y2z2_exact(a, b, c, d) + &f!(2.0 * u) * &orient_3d_exact(i, j, k, l))
    }

    fn check_sign_det_x_y_x2y2z2_plus_2x_det_x_y_z(a: Vec3, b: Vec3, c: Vec3, d: Vec3, u: f64, i: Vec3, j: Vec3, k: Vec3, l: Vec3) {
        let adapt = sign_det_x_y_x2y2z2_plus_2x_det_x_y_z(a, b, c, d, u, i, j, k, l);
        let exact = sign_det_x_y_x2y2z2_plus_2x_det_x_y_z_exact(a, b, c, d, u, i, j, k, l);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}, {}, {}, {}, {}, {}) gave wrong result: {} vs {}",
            a,
            b,
            c,
            d,
            u,
            i,
            j,
            k,
            l,
            adapt,
            exact
        );

        let adapt = sign_det_x_y_x2y2z2_plus_2x_det_x_y_z_adapt(a, b, c, d, u, i, j, k, l);
        let exact = sign_det_x_y_x2y2z2_plus_2x_det_x_y_z_exact(a, b, c, d, u, i, j, k, l);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}, {}, {}, {}, {}, {}) adapted gave wrong result: {} vs {}",
            a,
            b,
            c,
            d,
            u,
            i,
            j,
            k,
            l,
            adapt,
            exact
        );
    }

    #[test]
    fn test_sign_det_x_y_x2y2z2_plus_2x_det_x_y_z_uniform_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(25).collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            let d = Vec3::new(vals[9], vals[10], vals[11]);
            let u = vals[12];
            let i = Vec3::new(vals[13], vals[14], vals[15]);
            let j = Vec3::new(vals[16], vals[17], vals[18]);
            let k = Vec3::new(vals[19], vals[20], vals[21]);
            let l = Vec3::new(vals[22], vals[23], vals[24]);
            check_sign_det_x_y_x2y2z2_plus_2x_det_x_y_z(a, b, c, d, u, i, j, k, l);
        }
    }

    #[test]
    fn test_sign_det_x_y_x2y2z2_plus_2x_det_x_y_z_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(25)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            let d = Vec3::new(vals[9], vals[10], vals[11]);
            let u = vals[12];
            let i = Vec3::new(vals[13], vals[14], vals[15]);
            let j = Vec3::new(vals[16], vals[17], vals[18]);
            let k = Vec3::new(vals[19], vals[20], vals[21]);
            let l = Vec3::new(vals[22], vals[23], vals[24]);
            check_sign_det_x_y_x2y2z2_plus_2x_det_x_y_z(a, b, c, d, u, i, j, k, l);
        }
    }

    fn sign_det_x_x2y2z2_plus_2x_det_x_y_exact(a: Vec3, b: Vec3, c: Vec3, u: f64, i: Vec2, j: Vec2, k: Vec2) -> Float {
        const PREC: u32 = 1000; // lazy; this should be enough
        macro_rules! f {
            ($expr:expr) => {
                Float::with_val(PREC, $expr)
            };
        }
        f!(&sign_det_x_x2y2z2_exact(a, b, c) + &f!(2.0 * u) * &orient_2d_exact(i, j, k))
    }

    fn check_sign_det_x_x2y2z2_plus_2x_det_x_y(a: Vec3, b: Vec3, c: Vec3, u: f64, i: Vec2, j: Vec2, k: Vec2) {
        let adapt = sign_det_x_x2y2z2_plus_2x_det_x_y(a, b, c, u, i, j, k);
        let exact = sign_det_x_x2y2z2_plus_2x_det_x_y_exact(a, b, c, u, i, j, k);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}, {}, {}, {}) gave wrong result: {} vs {}",
            a,
            b,
            c,
            u,
            i,
            j,
            k,
            adapt,
            exact
        );

        let adapt = sign_det_x_x2y2z2_plus_2x_det_x_y_adapt(a, b, c, u, i, j, k);
        let exact = sign_det_x_x2y2z2_plus_2x_det_x_y_exact(a, b, c, u, i, j, k);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}, {}, {}, {}) adapted gave wrong result: {} vs {}",
            a,
            b,
            c,
            u,
            i,
            j,
            k,
            adapt,
            exact
        );
    }

    #[test]
    fn test_sign_det_x_x2y2z2_plus_2x_det_x_y_uniform_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(16).collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            let u = vals[9];
            let i = Vec2::new(vals[10], vals[11]);
            let j = Vec2::new(vals[12], vals[13]);
            let k = Vec2::new(vals[14], vals[15]);
            check_sign_det_x_x2y2z2_plus_2x_det_x_y(a, b, c, u, i, j, k);
        }
    }

    #[test]
    fn test_sign_det_x_x2y2z2_plus_2x_det_x_y_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(16)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let c = Vec3::new(vals[6], vals[7], vals[8]);
            let u = vals[9];
            let i = Vec2::new(vals[10], vals[11]);
            let j = Vec2::new(vals[12], vals[13]);
            let k = Vec2::new(vals[14], vals[15]);
            check_sign_det_x_x2y2z2_plus_2x_det_x_y(a, b, c, u, i, j, k);
        }
    }

    fn sign_det_x2y2z2_plus_2x_det_x_exact(a: Vec3, b: Vec3, u: f64, i: f64, j: f64) -> Float {
        const PREC: u32 = 1000; // lazy; this should be enough
        macro_rules! f {
            ($expr:expr) => {
                Float::with_val(PREC, $expr)
            };
        }
        f!(&magnitude_cmp_3d_exact(a, b) + &f!(2.0 * u) * &f!(&f!(i) - &f!(j)))
    }

    fn check_sign_det_x2y2z2_plus_2x_det_x(a: Vec3, b: Vec3, u: f64, i: f64, j: f64) {
        let adapt = sign_det_x2y2z2_plus_2x_det_x(a, b, u, i, j);
        let exact = sign_det_x2y2z2_plus_2x_det_x_exact(a, b, u, i, j);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}, {}) gave wrong result: {} vs {}",
            a,
            b,
            u,
            i,
            j,
            adapt,
            exact
        );

        let adapt = sign_det_x2y2z2_plus_2x_det_x_adapt(a, b, u, i, j);
        let exact = sign_det_x2y2z2_plus_2x_det_x_exact(a, b, u, i, j);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}, {}) adapted gave wrong result: {} vs {}",
            a,
            b,
            u,
            i,
            j,
            adapt,
            exact
        );
    }

    #[test]
    fn test_sign_det_x2y2z2_plus_2x_det_x_uniform_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(9).collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let u = vals[6];
            let i = vals[7];
            let j = vals[8];
            check_sign_det_x2y2z2_plus_2x_det_x(a, b, u, i, j);
        }
    }

    #[test]
    fn test_sign_det_x2y2z2_plus_2x_det_x_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(16)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec3::new(vals[0], vals[1], vals[2]);
            let b = Vec3::new(vals[3], vals[4], vals[5]);
            let u = vals[6];
            let i = vals[7];
            let j = vals[8];
            check_sign_det_x2y2z2_plus_2x_det_x(a, b, u, i, j);
        }
    }

    fn sign_det_x_x2y2_plus_2x_det_x_y_exact(a: Vec2, b: Vec2, c: Vec2, u: f64, i: Vec2, j: Vec2, k: Vec2) -> Float {
        const PREC: u32 = 1000; // lazy; this should be enough
        macro_rules! f {
            ($expr:expr) => {
                Float::with_val(PREC, $expr)
            };
        }
        f!(&sign_det_x_x2y2_exact(a, b, c) + &f!(2.0 * u) * &orient_2d_exact(i, j, k))
    }

    fn check_sign_det_x_x2y2_plus_2x_det_x_y(a: Vec2, b: Vec2, c: Vec2, u: f64, i: Vec2, j: Vec2, k: Vec2) {
        let adapt = sign_det_x_x2y2_plus_2x_det_x_y(a, b, c, u, i, j, k);
        let exact = sign_det_x_x2y2_plus_2x_det_x_y_exact(a, b, c, u, i, j, k);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}, {}, {}, {}) gave wrong result: {} vs {}",
            a,
            b,
            c,
            u,
            i,
            j,
            k,
            adapt,
            exact
        );

        let adapt = sign_det_x_x2y2_plus_2x_det_x_y_adapt(a, b, c, u, i, j, k);
        let exact = sign_det_x_x2y2_plus_2x_det_x_y_exact(a, b, c, u, i, j, k);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}, {}, {}, {}) adapted gave wrong result: {} vs {}",
            a,
            b,
            c,
            u,
            i,
            j,
            k,
            adapt,
            exact
        );
    }

    #[test]
    fn test_sign_det_x_x2y2_plus_2x_det_x_y_uniform_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(13).collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            let c = Vec2::new(vals[4], vals[5]);
            let u = vals[6];
            let i = Vec2::new(vals[7], vals[8]);
            let j = Vec2::new(vals[9], vals[10]);
            let k = Vec2::new(vals[11], vals[12]);
            check_sign_det_x_x2y2_plus_2x_det_x_y(a, b, c, u, i, j, k);
        }
    }

    #[test]
    fn test_sign_det_x_x2y2_plus_2x_det_x_y_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(16)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            let c = Vec2::new(vals[4], vals[5]);
            let u = vals[6];
            let i = Vec2::new(vals[7], vals[8]);
            let j = Vec2::new(vals[9], vals[10]);
            let k = Vec2::new(vals[11], vals[12]);
            check_sign_det_x_x2y2_plus_2x_det_x_y(a, b, c, u, i, j, k);
        }
    }

    fn sign_det_x2y2_plus_2x_det_x_exact(a: Vec2, b: Vec2, u: f64, i: f64, j: f64) -> Float {
        const PREC: u32 = 1000; // lazy; this should be enough
        macro_rules! f {
            ($expr:expr) => {
                Float::with_val(PREC, $expr)
            };
        }
        f!(&magnitude_cmp_2d_exact(a, b) + &f!(2.0 * u) * &f!(&f!(i) - &f!(j)))
    }

    fn check_sign_det_x2y2_plus_2x_det_x(a: Vec2, b: Vec2, u: f64, i: f64, j: f64) {
        let adapt = sign_det_x2y2_plus_2x_det_x(a, b, u, i, j);
        let exact = sign_det_x2y2_plus_2x_det_x_exact(a, b, u, i, j);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}, {}) gave wrong result: {} vs {}",
            a,
            b,
            u,
            i,
            j,
            adapt,
            exact
        );

        let adapt = sign_det_x2y2_plus_2x_det_x_adapt(a, b, u, i, j);
        let exact = sign_det_x2y2_plus_2x_det_x_exact(a, b, u, i, j);
        assert_eq!(
            adapt.partial_cmp(&0.0),
            exact.partial_cmp(&0.0),
            "({}, {}, {}, {}, {}) adapted gave wrong result: {} vs {}",
            a,
            b,
            u,
            i,
            j,
            adapt,
            exact
        );
    }

    #[test]
    fn test_sign_det_x2y2_plus_2x_det_x_uniform_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-10.0, 10.0).unwrap();

        for _ in 0..1000 {
            let vals = dist.sample_iter(&mut rng).take(7).collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            let u = vals[4];
            let i = vals[5];
            let j = vals[6];
            check_sign_det_x2y2_plus_2x_det_x(a, b, u, i, j);
        }
    }

    #[test]
    fn test_sign_det_x2y2_plus_2x_det_x_geometric_random() {
        let mut rng = Pcg64::new(PCG_STATE, PCG_STREAM);
        let mut rng2 = Pcg64::new(PCG_STATE, PCG_STREAM);
        let dist = Uniform::new(-30.0, 30.0).unwrap();

        for _ in 0..1000 {
            let vals = dist
                .sample_iter(&mut rng)
                .take(16)
                .map(|x: f64| if rng2.random() { -1.0 } else { 1.0 } * x.exp2())
                .collect::<Vec<_>>();
            let a = Vec2::new(vals[0], vals[1]);
            let b = Vec2::new(vals[2], vals[3]);
            let u = vals[4];
            let i = vals[5];
            let j = vals[6];
            check_sign_det_x2y2_plus_2x_det_x(a, b, u, i, j);
        }
    }
}