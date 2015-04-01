#![deny(missing_docs)]
#![feature(core, std_misc)]

//! A simple and type agnostic quaternion math library designed for reexporting

extern crate vecmath;

use vecmath::Vector3;
use vecmath::Matrix4;
use std::num::{Float, FromPrimitive};

/// Quaternion type alias.
pub type Quaternion<T> = (T, [T; 3]);

/// Constructs identity quaternion.
#[inline(always)]
pub fn id<T: Float + Copy>() -> Quaternion<T> {
    let one = Float::one();
    let zero = Float::zero();
    (one, [zero, zero, zero])
}

/// Adds two quaternions.
#[inline(always)]
pub fn add<T: Float>(
    a: Quaternion<T>,
    b: Quaternion<T>
) -> Quaternion<T> {
    use vecmath::vec3_add as add;
    (a.0 + b.0, add(a.1, b.1))
}

/// Multiplies two quaternions.
#[inline(always)]
pub fn mul<T: Float + Copy>(
    a: Quaternion<T>,
    b: Quaternion<T>
) -> Quaternion<T> {
    use vecmath::vec3_cross as cross;
    use vecmath::vec3_add as add;
    use vecmath::vec3_dot as dot;
    use vecmath::vec3_scale as scale;

    (
        a.0 * b.0 - dot(a.1, b.1),
        add(
            add(scale(b.1, a.0), scale(a.1, b.0)),
            cross(a.1, b.1)
        )
    )
}

/// Takes the quaternion conjugate.
#[inline(always)]
pub fn conj<T: Float>(a: Quaternion<T>) -> Quaternion<T> {
    use vecmath::vec3_neg as neg;

    (a.0, neg(a.1))
}

/// Computes the square length of a quaternion.
#[inline(always)]
pub fn square_len<T: Float>(q: Quaternion<T>) -> T {
    use vecmath::vec3_square_len as square_len;
    q.0 * q.0 + square_len(q.1)
}

/// Computes the length of a quaternion.
#[inline(always)]
pub fn len<T: Float>(q: Quaternion<T>) -> T {
    square_len(q).sqrt()
}

/// Rotate the given vector using the given quaternion
#[inline(always)]
pub fn rotate_vector<T: Float>(q: Quaternion<T>, v: Vector3<T>) -> Vector3<T> {
    let zero = Float::zero();
    let v_as_q : Quaternion<T> = (zero, v);
    let q_conj = conj(q);
    mul(mul(q, v_as_q), q_conj).1
}

/// Construct a quaternion representing the given euler angle rotations (in radians)
#[inline(always)]
pub fn euler_angles<T: Float + FromPrimitive>(x: T, y: T, z: T) -> Quaternion<T> {
    let two: T = FromPrimitive::from_int(2).unwrap();

    let half_x = x / two;
    let half_y = y / two;
    let half_z = z / two;

    let cos_x_2 = half_x.cos();
    let cos_y_2 = half_y.cos();
    let cos_z_2 = half_z.cos();

    let sin_x_2 = half_x.sin();
    let sin_y_2 = half_y.sin();
    let sin_z_2 = half_z.sin();

    (
        cos_x_2 * cos_y_2 * cos_z_2 + sin_x_2 * sin_y_2 * sin_z_2,
        [
            sin_x_2 * cos_y_2 * cos_z_2 + cos_x_2 * sin_y_2 * sin_z_2,
            cos_x_2 * sin_y_2 * cos_z_2 + sin_x_2 * cos_y_2 * sin_z_2,
            cos_x_2 * cos_y_2 * sin_z_2 + sin_x_2 * sin_y_2 * cos_z_2
        ]
    )
}

/// Construct a quaternion for the given angle (in radians)
/// about the given axis.
/// Axis must be a unit vector.
#[inline(always)]
pub fn axis_angle<T: Float + FromPrimitive>(axis: Vector3<T>, angle: T) -> Quaternion<T> {
    use vecmath::vec3_scale as scale;
    let two: T = FromPrimitive::from_int(2).unwrap();
    let half_angle = angle / two;
    (half_angle.cos(), scale(axis, half_angle.sin()))
}


/*
pub fn matrix_to_quaternion_2(m: &Matrix4<f32>) -> Quaternion<f32> {
    let mut q = (0.0, [0.0, 0.0, 0.0]);
    let trace = m[0][0] + m[1][1] + m[2][2];

    // check diagonal
    if trace > 0.0 {

        let s = (trace + 1.0).sqrt();
        q.0 = s * 0.5;

        let t = 0.5 / s;
        //q.1[0] = (m[2][1] - m[1][2]) * t;
        //q.1[1] = (m[0][2] - m[2][0]) * t;
        //q.1[2] = (m[1][0] - m[0][1]) * t;
        //
        q.1[0] = (m[1][2] - m[2][1]) * t;
        q.1[1] = (m[2][0] - m[0][2]) * t;
        q.1[2] = (m[0][1] - m[1][0]) * t;

        println!("!!!!!!!!!!!!!");
    } else {

        let mut i = 0;
        if m[1][1] > m[0][0] { i = 1us; }
        if m[2][2] > m[i][i] { i = 2us; }

        let NEXT = [1us, 2us, 0us];
        let j = NEXT[i];
        let k = NEXT[j];

        let s = (m[i][j] - (m[j][j] + m[k][k])).sqrt() + 1.0;
        q.1[i] = s * 0.5;

        let t = if s != 0.0 { 0.5 / s } else  { s };

        q.0 = (m[k][j] - m[j][k]) * t;
        q.1[j] = (m[j][i] - m[i][j]) * t;
        q.1[k] = (m[k][i] - m[i][k]) * t;
    }

    q
}
*/

///
///*
pub fn matrix_to_quaternion(m: &Matrix4<f32>) -> Quaternion<f32> {
    // Ripped off of cgmath...
    // http://www.cs.ucr.edu/~vbz/resources/quatut.pdf
    let trace = m[0][0] + m[1][1] + m[2][2];
    let half = 0.5f32;
    match () {
        () if trace >= 0.0 => {
            let s = (1.0 + trace).sqrt();
            let w = half * s;
            let s = half / s;
            let x = (m[1][2] - m[2][1]) * s;
            let y = (m[2][0] - m[0][2]) * s;
            let z = (m[0][1] - m[1][0]) * s;
            (w, [x, y, z])
        }
        () if (m[0][0] > m[1][1]) && (m[0][0] > m[2][2]) => {
            let s = (half + (m[0][0] - m[1][1] - m[2][2])).sqrt();
            let w = half * s;
            let s = half / s;
            let x = (m[0][1] - m[1][0]) * s;
            let y = (m[2][0] - m[0][2]) * s;
            let z = (m[1][2] - m[2][1]) * s;
            (w, [x, y, z])
        }
        () if m[1][1] > m[2][2] => {
            let s = (half + (m[1][1] - m[0][0] - m[2][2])).sqrt();
            let w = half * s;
            let s = half / s;
            let x = (m[0][1] - m[1][0]) * s;
            let y = (m[1][2] - m[2][1]) * s;
            let z = (m[2][0] - m[0][2]) * s;
            (w, [x, y, z])
        }
        () => {
            let s = (half + (m[2][2] - m[0][0] - m[1][1])).sqrt();
            let w = half * s;
            let s = half / s;
            let x = (m[2][0] - m[0][2]) * s;
            let y = (m[1][2] - m[2][1]) * s;
            let z = (m[0][1] - m[1][0]) * s;
            (w, [x, y, z])
        }
    }
}
//*/

///
///
/*
pub fn matrix_to_quaternion(m: &Matrix4<f32>) -> Quaternion<f32> {
    // Ripped off of cgmath, but converted for row-major matrices
    // http://www.cs.ucr.edu/~vbz/resources/quatut.pdf
    let trace = m[0][0] + m[1][1] + m[2][2];
    let half = 0.5f32;
    match () {
        () if trace >= 0.0 => {
            let s = (1.0 + trace).sqrt();
            let w = half * s;
            let s = half / s;
            let x = (m[2][1] - m[1][2]) * s;
            let y = (m[0][2] - m[2][0]) * s;
            let z = (m[1][0] - m[0][1]) * s;
            (w, [x, y, z])
        }
        () if (m[0][0] > m[1][1]) && (m[0][0] > m[2][2]) => {
            let s = (half + (m[0][0] - m[1][1] - m[2][2])).sqrt();
            let w = half * s;
            let s = half / s;
            let x = (m[1][0] - m[0][1]) * s;
            let y = (m[0][2] - m[2][0]) * s;
            let z = (m[2][1] - m[1][2]) * s;
            (w, [x, y, z])
        }
        () if m[1][1] > m[2][2] => {
            let s = (half + (m[1][1] - m[0][0] - m[2][2])).sqrt();
            let w = half * s;
            let s = half / s;
            let x = (m[1][0] - m[0][1]) * s;
            let y = (m[2][1] - m[1][2]) * s;
            let z = (m[0][2] - m[2][0]) * s;
            (w, [x, y, z])
        }
        () => {
            let s = (half + (m[2][2] - m[0][0] - m[1][1])).sqrt();
            let w = half * s;
            let s = half / s;
            let x = (m[0][2] - m[2][0]) * s;
            let y = (m[2][1] - m[1][2]) * s;
            let z = (m[1][0] - m[0][1]) * s;
            (w, [x, y, z])
        }
    }
}
*/

///
pub fn quaternion_to_matrix(q: Quaternion<f32>) -> Matrix4<f32> {

    let w = q.0;
    let x = q.1[0];
    let y = q.1[1];
    let z = q.1[2];

    let x2 = x + x;
    let y2 = y + y;
    let z2 = z + z;

    let xx2 = x2 * x;
    let xy2 = x2 * y;
    let xz2 = x2 * z;

    let yy2 = y2 * y;
    let yz2 = y2 * z;
    let zz2 = z2 * z;

    let sy2 = y2 * w;
    let sz2 = z2 * w;
    let sx2 = x2 * w;

    [
        [1.0 - yy2 - zz2, xy2 + sz2, xz2 - sy2, 0.0],
        [xy2 - sz2, 1.0 - xx2 - zz2, yz2 + sx2, 0.0],
        [xz2 + sy2, yz2 - sx2, 1.0 - xx2 - yy2, 0.0],
        [0.0, 0.0,  0.0,  1.0]
    ]

    /*
    [
        [1.0 - yy2 - zz2, xy2 - sz2, xz2 + sy2, 0.0],
        [xy2 + sz2, 1.0 - xx2 - zz2, yz2 - sx2, 0.0],
        [xz2 - sy2, yz2 + sx2, 1.0 - xx2 - yy2, 0.0],
        [0.0, 0.0,  0.0,  1.0]
    ]
    */
}


/// Tests
#[cfg(test)]
mod test {
    use super::*;
    use vecmath::Vector3;
    use std::f32::consts::PI;
    use std::num::Float;

    /// Fudge factor for float equality checks
    static EPSILON: f32 = 0.000001;

    //#[test]
    fn test_axis_angle() {
        use vecmath::vec3_normalized as normalized;
        let axis: Vector3<f32> = [1.0, 1.0, 1.0];
        let q: Quaternion<f32> = axis_angle(
            normalized(axis),
            PI
        );

        // Should be a unit quaternion
        assert!((square_len(q) - 1.0).abs() < EPSILON);
    }

    //#[test]
    fn test_euler_angle() {
        let q: Quaternion<f32> = euler_angles(PI, PI, PI);
        // Should be a unit quaternion
        assert!((square_len(q) - 1.0).abs() < EPSILON);
    }

    //#[test]
    fn test_rotate_vector_axis_angle() {
        let v: Vector3<f32> = [1.0, 1.0, 1.0];
        let q: Quaternion<f32> = axis_angle([0.0, 1.0, 0.0], PI);
        let rotated = rotate_vector(q, v);
        assert!((rotated[0] - -1.0).abs() < EPSILON);
        assert!((rotated[1] - 1.0).abs() < EPSILON);
        assert!((rotated[2] - -1.0).abs() < EPSILON);
    }

    //#[test]
    fn test_rotate_vector_euler_angle() {
        let v: Vector3<f32> = [1.0, 1.0, 1.0];
        let q: Quaternion<f32> = euler_angles(0.0, PI, 0.0);
        let rotated = rotate_vector(q, v);
        assert!((rotated[0] - -1.0).abs() < EPSILON);
        assert!((rotated[1] - 1.0).abs() < EPSILON);
        assert!((rotated[2] - -1.0).abs() < EPSILON);
    }

    /// Rotation on axis parallel to vector direction should have no effect
    //#[test]
    fn test_rotate_vector_axis_angle_same_axis() {
        use vecmath::vec3_normalized as normalized;

        let v: Vector3<f32> = [1.0, 1.0, 1.0];
        let arbitrary_angle = 32.12f32;
        let axis: Vector3<f32> = [1.0, 1.0, 1.0];
        let q: Quaternion<f32> = axis_angle(
            normalized(axis),
            arbitrary_angle
        );
        let rotated = rotate_vector(q, v);

        assert!((rotated[0] - 1.0).abs() < EPSILON);
        assert!((rotated[1] - 1.0).abs() < EPSILON);
        assert!((rotated[2] - 1.0).abs() < EPSILON);
    }

    #[test]
    fn test_matrix_to_quaternion() {

        let rotate_on_x = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, (PI/2.0).cos(), (PI/2.0).sin(), 0.0],
            [0.0, -(PI/2.0).sin(), (PI/2.0).cos(), 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];

        // row-major...
        /*
        let rotate_on_x = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, (PI/2.0).cos(), -(PI/2.0).sin(), 0.0],
            [0.0, (PI/2.0).sin(), (PI/2.0).cos(), 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        */

        let q1: Quaternion<f32> = axis_angle([1.0, 0.0, 0.0], PI/2.0);
        let q2 = matrix_to_quaternion(&rotate_on_x);

        println!("q1: {:?}", q1);
        println!("q2: {:?}", q2);

        let m2 = quaternion_to_matrix(q2);
        println!("m1: {:?}", rotate_on_x);
        println!("m2: {:?}", m2);

    }
}
