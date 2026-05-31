
/// Type alias for 3D vectors used throughout the crate.
pub type Vec3 = nalgebra::Vector3<f64>;

/// Compute the unit vector from `from` to `to`.
pub fn unit_vector(from: Vec3, to: Vec3) -> Vec3 {
    let diff = to - from;
    let norm = diff.norm();
    if norm < f64::EPSILON {
        Vec3::zeros()
    } else {
        diff.normalize()
    }
}

/// Compute the distance between two points.
pub fn distance(a: Vec3, b: Vec3) -> f64 {
    (b - a).norm()
}

/// Cross product convenience wrapper.
pub fn cross(a: &Vec3, b: &Vec3) -> Vec3 {
    a.cross(b)
}

/// Dot product convenience wrapper.
pub fn dot(a: &Vec3, b: &Vec3) -> f64 {
    a.dot(b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::vector;

    #[test]
    fn test_distance() {
        let a = vector![0.0, 0.0, 0.0];
        let b = vector![3.0, 4.0, 0.0];
        assert!((distance(a, b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_unit_vector() {
        let a = vector![0.0, 0.0, 0.0];
        let b = vector![3.0, 0.0, 4.0];
        let u = unit_vector(a, b);
        assert!((u.norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cross_product() {
        let a = vector![1.0, 0.0, 0.0];
        let b = vector![0.0, 1.0, 0.0];
        let c = cross(&a, &b);
        assert!(c.x.abs() < 1e-10);
        assert!(c.y.abs() < 1e-10);
        assert!((c.z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_zero_distance() {
        let a = vector![1.0, 2.0, 3.0];
        assert!(distance(a, a).abs() < 1e-10);
    }
}
