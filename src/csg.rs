/// Constructive Solid Geometry and some math
use std::convert::From;

const EPSILON: f64 = 1e-10; // [cm], verrrry small number for boundaries

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Point {
    x: f64,
    y: f64,
    z: f64
}

impl Point {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self {x, y, z}
    }
    fn norm2(&self) -> f64 {
        self.x*self.x + self.y*self.y + self.z*self.z
    }
    fn norm(&self) -> f64 {
        self.norm2().sqrt()
    }
    fn dot(&self, other: &Point) -> f64 {
        self.x*other.x + self.y*other.y + self.z*other.z
    }
    fn cross(&self, other: &Point) -> Self {
        Point {
            x: self.y*other.z - self.z*other.y,
            y: self.z*other.x - self.x*other.z,
            z: self.x*other.y - self.y*other.x
        }
    }
}

impl std::ops::Add for Point {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Point {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z
        }
    }
}

impl std::ops::Sub for Point {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Point {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z
        }
    }
}

impl std::ops::Neg for Point {
    type Output = Self;
    fn neg(self) -> Self {
        Point {
            x: - self.x,
            y: - self.y,
            z: - self.z,
        }
    }
}

impl std::ops::Mul<f64> for Point {
    type Output = Self;
    fn mul(self, k: f64) -> Self {
        Self {
            x: k*self.x,
            y: k*self.y,
            z: k*self.z,
        }
    }
}

impl std::ops::Mul<Point> for f64 {
    type Output = Point;
    fn mul(self, p: Point) -> Point {
        Point {
            x: self*p.x,
            y: self*p.y,
            z: self*p.z,
        }
    }
}

#[derive(Debug)]
pub struct LineSegment {
    pub p: Point,
    pub v: Point,
}

impl LineSegment {
    pub fn new(px: f64, py: f64, pz: f64, vx: f64, vy: f64, vz: f64) -> Self {
        Self {
            p: Point {x: px, y: py, z: pz},
            v: Point {x: vx, y: vy, z: vz},
        }
    }
    pub fn k(&self, k: f64) -> Point {
        self.p + self.v*k
    }
    pub fn length(&self) -> f64 {
        self.v.norm()
    }
}

/// infinte plane
#[derive(Debug)]
struct Plane {
    p0: Point,
    v1: Point,
    v2: Point,
}

impl Plane {
    fn intersection(&self, ls: &LineSegment) -> Option<f64> {
        let a = self.v1.x;
        let b = self.v2.x;
        let c = -ls.v.x;
        let d = self.v1.y;
        let e = self.v2.y;
        let f = -ls.v.y;
        let g = self.v1.z;
        let h = self.v2.z;
        let i = -ls.v.z;
        let j = ls.p.x - self.p0.x;
        let k = ls.p.y - self.p0.y;
        let l = ls.p.z - self.p0.z;
        let upper = j*(d*h-e*g) - k*(a*h-b*g) + l*(a*e-b*d);
        let lower = a*(e*i-f*h) - b*(d*i-f*g) + c*(d*h-e*g);
        if lower == 0.0 {
            None
        } else {
            let k = upper/lower;
            if 0.0 <= k && k <= 1.0 {
                Some(k)
            } else {
                None
            }
        }
    }
}

/// infinte cylinder
#[derive(Debug)]
struct Cylinder {
    p0: Point,
    axis: Point,
    radius: f64,
}

impl Cylinder {
    fn intersection(&self, ls: &LineSegment) -> Vec<f64>{
        // ( (P+kV) - (P0+mAxis) ).dot(Axis) = 0   [1]
        // | (P+kV) - (P0+mAxis) | = radius        [2]
        let pp0 = ls.p - self.p0;
        // m = m0 + m1*k                           [1']
        let m0 = pp0.dot(&self.axis)/self.axis.norm2();
        let m1 = ls.v.dot(&self.axis)/self.axis.norm2();
        // | k(V-m1Axis) + P-P0-m0Axis | = radius  [2']
        // | kX + Y | = radius                     [2'']
        let x = ls.v - self.axis*m1;
        let y = pp0 - self.axis* m0;
        // ----------------------------------
        // X^2*k^2 + 2XY*k + Y^2-radius^2 = 0      [2''']
        let a = x.norm2();
        let b = x.dot(&y)*2.0;
        let c = y.norm2() - self.radius*self.radius;
        let d = b*b - 4.0*a*c;

        if a == 0.0 {
            vec![]
        } else if d.abs() < EPSILON {
            let k0 = -b/2.0/a;
            if 0. <= k0 && k0 <= 1. {
                vec![k0]
            } else{
                vec![]
            }
        } else if d > 0.0 {
            let k0 = (-b+d.powf(0.5))/2.0/a;
            let k1 = (-b-d.powf(0.5))/2.0/a;
            if 0. <= k0 && k0 <= 1.0 && 0. <= k1 && k1 <= 1. {
                vec![k0, k1]
            } else if 0. <= k0 && k0 <= 1.0 {
                vec![k0]
            } else if 0. <= k1 && k1 <= 1. {
                vec![k1]
            } else {
                vec![]
            }
        } else {
            vec![]
        }
    }
}


pub trait Primitive : Send + Sync {
    fn has(&self, _v: &Point) -> bool { true }
    fn intersection(&self, ls: &LineSegment) -> Vec<f64>;
}

#[derive(Debug)]
pub struct Rpp {
    dps: [Plane; 6],
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    zmin: f64,
    zmax: f64,
}

impl Rpp {
    pub fn new(xmin: f64, xmax: f64, ymin: f64, ymax: f64, zmin: f64, zmax: f64) -> Rpp {
        Rpp {
            dps: [
                Plane {p0: Point::new(xmin, ymin, zmin), v1: Point::new(1.0, 0.0, 0.0,),  v2: Point::new(0.0, 1.0, 0.0)},
                Plane {p0: Point::new(xmin, ymin, zmin), v1: Point::new(0.0, 1.0, 0.0,),  v2: Point::new(0.0, 0.0, 1.0)},
                Plane {p0: Point::new(xmin, ymin, zmin), v1: Point::new(0.0, 0.0, 1.0,),  v2: Point::new(1.0, 0.0, 0.0)},
                Plane {p0: Point::new(xmax, ymax, zmax), v1: Point::new(-1.0, 0.0, 0.0,), v2: Point::new(0.0, -1.0, 0.0)},
                Plane {p0: Point::new(xmax, ymax, zmax), v1: Point::new(0.0, -1.0, 0.0,), v2: Point::new(0.0, 0.0, -1.0)},
                Plane {p0: Point::new(xmax, ymax, zmax), v1: Point::new(0.0, 0.0, -1.0,), v2: Point::new(-1.0, 0.0, 0.0)},
            ],
            xmin, xmax, ymin, ymax, zmin, zmax
        }
    }
}

impl Primitive for Rpp {
    fn has(&self, p: &Point) -> bool {
        self.xmin <= p.x+EPSILON &&
            p.x <= self.xmax+EPSILON &&
            self.ymin <= p.y+EPSILON &&
            p.y <= self.ymax+EPSILON &&
            self.zmin <= p.z+EPSILON &&
            p.z <= self.zmax+EPSILON
    }
    fn intersection(&self, ls: &LineSegment) -> Vec<f64>{
        let mut ks = vec![];
        for dp in self.dps.iter() {
            if let Some(k) = dp.intersection(&ls) {
                ks.push(k);
            }
        }
        ks
    }
}

impl From<&Vec<f64>> for Rpp {
    fn from(v: &Vec<f64>) -> Self {
        Rpp::new(v[0], v[1], v[2], v[3], v[4], v[5])
    }
}


#[derive(Debug)]
pub struct Rcc {
    dp0: Plane,
    dp1: Plane,
    dcyl: Cylinder
}

impl Rcc {
    pub fn new(x0: f64, y0: f64, z0: f64, ax: f64, ay: f64, az: f64, radius:f64) -> Self {
        Self {
            dp0: Plane {
                p0: Point::new(x0, y0, z0),
                v1: if ax == 0.0 && ay == 0.0 { Point::new(1.0, 0.0, 0.0) } else { Point::new(-ay, ax, 0.0)},
                v2: if ax == 0.0 && ay == 0.0 { Point::new(0.0, 1.0, 0.0) } else { Point::new(ay, -ax, 0.0)},
            },
            dp1: Plane {
                p0: Point::new(x0+ax, y0+ay, z0+az),
                v1: if ax == 0.0 && ay == 0.0 { Point::new(1.0, 0.0, 0.0) } else { Point::new(-ay, ax, 0.0)},
                v2: if ax == 0.0 && ay == 0.0 { Point::new(0.0, 1.0, 0.0) } else { Point::new(ay, -ax, 0.0)},
            },
            dcyl: Cylinder {
                p0: Point::new(x0, y0, z0),
                axis: Point::new(ax, ay, az),
                radius
            }
        }
    }
}

impl Primitive for Rcc {
    fn has(&self, p: &Point) -> bool {
        let axislen = self.dcyl.axis.norm();
        let sin = (*p-self.dcyl.p0).cross(&self.dcyl.axis).norm() / axislen;
        if sin > self.dcyl.radius {
            return false;
        }
        let cos = (*p-self.dcyl.p0).dot(&self.dcyl.axis) / axislen;
        -EPSILON <= cos && cos <= self.dcyl.axis.norm()+EPSILON
    }
    fn intersection(&self, ls: &LineSegment) -> Vec<f64>{
        let mut ks = vec![];
        if let Some(k) = self.dp0.intersection(&ls) {
            if self.has(&ls.k(k)) {
                ks.push(k);
            }
        }
        if let Some(k) = self.dp1.intersection(&ls) {
            if self.has(&ls.k(k)) {
                ks.push(k);
            }
        }
        for k in self.dcyl.intersection(&ls) {
            if self.has(&ls.k(k)) {
                ks.push(k);
            }
        }
        ks
    }
}

impl From<&Vec<f64>> for Rcc {
    fn from(v: &Vec<f64>) -> Self {
        Rcc::new(v[0], v[1], v[2], v[3], v[4], v[5], v[6])
    }
}

pub struct Sph {
    p0: Point,
    radius: f64
}

impl Sph {
    pub fn new(p0x: f64, p0y: f64, p0z: f64, radius: f64) -> Sph {
        Sph {
            p0: Point::new(p0x, p0y, p0z),
            radius
        }
    }
}

impl Primitive for Sph {
    fn has(&self, p: &Point) -> bool {
        (*p - self.p0).norm() <= self.radius
    }
    fn intersection(&self, ls: &LineSegment) -> Vec<f64> {
        // |P + kV - P0| = r
        // k**2*|V|**2 + k*2V.dot(P-P0) + (|P-P0|**2 - r**2) = 0
        let pp0 = ls.p - self.p0;
        let a = ls.v.norm2();
        let b = ls.v.dot(&pp0)*2.;
        let c = pp0.norm2() - self.radius*self.radius;
        let d = b*b - 4.*a*c;

        if a == 0. {
            println!("a was 0.");
            vec![]
        } else if d == 0. {
            let k = -b/2./a;
            if 0. < k && k < 1. {
                vec![k]
            } else {
                vec![]
            }
        } else if d > 0. {
            let mut ks = vec![];
            let k0 = (-b-d.sqrt())/2./a;
            let k1 = (-b+d.sqrt())/2./a;
            if 0. < k0 && k0 < 1. {ks.push(k0)}
            if 0. < k1 && k1 < 1. {ks.push(k1)}
            ks
        } else {
            vec![]
        }
    }
}

impl From<&Vec<f64>> for Sph {
    fn from(v: &Vec<f64>) -> Self {
        Sph::new(v[0], v[1], v[2], v[3])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn csg_points() {
        let p = Point::new(1., 2., 2.);
        assert_eq!(p.norm2(), 9.0);
        assert_eq!(p.norm(), 3.0);
        let pk = 2.*p;
        assert_eq!(pk, Point::new(2., 4., 4.));
        assert_eq!(p + pk, Point::new(3., 6., 6.));
        assert_eq!(p - pk, -p);
    }

    #[test]
    fn csg_line_segment() {
        let ls = LineSegment::new(0., 0., 0., 2., 0., 0.);
        assert_eq!(ls.length(), 2.);
    }

    #[test]
    fn csg_plane_intersection() {
        let p = Plane {  // XY plane
            p0: Point::new(0., 0., 0.),
            v1: Point::new(1., 0., 0.),
            v2: Point::new(0., 1., 0.),
        };
        assert_eq!(p.intersection(&LineSegment {
            p: Point::new(-1., -1., -1.),
            v: Point::new(0., 0., 1.),
        }), Some(1.0));
        assert_eq!(p.intersection(&LineSegment {
            p: Point::new(-1., -1., -1.),
            v: Point::new(0., 0., 2.),
        }), Some(0.5));
        assert_eq!(p.intersection(&LineSegment {
            p: Point::new(-1., -1., -1.),
            v: Point::new(0., 0., -1.),
        }), None);
    }

    #[test]
    fn csg_cylinder_intersection() {
        let cyl = Cylinder {
            p0: Point::new(0., 0., 0.),
            axis: Point::new(0., 0., 1.),
            radius: 5.
        };
        let line0 = LineSegment {
            p: Point::new(-10., 0., 0.5),
            v: Point::new(20., 0., 0.),
        };
        let ks = cyl.intersection(&line0);
        assert_eq!(&ks, &[0.75, 0.25]);
    }

    #[test]
    fn csg_primitive_rpp() {
        let rpp = Rpp::new(-1., 1., -2., 3., -5., 10.);
        assert_eq!(rpp.intersection(&LineSegment {
            p: Point::new(-2., 0., 0.),
            v: Point::new(2., 0., 0.)
        }), vec![0.5]);
        assert_eq!(rpp.intersection(&LineSegment {
            p: Point::new(-1., 0., 0.),
            v: Point::new(1., 0., 0.)
        }), vec![0.]);
        assert_eq!(rpp.intersection(&LineSegment {
            p: Point::new(0., 4., 0.),
            v: Point::new(0., -10., 0.)
        }), vec![0.6, 0.1]);
        assert_eq!(rpp.intersection(&LineSegment {
            p: Point::new(0., 4., 0.),
            v: Point::new(0., 10., 0.)
        }), vec![]);
    }

    #[test]
    fn csg_primitive_rcc() {
        let rcc = Rcc::new(0., 0., 0., 0., 0., 3., 2.);
        assert_eq!(rcc.intersection(&LineSegment {
            p: Point::new(-5., 0., 1.),
            v: Point::new(10., 0., 0.)
        }), vec![0.7, 0.3]);
        assert_eq!(rcc.intersection(&LineSegment {
            p: Point::new(0., 0., 4.),
            v: Point::new(0., 0., -10.)
        }), vec![0.4, 0.1]);
    }

    #[test]
    fn csg_primitive_sph() {
        let sph = Sph::new(0., 0., 0., 1.0);
        let line = LineSegment { p: Point::new(-2., 0., 0.,), v: Point::new(4., 0., 0.)};
        assert_eq!(sph.intersection(&line), &[0.25, 0.75]);
        let line = LineSegment { p: Point::new(-2., 0., 0.,), v: Point::new(-4., 0., 0.)};
        assert_eq!(sph.intersection(&line), &[]);
        let line = LineSegment { p: Point::new(-2., 0., 2.,), v: Point::new(4., 0., 0.)};
        assert_eq!(sph.intersection(&line), &[]);
        let line = LineSegment { p: Point::new(-2., 0., 1.,), v: Point::new(4., 0., 0.)};
        assert_eq!(sph.intersection(&line), &[0.5]);
    }
}