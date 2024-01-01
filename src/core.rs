/// core.rs
use std::collections::HashSet;
use crate::csg;
use crate::buildup;
use crate::source;

pub struct Cell {
    pub name: String,
    pub material_index: usize,
    pub csg_operations: HashSet::<i64>,
}

pub struct Material {
    pub name: String,
    pub linear_attenuation_coefficient: Vec<f64>,
}

pub struct Geometry {
    pub primitives: Vec<Box<dyn csg::Primitive>>,
    pub cells: Vec<Cell>,
    pub materials: Vec<Material>,
}

impl Geometry {
    /// This method will remove duplicated primitives
    fn dedup_primitives(&self) {
        todo!("")
    }

    fn intersections(&self, ls: &csg::LineSegment) -> Vec<f64> {
        let mut ks: Vec<f64> = vec![0.0, 1.0];
        for s in self.primitives.iter() {
            ks.append(&mut s.intersection(&ls));
        }
        ks.sort_by(|a, b| a.partial_cmp(b).expect("NaN in ks"));
        ks.dedup();
        ks
    }

    fn attenuation(&self, ls: &csg::LineSegment, energy_group: usize) -> f64 {
        let mut atten = 0.;
        let ray_length = ls.length();

        let ks = self.intersections(&ls);
        // println!("ks: {:?}", ks);
        for kindex in 0..ks.len()-1 {
            // get the midpoint
            let kcenter = (ks[kindex]+ks[kindex+1])/2.;
            let midpoint = ls.k(kcenter);
            
            let mut cell_ids = HashSet::new();
            for (i, pr) in self.primitives.iter().enumerate() {
                let key = (i+1) as i64;
                let key = if pr.has(&midpoint) {-key} else {key};
                cell_ids.insert(key);
            }
            // println!("cell_ids for k={}: {:?}", kcenter, cell_ids);
            for c in self.cells.iter() {
                let count = c.csg_operations.difference(&cell_ids).count();
                if count != 0 {
                    continue;
                }
                let lac = self.materials[c.material_index].linear_attenuation_coefficient[energy_group];
                let flight_length = (ks[kindex+1] - ks[kindex])*ray_length;
                atten -= flight_length*lac;
                // println!("atten: {}, lac: {}", flight_length*lac, lac);
                break;
            }
        }
        atten
    }
}

pub fn run(geom: &Geometry, srcs: &(impl Iterator<Item = source::Source> + Clone), tally_point: &csg::Point, factor_by_group: &Vec<f64>, bu: &buildup::BuildUpFactorUsed) -> Vec<(f64, f64)> {
    factor_by_group.iter().enumerate().map(|(energy_group_index, factor)| {
        if *factor == 0. {
            return (0., 0.)
        }
        let cloned_srcs = srcs.clone();
        let mut total_k0 = 0.;
        let mut total_k1 = 0.;
        for src in cloned_srcs {
            let ray = csg::LineSegment { p: src.p0, v: *tally_point - src.p0};
            let atten = geom.attenuation(&ray, energy_group_index);
            let k0 = src.volume * factor * atten.exp() / (4.*std::f64::consts::PI * ray.length().powf(2.));
            let buf = bu.get_buf(energy_group_index, -atten);
            total_k0 += k0;
            total_k1 += k0*buf;
        }
        (total_k0, total_k1)
    }).collect()
}



#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_intersection1() {
        let g = Geometry {
            primitives: vec![
                Box::new(csg::Sph::new(0., 0., 0., 1.)),
                Box::new(csg::Sph::new(0., 0., 0., 2.)),
            ],
            cells: vec![
                Cell {
                    name: String::from("inside"),
                    material_index: 0,
                    csg_operations: vec![-1].into_iter().collect(),
                },
                Cell {
                    name: String::from("middle"),
                    material_index: 1,
                    csg_operations: vec![1, -2].into_iter().collect(),
                },
                Cell {
                    name: String::from("outside"),
                    material_index: 2,
                    csg_operations: vec![2].into_iter().collect(),
                },
            ],
            materials: vec![
                Material {
                    name: String::from("fuel"),
                    linear_attenuation_coefficient: vec![0.1],
                },
                Material {
                    name: String::from("iron"),
                    linear_attenuation_coefficient: vec![0.3],
                },
                Material {
                    name: String::from("air"),
                    linear_attenuation_coefficient: vec![0.01],
                },
            ]
        };
        assert_eq!(g.intersections(&csg::LineSegment { p: csg::Point::new(0., 0., 0.), v: csg::Point::new(2., 0., 0.) }), vec![0., 0.5, 1.]);
        assert_eq!(g.intersections(&csg::LineSegment { p: csg::Point::new(0., 0., 0.), v: csg::Point::new(4., 0., 0.) }), vec![0., 0.25, 0.5, 1.]);
        assert_eq!(g.intersections(&csg::LineSegment { p: csg::Point::new(0., 0., 0.), v: csg::Point::new(0., 4., 0.) }), vec![0., 0.25, 0.5, 1.]);
        assert_eq!(g.intersections(&csg::LineSegment { p: csg::Point::new(0., 0., 0.), v: csg::Point::new(0., 0., -4.)}), vec![0., 0.25, 0.5, 1.]);
        assert_eq!(
            g.attenuation(&csg::LineSegment { p: csg::Point::new(0., 0., 0.), v: csg::Point::new(0., 0., -4.)}, 0),
            -0.1*1.0 - 0.3*1.0 - 0.01*2.0
        );
        assert_eq!(
            g.attenuation(&csg::LineSegment { p: csg::Point::new(0., 0., 0.), v: csg::Point::new(0., 0., 4.)}, 0),
            -0.1*1.0 - 0.3*1.0 - 0.01*2.0
        );
        assert_eq!(
            g.attenuation(&csg::LineSegment { p: csg::Point::new(0., 0., 0.), v: csg::Point::new(0., 0., 1.5)}, 0),
            -0.1*1.0 - 0.3*0.5
        );
        assert_eq!(
            g.attenuation(&csg::LineSegment { p: csg::Point::new(3., 0., 0.), v: csg::Point::new(1., 0., 0.)}, 0),
            -0.01*1.0
        );
        assert_eq!(
            g.attenuation(&csg::LineSegment { p: csg::Point::new(3., 0., 0.), v: csg::Point::new(-6., 0., 0.)}, 0),
            -0.01*1.0 - 0.3*1.0 - 0.1*2.0 - 0.3*1.0 - 0.01*1.0
        );
    }


    #[test]
    fn test_run() {
        let g = Geometry {
            primitives: vec![
                Box::new(csg::Sph::new(0., 0., 0., 1.)),
                Box::new(csg::Sph::new(0., 0., 0., 2.)),
            ],
            cells: vec![
                Cell {
                    name: String::from("inside"),
                    material_index: 0,
                    csg_operations: vec![-1].into_iter().collect(),
                },
                Cell {
                    name: String::from("middle"),
                    material_index: 1,
                    csg_operations: vec![1, -2].into_iter().collect(),
                },
                Cell {
                    name: String::from("outside"),
                    material_index: 2,
                    csg_operations: vec![2].into_iter().collect(),
                },
            ],
            materials: vec![
                Material {
                    name: String::from("fuel"),
                    linear_attenuation_coefficient: vec![0.1, 0.1],
                },
                Material {
                    name: String::from("iron"),
                    linear_attenuation_coefficient: vec![0.3, 0.3],
                },
                Material {
                    name: String::from("air"),
                    linear_attenuation_coefficient: vec![0.01, 0.01],
                },
            ]
        };
        let mut source = source::TestSource::new(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, 1, 1, 1);
        let p = csg::Point::new(2., 0., 0.);
        let bf = 2.0;
        let bu = buildup::BuildUpFactorUsed {
            d: vec![
                buildup::BuildUpFactor::new_testing(bf),
                buildup::BuildUpFactor::new_testing(bf),
            ]
        };
        let factor_by_group = vec![1.0, 2.0];
        let ans = (-0.1_f64*1. - 0.3*1.).exp() / (4.*std::f64::consts::PI * 2.*2.);
        assert_eq!(
            run(&g, &mut source, &p, &factor_by_group, &bu),
            factor_by_group.into_iter().map(|factor| (factor*ans, factor*ans*bf)).collect::<Vec<(f64, f64)>>()
        );
    }
}
