/// Photon source
use crate::csg;

#[derive(Debug)]
pub struct Source {
    /// source volume, in the unit of cm3
    pub volume: f64,
    pub p0: csg::Point,
}

#[derive(Clone)]
pub struct TestSource {
    xmin: f64,
    ymin: f64,
    zmin: f64,
    nx: usize,
    ny: usize,
    nz: usize,
    current_nx: usize,
    current_ny: usize,
    current_nz: usize,
    dx: f64,
    dy: f64,
    dz: f64,
    dv: f64,
}

impl TestSource {
    pub fn new(xmin: f64, xmax: f64, ymin: f64, ymax: f64, zmin: f64, zmax: f64, nx: usize, ny: usize, nz: usize) -> Self {
        let dx = (xmax-xmin) / nx as f64;
        let dy = (ymax-ymin) / ny as f64;
        let dz = (zmax-zmin) / nz as f64;
        let dv = dx*dy*dz;
        TestSource {xmin, ymin, zmin, nx, ny, nz, current_nx: 0, current_ny: 0, current_nz: 0, dx, dy, dz, dv}
    }
}

impl Iterator for TestSource {
    type Item = Source;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_nz == self.nz {
            return None;
        }
        let x = self.dx*(self.current_nx as f64+0.5) + self.xmin;
        let y = self.dy*(self.current_ny as f64+0.5) + self.ymin;
        let z = self.dz*(self.current_nz as f64+0.5) + self.zmin;
        
        self.current_nx += 1;
        if self.current_nx == self.nx {
            self.current_nx = 0;
            self.current_ny += 1;
        }
        if self.current_ny == self.ny {
            self.current_ny = 0;
            self.current_nz += 1;
        }
        Some(Source {volume: self.dv, p0: csg::Point::new(x, y, z)})
    }
}


#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_iter() {
        let ts = TestSource::new(0., 1., 0., 1., 0., 1., 3, 4, 5);
        assert_eq!(ts.into_iter().count(), 3*4*5);
    }

    #[test]
    fn test_by_ref() {
        let factor_by_group = vec![1.0, 2.0];
        let srcs = TestSource::new(0., 1., 0., 1., 0., 1., 4, 4, 4);
        let _result = factor_by_group.iter().enumerate().map(|(_energy_group_index, factor)| -> f64 {
            let csrc = srcs.clone();
            for src in csrc {
                println!("{:?}, factor: {}", src.p0, factor);
            }
            *factor
        }).collect::<Vec<f64>>();
    }
}