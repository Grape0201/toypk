/// Build Up Factor
/// currently implemented
/// - GP from
/// - Capo form
/// - Testing form (for testing)
/// 
/// to add new form, implement `BuildUpFactor` trait

const TANHM2: f64 = -0.9640275800758169;  // tanh(-2)

fn gpk(x: f64, c: f64, a: f64, xk: f64, d: f64) -> f64 {
    //! calculate K,
    //! K = c*x^a + d*...
    c * x.powf(a) + d * (
        (x/xk - 2.).tanh() - TANHM2
    ) / (1. - TANHM2)
}

pub enum BuildUpFactor {
    GpForm(f64, f64, f64, f64, f64, f64, f64, f64),
    CapoForm(f64, f64, f64, f64),
    TestingForm(f64),
}

impl BuildUpFactor {
    /// # Arguments
    /// * `x` - the length the photon penetrated, in the unit of mfp
    /// # Returns
    /// build up factor
    fn interpolate(&self, x: f64) -> f64 {
        match self {
            &BuildUpFactor::GpForm(b, c, a, xk, d, k35, k40k35, k401k351) => {
                let k = if x <= 40. {
                    gpk(x, c, a, xk, d)
                } else {
                    let zeta = ((x/35.).powf(0.1) - 1.) / ((x/40.).powf(0.1) - 1.);
                    if k401k351 < 1.0 {
                        1. + (k35 - 1.) * k401k351.powf(zeta)
                    } else {
                        k35 * k40k35.powf(zeta.powf(0.8))
                    }
                };
                let buf = if k == 1.0 {
                    1. + (b-1.)*x
                } else {
                    1. + (b-1.) * (k.powf(x)-1.) / (k - 1.)
                };
                buf
            },
            &BuildUpFactor::TestingForm(bf) => bf,
            &BuildUpFactor::CapoForm(b0, b1, b2, b3) => b0 + b1*x + b2*x*x + b3*x*x*x
        }
    }

    pub fn new_gp(b: f64, c: f64, a: f64, xk: f64, d: f64) -> BuildUpFactor {
        let k35 = gpk(35.0, c, a, xk, d);
        let k40 = gpk(40.0, c, a, xk, d);
        let k40k35 = (k40/k35).abs();
        let k401k351 = ((k40 - 1.) / (k35 - 1.)).abs();
        BuildUpFactor::GpForm(b, c, a, xk, d, k35, k40k35, k401k351)
    }

    pub fn new_capo(b0: f64, b1: f64, b2: f64, b3: f64) -> BuildUpFactor {
        BuildUpFactor::CapoForm(b0, b1, b2, b3)
    }

    pub fn new_testing(b0: f64) -> BuildUpFactor { BuildUpFactor::TestingForm(b0) }
}

// -----------------------------------------------------------------------------
pub struct BuildUpFactorUsed {
    // BuidUpFactor by energy group
    pub d: Vec<BuildUpFactor>
}

impl BuildUpFactorUsed {
    pub fn get_buf(&self, energy_group: usize, x: f64) -> f64 {
        self.d[energy_group].interpolate(x)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bf_tanh2() {
        assert_eq!(TANHM2, (-2. as f64).tanh());
    }

    #[test]
    fn bf_gp_k() {
        let c = 1.0;
        let a = 1.0;
        let xk = 1.0;
        let d = 1.0;

        let x = 0.0;
        assert_eq!(gpk(x, c, a, xk, d), 0.0);
    }

    #[test]
    fn bf_gp_form() {
        // E(MEV)   B      C      A      XK      D
        // 0.015  1.004  1.561 -0.554   5.60    0.3524
        let gpd = BuildUpFactor::new_gp(1.004, 1.561, -0.554, 5.6, 0.3524);
        assert_eq!(gpd.interpolate(0.0), 1.0);
        assert_eq!(gpd.interpolate(1.0), 1.004);
        assert_eq!(gpd.interpolate(2.0), 1.0082789195367863);
        assert_eq!(gpd.interpolate(35.0), 1.0093043284601377);
        assert_eq!(gpd.interpolate(40.0), 1.0089812006220353);
        assert_eq!(gpd.interpolate(1000.0), 1.0089692597254705);

        let gpu = BuildUpFactorUsed { d: vec![gpd]};
        assert_eq!(gpu.get_buf(0, 0.0), 1.0);
        assert_eq!(gpu.get_buf(0, 1.0), 1.004);
        assert_eq!(gpu.get_buf(0, 2.0), 1.0082789195367863);
        assert_eq!(gpu.get_buf(0, 35.0), 1.0093043284601377);
        assert_eq!(gpu.get_buf(0, 40.0), 1.0089812006220353);
        assert_eq!(gpu.get_buf(0, 1000.0), 1.0089692597254705);
    }

    #[test]
    fn bf_capo() {
        let capo = BuildUpFactor::new_capo(1., 1., 1., 1.);
        assert_eq!(capo.interpolate(0.), 1.);
    }
}
