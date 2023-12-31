/// Build Up Factor
/// currently implemented
/// - GP from
/// - Capo form
/// - Testing form (for testing)
/// 
/// to add new form, implement `BuildUpFactor` trait

const TANHM2: f64 = -0.9640275800758169;  // tanh(-2)

pub trait BuildUpFactor : Send + Sync {
    /// # Arguments
    /// * `x` - the length the photon penetrated, in the unit of mfp
    /// # Returns
    /// build up factor
    fn interpolate(&self, x: f64) -> f64;
}

// -----------------------------------------------------------------------------
// Testing form
pub struct TestingForm {
    pub bf: f64,
}

impl BuildUpFactor for TestingForm {
    fn interpolate(&self, _x: f64) -> f64 {
        self.bf  // returns value unrelated to mean free path
    }
}

// testing form ends here

// -----------------------------------------------------------------------------
// GP form
#[derive(Debug)]
pub struct GpForm {
    b: f64,
    c: f64,
    a: f64,
    xk: f64,
    d: f64,
    k35: f64,       // K(x=35mfp)
    k40k35: f64,    // |K(x=40mfp)/K(x=35mfp)|
    k401k351: f64,  // |K(x=40mfp)-1| / |K(x=35mfp)-1|
}


fn gpk(x: f64, c: f64, a: f64, xk: f64, d: f64) -> f64 {
    //! calculate K,
    //! K = c*x^a + d*...
    c * x.powf(a) + d * (
        (x/xk - 2.).tanh() - TANHM2
    ) / (1. - TANHM2)
}

impl GpForm {
    pub fn new(b: f64, c: f64, a: f64, xk: f64, d: f64) -> GpForm {
        let k35 = gpk(35.0, c, a, xk, d);
        let k40 = gpk(40.0, c, a, xk, d);
        let k40k35 = (k40/k35).abs();
        let k401k351 = ((k40 - 1.) / (k35 - 1.)).abs();
        GpForm { b, c, a, xk, d, k35, k40k35, k401k351 }
    }
}

impl BuildUpFactor for GpForm {
    fn interpolate(&self, x: f64) -> f64 {
        //! does not check if x >= 0
        let k = if x <= 40. {
            gpk(x, self.c, self.a, self.xk, self.d)
        } else {
            let zeta = ((x/35.).powf(0.1) - 1.) / ((x/40.).powf(0.1) - 1.);
            if self.k401k351 < 1.0 {
                1. + (self.k35 - 1.) * self.k401k351.powf(zeta)
            } else {
                self.k35 * self.k40k35.powf(zeta.powf(0.8))
            }
        };
        let buf = if k == 1.0 {
            1. + (self.b-1.)*x
        } else {
            1. + (self.b-1.) * (k.powf(x)-1.) / (k - 1.)
        };
        buf
    }
}

// GP form ends here

// -----------------------------------------------------------------------------
// Capo form
pub struct CapoForm {
    b: [f64; 4]  // Sum_{j=0}^{4} CijE^{-J}
}

impl CapoForm {
    pub fn new(b0: f64, b1: f64, b2: f64, b3: f64) -> Self {
        CapoForm {b: [b0, b1, b2, b3]}
    }
}

impl BuildUpFactor for CapoForm {
    fn interpolate(&self, x: f64) -> f64 {
        (0..4).map(|i| self.b[i]*x.powi(i as i32)).sum()
    }
}
// Capo form ends here

// -----------------------------------------------------------------------------
pub struct BuildUpFactorUsed {
    // BuidUpFactor by energy group
    pub d: Vec<Box<dyn BuildUpFactor>>
}

impl BuildUpFactorUsed {
    pub fn get_buf(&self, energy_group: usize, x: f64) -> f64 {
        let d = &self.d[energy_group];
        d.interpolate(x)
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
        let gpd = GpForm::new(1.004, 1.561, -0.554, 5.6, 0.3524);
        assert_eq!(gpd.interpolate(0.0), 1.0);
        assert_eq!(gpd.interpolate(1.0), 1.004);
        assert_eq!(gpd.interpolate(2.0), 1.0082789195367863);
        assert_eq!(gpd.interpolate(35.0), 1.0093043284601377);
        assert_eq!(gpd.interpolate(40.0), 1.0089812006220353);
        assert_eq!(gpd.interpolate(1000.0), 1.0089692597254705);

        let gpu = BuildUpFactorUsed { d: vec![Box::new(gpd)]};
        assert_eq!(gpu.get_buf(0, 0.0), 1.0);
        assert_eq!(gpu.get_buf(0, 1.0), 1.004);
        assert_eq!(gpu.get_buf(0, 2.0), 1.0082789195367863);
        assert_eq!(gpu.get_buf(0, 35.0), 1.0093043284601377);
        assert_eq!(gpu.get_buf(0, 40.0), 1.0089812006220353);
        assert_eq!(gpu.get_buf(0, 1000.0), 1.0089692597254705);
    }

    #[test]
    fn bf_capo() {
        let capo = CapoForm::new(1., 1., 1., 1.);
        assert_eq!(capo.interpolate(0.), 1.);
    }
}
