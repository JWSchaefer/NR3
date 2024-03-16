use crate::interp::interpolator::Interpolate1D;
use crate::table::bisect_hunt::BisectHunt1D;
use crate::table::search::Search;
use ndarray::prelude::*;
use plotters::prelude::*;

/// 1 Dimensional Cubic Spine Interpolator
pub struct Spline1D {
    x: Array1<f64>,       // x   table
    y: Array1<f64>,       // y   table
    y2: Array1<f64>,      // y'' table
    search: BisectHunt1D, // Search Algorithm
}

impl Spline1D {
    /// Constructor
    /// #  Arguments
    /// * `x` - A monotonicly ascending or decending table of `x\[0\], ..., x\[n-1\]`
    /// * `y` - A table of `f(x\[0\]), ..., f(x\[n-1\])`
    /// FINISH DOCS
    /// # Returns
    /// * `Self` - Linear1D Instaciate
    pub fn new(x: Array1<f64>, y: Array1<f64>, yp1: Option<f64>, ypn: Option<f64>) -> Self {
        let search = BisectHunt1D::new(&x, 2);
        let y2 = Self::set_y2(&x, &y, &yp1, &ypn);

        Self {
            x: x,
            y: y,
            y2: y2,
            search: search,
        }
    }

    fn set_y2(
        x: &Array1<f64>,
        y: &Array1<f64>,
        yp1: &Option<f64>,
        ypn: &Option<f64>,
    ) -> Array1<f64> {
        let n = y.len();

        let mut y2 = Array1::<f64>::zeros([n]);
        let mut u = Array1::<f64>::zeros([n - 1]);

        if yp1.is_some() {
            y2[0] = -0.5;
            u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1.unwrap());
        }

        let a = &x.slice(s![0..n - 2]).to_owned(); //  i - 1
        let b = &x.slice(s![1..n - 1]).to_owned(); //  i
        let c = &x.slice(s![2..n]).to_owned(); //      i + 1

        let d = y2.slice(s![0..n - 2]).to_owned(); // i - 1

        let sig = &((b - a) / (c - a));
        let p = &((sig * d) + 2.);

        y2.slice_mut(s![1..n - 1]).assign(&((sig - 1.) / p));

        let mut un = 0.;
        let mut qn = 0.;

        if ypn.is_some() {
            qn = 0.5;
            un = (3.0 / (x[n - 1] - x[n - 2]))
                * (ypn.unwrap() - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
        }

        y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.);

        let a = &y2.slice(s![0..n-2;-1]).to_owned(); //  k
        let b = &y2.slice(s![1..n-1;-1]).to_owned(); // k + 1
        let c = &u.slice(s![0..n-2;-1]).to_owned(); //   k

        y2.slice_mut(s![..n-2;-1]).assign(&(a * b + c));

        y2
    }

    fn _interpolate(&mut self, x: f64, i: usize) -> f64 {
        //
        let h = self.x[i + 1] - self.x[i];

        // x's must be unique
        if h == 0.0 {
            panic!("Bad input to routine spline interpolation.")
        };
        let a = (self.x[i + 1] - x) / h;
        let b = (x - self.x[i]) / h;

        let p = a * self.y[i];
        let q = b * self.y[i + 1];
        let r = (a.powi(3) - a) * self.y2[i];
        let s = (b.powi(3) - b) * self.y2[i + 1];

        let y = p + q + (r + s) * (h * h) / 6.0;

        y
    }
}

impl Interpolate1D for Spline1D {
    type Dtype = f64;
    type Data = Array1<Self::Dtype>;
    type Index = usize;
    /// Interpolation
    /// #  Arguments
    /// * `x` - The x value for which `f(x)` is being approximated
    /// # Returns
    /// * `y : f64` - An approximation of `f(x)`
    fn interpolate(&mut self, x: Self::Dtype) -> Self::Dtype {
        // Get indecies
        let i = self.search.locate(&self.x, x);

        // Evaluate
        self._interpolate(x, i)
    }
}

pub fn proof() {
    // Revisit this
    println!("Starting proof...");

    let x_gt = Array1::linspace(0., 3.14159265 * 2., 1000);
    let y_gt = x_gt.clone().map(|&x| f64::sin(x));

    // let iter_gt = std::iter::zip(x_gt, y_gt);

    let x = Array1::linspace(0., 3.14159265 * 2., 7);

    let mut my_interp = Spline1D::new(x_gt.clone(), y_gt.clone(), None, None);

    let y = x.clone().map(|&x| my_interp.interpolate(x));

    let root_area =
        BitMapBackend::new("images/Spline1D_proof.png", (1920, 1040)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut ctx = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption("Linear 1D Interpolation Proof", ("sans-serif", 40))
        .build_cartesian_2d(0f64..(3.14159265f64 * 2.), -1.0f64..1.0f64)
        .unwrap();

    ctx.configure_mesh().draw().unwrap();

    ctx.draw_series(LineSeries::new(
        x_gt.iter().zip(y_gt.iter()).map(|(&x, &y)| (x, y)),
        &RED,
    ))
    .unwrap();

    ctx.draw_series(LineSeries::new(
        x.iter().zip(y.iter()).map(|(&x, &y)| (x, y)),
        &GREEN,
    ))
    .unwrap();

    println!("Proof complete.");
}
