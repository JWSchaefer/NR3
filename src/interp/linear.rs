use crate::interp::interpolator::Interpolate1D;
use crate::table::bisect_hunt::BisectHunt1D;
use crate::table::search::Search;
use ndarray::prelude::*;
use plotters::prelude::*;

/// 1 Dimensional linear Instaciatenterpolator
pub struct Linear1D {
    x: Array1<f64>,       // x table
    y: Array1<f64>,       // y table
    search: BisectHunt1D, // Search Algorithm
}

impl Linear1D {
    /// Constructor
    /// #  Arguments
    /// * `x` - A monotonicly ascending or decending table of `x\[0\], ..., x\[n-1\]`
    /// * `y` - A table of `f(x\[0\]), ..., f(x\[n-1\])`
    /// # Returns
    /// * `Self` - Linear1D Instaciate
    pub fn new(x: Array1<f64>, y: Array1<f64>) -> Self {
        // Instaciate search algorithm
        let search = BisectHunt1D::new(&x, 2);
        Self {
            x: x,
            y: y,
            search: search,
        }
    }

    /// Raw Interpolator
    /// #  Arguments
    /// * `x` - The x value for which `f(x)` is being approximated
    /// * `i` - The index for which `x\[i\] <= x <= x\[i + 1\]` is guarenteed
    /// # Returns
    /// * `y : f64` - An approximation of `f(x)`
    fn _interpolate(&self, x: f64, i: usize) -> f64 {
        let x_0 = self.x[i];
        let x_1 = self.x[i + 1];
        let y_0 = self.y[i];
        let y_1 = self.y[i + 1];

        return y_0 + ((y_1 - y_0) * ((x - x_0) / (x_1 - x_0)));
    }
}

impl Interpolate1D for Linear1D {
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

    let mut my_interp = Linear1D::new(x_gt.clone(), y_gt.clone());

    let y = x.clone().map(|&x| my_interp.interpolate(x));

    let root_area =
        BitMapBackend::new("images/Linear1D_proof.png", (1920, 1040)).into_drawing_area();
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
