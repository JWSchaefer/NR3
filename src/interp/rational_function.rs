/// Come back to this ->
use crate::interp::interpolator::Interpolate;
use crate::table::bisect_hunt::BisectHunt1D;
use crate::table::search::Search;
use ndarray::prelude::*;
use plotters::prelude::*;

/// 1 Dimensional polynomial interpolator
/// Based on Neville's Algorithm
pub struct Poly1D {
    x: Array1<f64>,       // x table
    y: Array1<f64>,       // y table
    m: usize,             // Polynomial degree / bracket size
    search: BisectHunt1D, // Search Algorithm
}

impl Poly1D {
    /// Constructor
    /// #  Arguments
    /// * `x` - A monotonicly ascending or decending table of `x\[0\], ..., x\[n-1\]`
    /// * `y` - A table of `f(x\[0\]), ..., f(x\[n-1\])`
    /// * `m` - Polynomial deree. Must be lesser than the length of `x`
    /// # Returns
    /// * `Self` - Poly1D Instaciate
    pub fn new(x: Array1<f64>, y: Array1<f64>, m: usize) -> Self {
        // Instaciate search algorithm
        let search = BisectHunt1D::new(&x, m);
        Self {
            x: x,
            y: y,
            m: m,
            search: search,
        }
    }

    /// Raw Interpolator
    /// #  Arguments
    /// * `x` - The x value for which `f(x)` is being approximated
    /// * `i` - The index for which `x\[i\] <= x <= x\[i + 1\]` is guarenteed
    /// # Returns
    /// * `y : f64` - An approximation of `f(x)`
    fn _interpolate(&mut self, x: f64, i: usize) -> f64 {
        let xa = self.x.slice(s![i..i + self.m + 1]);
        let mut ya = self.y.slice_mut(s![i..i + self.m + 1]).to_owned();

        for m in 1..self.m + 1 {
            for n in 0..self.m - m {
                let p_ip1_j = ya[n + 1];
                let p_i_jm1 = ya[n];
                let x_i = xa[n];
                let x_j = xa[n + m];
                ya[n] = (((x - x_i) * p_ip1_j) - ((x - x_j) * p_i_jm1)) / (x_j - x_i);
            }
        }

        return ya[0];
    }
}

impl Interpolate for Poly1D {
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

    let mut my_interp = Poly1D::new(x_gt.clone(), y_gt.clone(), 4);

    let y = x.clone().map(|&x| my_interp.interpolate(x));

    let root_area = BitMapBackend::new("images/Poly1D_proof.png", (1920, 1040)).into_drawing_area();
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
