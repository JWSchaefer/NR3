use crate::table::bisect_hunt::BisectHunt1D;
use crate::table::search::Search;
use crate::{interp::interpolator::Interpolate, table::search};
use ndarray::{prelude::*, stack, Axis};
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

    let x = Array1::linspace(-1., 1., 4);
    // let y = x.map(|i| i * i);
    let y = Array1::linspace(-1., 1., 4);

    let coefs1 = ceoficients_2(&x.clone().to_owned(), &y.clone().to_owned());

    let _x = stack![
        Axis(1),
        x.clone().to_owned().mapv(|a| a.powi(0)),
        x.clone().to_owned().mapv(|a| a.powi(1)),
        x.clone().to_owned().mapv(|a| a.powi(2)),
        x.clone().to_owned().mapv(|a| a.powi(3))
    ];

    println!("Coeficients: {}", coefs1);
    println!("True y {}", y);
    println!("X . ceofs {}", _x.dot(&coefs1));

    // let x_gt = Array1::linspace(0., 3.14159265 * 2., 1000);
    // let y_gt = x_gt.clone().map(|&x| f64::sin(x));

    // // let iter_gt = std::iter::zip(x_gt, y_gt);

    // let x = Array1::linspace(0., 3.14159265 * 2., 7);

    // let mut my_interp = Poly1D::new(x_gt.clone(), y_gt.clone(), 4);

    // let y = x.clone().map(|&x| my_interp.interpolate(x));

    // let root_area = BitMapBackend::new("images/Poly1D_proof.png", (1920, 1040)).into_drawing_area();
    // root_area.fill(&WHITE).unwrap();

    // let mut ctx = ChartBuilder::on(&root_area)
    //     .set_label_area_size(LabelAreaPosition::Left, 40)
    //     .set_label_area_size(LabelAreaPosition::Bottom, 40)
    //     .caption("Linear 1D Interpolation Proof", ("sans-serif", 40))
    //     .build_cartesian_2d(0f64..(3.14159265f64 * 2.), -1.0f64..1.0f64)
    //     .unwrap();

    // ctx.configure_mesh().draw().unwrap();

    // ctx.draw_series(LineSeries::new(
    //     x_gt.iter().zip(y_gt.iter()).map(|(&x, &y)| (x, y)),
    //     &RED,
    // ))
    // .unwrap();

    // ctx.draw_series(LineSeries::new(
    //     x.iter().zip(y.iter()).map(|(&x, &y)| (x, y)),
    //     &GREEN,
    // ))
    // .unwrap();

    // println!("Proof complete.");
}

fn ceoficients_1(x: Array1<f64>, y: Array1<f64>) -> Array1<f64> {
    let n = x.len();
    let mut s = Array1::<f64>::from_elem([n], 0.);
    let mut coef = Array1::<f64>::from_elem([n], 0.);

    s[n - 1] = -x[0];

    for i in 1..n {
        for j in n - 1 - i..n - 1 {
            s[j] -= x[i] * s[j + 1];
        }
        s[n - 1] -= x[i];
    }

    for j in 0..n {
        let mut phi = n as f64;
        for k in (1..n).rev() {
            phi = (k as f64) * s[k] + x[j] * phi
        }
        let ff = y[j] / phi;
        let mut b = 1.0;
        for k in (0..n).rev() {
            coef[k] += b * ff;
            b = s[k] + x[j] * b;
        }
    }

    coef
}

fn ceoficients_2(x: &Array1<f64>, y: &Array1<f64>) -> Array1<f64> {
    let n = x.len();
    let mut coef = Array1::<f64>::from_elem([n], 0.);

    let mut xa = x.clone().to_owned();
    let mut ya = y.clone().to_owned();

    for j in 0..n {
        println!("{n}-{j} : {}", n - j);

        let mut interp = Poly1D::new(xa.clone().to_owned(), ya.clone().to_owned(), n - j - 1);

        coef[j] = interp._interpolate(0., 0);

        let mut k = None;
        let mut xmin = 1.0e99;

        for i in 0..n - j {
            if xa[i].abs() < xmin {
                xmin = xa[i].abs();
                k = Some(i);
            }
            if xa[i] != 0.0 {
                ya[i] = (ya[i] - coef[j]) / xa[i]
            }
        }

        for i in (k.unwrap() + 1)..n - j {
            ya[i - 1] = ya[i];
            xa[i - 1] = xa[i];
        }
    }

    coef
}
