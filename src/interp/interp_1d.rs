use crate::interp::interpolator::Interpolate;
use crate::table::bisect_hunt::BisectHunt1D;
use crate::table::search::Search;
use ndarray::prelude::*;
use plotters::prelude::*;

pub struct Linear1D {
    x: Array1<f64>,
    y: Array1<f64>,
    search: BisectHunt1D,
}

impl Linear1D {
    fn lerp(x: f64, x_0: f64, x_1: f64, y_0: f64, y_1: f64) -> f64 {
        return y_0 + ((y_1 - y_0) * ((x - x_0) / (x_1 - x_0)));
    }
}

impl Interpolate for Linear1D {
    type Dtype = f64;
    type Input = Array1<Self::Dtype>;
    type Index = usize;

    fn new(x: Self::Input, y: Self::Input) -> Self {
        let mut search = BisectHunt1D::new(&x);
        search.set_m(2);
        Self {
            x: x,
            y: y,
            search: search,
        }
    }

    fn interpolate(&mut self, x: Self::Dtype) -> Self::Dtype {
        // Get indecies for interpolation
        let (i_lower, i_upper) = self.search.locate(&self.x, x);
        // Evaluate
        Self::lerp(
            x,
            self.x[i_lower],
            self.x[i_upper],
            self.y[i_lower],
            self.y[i_upper],
        )
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
        .build_cartesian_2d(0f64..(3.14159265f64 * 4.), -1.0f64..1.0f64)
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
