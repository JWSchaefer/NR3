// use crate::interp::interpolator::Interpolate;
// use crate::table::bisect_hunt::BisectHunt1D;
// use crate::table::search::Search;
// use ndarray::{prelude::*, stack, Axis};
// use plotters::prelude::*;

// /// 1 Dimensional linear Instaciatenterpolator
// pub struct Linear2D {
//     x1: Array1<f64>,       // x1 table
//     x2: Array1<f64>,       // x2 table
//     y: Array2<f64>,        // y table
//     search1: BisectHunt1D, // Search Algorithm x1
//     search2: BisectHunt1D, // Search Algorithm x2
// }

// impl Linear2D {
//     /// Constructor
//     /// #  Arguments
//     /// * `x1` - A monotonicly ascending or decending table of `x\[0\], ..., x\[n-1\]`
//     /// * `x2` - A monotonicly ascending or decending table of `x\[0\], ..., x\[n-1\]`
//     /// * `y` - A table of `f(x\[0\]), ..., f(x\[n-1\])`
//     /// # Returns
//     /// * `Self` - Linear2D Instance
//     pub fn new(x1: Array1<f64>, x2: Array1<f64>, y: Array2<f64>) -> Self {
//         // Instaciate search algorithm
//         let search1 = BisectHunt1D::new(&x1, 2);
//         let search2 = BisectHunt1D::new(&x2, 2);
//         Self {
//             x1: x1,
//             x2: x2,
//             y: y,
//             search1: search1,
//             search2: search2,
//         }
//     }

//     fn interpolate(&mut self, x: [f64; 2]) -> f64 {
//         let [x1, x2] = x;

//         // Get indicies
//         let i = self.search1.locate(&self.x1, x1);
//         let j = self.search2.locate(&self.x2, x2);

//         // Evaluate
//         let t = (x1 - self.x1[i]) / (self.x1[i + 1] - self.x1[i]);
//         let u = (x2 - self.x2[j]) / (self.x2[j + 1] - self.x2[j]);
//         let y = (1. - t) * (1. - u) * self.y[[i, j]]
//             + t * (1. - u) * self.y[[i + 1, j]]
//             + (1. - t) * u * self.y[[i, j + 1]]
//             + t * u * self.y[[i + 1, j + 1]];

//         return y;
//     }
// }

// pub fn proof() {
//     // Revisit this
//     println!("Starting proof...");

//     fn f(x1: &f64, x2: &f64) -> f64 {
//         x1 - x2
//     }

//     let x1 = Array1::linspace(0., 1., 5);
//     let x2 = Array1::linspace(0., 2., 5);

//     let mut y = Array2::<f64>::zeros([x1.len(), x2.len()]);

//     for i in 0..x1.len() {
//         for j in 0..x2.len() {
//             y[[i, j]] = f(&x1[i], &x2[j]);
//         }
//     }

//     for row in y.iter() {
//         println!("{}", row)
//     }
// }
