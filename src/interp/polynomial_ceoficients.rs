use crate::interp::polynomial::Poly1D;
use ndarray::{prelude::*, Zip};
use std::ops::SubAssign;

fn polcoe_0(x: &Array1<f64>, y: &Array1<f64>) -> Array1<f64> {
    let n = x.len();
    let mut s = Array1::<f64>::zeros([n]);
    let mut coef = Array1::<f64>::zeros([n]);

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

fn polcoe_1(x: &Array1<f64>, y: &Array1<f64>) -> Array1<f64> {
    let shape = x.raw_dim();
    let n = x.len();

    let mut s = Array1::<f64>::zeros(shape);
    let mut b = Array1::<f64>::ones(shape);
    let mut coef = Array1::<f64>::zeros(shape);

    s[n - 1] = -x[0];

    for i in 1..n {
        let s_j_1 = s.slice(s![n - i..n]).to_owned();
        let mut s_j = s.slice_mut(s![n - 1 - i..n - 1]);

        s_j.sub_assign(&(x[i] * s_j_1));

        s[n - 1] -= x[i];
    }

    let mut phi = Array1::<f64>::from_elem(x.raw_dim(), n as f64);
    let k = Array1::linspace((n - 1) as f64, 1., n - 1);
    let s_k = s.slice(s![1..n;-1]);

    phi.zip_mut_with(&x, |p, x_j| {
        *p = (&k * &s_k).iter().fold(*p, |acc, _a| _a + x_j * acc)
    });

    phi = y / phi;

    let k = s![0..n;-1];

    let mut coef_k = coef.slice_mut(k);
    let s_k = s.slice(k);

    Zip::from(&mut b).and(&phi).and(x).for_each(|_b, &_p, &_x| {
        Zip::from(&mut coef_k).and(&s_k).for_each(|_c, &_s| {
            *_c += *_b * _p;
            *_b = _s + _x * *_b;
        })
    });

    coef
}

fn polcof(x: &Array1<f64>, y: &Array1<f64>) -> Array1<f64> {
    let n = x.len();
    let mut coef = Array1::<f64>::from_elem([n], 0.);

    let mut xa = x.clone().to_owned();
    let mut ya = y.clone().to_owned();

    for j in 0..n {
        let mut interp = Poly1D::new(xa.clone().to_owned(), ya.clone().to_owned(), n - j - 1);

        coef[j] = interp.raw_interpolate(0., 0);

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

pub fn proof() {
    println!("Starting proof...");

    // Test data
    let x = Array1::linspace(0., 1., 10);
    let y = Array1::linspace(400., 500., 10);

    // Calculate coeficients
    let coefs_1_0 = polcoe_0(&x.clone().to_owned(), &y.clone().to_owned());
    let coefs_1_1 = polcoe_1(&x.clone().to_owned(), &y.clone().to_owned());
    let coefs_2 = polcof(&x.clone().to_owned(), &y.clone().to_owned());

    // Define test routine
    fn test(coefs: Array1<f64>, x: &Array1<f64>, y: &Array1<f64>) {
        let _x = Array2::<f64>::from_shape_fn([x.len(), x.len()], |(i, j)| x[i].powi(j as i32));

        let err = 100. * (y - _x.dot(&coefs)).dot(&(y - _x.dot(&coefs))) / (x.len() as f64);
        println!("Error:\t{err:.3} %");
    }

    // Run tests
    println!("\npolcoe_0");
    test(coefs_1_0, &x, &y);

    println!("\npolcoe_0");
    test(coefs_1_1, &x, &y);

    println!("\npolcof");
    test(coefs_2, &x, &y);
}
