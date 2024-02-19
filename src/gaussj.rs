use ndarray::prelude::*;
use ndarray::{Data, DataMut, RawData};

fn swap_row<S, D>(a: &mut ArrayBase<S, D>, b: &mut ArrayBase<S, D>)
where
    S: Data + DataMut + RawData,
    <S as RawData>::Elem: Clone,
    D: Dimension,
{
    let temp = a.to_owned().clone();
    a.assign(&b);
    b.assign(&temp);
}

fn gaussj() {
    let mut a: Array2<f64> = arr2(&[
        [1., 2., 3., 4.],
        [4., 5., 6., 7.],
        [8., 9., 10., 11.],
        [12., 13., 14., 15.],
    ]); // n x n -> 4 x 4

    let mut b: Array2<f64> = arr2(&[[0., 1.], [2., 3.], [4., 5.], [6., 7.]]); // n x m -> 4 x 2

    println!("a\n{a}\n");
    println!("b\n{b}\n");

    let n: usize = a.nrows();
    let m: usize = b.ncols();

    let (mut irow, mut icol): (usize, usize) = (0, 0);

    let (mut big, mut dum, mut pivinv): (f64, f64, f64);

    let mut indxc = vec![0; n];
    let mut indxr = vec![0; n];
    let mut ipiv = vec![0; n];

    for i in 0..n {
        big = 0.;

        for j in 0..n {
            if ipiv[j] != 1 {
                for k in 0..n {
                    if ipiv[k] == 0 {
                        if a[[j, k]] > big {
                            big = a[[j, k]];
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }

            ipiv[icol] += 1;

            if irow != icol {
                swap_row(
                    &mut a.index_axis_mut(Axis(0), irow),
                    &mut a.index_axis_mut(Axis(0), icol),
                );

                swap_row(
                    &mut b.index_axis_mut(Axis(0), irow),
                    &mut b.index_axis_mut(Axis(0), icol),
                );

                // for l in 0..n { swap(&a[[irow, l]], &a[[icol, l]]) };
                // for l in 0..m { swap(&b[[irow, l]], &b[[icol, l]])  };
            }

            indxr[i] = irow;
            indxc[i] = icol;

            if a[[icol, icol]] == 0.0 {
                panic!("gaussj: Singular Matrix");
            }

            pivinv = 1.0 / a[[icol, icol]];

            a[[icol, icol]] = 1.0;

            for l in 0..n {
                a[[icol, l]] *= pivinv;
            }
            for l in 0..m {
                b[[icol, l]] *= pivinv;
            }

            for ll in 0..n {
                if ll != icol {
                    dum = a[[ll, icol]];
                    a[[ll, icol]] = 0.0;
                    for l in 0..n {
                        a[[ll, l]] -= a[[icol, l]] * dum;
                    }
                    for l in 0..m {
                        b[[ll, l]] -= b[[icol, l]] * dum;
                    }
                }
            }
        }
        for l in n - 1..0 {
            if indxr[l] != indxc[l] {
                swap_row(
                    &mut b.index_axis_mut(Axis(1), indxr[l]),
                    &mut b.index_axis_mut(Axis(1), indxc[l]),
                );

                // for k in 0..n {
                //     swap(&a[[k,indxr[l]]], &a[[k, indxc[l]]]);
                //     }
            }
        }
    }

    println!("a\n{a}\n");
    println!("b\n{b}\n");
}
