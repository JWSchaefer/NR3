use crate::table::search::Search;

use ndarray::Array1;
use std::cmp;

/// Bookmarking variables
pub struct BisectHunt1D {
    m: usize,              // Interval size : x centered in x[j_lo], ..., x[j_lo + m]
    corr: bool,            // Previous searches coreelated?
    ascend: bool,          // x ascending?
    dj: usize,             // Determines when worth hunting
    i_save: Option<usize>, // Previous x : Used if correlated
}

/// Defines functions unique to this type of table search
impl BisectHunt1D {
    pub fn set_m(&mut self, m: usize) {
        if m < 2 {
            panic!("Invalid search parameter m")
        }
        self.m = m;
    }
    fn hunt_bisect(&mut self, table: &Array1<f64>, x: f64) -> usize {
        // Check search value exsits
        // if match self.ascend {
        //     true => (x < table[0]) || (x > table[table.len() - 1]),
        //     false => (x > table[0]) || (x < table[table.len() - 1]),
        // } {
        //     panic!("Search value outside table bounds ")
        // }

        // Hunt step size, doubles each iteration
        let mut inc: usize = 1;

        let mut lower: usize;
        let mut upper: usize;

        // If no saved value, defaults to a normal bisection & full table bracket
        if !self.corr {
            lower = 0;
            upper = table.len() - 1;
        }
        // If correlated value
        else {
            lower = self.i_save.unwrap();
            upper = self.i_save.unwrap();
            // Hunt up
            if (x >= table[lower]) == self.ascend {
                loop {
                    // upper index bound to be exceded -> Sucess! - take limit
                    if inc + upper > table.len() - 1 {
                        upper = table.len() - 1;
                        break;
                    }

                    // Update bracket
                    upper = lower + inc;

                    // Value bracketed -> Sucess!
                    if (x < table[upper]) == self.ascend {
                        break;
                    }
                    // Failure... -> Double bracket size
                    else {
                        upper = lower;
                        inc <<= 1;
                    }
                }
            }
            // Hunt down
            else {
                loop {
                    // Lower index bound to be exceded -> Sucess! - take limit
                    if inc > lower {
                        lower = 0;
                        break;
                    }
                    // Update bracket
                    lower = lower - inc;

                    // Value bracketed -> Sucess!
                    if (x >= table[lower]) == self.ascend {
                        break;
                    }
                    // Failure... -> Double bracket size
                    else {
                        upper = lower;
                        inc <<= 1
                    };
                }
            }
        }

        // Bisect on identified bracket
        while upper - lower > 1 {
            let mid = upper + lower >> 1; // Calculate midpoint
            match (x >= table[mid]) == self.ascend {
                true => lower = mid,
                false => upper = mid,
            }
        }

        // If previous value was saved, determine if current point is close to previosu value
        if self.i_save.is_some() {
            self.corr = lower.abs_diff(self.i_save.unwrap()) < self.dj
        }

        // Save location
        self.i_save = Some(lower);

        // Return j_lo, while respecting index bounds
        cmp::max(
            0,
            cmp::min(table.len() - self.m, lower - ((self.m - 2) >> 2)),
        )
    }
}

impl Search for BisectHunt1D {
    type Dtype = f64;
    type Input = Array1<f64>;
    type Index = usize;

    fn new(table: &Self::Input) -> Self {
        let dj = cmp::min(1, f64::powf(table.len() as f64, 0.25).trunc() as usize);
        let ascend = table[table.len() - 1] > table[0];

        Self {
            m: 2,
            corr: false,
            ascend: ascend,
            dj: dj,
            i_save: None,
        }
    }

    fn locate(&mut self, table: &Self::Input, x: f64) -> (Self::Index, Self::Index) {
        let res = self.hunt_bisect(table, x);
        return (res, res + (self.m - 1));
    }
}
