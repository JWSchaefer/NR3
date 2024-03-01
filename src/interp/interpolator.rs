/// Interpolator Trait
pub trait Interpolate {
    type Dtype; // Base datatype, floating point
    type Data; // N dimesional array of the base datatype
    type Index; // Index to retirieve a single data from the input/output array
                // Aproximates f(x) from x given a monotonicly increasing or decending
                // tables x[0], ..., x[n-1], and y[0],  ..., y[n-1]
    fn interpolate(&mut self, x: Self::Dtype) -> Self::Dtype;
}
