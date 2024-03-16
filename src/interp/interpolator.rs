/// Interpolator Trait
pub trait Interpolate1D {
    type Dtype; // Base datatype
    type Data; // N dimesional array of the base datatype
    type Index; // Index to retirieve a single data from the input/output array
                // Aproximates f(x) from x given a monotonicly increasing or decending
                // tables x[0], ..., x[n-1], and y[0],  ..., y[n-1]
    fn interpolate(&mut self, x: Self::Dtype) -> Self::Dtype;
}

pub trait InterpolateND {
    type YDim;

    type Dtype;

    type X;
    type IndX;

    type Y;
    type IndY;

    fn interpolate(&mut self, x: Self::X) -> Self::Y;
}
