pub trait Interpolate {
    type Dtype;
    type Input;
    type Index;

    fn new(x: Self::Input, y: Self::Input) -> Self;

    fn interpolate(&mut self, x: Self::Dtype) -> Self::Dtype;
}
