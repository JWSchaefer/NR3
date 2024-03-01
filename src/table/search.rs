pub trait Search {
    type Dtype;
    type Input;
    type Index;

    fn locate(&mut self, table: &Self::Input, x: Self::Dtype) -> Self::Index;
}
