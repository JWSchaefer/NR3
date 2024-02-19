pub trait Search {
    type Dtype;
    type Input;
    type Index;

    fn new(table: &Self::Input) -> Self;
    fn locate(&mut self, table: &Self::Input, x: Self::Dtype) -> (Self::Index, Self::Index);
}
