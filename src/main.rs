use recipies::interp::linear;
use recipies::interp::polynomial;
use recipies::interp::spline;

fn main() {
    linear::proof();
    polynomial::proof();
    spline::proof();
}
