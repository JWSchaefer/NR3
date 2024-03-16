// use recipies::interp::biliniar;
use recipies::interp::linear;
use recipies::interp::polynomial_ceoficients;
use recipies::interp::spline;

fn main() {
    linear::proof();
    polynomial_ceoficients::proof();
    spline::proof();
    // biliniar::proof();
}
