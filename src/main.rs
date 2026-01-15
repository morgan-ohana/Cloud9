use crate::hydrostatics::isothermal_abg_background;

mod hydrostatics;
mod plotting;

fn main() {
    isothermal_abg_background(300.0, 1.0, 3.0, 1.0, 0.0e7, 3.0, 1.0e3);
}
