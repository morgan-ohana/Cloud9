use crate::hydrostatics::isothermal_abg_background;

mod hydrostatics;
mod plotting;

fn main() {
    let init_num_density = 2e18;

    isothermal_abg_background(
        300.0,
        1.0,
        3.0,
        1.0,
        1.0e7,
        3.0,
        init_num_density * hydrostatics::M_P,
    );
}
