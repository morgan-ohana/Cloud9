use crate::fitting;

pub fn make_file_name(
    num_walkers: usize,
    steps: usize,
    prior: &fitting::Prior,
    fixed_cross_section: &Option<f64>,
) -> String {
    let mut file_name = format!("{}_x_{}k", num_walkers, steps / 1000);

    match prior {
        fitting::Prior::MassConcentrationRelation(_) => {
            file_name.push_str("_MCR");
        }
        fitting::Prior::None => {}
    }

    if let Some(cross_section) = fixed_cross_section {
        file_name.push_str(&format!("_sigma={cross_section}"));
    }

    file_name
}
