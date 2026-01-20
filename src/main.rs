use std::f64::consts::PI;

use crate::hydrostatics::{isothermal_abg_background, isothermal_core_collapse_background};
use crate::fitting::find_parameters;
use crate::plotting::plot_function;

mod hydrostatics;
mod plotting;
mod fitting;

fn main() {
    let data = vec![
        (0.15649748494408605,	15595734576138530000.0),
        (0.18722388772639217,	15290963190523656000.0),
        (0.2252243432550796,	14924613361774703000.0),
        (0.26936457583283363,	14424784086219026000.0),
        (0.32301280500582225,	13894659794863262000.0),
        (0.38799031821006913,	13091285860285037000.0),
        (0.4647330501352986,	12151824423121148000.0),
        (0.5571917155289376,	11103690653134041000.0),
        (0.6692414165028057,	9918955303017603000.0),
        (0.801419168108796, 	8606124456897346000.0),
        (0.9595907722097119,	7014426686842890000.0),
        (1.1500508076747187,	5250966932011177000.0),
        (1.3771319397378698,	3511010519320841000.0),
        (1.6463993031339876,	2293590617222256000.0),
        (1.983419488324015,	    1615272370298475800.0),
    ];

    let sidm_fit_params = [
        6282772.676310997,
        2.8897601910357062,
        0.18438405302532834,
        54574.17044567845,
    ];

    let cdm_fit_params = [
        7585724.648997071,
        2.195234319751577,
        0.0,
        169116.20672504138,
    ];
    
    //let params = find_parameters(&data, [2e7, 4.0, 0.2, 5e5], Some(0.0));
    let params = cdm_fit_params;

    dbg!(params);
    let fit = isothermal_core_collapse_background(
        1e4,
        params[0],
        params[1],
        params[2],
        (0.1*data[0].0, 10.0*data.last().unwrap().0),
        params[3],
        false
    );

    let legend_text = format!("rho_s_0: {}\nr_s_0: {}\ntau: {}\nrho_c: {}", params[0], params[1], params[2], params[3]); 
    println!("{}", &legend_text);

    plot_function(
        &fit.0,
        &fit.1,
        "fit.png",
        "Fit To Observed Profile",
        "r (arcmin)",
        "n_H (num / cm^2)",
        Some(legend_text),
        Some(&data)
    ).unwrap();
}

fn plot_distribution() {
    let target_init_num_col_density = 1e19;
    let target_init_col_density: f64 = target_init_num_col_density * hydrostatics::M_PROTON * hydrostatics::CM_IN_KPC.powi(2);
    let temperature = 1e4;
    let mut sound_speed_squared = temperature / (hydrostatics::MP_OVER_KB); // kpc^2 / s^2
    sound_speed_squared *= hydrostatics::KM_IN_KPC.powi(2); // km^2 / s^2
    
    // r_s = c_s / sqrt(4pi G rho_c)
    // uniform sphere approx:
    // Sigma ~ 2*r_s*rho_c = 2 c_s sqrt(rho_c / 4pi G) => rho_c ~ 4pi G * (Sigma / 2c_s)^2 = pi G Sigma^2 / c_s^2 
    // rho_c (kpc^-3) ~ pi G * (n_sigma (cm^-2) * MP * CMinKPC^2)^2 / c_s^2 
    // Units: [rho_c] = M_sun kpc^-3 = [G * sigma^2 / c_s^2] = km^2 kpc M_sun^-1 s^-2 * M_sun^2 kpc^-4 * [c_s]^-2 = km^2 kpc^-3 M_sun s^-2 * [c_s]^-2
    // => [c_s^2] = km^2 kpc^-3 M_sun s^-2 / M_sun kpc^-3 = km^2 / s^2
    
    let rho_center_approx = PI * hydrostatics::GG * (target_init_col_density).powi(2) / sound_speed_squared;
    let scale_radius = sound_speed_squared.sqrt() / (4.0 * PI * hydrostatics::GG * rho_center_approx);
    println!("Center rho set to: {rho_center_approx}");

    // isothermal_abg_background(
    //     temperature,
    //     1.0,
    //     3.0,
    //     1.0,
    //     0.0e7,
    //     3.0,
    //     (1e-1, 10.0*scale_radius),
    //     rho_center_approx,
    // );

    isothermal_core_collapse_background(
        temperature,
        3e7,
        3.0,
        0.2,
        (1e-1, 10.0*scale_radius),
        rho_center_approx,
        true
    );
}
