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

pub fn svg_to_pdf(output_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let output_path = String::from(output_path);
    let svg_path = &(output_path.clone() + ".svg");
    let pdf_path = &(output_path.clone() + ".pdf");

    let status = std::process::Command::new("inkscape")
        .args(&["--export-type=pdf", "--export-filename", pdf_path, svg_path])
        .status()?;

    if !status.success() {
        return Err(format!("Inkscape failed with status: {}", status).into());
    }

    let compressed_path = &(output_path.clone() + "_compressed.pdf");
    let status = std::process::Command::new("gs")
        .args(&[
            "-dBATCH",
            "-dNOPAUSE",
            "-dQUIET",
            "-sDEVICE=pdfwrite",
            "-dCompatibilityLevel=1.5",
            &format!("-sOutputFile={}", compressed_path),
            pdf_path,
        ])
        .status()?;
    if status.success() {
        std::fs::rename(compressed_path, pdf_path)?;
        Ok(())
    } else {
        Err(format!("Ghostscript failed with status: {}", status).into())
    }
}
