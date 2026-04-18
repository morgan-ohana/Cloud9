use rkyv::{deserialize, rancor::Error};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};

pub fn write_cells_to_csv(
    cells: &Vec<(f64, (f64, f64, f64), (f64, f64, f64), f64)>,
    path: &str,
) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);

    // Header
    writeln!(
        w,
        "m200_min,c200_min,tau_min,m200_max,c200_max,tau_max,density,level"
    )?;

    for (density, (x0, y0, z0), (x1, y1, z1), level) in cells {
        writeln!(
            w,
            "{},{},{},{},{},{},{},{}",
            x0, y0, z0, x1, y1, z1, density, level
        )?;
    }

    Ok(())
}
