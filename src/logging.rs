use rkyv::{Archive, Deserialize, Serialize, deserialize, rancor::Error};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};

#[derive(Archive, Deserialize, Serialize, Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct MCMCOutput {
    pub best_params: [f64; 4],
    pub chain: Vec<[f64; 4]>,
    pub likelihoods: Vec<f64>,
}

pub fn load_file(file_name: String) -> anyhow::Result<MCMCOutput> {
    let file = File::open(&file_name)?;
    let mut reader = BufReader::new(file);

    let mut buffer = Vec::new();
    reader.read_to_end(&mut buffer)?;

    let archived = rkyv::access::<ArchivedMCMCOutput, Error>(&buffer)
        .map_err(|e| anyhow::anyhow!("Data verification failed: {}", e))?;

    let output: MCMCOutput = deserialize::<MCMCOutput, Error>(archived)?;

    println!("File loaded from {file_name}");
    Ok(output)
}

pub fn save_output(
    file_name: String,
    best_params: [f64; 4],
    chain: Vec<[f64; 4]>,
    likelihoods: Vec<f64>,
) -> anyhow::Result<()> {
    let output = MCMCOutput {
        best_params,
        chain,
        likelihoods,
    };

    let file = File::create(&file_name)?;
    let mut writer = BufWriter::new(file);

    let bytes = rkyv::to_bytes::<Error>(&output)
        .map_err(|e| anyhow::anyhow!("Serialization failed: {}", e))?;

    writer.write_all(&bytes)?;

    println!("Saved output at {file_name}");
    Ok(())
}

pub fn save_output_json(
    file_name: String,
    best_params: [f64; 4],
    chain: Vec<[f64; 4]>,
    likelihoods: Vec<f64>,
) -> anyhow::Result<()> {
    let output = MCMCOutput {
        best_params,
        chain,
        likelihoods,
    };

    let file = File::create(&file_name)?;
    let writer = BufWriter::new(file);

    // Serialize to JSON
    serde_json::to_writer_pretty(writer, &output)
        .map_err(|e| anyhow::anyhow!("JSON serialization failed: {}", e))?;

    println!("Saved json data at {file_name}");
    Ok(())
}

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
