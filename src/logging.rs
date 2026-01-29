use rkyv::{Archive, Deserialize, Serialize, deserialize, rancor::Error};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};

#[derive(Archive, Deserialize, Serialize, Clone, Debug)]
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

    println!("Saved Initial Conditions at {file_name}");
    Ok(())
}
