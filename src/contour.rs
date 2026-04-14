use crate::logging::write_cells_to_csv;

pub fn get_3d_contour(chain: &Vec<Vec<f64>>, bounds: &[[f64; 2]; 4]) {
    let [[x_min, x_max], [y_min, y_max], [z_min, z_max], [_, _]] = bounds;
    // Simple 3D density estimation
    let grid_size = 100;

    let mut edges = vec![vec![vec![(0.0, 0.0, 0.0); grid_size + 1]; grid_size + 1]; grid_size + 1];
    for i in 0..=grid_size {
        for j in 0..=grid_size {
            for k in 0..=grid_size {
                // log scale m200
                let x_edge =
                    (x_min.ln() + (x_max.ln() - x_min.ln()) * i as f64 / grid_size as f64).exp();
                let y_edge = y_min + (y_max - y_min) * j as f64 / grid_size as f64;
                let z_edge = z_min + (z_max - z_min) * k as f64 / grid_size as f64;
                edges[i][j][k] = (x_edge, y_edge, z_edge);
            }
        }
    }

    let mut density = vec![vec![vec![0.0; grid_size]; grid_size]; grid_size];
    let mut count = 0;
    for i in 0..chain.len() {
        let (x, y, z) = (chain[i][0], chain[i][1], chain[i][2]);

        if &x >= x_max || &x < x_min || &y >= y_max || &y < y_min || &z >= z_max || &z < z_min {
            continue;
        }

        let x_bin =
            ((x.ln() - x_min.ln()) / (x_max.ln() - x_min.ln()) * grid_size as f64).floor() as usize;
        let y_bin = ((y - y_min) / (y_max - y_min) * grid_size as f64).floor() as usize;
        let z_bin = ((z - z_min) / (z_max - z_min) * grid_size as f64).floor() as usize;

        count += 1;
        density[x_bin][y_bin][z_bin] += 1.0;
    }

    let mut density = gaussian_smooth_3d(&density, 1.0);

    //Normalize
    for i in 0..grid_size {
        for j in 0..grid_size {
            for k in 0..grid_size {
                let vol = (edges[i + 1][j][k].0 - edges[i][j][k].0)
                    * (edges[i][j + 1][k].1 - edges[i][j][k].1)
                    * (edges[i][j][k + 1].2 - edges[i][j][k].2);
                density[i][j][k] /= vol * count as f64
            }
        }
    }

    let mut cells = Vec::with_capacity(grid_size.pow(3));
    for i in 0..grid_size {
        for j in 0..grid_size {
            for k in 0..grid_size {
                let x_start = edges[i][j][k].0;
                let x_end = edges[i + 1][j][k].0;
                let y_start = edges[i][j][k].1;
                let y_end = edges[i][j + 1][k].1;
                let z_start = edges[i][j][k].2;
                let z_end = edges[i][j][k + 1].2;

                cells.push((
                    density[i][j][k],
                    (x_start, y_start, z_start),
                    (x_end, y_end, z_end),
                    3.5,
                ));
            }
        }
    }

    cells.sort_by(|cell_a, cell_b| cell_b.0.partial_cmp(&cell_a.0).unwrap());

    const CONTOUR_LEVELS: [f64; 6] = [0.031, 0.199, 0.478, 0.739, 0.900, 0.971]; //[0.5, 1, 1.5, 2, 2.5, 3] sigma for 3D gaussian

    let mut cumulative_probability = 0.0;

    for cell in &mut cells {
        let vol = (cell.2.0 - cell.1.0) * (cell.2.1 - cell.1.1) * (cell.2.2 - cell.1.2);
        cumulative_probability += cell.0 * vol;

        for i in 0..CONTOUR_LEVELS.len() {
            if cumulative_probability <= CONTOUR_LEVELS[i] {
                cell.3 = (i as f64 / 2.0) + 0.5;
                break;
            }
        }
    }

    write_cells_to_csv(&cells, "data/cells.csv").unwrap();
}

fn gaussian_kernel(sigma: f64) -> Vec<f64> {
    let radius = (3.0 * sigma).ceil() as isize;
    let mut kernel = Vec::new();
    let mut sum = 0.0;

    for i in -radius..=radius {
        let x = i as f64;
        let v = (-x * x / (2.0 * sigma * sigma)).exp();
        kernel.push(v);
        sum += v;
    }

    for v in kernel.iter_mut() {
        *v /= sum;
    }

    kernel
}

fn blur_x_3d(field: &Vec<Vec<Vec<f64>>>, sigma: f64) -> Vec<Vec<Vec<f64>>> {
    let kernel = gaussian_kernel(sigma);
    let r = (kernel.len() / 2) as isize;

    let nx = field.len();
    let ny = field[0].len();
    let nz = field[0][0].len();

    let mut out = vec![vec![vec![0.0; nz]; ny]; nx];

    for x in 0..nx {
        for y in 0..ny {
            for z in 0..nz {
                let mut sum = 0.0;
                for k in -r..=r {
                    let xx = (x as isize + k).clamp(0, (nx - 1) as isize) as usize;
                    sum += field[xx][y][z] * kernel[(k + r) as usize];
                }
                out[x][y][z] = sum;
            }
        }
    }

    out
}

fn blur_y_3d(field: &Vec<Vec<Vec<f64>>>, sigma: f64) -> Vec<Vec<Vec<f64>>> {
    let kernel = gaussian_kernel(sigma);
    let r = (kernel.len() / 2) as isize;

    let nx = field.len();
    let ny = field[0].len();
    let nz = field[0][0].len();

    let mut out = vec![vec![vec![0.0; nz]; ny]; nx];

    for x in 0..nx {
        for y in 0..ny {
            for z in 0..nz {
                let mut sum = 0.0;
                for k in -r..=r {
                    let yy = (y as isize + k).clamp(0, (ny - 1) as isize) as usize;
                    sum += field[x][yy][z] * kernel[(k + r) as usize];
                }
                out[x][y][z] = sum;
            }
        }
    }

    out
}

fn blur_z_3d(field: &Vec<Vec<Vec<f64>>>, sigma: f64) -> Vec<Vec<Vec<f64>>> {
    let kernel = gaussian_kernel(sigma);
    let r = (kernel.len() / 2) as isize;

    let nx = field.len();
    let ny = field[0].len();
    let nz = field[0][0].len();

    let mut out = vec![vec![vec![0.0; nz]; ny]; nx];

    for x in 0..nx {
        for y in 0..ny {
            for z in 0..nz {
                let mut sum = 0.0;
                for k in -r..=r {
                    let zz = (z as isize + k).clamp(0, (nz - 1) as isize) as usize;
                    sum += field[x][y][zz] * kernel[(k + r) as usize];
                }
                out[x][y][z] = sum;
            }
        }
    }

    out
}

fn gaussian_smooth_3d(field: &Vec<Vec<Vec<f64>>>, sigma: f64) -> Vec<Vec<Vec<f64>>> {
    let tmp = blur_x_3d(field, sigma);
    let tmp = blur_y_3d(&tmp, sigma);
    blur_z_3d(&tmp, sigma)
}
