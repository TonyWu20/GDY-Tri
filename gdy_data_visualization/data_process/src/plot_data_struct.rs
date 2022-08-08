use std::{
    fmt::Display,
    fs,
    io::{Error, Write},
};

use ndarray::Array1;

use super::misc_methods::meshgrid;

use super::export_format::{ExportGnuData, MyOutput};
/// Preparing data for plotti&ng in various form.

/// Mesh-like 3d data in cartesian coordinate
pub struct Cart3dMeshData<T>
where
    T: Copy + Display + MyOutput,
{
    num_cols: usize,
    num_rows: usize,
    raw_data: Vec<T>,
}

impl<T> Cart3dMeshData<T>
where
    T: Copy + Display + MyOutput,
{
    pub fn new(num_cols: usize, num_rows: usize, raw_data: Vec<T>) -> Self {
        Self {
            num_cols,
            num_rows,
            raw_data,
        }
    }
}

/// Mesh-like 3d data in polar coordinate.
pub struct Polar3dMeshData<T>
where
    T: Copy + Display + MyOutput,
{
    /// Range of radial coordinate.
    rad_range: (f64, f64),
    /// Radial length of each grid.
    rad_steps: usize,
    /// Range of theta in polar plot.
    theta_range: (f64, f64),
    /// Total number of sectors in the polar mesh.
    num_sectors: usize,
    /// Number of columns of shaped-data. Equal to rad_range / rad_step.
    raw_data: Vec<T>,
}

/// Shared behaviour of `Cart3dMeshData` and `Polar3dMeshData`
trait MeshGrid {
    fn get_meshgrid(&self) -> (Vec<f64>, Vec<f64>);
}

impl<T> MeshGrid for Cart3dMeshData<T>
where
    T: Copy + Display + MyOutput,
{
    fn get_meshgrid(&self) -> (Vec<f64>, Vec<f64>) {
        let col_space = Array1::linspace(0.0, (self.num_cols - 1) as f64, self.num_cols).to_vec();
        let row_space = Array1::linspace(0.0, (self.num_rows - 1) as f64, self.num_rows).to_vec();
        meshgrid(col_space, row_space)
    }
}

impl<T> MeshGrid for Polar3dMeshData<T>
where
    T: Copy + Display + MyOutput,
{
    /// Returns (rads, thetas)
    fn get_meshgrid(&self) -> (Vec<f64>, Vec<f64>) {
        let (rad_start, rad_end) = self.rad_range;
        let rad_space: Vec<f64> = Array1::linspace(rad_start, rad_end, self.rad_steps).to_vec();
        let (theta_start, theta_end) = self.theta_range;
        let theta_space: Vec<f64> =
            Array1::linspace(theta_start, theta_end, self.num_sectors).to_vec();
        meshgrid(rad_space, theta_space)
    }
}

impl<T> Polar3dMeshData<T>
where
    T: Copy + Display + MyOutput,
{
    pub fn new(
        rad_range: (f64, f64),
        rad_steps: usize,
        theta_range: (f64, f64),
        num_sectors: usize,
        raw_data: Vec<T>,
    ) -> Self {
        Self {
            rad_range,
            rad_steps,
            theta_range,
            num_sectors,
            raw_data,
        }
    }
}

impl<T> ExportGnuData for Cart3dMeshData<T>
where
    T: Display + Copy + MyOutput + PartialOrd,
{
    fn to_gnu_data<P: AsRef<std::path::Path>>(&self, filename: P, mode: &str) -> Result<(), Error> {
        assert_eq!(
            self.num_cols * self.num_rows,
            self.raw_data.len(),
            "In consistent shape of data: cols * rows = {}, data.len() = {}",
            self.num_cols * self.num_rows,
            self.raw_data.len()
        );
        let accepted_modes = vec!["w", "write", "a", "append"];
        assert!(
            accepted_modes.contains(&mode),
            "Invalid mode parameter: {}",
            mode
        );
        let (xv, yv) = self.get_meshgrid();
        let min_value = self
            .raw_data
            .iter()
            .min_by(|&a, &b| a.partial_cmp(b).unwrap());
        let max_value = self
            .raw_data
            .iter()
            .max_by(|&a, &b| a.partial_cmp(b).unwrap());
        let min_max_line = format!(
            "# max: {:.8}, min: {:.8}",
            max_value.unwrap(),
            min_value.unwrap()
        );
        let mut lines = vec![min_max_line];
        let mut data_text: Vec<String> = xv
            .iter()
            .zip(yv.iter())
            .zip(self.raw_data.iter())
            .enumerate()
            .map(|(i, data)| {
                let ((x, y), val): ((&f64, &f64), &T) = data;
                format!(
                    "{:.8} {:.8} {}{}",
                    x,
                    y,
                    val.formatted_output(),
                    if (i + 1) % self.num_cols == 0 {
                        "\n"
                    } else {
                        ""
                    }
                )
            })
            .collect();
        lines.append(&mut data_text);
        if mode == "w" || mode == "write" {
            fs::write(filename, lines.join("\n"))
        } else {
            let mut file = fs::OpenOptions::new()
                .append(true)
                .create(true)
                .open(filename)
                .unwrap();
            file.write_all(lines.join("\n").as_bytes())?;
            file.write_all("\n\n".as_bytes())
        }
    }
}

impl<T> ExportGnuData for Polar3dMeshData<T>
where
    T: Display + Copy + MyOutput,
{
    fn to_gnu_data<P: AsRef<std::path::Path>>(&self, filename: P, mode: &str) -> Result<(), Error> {
        assert_eq!(
            self.rad_steps * self.num_sectors,
            self.raw_data.len(),
            "In consistent shape of data: rad * theta = {}, data.len() = {}",
            self.rad_steps * self.num_sectors,
            self.raw_data.len()
        );
        let accepted_modes = vec!["w", "write", "a", "append"];
        assert!(
            accepted_modes.contains(&mode),
            "Invalid mode parameter: {}",
            mode
        );
        let (rad_col, theta_col) = self.get_meshgrid();
        let data_text: Vec<String> = rad_col
            .iter()
            .zip(theta_col.iter())
            .zip(self.raw_data.iter())
            .enumerate()
            .map(|(i, ((&rad, &theta), &val))| {
                format!(
                    "{:.8} {:.8} {}{}",
                    theta,
                    rad,
                    val.formatted_output(),
                    if (i + 1) % 8 == 0 { "\n" } else { "" }
                )
            })
            .collect();
        if mode == "w" || mode == "write" {
            fs::write(filename, data_text.join("\n"))
        } else {
            let mut file = fs::OpenOptions::new()
                .append(true)
                .create(true)
                .open(filename)
                .unwrap();
            file.write_all(data_text.join("\n").as_bytes())?;
            file.write_all("\n\n".as_bytes())
        }
    }
}
