use std::{fmt::Display, fs, io::Error};

use ndarray::Array1;

use crate::misc_methods::meshgrid;

use super::export_format::{ExportGnuData, MyOutput};
/// Preparing data for plotti&ng in various form.

/// Mesh-like data in polar coordinate.
pub struct PolarMeshData<T>
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

impl<T> PolarMeshData<T>
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
    /// Returns (rads, thetas)
    fn get_mesh(&self) -> (Vec<f64>, Vec<f64>) {
        let (rad_start, rad_end) = self.rad_range;
        let rad_space: Vec<f64> = Array1::linspace(rad_start, rad_end, self.rad_steps).to_vec();
        let (theta_start, theta_end) = self.theta_range;
        let theta_space: Vec<f64> =
            Array1::linspace(theta_start, theta_end, self.num_sectors).to_vec();
        meshgrid(theta_space, rad_space)
    }
}

impl<T> ExportGnuData for PolarMeshData<T>
where
    T: Display + Copy + MyOutput,
{
    fn to_gnu_data<P: AsRef<std::path::Path>>(&self, filename: P) -> Result<(), Error> {
        assert_eq!(
            self.rad_steps * self.num_sectors,
            self.raw_data.len(),
            "In consistent shape of data: rad * theta = {}, data.len() = {}",
            self.rad_steps * self.num_sectors,
            self.raw_data.len()
        );
        let (theta_col, rad_col) = self.get_mesh();
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
        fs::write(filename, data_text.join("\n"))
    }
}
