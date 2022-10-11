pub mod msi_parser {
    use std::{error::Error, fs::read_to_string, path::Path};

    use nalgebra::{Matrix3, Point3, Vector3};
    use regex::{Captures, Regex};

    use crate::{atom::Atom, lattice::Lattice};

    pub fn parse_atom(text: &str) -> Result<Vec<Atom>, Box<dyn Error>> {
        let atom_re = Regex::new(
            r#"ACL "([0-9]+) ([a-zA-Z]+).*
.*Label ".*
.*XYZ \(([0-9.e-]+) ([0-9.e-]+) ([0-9.e-]+).*
.*Id ([0-9]+)"#,
        )?;
        assert!(atom_re.is_match(text));
        let mut atom_struct_vec: Vec<Atom> = vec![];
        for cap in atom_re.captures_iter(text) {
            let element: String = cap[2].to_string();
            let element_id: u32 = cap[1].to_string().parse::<u32>()?;
            let point: Point3<f64> = na::point![
                cap[3].to_string().parse::<f64>()?,
                cap[4].to_string().parse::<f64>()?,
                cap[5].to_string().parse::<f64>()?
            ];
            let atom_id: u32 = cap[6].to_string().parse::<u32>()?;
            atom_struct_vec.push(Atom::new(element, element_id, point, atom_id));
        }
        Ok(atom_struct_vec)
    }

    pub fn parse_lattice_vectors(text: &str) -> Result<Matrix3<f64>, Box<dyn Error>> {
        let lattice_vec_re = Regex::new(
            r#".*A3 \(([0-9e. -]+)\)\)
.*B3 \(([0-9e. -]+)\)\)
.*C3 \(([0-9e. -]+)\)\)"#,
        )?;
        let match_result: Captures = lattice_vec_re
            .captures(text)
            .unwrap_or_else(|| panic!("Failed to match lattice vectors with regex."));
        let mut lattice_vectors: Vec<Vector3<f64>> = vec![];
        for i in 1..4 {
            let vector = Vector3::from_iterator(
                match_result[i]
                    .to_string()
                    .split_whitespace()
                    .flat_map(str::parse::<f64>)
                    .collect::<Vec<f64>>(),
            );
            lattice_vectors.push(vector);
        }
        Ok(Matrix3::from_columns(&lattice_vectors))
    }
    pub fn parse_lattice(filename: &str) -> Result<Lattice, Box<dyn Error>> {
        let contents = read_to_string(filename)?;
        let lat_name: String = Path::new(filename)
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap()
            .to_string();
        let atoms: Vec<Atom> = parse_atom(&contents)?;
        let lat_vectors: Matrix3<f64> = parse_lattice_vectors(&contents)?;
        let lattice: Lattice =
            Lattice::new(lat_name, atoms, lat_vectors, vec![73, 74, 75], None, false);
        Ok(lattice)
    }
}
