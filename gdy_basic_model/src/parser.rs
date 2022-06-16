pub mod msi_parser {
    use std::{fs::read_to_string, path::Path};

    use nalgebra::{Matrix3, Point3, Vector3};
    use regex::{Captures, Regex};

    use crate::{Atom, Lattice, Molecule};

    pub fn parse_atom(text: &str) -> Vec<Atom> {
        let atom_re = Regex::new(
            r#"ACL "([0-9]+) ([a-zA-Z]+).*
.*Label ".*
.*XYZ \(([0-9.e-]+) ([0-9.e-]+) ([0-9.e-]+).*
.*Id ([0-9]+)"#,
        )
        .unwrap();
        assert!(atom_re.is_match(text));
        let mut atom_struct_vec: Vec<Atom> = vec![];
        for cap in atom_re.captures_iter(text) {
            let element: String = cap[2].to_string();
            let element_id: u8 = cap[1].to_string().parse::<u8>().unwrap();
            let point: Point3<f64> = na::point![
                cap[3].to_string().parse::<f64>().unwrap(),
                cap[4].to_string().parse::<f64>().unwrap(),
                cap[5].to_string().parse::<f64>().unwrap()
            ];
            let atom_id: u8 = cap[6].to_string().parse::<u8>().unwrap();
            atom_struct_vec.push(Atom::new(element, element_id, point, atom_id));
        }
        atom_struct_vec
    }

    pub fn parse_lattice_vectors(text: &str) -> Matrix3<f64> {
        let lattice_vec_re = Regex::new(
            r#".*A3 \(([0-9e. -]+)\)\)
.*B3 \(([0-9e. -]+)\)\)
.*C3 \(([0-9e. -]+)\)\)"#,
        )
        .unwrap();
        assert!(lattice_vec_re.is_match(text));
        let match_result: Captures = lattice_vec_re.captures(text).unwrap();
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
        Matrix3::from_columns(&lattice_vectors)
    }
    pub fn parse_lattice(filename: &str) -> Lattice {
        let contents = read_to_string(filename).expect(&format!(
            "Something went wrong in reading file {}",
            filename
        ));
        let mol_name: String = Path::new(filename)
            .file_stem()
            .expect("Failed getting file_stem")
            .to_str()
            .unwrap()
            .to_string();
        let atoms: Vec<Atom> = parse_atom(&contents);
        let mol: Molecule = Molecule::new(mol_name, atoms);
        let lat_vectors: Matrix3<f64> = parse_lattice_vectors(&contents);
        let lattice: Lattice = Lattice::new(mol, lat_vectors, vec![73, 74, 75], None);
        lattice
    }
}
