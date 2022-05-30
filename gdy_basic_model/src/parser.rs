pub mod msi_parser {
    use gdy_model::Atom;
    use nalgebra::{Matrix3, Point3, Vector3};
    use regex::{Captures, Regex};

    pub fn parse_atom(text: &str) -> Vec<Atom> {
        let atom_re = Regex::new(
        r#"ACL "([0-9]+) ([a-zA-Z]+).*\r\n.*Label ".*\r\n.*XYZ \(([0-9.e-]+) ([0-9.e-]+) ([0-9.e-]+).*\r\n.*Id ([0-9]+)"#,
    )
    .unwrap();
        assert!(atom_re.is_match(text));
        let mut atom_struct_vec: Vec<Atom> = vec![];
        for cap in atom_re.captures_iter(text) {
            let element: String = cap[2].to_string();
            let element_id: u32 = cap[1].to_string().parse::<u32>().unwrap();
            let point: Point3<f64> = na::point![
                cap[3].to_string().parse::<f64>().unwrap(),
                cap[4].to_string().parse::<f64>().unwrap(),
                cap[5].to_string().parse::<f64>().unwrap()
            ];
            let atom_id: u32 = cap[6].to_string().parse::<u32>().unwrap();
            atom_struct_vec.push(Atom::new(element, element_id, point, atom_id));
        }
        atom_struct_vec
    }

    pub fn parse_lattice_vectors(text: &str) -> Matrix3<f64> {
        let lattice_vec_re = Regex::new(
            r#".*A3 \(([0-9e. -]+)\)\)\r\n.*B3 \(([0-9e. -]+)\)\)\r\n.*C3 \(([0-9e. -]+)\)\)\r\n"#,
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
}
