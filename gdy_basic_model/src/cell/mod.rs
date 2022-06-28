use std::collections::HashMap;

use na::Matrix3;

use crate::{external_info::element_table::Element, lattice::Lattice};

pub struct Cell<'a> {
    pub lattice: &'a mut Lattice,
    sorted: bool,
}

impl<'a> Cell<'a> {
    pub fn new(lattice: &'a mut Lattice, sorted: bool) -> Self {
        lattice.rotate_to_standard_orientation();
        Self { lattice, sorted }
    }
    pub fn sort_atoms_by_elements(&mut self) {
        self.lattice.molecule.vector_atoms.sort();
        self.sorted = true;
    }
    pub fn get_cell_name(&self) -> String {
        self.lattice.molecule.mol_name.to_string()
    }
    // Accept tuple which has name + content
    fn write_block(&self, block: (String, String)) -> String {
        let (block_name, content) = block;
        format!(
            "%BlOCK {}\n{}%ENDBLOCK {}\n\n",
            block_name, content, block_name
        )
    }
    fn lattice_vector_str(&self) -> (String, String) {
        let vectors = self.lattice.get_lattice_vectors();
        let mut vectors_string = String::new();
        vectors.column_iter().for_each(|col| {
            vectors_string.push_str(&format!("{:24.18}{:24.18}{:24.18}\n", col.x, col.y, col.z));
        });
        ("LATTICE_CART".to_string(), vectors_string)
    }
    fn positions_str(&self, element_info: &HashMap<String, Element>) -> (String, String) {
        assert!(self.sorted == true);
        let mut pos_strings = String::new();
        self.lattice.molecule.vector_atoms.iter().for_each(|atom| {
            let frac_coord = fractional_coord_matrix(&self.lattice) * atom.xyz();
            let atom_info = element_info.get(atom.element_name()).expect(&format!(
                "Element {} not in element hash table!",
                atom.element_name()
            ));
            if atom_info.spin > 0 {
                let line = format!(
                    "{:>3}{:20.16}{:20.16}{:20.16} SPIN={:14.10}\n",
                    atom.element_name(),
                    frac_coord.x,
                    frac_coord.y,
                    frac_coord.z,
                    atom_info.spin as f64
                );
                pos_strings.push_str(&line);
            } else {
                let line = format!(
                    "{:>3}{:20.16}{:20.16}{:20.16}\n",
                    atom.element_name(),
                    frac_coord.x,
                    frac_coord.y,
                    frac_coord.z,
                );
                pos_strings.push_str(&line);
            }
        });
        ("POSITIONS_FRAC".to_string(), pos_strings)
    }
    fn kpoints_list_str(&self) -> (String, String) {
        ("KPOINTS_LIST".to_string(), "   0.0000000000000000   0.0000000000000000   0.0000000000000000       1.000000000000000
".to_string())
    }
    fn misc_str(&self) -> String {
        let options_1: String = format!(
            "FIX_ALL_CELL : true\n\nFIX_COM : false\n{}",
            self.write_block(("IONIC_CONSTRAINTS".to_string(), "".to_string()))
        );
        let external_efield = self.write_block((
            "EXTERNAL_EFIELD".to_string(),
            "    0.0000000000     0.0000000000     0.0000000000\n".to_string(),
        ));
        let external_pressure = self.write_block((
            "EXTERNAL_PRESSURE".to_string(),
            r#"    0.0000000000    0.0000000000    0.0000000000
                    0.0000000000    0.0000000000
                                    0.0000000000
"#
            .to_string(),
        ));
        let mut misc = String::new();
        misc.push_str(&options_1);
        misc.push_str(&external_efield);
        misc.push_str(&external_pressure);
        misc
    }
    fn species_mass_str(&self, element_info: &HashMap<String, Element>) -> (String, String) {
        let element_list = self.lattice.get_element_list();
        let mut mass_strings = String::new();
        element_list.iter().for_each(|elm| {
            let mass: f64 = element_info.get(elm).unwrap().mass;
            let mass_line: String = format!("{:>8}{:17.10}\n", elm, mass);
            mass_strings.push_str(&mass_line);
        });
        ("SPECIES_MASS".to_string(), mass_strings)
    }
    fn species_pot_str(&self, element_info: &HashMap<String, Element>) -> (String, String) {
        let element_list = self.lattice.get_element_list();
        let mut pot_strings = String::new();
        element_list.iter().for_each(|elm| {
            let pot_file: &String = &element_info.get(elm).unwrap().pot;
            let pot_line: String = format!("{:>8}  {}\n", elm, pot_file);
            pot_strings.push_str(&pot_line);
        });
        ("SPECIES_POT".to_string(), pot_strings)
    }
    fn species_lcao_str(&self, element_info: &HashMap<String, Element>) -> (String, String) {
        let element_list = self.lattice.get_element_list();
        let mut lcao_strings = String::new();
        element_list.iter().for_each(|elm| {
            let lcao_state = &element_info.get(elm).unwrap().lcao;
            let lcao_line: String = format!("{:>8}{:9}\n", elm, lcao_state);
            lcao_strings.push_str(&lcao_line);
        });
        ("SPECIES_LCAO_STATES".to_string(), lcao_strings)
    }
    pub fn format_output(&self, element_info: &HashMap<String, Element>) -> String {
        let mut content = String::new();
        let block_lat_vec = self.write_block(self.lattice_vector_str());
        content.push_str(&block_lat_vec);
        let block_pos = self.write_block(self.positions_str(&element_info));
        content.push_str(&block_pos);
        let block_kpoints_list = self.write_block(self.kpoints_list_str());
        content.push_str(&block_kpoints_list);
        let block_misc = self.misc_str();
        content.push_str(&block_misc);
        let block_mass = self.write_block(self.species_mass_str(&element_info));
        content.push_str(&block_mass);
        let block_pot = self.write_block(self.species_pot_str(&element_info));
        content.push_str(&block_pot);
        let block_lcao = self.write_block(self.species_lcao_str(&element_info));
        content.push_str(&block_lcao);
        content
    }
}
pub fn fractional_coord_matrix(lattice: &Lattice) -> Matrix3<f64> {
    let lattice_vectors = lattice.get_lattice_vectors();
    let vec_a = lattice_vectors.column(0);
    let vec_b = lattice_vectors.column(1);
    let vec_c = lattice_vectors.column(2);
    let len_a: f64 = vec_a.norm();
    let len_b: f64 = vec_b.norm();
    let len_c: f64 = vec_c.norm();
    let (alpha, beta, gamma) = (
        vec_b.angle(&vec_c),
        vec_a.angle(&vec_c),
        vec_a.angle(&vec_b),
    );
    let vol = vec_a.dot(&vec_b.cross(&vec_c));
    let to_cart = Matrix3::new(
        len_a,
        len_b * gamma.cos(),
        len_c * beta.cos(),
        0.0,
        len_b * gamma.sin(),
        len_c * (alpha.cos() - beta.cos() * gamma.cos()) / gamma.sin(),
        0.0,
        0.0,
        vol / (len_a * len_b * gamma.sin()),
    );
    let to_frac = to_cart.try_inverse().unwrap();
    to_frac
}
