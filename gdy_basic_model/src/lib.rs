#![allow(dead_code)]
pub mod editor;
pub mod external_info;
pub mod param_writer;
pub mod parser;
pub mod test;

extern crate nalgebra as na;

use ::core::fmt;

use external_info::element_table::Element;
use periodic_table_on_an_enum as pt_enum;
use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet},
};

use na::*;

// Shared behaviour
pub trait Export {
    fn format_output(&self) -> String;
}
trait Transformation {
    fn rotate(&mut self, rotate_quatd: na::UnitQuaternion<f64>);
    fn translate(&mut self, translate_matrix: na::Translation<f64, 3>);
}
// trait ends

// Atom
#[derive(Debug, Clone)]
pub struct Atom {
    element_name: String,
    element_id: u8,
    xyz: Point3<f64>,
    atom_id: u8,
    // lcao: u8,
    // mass: f64,
    // pot: String,
    // spin: u8,
}

impl Atom {
    pub fn new(
        element_name: String,
        element_id: u8,
        xyz: Point3<f64>,
        atom_id: u8,
        // lcao: u8,
        // mass: f64,
        // pot: String,
        // spin: u8,
    ) -> Self {
        Self {
            element_name,
            element_id,
            xyz,
            atom_id,
            // lcao,
            // mass,
            // pot,
            // spin,
        }
    }

    pub fn element_name(&self) -> &str {
        &self.element_name
    }
    pub fn set_element_name(&mut self, new_name: &str) {
        self.element_name = new_name.to_string();
    }
    pub fn element_id(&self) -> u8 {
        self.element_id
    }
    pub fn set_element_id(&mut self, new_id: u8) {
        self.element_id = new_id;
    }
    pub fn xyz(&self) -> &Point3<f64> {
        &self.xyz
    }
    pub fn set_xyz(&mut self, new_xyz: Point3<f64>) {
        self.xyz = new_xyz;
    }
    pub fn atom_id(&self) -> u8 {
        self.atom_id
    }
    // pub fn lcao(&self) -> u8 {
    //     self.lcao
    // }
    // pub fn spin(&self) -> u8 {
    //     self.spin
    // }
    // pub fn potential_file(&self) -> String {
    //     self.pot.to_string()
    // }
    // pub fn mass(&self) -> f64 {
    //     self.mass
    // }
}

impl Export for Atom {
    fn format_output(&self) -> String {
        let msi_output: String = format!(
            r#"  ({item_id} Atom
    (A C ACL "{elm_id} {elm}")
    (A C Label "{elm}")
    (A D XYZ ({x:.12} {y:.12} {z:.12}))
    (A I Id {atom_id})
  )
"#,
            item_id = self.atom_id() + 1,
            elm_id = self.element_id(),
            elm = self.element_name(),
            x = self.xyz().x,
            y = self.xyz().y,
            z = self.xyz().z,
            atom_id = self.atom_id(),
        );
        msi_output
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Element: {}\nElement ID: {}\ncoord: {}\nAtom ID: {}",
            self.element_name, self.element_id, self.xyz, self.atom_id
        )
    }
}

impl Ord for Atom {
    fn cmp(&self, other: &Self) -> Ordering {
        self.atom_id.cmp(&other.atom_id)
    }
}

impl PartialOrd for Atom {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Atom {
    fn eq(&self, other: &Self) -> bool {
        self.atom_id == other.atom_id
    }
}

impl Eq for Atom {}
// End Atom

#[derive(Debug, Clone)]
pub struct Molecule {
    pub mol_name: String,
    vector_atoms: Vec<Atom>,
}

impl Molecule {
    pub fn new(mol_name: String, vector_atoms: Vec<Atom>) -> Self {
        Self {
            mol_name,
            vector_atoms,
        }
    }
    /// Return an atom with the given id, starting from 1 (caution!)
    pub fn get_atom_by_id(&self, atom_id: u8) -> Option<&Atom> {
        self.vector_atoms.get(atom_id as usize - 1)
    }
    pub fn get_mut_atom_by_id(&mut self, atom_id: u8) -> Option<&mut Atom> {
        self.vector_atoms.get_mut(atom_id as usize - 1)
    }
    /// Push a new atom to the atom vectors.
    pub fn append_atom(&mut self, new_atom: Atom) {
        self.vector_atoms.push(new_atom);
    }
    /// Return the current number of atoms of the molecule.
    pub fn number_of_atoms(&self) -> usize {
        self.vector_atoms.len()
    }
    pub fn atoms_iterator(&self) -> std::slice::Iter<Atom> {
        self.vector_atoms.iter()
    }
    pub fn get_vector_ab(&self, a_id: u8, b_id: u8) -> Vector3<f64> {
        let atom_a: &Atom = self.get_atom_by_id(a_id).unwrap();
        let atom_a_xyz = atom_a.xyz();
        let atom_b: &Atom = self.get_atom_by_id(b_id).unwrap();
        let atom_b_xyz = atom_b.xyz();
        atom_b_xyz - atom_a_xyz
    }
    pub fn set_molecule_name(&mut self, new_name: &str) {
        self.mol_name = new_name.to_string();
    }
}

impl Export for Molecule {
    fn format_output(&self) -> String {
        self.atoms_iterator()
            .map(|x| x.format_output())
            .collect::<Vec<String>>()
            .join("")
    }
}

#[derive(Debug, Clone)]
pub struct Adsorbate<'a> {
    pub molecule: &'a Molecule,
    coord_atom_nums: u8,
    coord_atom_ids: Vec<u8>,
    stem_atom_ids: [u8; 2],
    plane_atom_ids: [u8; 3],
    vertical: bool,
    symmetric: bool,
    upper_atom_id: u8,
    pathway_name: String,
}

impl<'a> Adsorbate<'a> {
    pub fn new(
        molecule: &'a Molecule,
        coord_atom_nums: u8,
        coord_atom_ids: Vec<u8>,
        stem_atom_ids: [u8; 2],
        plane_atom_ids: [u8; 3],
        vertical: bool,
        symmetric: bool,
        upper_atom_id: u8,
        pathway_name: String,
    ) -> Self {
        Self {
            molecule,
            coord_atom_nums,
            coord_atom_ids,
            stem_atom_ids,
            plane_atom_ids,
            vertical,
            symmetric,
            upper_atom_id,
            pathway_name,
        }
    }
    pub fn get_stem_vector(&self) -> Vector3<f64> {
        self.molecule
            .get_vector_ab(self.stem_atom_ids[0], self.stem_atom_ids[1])
    }
    pub fn get_plane_normal(&self) -> Vector3<f64> {
        let ba = self
            .molecule
            .get_vector_ab(self.plane_atom_ids[0], self.plane_atom_ids[1]);
        let ca = self
            .molecule
            .get_vector_ab(self.plane_atom_ids[0], self.plane_atom_ids[2]);
        let plane_normal = ba.cross(&ca).normalize();
        plane_normal
    }
    pub fn make_upright(&mut self) {
        todo!();
    }
}

pub struct Lattice {
    pub molecule: Molecule,
    lattice_vectors: Matrix3<f64>,
    metal_sites: Vec<u32>,
    adsorbate: Option<String>,
}

impl Lattice {
    pub fn new(
        molecule: Molecule,
        lattice_vectors: Matrix3<f64>,
        metal_sites: Vec<u32>,
        adsorbate: Option<String>,
    ) -> Self {
        Self {
            molecule,
            lattice_vectors,
            metal_sites,
            adsorbate,
        }
    }
    pub fn get_mut_molecule(&mut self) -> &mut Molecule {
        &mut self.molecule
    }
    pub fn get_lattice_vectors(&self) -> &Matrix3<f64> {
        &self.lattice_vectors
    }
    pub fn set_lattice_vectors(&mut self, new_lattice_vec: Matrix3<f64>) {
        self.lattice_vectors = new_lattice_vec;
    }
    pub fn get_metal_sites(&self) -> &Vec<u32> {
        &self.metal_sites
    }
    // Get element list sorted by atomic number
    pub fn get_element_list(&self) -> Vec<String> {
        let mut elm_list: Vec<String> = vec![];
        elm_list.extend(
            self.molecule
                .atoms_iterator()
                .map(|atom| atom.element_name().to_string())
                .collect::<Vec<String>>()
                .drain(..)
                .collect::<HashSet<String>>()
                .into_iter(),
        );
        elm_list.sort_unstable_by(|a, b| {
            pt_enum::Element::from_symbol(&a)
                .unwrap()
                .get_atomic_number()
                .cmp(
                    &pt_enum::Element::from_symbol(&b)
                        .unwrap()
                        .get_atomic_number(),
                )
        });
        elm_list
    }
    pub fn get_adsorbate_name(&self) -> Option<&str> {
        match &self.adsorbate {
            Some(x) => Some(x.as_str()),
            None => None,
        }
    }
    pub fn set_adsorbate_name(&mut self, ads_name: String) {
        self.adsorbate = Some(ads_name);
    }
    pub fn rotate_to_standard_orientation(&mut self) {
        let x_axis: Vector3<f64> = Vector::x();
        let a_vec = &self.get_lattice_vectors().column(0);
        let a_to_x_angle: f64 = a_vec.angle(&x_axis);
        if a_to_x_angle == 0.0 {
            return;
        }
        let rot_axis = a_vec.cross(&x_axis).normalize();
        let rot_quatd: UnitQuaternion<f64> = UnitQuaternion::new(rot_axis * a_to_x_angle);
        self.rotate(rot_quatd);
    }
    pub fn update_base_name(&mut self) {
        // Collect all metal's symbols
        let metal_names: Vec<String> = self
            .metal_sites
            .iter()
            .map(|metal_id| -> String {
                self.molecule
                    .get_atom_by_id(*metal_id as u8)
                    .unwrap()
                    .element_name
                    .to_string()
            })
            .collect::<Vec<String>>();
        // Because we have only 3 metal elements
        let new_name = format!(
            "GDY_{}_{}_{}",
            metal_names[0], metal_names[1], metal_names[2]
        );
        self.molecule.set_molecule_name(&new_name);
    }
}

impl Export for Lattice {
    fn format_output(&self) -> String {
        let headers: String = concat!(
            "# MSI CERIUS2 DataModel File Version 4 0\n",
            "(1 Model\n",
            "  (A I CRY/DISPLAY (192 256))\n",
            "  (A I PeriodicType 100)\n",
            "  (A C SpaceGroup \"1 1\")\n",
            "  (A D A3 (16.39518593025 -9.465765010246 0))\n",
            "  (A D B3 (0 18.93153002049 0))\n",
            "  (A D C3 (0 0 9.999213039981))\n",
            "  (A D CRY/TOLERANCE 0.05)\n"
        )
        .to_string();
        let atom_strings: String = self.molecule.format_output();
        let contents: String = format!("{}{})", headers, atom_strings);
        contents
    }
}

pub struct Cell<'a> {
    pub lattice: &'a mut Lattice,
    sorted: bool,
}

impl<'a> Cell<'a> {
    pub fn new(lattice: &'a mut Lattice, sorted: bool) -> Self {
        Self { lattice, sorted }
    }
    pub fn sort_atoms_by_elements(&mut self) {
        self.lattice.molecule.vector_atoms.sort();
        self.sorted = true;
    }
    // Accept tuple which has name + content
    fn write_block(&self, block: (String, String)) -> String {
        let (block_name, content) = block;
        format!(
            "%%BlOCK {}\n{}%%ENDBLOCK {}\n\n",
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
        ("SPECIES_LCAO".to_string(), lcao_strings)
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

impl Transformation for Molecule {
    fn rotate(&mut self, rotate_quatd: na::UnitQuaternion<f64>) {
        self.vector_atoms
            .iter_mut()
            .for_each(|atom: &mut Atom| atom.set_xyz(rotate_quatd.transform_point(atom.xyz())));
    }
    fn translate(&mut self, translate_matrix: na::Translation<f64, 3>) {
        self.vector_atoms
            .iter_mut()
            .for_each(|atom| atom.set_xyz(translate_matrix.transform_point(atom.xyz())));
    }
}

impl Transformation for Lattice {
    fn rotate(&mut self, rotate_quatd: na::UnitQuaternion<f64>) {
        self.molecule.rotate(rotate_quatd);
        let rotation_matrix = rotate_quatd.to_rotation_matrix();
        let new_lat_vec: Matrix3<f64> = rotation_matrix * self.get_lattice_vectors();
        self.set_lattice_vectors(new_lat_vec);
    }
    fn translate(&mut self, translate_matrix: na::Translation<f64, 3>) {
        self.molecule.translate(translate_matrix);
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
