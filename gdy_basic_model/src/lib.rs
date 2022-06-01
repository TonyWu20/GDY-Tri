extern crate nalgebra as na;

use ::core::fmt;

use na::*;
pub struct Atom {
    element_name: String,
    element_id: u32,
    xyz: Point3<f64>,
    atom_id: u32,
}

impl Atom {
    pub fn new(element_name: String, element_id: u32, xyz: Point3<f64>, atom_id: u32) -> Self {
        Self {
            element_name,
            element_id,
            xyz,
            atom_id,
        }
    }
    pub fn element_name(&self) -> &str {
        &self.element_name
    }
    pub fn set_element_name(&mut self, new_name: &str) {
        self.element_name = new_name.to_string();
    }
    pub fn element_id(&self) -> u32 {
        self.element_id
    }
    pub fn set_element_id(&mut self, new_id: u32) {
        self.element_id = new_id;
    }
    pub fn xyz(&self) -> &Point3<f64> {
        &self.xyz
    }
    pub fn set_xyz(&mut self, new_xyz: Point3<f64>) {
        self.xyz = new_xyz;
    }
    pub fn atom_id(&self) -> u32 {
        self.atom_id
    }
    pub fn text(&self) -> String {
        let msi_output: String = format!(
            "  ({item_id} Atom\n    (A C ACL \"{elm_id} {elm}\")\n    (A C Label \"{elm}\")\n    (A D XYZ ({x:.12} {y:.12} {z:.12}))\n    (A I Id {atom_id})\n  )\n",
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
    pub fn get_atom_by_id(&self, atom_id: u32) -> Option<&Atom> {
        self.vector_atoms.get(atom_id as usize - 1)
    }
    pub fn get_mut_atom_by_id(&mut self, atom_id: u32) -> Option<&mut Atom> {
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
    pub fn get_vector_ab(&self, atom_a_id: u32, atom_b_id: u32) -> Vector3<f64> {
        let atom_a: &Atom = self.molecule.get_atom_by_id(atom_a_id).unwrap();
        let atom_a_xyz = atom_a.xyz();
        let atom_b: &Atom = self.molecule.get_atom_by_id(atom_b_id).unwrap();
        let atom_b_xyz = atom_b.xyz();
        atom_b_xyz - atom_a_xyz
    }
    pub fn export_msi(&self) -> String {
        let headers: String = concat!(
            "#MSI CERIUS2 DataModel File Version 4 0\n",
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
        let mut atom_strings: String = "".to_string();
        for atom in self.molecule.atoms_iterator() {
            atom_strings.push_str(&atom.text());
        }
        let contents: String = format!("{}{})", headers, atom_strings);
        contents
    }
}

trait Transformation {
    fn rotate(&mut self, rotate_quatd: na::UnitQuaternion<f64>);
    fn translate(&mut self, translate_matrix: na::Translation<f64, 3>);
}

impl Transformation for Molecule {
    fn rotate(&mut self, rotate_quatd: na::UnitQuaternion<f64>) {
        for atom in self.vector_atoms.iter_mut() {
            let rotated_point: Point3<f64> = rotate_quatd.transform_point(atom.xyz());
            atom.set_xyz(rotated_point);
        }
    }
    fn translate(&mut self, translate_matrix: na::Translation<f64, 3>) {
        for atom in self.vector_atoms.iter_mut() {
            let translated_point: Point3<f64> = translate_matrix.transform_point(atom.xyz());
            atom.set_xyz(translated_point);
        }
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
