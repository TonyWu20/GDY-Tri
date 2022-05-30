extern crate nalgebra as na;

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
    pub fn element_id(&self) -> u32 {
        self.element_id
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
    pub fn get_atom_by_id(&self, atom_id: u32) -> &Atom {
        &self.vector_atoms[atom_id as usize]
    }
    pub fn append_atom(&mut self, new_atom: Atom) {
        self.vector_atoms.push(new_atom);
    }
    pub fn number_of_atoms(&self) -> usize {
        self.vector_atoms.len()
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
