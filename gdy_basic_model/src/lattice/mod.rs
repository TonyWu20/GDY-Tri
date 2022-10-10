use periodic_table_on_an_enum as pt_enum;
use std::collections::HashSet;

use na::{Matrix3, UnitQuaternion, Vector, Vector3};

use crate::{
    atom::{Atom, AtomArray},
    Export, Transformation,
};

pub struct Lattice {
    lattice_name: String,
    atoms_vec: Vec<Atom>,
    lattice_vectors: Matrix3<f64>,
    metal_sites: Vec<u32>,
    adsorbate: Option<String>,
    sorted: bool,
}

impl Lattice {
    pub fn new(
        lattice_name: String,
        atoms_vec: Vec<Atom>,
        lattice_vectors: Matrix3<f64>,
        metal_sites: Vec<u32>,
        adsorbate: Option<String>,
        sorted: bool,
    ) -> Self {
        Self {
            lattice_name,
            atoms_vec,
            lattice_vectors,
            metal_sites,
            adsorbate,
            sorted,
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
    // Get element list sorted by atomic number
    pub fn get_element_list(&self) -> Vec<String> {
        let mut elm_list: Vec<String> = vec![];
        elm_list.extend(
            self.atoms_vec
                .iter()
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
                self.atoms_vec
                    .get_atom_by_id(*metal_id)
                    .unwrap()
                    .element_name()
                    .to_string()
            })
            .collect::<Vec<String>>();
        // Because we have only 3 metal elements
        let new_name = format!(
            "GDY_{}_{}_{}",
            metal_names[0], metal_names[1], metal_names[2]
        );
        self.set_lattice_name(new_name);
    }

    pub fn lattice_name(&self) -> &str {
        self.lattice_name.as_ref()
    }

    pub fn set_lattice_name(&mut self, lattice_name: String) {
        self.lattice_name = lattice_name;
    }

    pub fn atoms_vec(&self) -> &[Atom] {
        self.atoms_vec.as_ref()
    }

    pub fn atoms_vec_mut(&mut self) -> &mut Vec<Atom> {
        &mut self.atoms_vec
    }

    pub fn set_sorted(&mut self, sorted: bool) {
        self.sorted = sorted;
    }
    pub fn sort_atoms_by_elements(&mut self) {
        self.atoms_vec_mut().sort();
        self.set_sorted(true);
    }
    pub fn lattice_vector_str(&self) -> (String, String) {
        let vectors = self.get_lattice_vectors();
        let mut vectors_string = String::new();
        vectors.column_iter().for_each(|col| {
            vectors_string.push_str(&format!("{:24.18}{:24.18}{:24.18}\n", col.x, col.y, col.z));
        });
        ("LATTICE_CART".to_string(), vectors_string)
    }
    pub fn fractional_coord_matrix(&self) -> Matrix3<f64> {
        let lattice_vectors = self.get_lattice_vectors();
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

    pub fn sorted(&self) -> bool {
        self.sorted
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
        let atom_strings: String = self.atoms_vec.format_output();
        let contents: String = format!("{}{})", headers, atom_strings);
        contents
    }
}
impl Transformation for Lattice {
    fn rotate(&mut self, rotate_quatd: na::UnitQuaternion<f64>) {
        self.atoms_vec.rotate(rotate_quatd);
        let rotation_matrix = rotate_quatd.to_rotation_matrix();
        let new_lat_vec: Matrix3<f64> = rotation_matrix * self.get_lattice_vectors();
        self.set_lattice_vectors(new_lat_vec);
    }
    fn translate(&mut self, translate_matrix: na::Translation<f64, 3>) {
        self.atoms_vec.translate(translate_matrix);
    }
}
