use std::f64::consts::PI;

use na::{vector, Unit, UnitQuaternion, Vector3};

use crate::{
    atom::{Atom, AtomArray},
    Export, Transformation,
};
#[derive(Debug, Clone)]
pub struct Adsorbate {
    mol_name: String,
    atoms_vec: Vec<Atom>,
    coord_atom_nums: u32,
    coord_atom_ids: Vec<u32>,
    stem_atom_ids: [u32; 2],
    plane_atom_ids: [u32; 3],
    vertical: bool,
    symmetric: bool,
    upper_atom_id: u32,
    pathway_name: String,
}

impl Adsorbate {
    pub fn new(
        mol_name: String,
        atoms_vec: Vec<Atom>,
        coord_atom_nums: u32,
        coord_atom_ids: Vec<u32>,
        stem_atom_ids: [u32; 2],
        plane_atom_ids: [u32; 3],
        vertical: bool,
        symmetric: bool,
        upper_atom_id: u32,
        pathway_name: String,
    ) -> Self {
        Self {
            mol_name,
            atoms_vec,
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

    pub fn set_adsorate_name(&mut self, new_name: &str) {
        self.mol_name = new_name.to_string();
    }
    pub fn get_stem_vector(&self) -> Result<Vector3<f64>, String> {
        self.atoms_vec
            .get_vector_ab(self.stem_atom_ids[0], self.stem_atom_ids[1])
    }
    pub fn get_plane_normal(&self) -> Result<Vector3<f64>, String> {
        let ba = self
            .atoms_vec
            .get_vector_ab(self.plane_atom_ids[0], self.plane_atom_ids[1])?;
        let ca = self
            .atoms_vec
            .get_vector_ab(self.plane_atom_ids[0], self.plane_atom_ids[2])?;
        let plane_normal = ba.cross(&ca).normalize();
        Ok(plane_normal)
    }
    pub fn make_upright(&mut self) {
        let stem_vector: Vector3<f64> = self
            .get_stem_vector()
            .unwrap_or_else(|_| panic!("Failed to get stem vector! Adsorbate: {}", self.mol_name));
        if self.vertical {
            let plane_normal: Vector3<f64> = self.get_plane_normal().unwrap_or_else(|_| {
                panic!("Failed to get plane normal! Adsorbate: {}", self.mol_name)
            });
            let plane_normal_xy_proj: Vector3<f64> = vector![plane_normal[0], plane_normal[1], 0.0];
            let rotate_angle = plane_normal.angle(&plane_normal_xy_proj);
            let rot_axis = plane_normal.cross(&plane_normal_xy_proj);
            let rot_axis_stem_angle = rot_axis.angle(&stem_vector);
            let rotation_quaternion = if rot_axis_stem_angle < PI / 2.0 {
                UnitQuaternion::from_axis_angle(&Unit::new_normalize(stem_vector), rotate_angle)
            } else {
                UnitQuaternion::from_axis_angle(
                    &Unit::new_normalize(stem_vector.scale(-1.0)),
                    rotate_angle,
                )
            };
            self.atoms_vec.rotate(rotation_quaternion);
        }
        todo!();
    }
}

impl Export for Adsorbate {
    fn format_output(&self) -> String {
        self.atoms_vec
            .iter()
            .map(|x| x.format_output())
            .collect::<Vec<String>>()
            .join("")
    }
}

impl Transformation for Adsorbate {
    fn rotate(&mut self, rotate_quatd: na::UnitQuaternion<f64>) {
        self.atoms_vec.rotate(rotate_quatd)
    }
    fn translate(&mut self, translate_matrix: na::Translation<f64, 3>) {
        self.atoms_vec.translate(translate_matrix)
    }
}
