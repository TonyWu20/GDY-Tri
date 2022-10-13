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
    coord_atom_nums: usize,
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
        coord_atom_nums: usize,
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
    fn is_face_up(&self) -> Result<bool, String> {
        let cd_z = self.atoms_vec.get_atom_by_id(self.coord_atom_ids[0])?.xyz()[2];
        let up_z = self.atoms_vec.get_atom_by_id(self.upper_atom_id)?.xyz()[2];
        if cd_z < up_z {
            Ok(true)
        } else {
            Ok(false)
        }
    }
    pub fn make_upright(&mut self) -> Result<(), String> {
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
            let rot_quatd = if rot_axis_stem_angle < PI / 2.0 {
                UnitQuaternion::from_axis_angle(&Unit::new_normalize(stem_vector), rotate_angle)
            } else {
                UnitQuaternion::from_axis_angle(
                    &Unit::new_normalize(stem_vector.scale(-1.0)),
                    rotate_angle,
                )
            };
            self.atoms_vec.rotate(rot_quatd);
        } else {
            let z_axis: Vector3<f64> = Vector3::from_vec(vec![0.0, 0.0, 1.0]);
            let angle = stem_vector.angle(&z_axis);
            let rot_axis = Unit::new_normalize(stem_vector.cross(&z_axis));
            let rot_quatd = UnitQuaternion::from_axis_angle(&rot_axis, angle);
            self.atoms_vec.rotate(rot_quatd);
        }
        if self.is_face_up()? == false {
            let stem_vector = self.get_stem_vector()?;
            let invert_quatd =
                UnitQuaternion::from_axis_angle(&Unit::new_normalize(stem_vector), PI);
            self.atoms_vec.rotate(invert_quatd);
            Ok(())
        } else {
            Ok(())
        }
    }

    pub fn atoms_vec(&self) -> &[Atom] {
        self.atoms_vec.as_ref()
    }

    pub fn mol_name(&self) -> &str {
        self.mol_name.as_ref()
    }

    pub fn coord_atom_nums(&self) -> usize {
        self.coord_atom_nums
    }

    pub fn coord_atom_ids(&self) -> &[u32] {
        self.coord_atom_ids.as_ref()
    }

    pub fn stem_atom_ids(&self) -> [u32; 2] {
        self.stem_atom_ids
    }

    pub fn plane_atom_ids(&self) -> [u32; 3] {
        self.plane_atom_ids
    }

    pub fn vertical(&self) -> bool {
        self.vertical
    }

    pub fn symmetric(&self) -> bool {
        self.symmetric
    }

    pub fn upper_atom_id(&self) -> u32 {
        self.upper_atom_id
    }

    pub fn pathway_name(&self) -> &str {
        self.pathway_name.as_ref()
    }

    pub fn atoms_vec_mut(&mut self) -> &mut Vec<Atom> {
        &mut self.atoms_vec
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
