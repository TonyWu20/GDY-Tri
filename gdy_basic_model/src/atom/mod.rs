use crate::{Export, Transformation};
use ::core::fmt;
use std::cmp::Ordering;

use na::{Point3, Vector3};
// Atom
#[derive(Debug, Clone)]
pub struct Atom {
    element_name: String,
    element_id: u32,
    xyz: Point3<f64>,
    atom_id: u32,
}

pub trait AtomArray {
    fn get_atom_by_id(&self, atom_id: u32) -> Result<&Atom, String>;
    fn get_mut_atom_by_id(&mut self, atom_id: u32) -> Result<&mut Atom, String>;
    fn append_atom(&mut self, new_atom: Atom);
    fn number_of_atoms(&self) -> usize;
    fn get_vector_ab(&self, a_id: u32, b_id: u32) -> Result<Vector3<f64>, String>;
}
pub trait AtomArrayRef {
    fn get_atom_by_id(&self, atom_id: u32) -> Result<&Atom, String>;
    fn number_of_atoms(&self) -> usize;
    fn get_vector_ab(&self, a_id: u32, b_id: u32) -> Result<Vector3<f64>, String>;
}

impl AtomArray for Vec<Atom> {
    /// Get atom reference by id. The atom id is one-based index
    fn get_atom_by_id(&self, atom_id: u32) -> Result<&Atom, String> {
        self.get(atom_id as usize - 1)
            .ok_or("Invalid atom index".to_string())
    }

    /// Get atom mutable reference by id. The atom id is one-based index
    fn get_mut_atom_by_id(&mut self, atom_id: u32) -> Result<&mut Atom, String> {
        self.get_mut(atom_id as usize - 1)
            .ok_or("Invalid atom index".to_string())
    }

    fn append_atom(&mut self, new_atom: Atom) {
        self.push(new_atom)
    }

    fn number_of_atoms(&self) -> usize {
        self.len()
    }

    fn get_vector_ab(&self, a_id: u32, b_id: u32) -> Result<Vector3<f64>, String> {
        let atom_a: &Atom = self.get_atom_by_id(a_id)?;
        let atom_b: &Atom = self.get_atom_by_id(b_id)?;
        let atom_a_xyz = atom_a.xyz();
        let atom_b_xyz = atom_b.xyz();
        Ok(atom_b_xyz - atom_a_xyz)
    }
}

impl AtomArray for &mut Vec<Atom> {
    fn get_atom_by_id(&self, atom_id: u32) -> Result<&Atom, String> {
        self.get(atom_id as usize - 1)
            .ok_or("Invalid atom index".to_string())
    }

    fn get_mut_atom_by_id(&mut self, atom_id: u32) -> Result<&mut Atom, String> {
        self.get_mut(atom_id as usize - 1)
            .ok_or("Invalid atom index".to_string())
    }

    fn append_atom(&mut self, new_atom: Atom) {
        self.push(new_atom)
    }

    fn number_of_atoms(&self) -> usize {
        self.len()
    }

    fn get_vector_ab(&self, a_id: u32, b_id: u32) -> Result<Vector3<f64>, String> {
        let atom_a: &Atom = self.get_atom_by_id(a_id)?;
        let atom_b: &Atom = self.get_atom_by_id(b_id)?;
        let atom_a_xyz = atom_a.xyz();
        let atom_b_xyz = atom_b.xyz();
        Ok(atom_b_xyz - atom_a_xyz)
    }
}

impl AtomArrayRef for &[Atom] {
    fn get_atom_by_id(&self, atom_id: u32) -> Result<&Atom, String> {
        self.get(atom_id as usize - 1)
            .ok_or("Invalid atom index".to_string())
    }

    fn number_of_atoms(&self) -> usize {
        self.len()
    }

    fn get_vector_ab(&self, a_id: u32, b_id: u32) -> Result<Vector3<f64>, String> {
        let atom_a: &Atom = self.get_atom_by_id(a_id)?;
        let atom_b: &Atom = self.get_atom_by_id(b_id)?;
        let atom_a_xyz = atom_a.xyz();
        let atom_b_xyz = atom_b.xyz();
        Ok(atom_b_xyz - atom_a_xyz)
    }
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

    pub fn set_atom_id(&mut self, atom_id: u32) {
        self.atom_id = atom_id;
    }
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

impl Export for Vec<Atom> {
    fn format_output(&self) -> String {
        let atom_strings: Vec<String> = self.iter().map(|atom| atom.format_output()).collect();
        atom_strings.concat()
    }
}

impl Transformation for Vec<Atom> {
    fn rotate(&mut self, rotate_quatd: na::UnitQuaternion<f64>) {
        self.iter_mut()
            .for_each(|atom: &mut Atom| atom.set_xyz(rotate_quatd.transform_point(atom.xyz())));
    }
    fn translate(&mut self, translate_matrix: na::Translation<f64, 3>) {
        self.iter_mut()
            .for_each(|atom: &mut Atom| atom.set_xyz(translate_matrix.transform_point(atom.xyz())));
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
