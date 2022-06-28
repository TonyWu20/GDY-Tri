use na::Vector3;

use crate::{atom::Atom, Export, Transformation};
#[derive(Debug, Clone)]
pub struct Molecule {
    pub mol_name: String,
    pub vector_atoms: Vec<Atom>,
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
