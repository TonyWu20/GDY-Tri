extern crate gdy_model;
extern crate nalgebra as na;
mod parser;
fn main() {
    println!("Hello World!");
}

#[test]
fn test_reg() {
    use gdy_model::*;
    use na::Matrix3;
    use std::fs;
    let filename = "SAC_GDY_V.msi";
    let contents = fs::read_to_string(filename).expect("Something went wrong reading the file");
    let atoms: Vec<Atom> = parser::msi_parser::parse_atom(&contents);
    let lattice_vectors: Matrix3<f64> = parser::msi_parser::parse_lattice_vectors(&contents);
    let filename_stem: &str = filename.split(".").collect::<Vec<&str>>()[0];
    let lattice_molecule: Molecule = Molecule::new(filename_stem.to_string(), atoms);
    let mut lat: Lattice = Lattice::new(lattice_molecule, lattice_vectors, vec![73], None);
    lat.rotate_to_standard_orientation();
    lat.set_adsorbate_name("C2H4".to_string());
    match lat.get_adsorbate_name() {
        Some(x) => println!("Adsorbate name: {}", x),
        None => println!("There is no adsorbate in the lattice now."),
    }
    println!(
        "Name: {}, number of atoms: {}",
        lat.molecule.mol_name,
        lat.molecule.number_of_atoms()
    );
}
