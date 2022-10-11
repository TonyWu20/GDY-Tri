#[cfg(test)]
mod test {
    use std::error::Error;

    use crate::atom::AtomArray;
    use crate::cell::CellOutput;
    use crate::lattice::Lattice;
    use crate::param_writer::param_writer::{export_destination, write_seed_files_for_cell};
    use crate::parser::msi_parser::parse_lattice;
    use crate::{external_info::element_table, *};
    use glob::glob;
    use rayon::iter::{IntoParallelRefIterator, ParallelBridge, ParallelIterator};

    use crate::{editor::msi_editor::change_atom_element, parser};

    // use crate::{editor, parser};

    // #[test]
    // fn iterate_elements() {
    //     use gdy_model::*;
    //     let filename = "GDY_tri.msi";
    //     let mut base_lat: Lattice = parser::msi_parser::parse_lattice(filename);
    //     editor::msi_editor::iterate_over_elements(&mut base_lat);
    // }

    #[test]
    #[ignore]
    fn read_and_write() -> Result<(), Box<dyn Error>> {
        let filename = "./resources/GDY_tri.msi";
        let base_lat: Lattice = parser::msi_parser::parse_lattice(filename)?;
        Ok(println!("{}", base_lat.format_output()))
    }
    #[test]
    #[ignore]
    fn sort_atoms_in_cell() -> Result<(), Box<dyn Error>> {
        let filename = "./resources/GDY_tri.msi";
        let mut base_lat: Lattice = parser::msi_parser::parse_lattice(filename)?;
        println!("{:#?}", &base_lat.get_element_list());
        base_lat.sort_atoms_by_elements();
        Ok(())
    }
    #[test]
    #[ignore]
    fn rotate_lattice() -> Result<(), Box<dyn Error>> {
        let filename = "./resources/GDY_tri.msi";
        let mut base_lat: Lattice = parser::msi_parser::parse_lattice(filename)?;
        base_lat.rotate_to_standard_orientation();
        println!("{:#.12?}", base_lat.get_lattice_vectors());
        println!("{}", base_lat.format_output());
        Ok(())
    }
    #[test]
    fn cell_test() -> Result<(), Box<dyn Error>> {
        let filename = "./resources/GDY_tri.msi";
        let mut base_lat: Lattice = parser::msi_parser::parse_lattice(filename)?;
        change_atom_element(
            base_lat.atoms_vec_mut().get_mut_atom_by_id(73).unwrap(),
            "Mn",
            25,
        );
        change_atom_element(
            base_lat.atoms_vec_mut().get_mut_atom_by_id(74).unwrap(),
            "Mn",
            25,
        );
        change_atom_element(
            base_lat.atoms_vec_mut().get_mut_atom_by_id(75).unwrap(),
            "Ni",
            28,
        );
        // let frac_mat = fractional_coord_matrix(&base_lat);
        base_lat.sort_atoms_by_elements();
        base_lat.rotate_to_standard_orientation();
        let element_info = element_table::hash_table();
        let cell_output = base_lat.cell_output(&element_info);
        let spin_total = base_lat
            .atoms_vec()
            .iter()
            .map(|atom| -> u8 { element_info.get(atom.element_name()).unwrap().spin })
            .reduce(|total, x| total + x)
            .unwrap();
        println!("{}", cell_output);
        println!("{spin_total}");
        Ok(())
    }
    #[test]
    fn test_write_param() -> Result<(), Box<dyn Error>> {
        let filename = "./resources/GDY_tri.msi";
        let mut base_lat: Lattice = parser::msi_parser::parse_lattice(filename)?;
        change_atom_element(
            base_lat.atoms_vec_mut().get_mut_atom_by_id(73).unwrap(),
            "Mn",
            25,
        );
        change_atom_element(
            base_lat.atoms_vec_mut().get_mut_atom_by_id(74).unwrap(),
            "Mn",
            25,
        );
        change_atom_element(
            base_lat.atoms_vec_mut().get_mut_atom_by_id(75).unwrap(),
            "Ni",
            28,
        );
        // let frac_mat = fractional_coord_matrix(&base_lat);
        base_lat.update_base_name();
        let element_info = element_table::hash_table();
        write_seed_files_for_cell(&mut base_lat, &element_info)?;
        Ok(())
    }
    #[test]
    fn test_glob() {
        let root_dir = "GDY_TAC_models";
        let msi_pattern = format!("{root_dir}/**/*.msi");
        glob(&msi_pattern)
            .expect("Failed to read glob pattern")
            .into_iter()
            .par_bridge()
            .collect::<Vec<_>>()
            .par_iter()
            .for_each(|entry| match entry {
                Ok(path) => {
                    println!("{}", path.to_str().unwrap());
                    let lattice = parse_lattice(path.to_str().unwrap()).unwrap();
                    println!(
                        "{}",
                        export_destination(&lattice)
                            .expect("Export destination creation failed")
                            .to_str()
                            .unwrap()
                    );
                }
                Err(e) => println!("{:?}", e),
            });
    }
}
