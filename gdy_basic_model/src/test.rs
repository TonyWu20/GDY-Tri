#[cfg(test)]
mod test {
    use crate::{external_info::element_table, param_writer::param_writer::write_param, *};

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
    fn read_and_write() {
        let filename = "./resources/GDY_tri.msi";
        let base_lat: Lattice = parser::msi_parser::parse_lattice(filename);
        println!("{}", base_lat.format_output());
    }
    #[test]
    #[ignore]
    fn sort_atoms_in_cell() {
        let filename = "./resources/GDY_tri.msi";
        let mut base_lat: Lattice = parser::msi_parser::parse_lattice(filename);
        println!("{:#?}", &base_lat.get_element_list());
        let mut cell: Cell = Cell::new(&mut base_lat, false);
        cell.sort_atoms_by_elements();
    }
    #[test]
    #[ignore]
    fn rotate_lattice() {
        let filename = "./resources/GDY_tri.msi";
        let mut base_lat: Lattice = parser::msi_parser::parse_lattice(filename);
        base_lat.rotate_to_standard_orientation();
        println!("{:#.12?}", base_lat.get_lattice_vectors());
        println!("{}", base_lat.format_output());
    }
    #[test]
    fn cell_test() {
        let filename = "./resources/GDY_tri.msi";
        let mut base_lat: Lattice = parser::msi_parser::parse_lattice(filename);
        change_atom_element(base_lat.molecule.get_mut_atom_by_id(73).unwrap(), "Mn", 25);
        change_atom_element(base_lat.molecule.get_mut_atom_by_id(74).unwrap(), "Mn", 25);
        change_atom_element(base_lat.molecule.get_mut_atom_by_id(75).unwrap(), "Ni", 28);
        // let frac_mat = fractional_coord_matrix(&base_lat);
        let mut cell = Cell::new(&mut base_lat, false);
        cell.sort_atoms_by_elements();
        cell.lattice.rotate_to_standard_orientation();
        let element_info = element_table::hash_table();
        let cell_output = cell.format_output(&element_info);
        let spin_total = cell
            .lattice
            .molecule
            .atoms_iterator()
            .map(|atom| -> u8 { element_info.get(atom.element_name()).unwrap().spin })
            .reduce(|total, x| total + x)
            .unwrap();
        println!("{}", cell_output);
        println!("{spin_total}");
    }
    #[test]
    fn test_write_param() {
        let filename = "./resources/GDY_tri.msi";
        let mut base_lat: Lattice = parser::msi_parser::parse_lattice(filename);
        change_atom_element(base_lat.molecule.get_mut_atom_by_id(73).unwrap(), "Mn", 25);
        change_atom_element(base_lat.molecule.get_mut_atom_by_id(74).unwrap(), "Mn", 25);
        change_atom_element(base_lat.molecule.get_mut_atom_by_id(75).unwrap(), "Ni", 28);
        // let frac_mat = fractional_coord_matrix(&base_lat);
        base_lat.update_base_name();
        let cell = Cell::new(&mut base_lat, false);
        let element_info = element_table::hash_table();
        write_param(&cell, &element_info);
    }
}
