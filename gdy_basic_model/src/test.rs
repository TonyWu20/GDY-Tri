#[cfg(test)]
mod test {
    use gdy_model::*;

    use crate::parser;

    // use crate::{editor, parser};

    // #[test]
    // fn iterate_elements() {
    //     use gdy_model::*;
    //     let filename = "GDY_tri.msi";
    //     let mut base_lat: Lattice = parser::msi_parser::parse_lattice(filename);
    //     editor::msi_editor::iterate_over_elements(&mut base_lat);
    // }

    #[test]
    fn read_and_write() {
        let filename = "./resources/GDY_tri.msi";
        let base_lat: Lattice = parser::msi_parser::parse_lattice(filename);
        println!("{}", base_lat.export_msi());
    }
    #[test]
    fn sort_atoms_in_cell() {
        let filename = "./resources/GDY_tri.msi";
        let base_lat: Lattice = parser::msi_parser::parse_lattice(filename);
        let mut cell: Cell = Cell::new(base_lat, false);
        cell.sort_atoms_by_elements();
        for atom in cell.lattice.molecule.atoms_iterator() {
            println!("{}", atom.element_name());
        }
    }
}
