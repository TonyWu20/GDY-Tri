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
        let filename = "GDY_tri.msi";
        let base_lat: Lattice = parser::msi_parser::parse_lattice(filename);
        println!("{}", base_lat.export_msi());
    }
}
