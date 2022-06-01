pub mod msi_editor {
    use std::{
        fs::{self, create_dir_all},
        path::{Path, PathBuf},
    };

    use gdy_model::{Atom, Lattice};
    use periodic_table as pt;
    use pt::Element;

    use crate::parser::msi_parser::parse_lattice;
    pub fn change_atom_element(
        target_atom: &mut Atom,
        new_element_name: &str,
        new_element_id: u32,
    ) {
        target_atom.set_element_name(new_element_name);
        target_atom.set_element_id(new_element_id);
    }
    pub fn iterate_over_elements(
        target_lattice: &mut Lattice,
        family: &str,
        to_use_metals: &Vec<&Element>,
    ) {
        let mut export_dirs: Vec<PathBuf> = vec![];
        for metal in to_use_metals.iter() {
            export_dirs.push(export_destination(family, metal.symbol));
        }
        let metals_dirs = std::iter::zip(to_use_metals, export_dirs);
        for (item, dir) in metals_dirs {
            for item_b in to_use_metals.iter() {
                let atom_1 = target_lattice.molecule.get_mut_atom_by_id(73).unwrap();
                change_atom_element(atom_1, item.symbol, item.atomic_number);
                let atom_2 = target_lattice.molecule.get_mut_atom_by_id(74).unwrap();
                change_atom_element(atom_2, item.symbol, item.atomic_number);
                let atom_3 = target_lattice.molecule.get_mut_atom_by_id(75).unwrap();
                change_atom_element(atom_3, item_b.symbol, item_b.atomic_number);
                let text = target_lattice.export_msi();
                let filepath = dir.join(format!(
                    "GDY_{}_{}_{}.msi",
                    item.symbol, item.symbol, item_b.symbol
                ));
                fs::write(filepath, text).expect("unable to write file");
            }
        }
    }
    fn export_destination(family: &str, main_element: &str) -> PathBuf {
        let dir_path = format!("./GDY_TAC_models/{}/{}", family, main_element);
        create_dir_all(&dir_path).unwrap_or_else(|why| {
            println!("! {:?}", why.kind());
        });
        Path::new(&dir_path).to_path_buf()
    }
    pub fn generate_all_base_models(src_filename: &str) {
        let mut src_lattice = parse_lattice(src_filename);
        let elements = pt::periodic_table();
        let metals_3d = elements[20..30].to_vec();
        let metals_4d = elements[38..48].to_vec();
        let metals_5d = elements[71..80].to_vec();
        let metals_rare_earth = elements[56..71].to_vec();
        iterate_over_elements(&mut src_lattice, "3d", &metals_3d);
        iterate_over_elements(&mut src_lattice, "4d", &metals_4d);
        iterate_over_elements(&mut src_lattice, "5d", &metals_5d);
        iterate_over_elements(&mut src_lattice, "rare_earth", &metals_rare_earth);
    }
}
