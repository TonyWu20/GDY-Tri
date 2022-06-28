pub mod msi_editor {
    use std::{
        fs::{self, create_dir_all},
        path::{Path, PathBuf},
    };

    use crate::lattice::Lattice;
    use crate::{atom::Atom, Export};
    use indicatif::ProgressBar;
    use periodic_table as pt;
    use pt::Element;

    use crate::parser::msi_parser::parse_lattice;
    pub fn change_atom_element(target_atom: &mut Atom, new_element_name: &str, new_element_id: u8) {
        target_atom.set_element_name(new_element_name);
        target_atom.set_element_id(new_element_id);
    }
    pub fn iterate_over_elements(target_lattice: &mut Lattice, to_use_metals: &Vec<&Element>) {
        let mut export_dirs: Vec<PathBuf> = vec![];
        for metal in to_use_metals.iter() {
            export_dirs.push(export_destination(metal));
        }
        let metals_dirs = to_use_metals.iter().zip(export_dirs);
        let bar = ProgressBar::new((to_use_metals.len().pow(2)) as u64);
        for (item, dir) in metals_dirs {
            for item_b in to_use_metals.iter() {
                let atom_1 = target_lattice.molecule.get_mut_atom_by_id(73).unwrap();
                change_atom_element(atom_1, item.symbol, item.atomic_number as u8);
                let atom_2 = target_lattice.molecule.get_mut_atom_by_id(74).unwrap();
                change_atom_element(atom_2, item.symbol, item.atomic_number as u8);
                let atom_3 = target_lattice.molecule.get_mut_atom_by_id(75).unwrap();
                change_atom_element(atom_3, item_b.symbol, item_b.atomic_number as u8);
                target_lattice.update_base_name();
                let text = target_lattice.format_output();
                let lat_name = &target_lattice.molecule.mol_name;
                let filepath = dir.join(format!("{}_opt/{}.msi", &lat_name, &lat_name));
                if !filepath.exists() {
                    let parent = filepath.parent().unwrap();
                    if !parent.exists() {
                        create_dir_all(&parent).unwrap_or_else(|why| {
                            println!("! {:?}", why.kind());
                        })
                    }
                    fs::write(filepath, text).expect("unable to write file");
                }
                bar.inc(1)
            }
        }
        bar.finish();
    }
    pub fn export_destination(element: &Element) -> PathBuf {
        let family: &str = match element.atomic_number {
            21..=30 => "3d",
            39..=48 => "4d",
            72..=80 => "5d",
            57..=71 => "rare_earth",
            _ => "else",
        };
        let dir_path = format!("./GDY_TAC_models/{}/{}", family, element.symbol);
        create_dir_all(&dir_path).unwrap_or_else(|why| {
            println!("! {:?}", why.kind());
        });
        Path::new(&dir_path).to_path_buf()
    }
    pub fn generate_all_base_models(src_filename: &str) {
        let mut src_lattice = parse_lattice(src_filename);
        let elements: &[&Element] = pt::periodic_table();
        let metals_3d: &[&Element] = &elements[20..30];
        let metals_4d: &[&Element] = &elements[38..48];
        let metals_5d: &[&Element] = &elements[71..80];
        let metals_rare_earth: &[&Element] = &elements[56..71];
        let mut total_elements: Vec<&Element> = vec![];
        total_elements.extend_from_slice(metals_3d);
        total_elements.extend_from_slice(metals_4d);
        total_elements.extend_from_slice(metals_5d);
        total_elements.extend_from_slice(metals_rare_earth);
        iterate_over_elements(&mut src_lattice, &total_elements);
    }
}
