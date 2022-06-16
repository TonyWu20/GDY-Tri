pub mod param_writer {
    use indicatif::ProgressBar;
    use std::fs::{self, create_dir_all};
    use std::path::{Path, PathBuf};
    use std::{collections::HashMap, fs::read_to_string};

    use crate::external_info::element_table::{self, Element};
    use crate::parser::msi_parser::parse_lattice;
    use crate::{Atom, Cell};
    use glob::glob;
    use rayon::prelude::*;
    use regex::Regex;

    pub fn generate_all_seed_files(root_dir: &str) {
        let element_infotab = element_table::hash_table();
        let msi_pattern = format!("{root_dir}/**/*.msi");
        let file_iter = glob(&msi_pattern)
            .expect("Failed to read glob pattern")
            .into_iter();
        let file_count = glob(&msi_pattern)
            .expect("Failed to read glob pattern")
            .into_iter()
            .count();
        let bar = ProgressBar::new(file_count as u64);
        file_iter
            .par_bridge()
            .into_par_iter()
            .for_each(|entry| match entry {
                Ok(path) => {
                    let mut lattice = parse_lattice(path.to_str().unwrap());
                    let mut cell = Cell::new(&mut lattice, false);
                    cell.sort_atoms_by_elements();
                    write_seed_files_for_cell(&cell, &element_infotab);
                    bar.inc(1);
                }
                Err(e) => println!("{:?}", e),
            });
        bar.finish();
    }

    pub fn write_seed_files_for_cell(cell: &Cell, element_infotab: &HashMap<String, Element>) {
        write_param(cell, element_infotab);
        write_kptaux(cell);
        write_trjaux(cell);
        copy_potentials(cell, element_infotab);
        copy_smcastep_extension(cell);
        let export_dir = export_destination(cell);
        let msi_path = export_dir
            .parent()
            .unwrap()
            .join(&format!("{}.msi", cell.lattice.molecule.mol_name));
        let moved_dest = export_dir.join(&msi_path.file_name().unwrap());
        if moved_dest.exists() == false {
            fs::rename(&msi_path, moved_dest).expect("Move msi file failed!");
        }
    }

    pub fn export_destination(cell: &Cell) -> PathBuf {
        let main_metal_element: &Atom = cell
            .lattice
            .molecule
            .get_atom_by_id(cell.lattice.get_metal_sites()[0 as usize] as u8)
            .unwrap();
        let family = match main_metal_element.element_id() {
            21..=30 => "3d",
            39..=48 => "4d",
            72..=80 => "5d",
            57..=71 => "rare_earth",
            _ => "else",
        };
        let dir_path = format!(
            "GDY_TAC_models/{}/{}/{}_opt",
            family,
            main_metal_element.element_name(),
            cell.lattice.molecule.mol_name
        );
        create_dir_all(&dir_path).unwrap_or_else(|why| {
            println!("! {:?}", why.kind());
        });
        Path::new(&dir_path).to_path_buf()
    }

    fn export_filepath(cell: &Cell, filename: &str) -> PathBuf {
        let export_dest = export_destination(cell);
        let export_filename = format!("{}{}", cell.lattice.molecule.mol_name, filename);
        export_dest.join(export_filename)
    }

    fn get_final_cutoff_energy(cell: &Cell, element_infotab: &HashMap<String, Element>) -> f64 {
        let mut energy: f64 = 0.0;
        let element_lists = cell.lattice.get_element_list();
        let fine_cutoff_energy_regex =
            Regex::new(r"([0-9]+).*FINE").expect("Error in compiling regex pattern");
        element_lists.iter().for_each(|elm| {
            let potential_file = &element_infotab.get(elm).unwrap().pot;
            let potential_file_contents =
                read_to_string(format!("./resources/Potentials/{potential_file}"))
                    .expect("Errors in opening potential file");
            let fine_cutoff_energy: u32 = fine_cutoff_energy_regex
                .captures(&potential_file_contents)
                .expect(&format!(
                    "Error in capturing fine cutoff energy for {}",
                    elm
                ))
                .get(1)
                .unwrap()
                .as_str()
                .parse::<u32>()
                .expect("Error in parsing fine cutoff energy as u32");
            let ultra_fine_energy = ((fine_cutoff_energy / 10 + 1) * 10) as f64;
            energy = if energy > ultra_fine_energy {
                energy
            } else {
                ultra_fine_energy
            };
        });
        energy
    }

    pub fn write_param(cell: &Cell, element_infotab: &HashMap<String, Element>) {
        let geom_param_path = export_filepath(cell, ".param");
        if !geom_param_path.exists() {
            let cutoff_energy = get_final_cutoff_energy(cell, element_infotab);
            let spin_total = cell
                .lattice
                .molecule
                .atoms_iterator()
                .map(|atom| -> u8 { element_infotab.get(atom.element_name()).unwrap().spin })
                .reduce(|total, i| total + i)
                .unwrap();
            let geom_param_content = format!(
                r#"task : BandStructure
continuation : default
comment : CASTEP calculation from Materials Studio
xc_functional : PBE
spin_polarized : true
spin :        {spin_total}
opt_strategy : Speed
page_wvfns :        0
cut_off_energy : {cutoff_energy:18.15}
grid_scale :        1.500000000000000
fine_grid_scale :        1.500000000000000
finite_basis_corr :        0
elec_energy_tol :   1.000000000000000e-005
max_scf_cycles :     6000
fix_occupancy : false
metals_method : dm
mixing_scheme : Pulay
mix_charge_amp :        0.500000000000000
mix_spin_amp :        2.000000000000000
mix_charge_gmax :        1.500000000000000
mix_spin_gmax :        1.500000000000000
mix_history_length :       20
perc_extra_bands :      72
smearing_width :        0.100000000000000
spin_fix :        6
num_dump_cycles : 0
bs_nextra_bands :       72
bs_xc_functional : PBE
bs_eigenvalue_tol :   1.000000000000000e-005
calculate_stress : false
calculate_ELF : false
popn_calculate : false
calculate_hirshfeld : false
calculate_densdiff : false
pdos_calculate_weights : true
bs_write_eigenvalues : true
"#
            );
            fs::write(&geom_param_path, geom_param_content).expect(&format!(
                "Unable to write geom param for {}",
                geom_param_path.to_str().unwrap()
            ));
        }
        let dos_param_path = export_filepath(cell, "_DOS.param");
        if !dos_param_path.exists() {
            let cutoff_energy = get_final_cutoff_energy(cell, element_infotab);
            let spin_total = cell
                .lattice
                .molecule
                .atoms_iterator()
                .map(|atom| -> u8 { element_infotab.get(atom.element_name()).unwrap().spin })
                .reduce(|total, i| total + i)
                .unwrap();
            let dos_param_content = format!(
                r#"task : BandStructure
continuation : default
comment : CASTEP calculation from Materials Studio
xc_functional : PBE
spin_polarized : true
spin :        {spin_total}
opt_strategy : Speed
page_wvfns :        0
cut_off_energy :      {cutoff_energy:.15}
grid_scale :        1.500000000000000
fine_grid_scale :        1.500000000000000
finite_basis_corr :        0
elec_energy_tol :   1.000000000000000e-005
max_scf_cycles :     6000
fix_occupancy : false
metals_method : dm
mixing_scheme : Pulay
mix_charge_amp :        0.500000000000000
mix_spin_amp :        2.000000000000000
mix_charge_gmax :        1.500000000000000
mix_spin_gmax :        1.500000000000000
mix_history_length :       20
perc_extra_bands :      72
smearing_width :        0.100000000000000
spin_fix :        6
num_dump_cycles : 0
bs_nextra_bands :       72
bs_xc_functional : PBE
bs_eigenvalue_tol :   1.000000000000000e-005
calculate_stress : false
calculate_ELF : false
popn_calculate : false
calculate_hirshfeld : false
calculate_densdiff : false
pdos_calculate_weights : true
bs_write_eigenvalues : true
"#
            );
            fs::write(dos_param_path, dos_param_content).expect(&format!(
                "Unable to write dos param for {}",
                cell.lattice.molecule.mol_name
            ));
        }
    }
    fn write_kptaux(cell: &Cell) {
        let kptaux_contents = r#"MP_GRID :        1       1       1
MP_OFFSET :   0.000000000000000e+000  
0.000000000000000e+000  0.000000000000000e+000
%BLOCK KPOINT_IMAGES
   1   1
%ENDBLOCK KPOINT_IMAGES"#
            .to_string();
        let kptaux_path = export_filepath(cell, ".kptaux");
        if !kptaux_path.exists() {
            fs::write(kptaux_path, &kptaux_contents).expect(&format!(
                "Unable to write kptaux for {}",
                cell.lattice.molecule.mol_name
            ));
        }
        let kptaux_dos_path = export_filepath(cell, "_DOS.kptaux");
        if !kptaux_dos_path.exists() {
            fs::write(kptaux_dos_path, &kptaux_contents).expect(&format!(
                "Unable to write dos_kptaux for {}",
                cell.lattice.molecule.mol_name
            ));
        }
    }
    fn write_trjaux(cell: &Cell) {
        let trjaux_path = export_filepath(cell, ".trjaux");
        if !trjaux_path.exists() {
            let mut trjaux_contents = String::new();
            let trjaux_header = r#"# Atom IDs to appear in any .trj file to be generated.
# Correspond to atom IDs which will be used in exported .msi file
# required for animation/analysis of trajectory within Cerius2.
"#;
            trjaux_contents.push_str(trjaux_header);
            cell.lattice.molecule.atoms_iterator().for_each(|atom| {
                trjaux_contents.push_str(&format!("{}\n", atom.atom_id()));
            });
            let trjaux_ending = r#"#Origin  0.000000000000000e+000  0.000000000000000e+000  0.000000000000000e+000"#;
            trjaux_contents.push_str(trjaux_ending);
            fs::write(trjaux_path, trjaux_contents).expect(&format!(
                "Unable to write trjaux for {}",
                cell.lattice.molecule.mol_name
            ));
        }
    }
    fn copy_potentials(cell: &Cell, element_infotab: &HashMap<String, Element>) {
        let target_dir = export_destination(cell);
        cell.lattice.get_element_list().iter().for_each(|elm| {
            let pot_file = &element_infotab.get(elm).unwrap().pot;
            let original_file = format!("./resources/Potentials/{}", pot_file);
            let original_path = Path::new(&original_file);
            let dest_path = target_dir.join(pot_file);
            if !dest_path.exists() {
                fs::copy(original_path, dest_path).expect("Error in copying potential file!");
            }
        });
    }
    fn copy_smcastep_extension(cell: &Cell) {
        let target_dir = export_destination(cell);
        let target_filename = format!("SMCastep_Extension_{}.xms", cell.lattice.molecule.mol_name);
        let target_path = target_dir.join(target_filename);
        if !target_path.exists() {
            fs::copy("./resources/SMCastep_Extension.xms", target_path)
                .expect("Error in copying SMCastep_Extension.xms!");
        }
    }
}
