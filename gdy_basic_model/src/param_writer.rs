pub mod param_writer {
    use std::fs::{self, create_dir_all};
    use std::path::{Path, PathBuf};
    use std::{collections::HashMap, fs::read_to_string};

    use crate::external_info::element_table::Element;
    use crate::{Atom, Cell};
    use regex::Regex;

    fn get_final_cutoff_energy(cell: &Cell, element_infotab: &HashMap<String, Element>) -> f64 {
        let mut energy: f64 = 0.0;
        let element_lists = cell.lattice.get_element_list();
        let fine_cutoff_energy_regex =
            Regex::new(r"([0-9]+) FINE").expect("Error in compiling regex pattern");
        element_lists.iter().for_each(|elm| {
            let potential_file = &element_infotab.get(elm).unwrap().pot;
            let potential_file_contents =
                read_to_string(format!("./resources/Potentials/{potential_file}"))
                    .expect("Errors in opening potential file");
            let fine_cutoff_energy: u32 = fine_cutoff_energy_regex
                .captures(&potential_file_contents)
                .expect("Error in capturing fine cutoff energy")
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
        let export_destination = |cell: &Cell| -> PathBuf {
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
                "./GDY_TAC_models/{}/{}",
                family,
                main_metal_element.element_name()
            );
            create_dir_all(&dir_path).unwrap_or_else(|why| {
                println!("! {:?}", why.kind());
            });
            Path::new(&dir_path).to_path_buf()
        };
        let param_write_dest = export_destination(cell);
        let geom_param_filename = format!("{}.param", cell.lattice.molecule.mol_name);
        let geom_param_path = param_write_dest.join(&geom_param_filename);
        fs::write(geom_param_path, geom_param_content).expect(&format!(
            "Unable to write geom param for {}",
            cell.lattice.molecule.mol_name
        ));
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
        let dos_param_filename = format!("{}_DOS.param", cell.lattice.molecule.mol_name);
        let dos_param_path = param_write_dest.join(&dos_param_filename);
        fs::write(dos_param_path, dos_param_content).expect(&format!(
            "Unable to write dos param for {}",
            cell.lattice.molecule.mol_name
        ));
    }
}
