use castep_cell_io::{cell_document::BSKpointPathSpacing, InvLengthUnit};
use std::{
    fs::{self, create_dir, File},
    io::{self, BufRead, BufReader},
    path::{Path, PathBuf},
};

use castep_cell_io::{
    cell_document::{
        CellEntries, ExtEFieldBlock, ExtPressureBlock, FixAllCell, FixCom, IonicConstraintsBlock,
        KpointMPSpacing, KpointQuality, KpointSettings, NCKpointSettings, SpeciesLCAOStatesBlock,
        SpeciesMassBlock, SpeciesPotBlock,
    },
    CastepParams, CastepTask, CellDocument,
};
use castep_periodic_table::{
    data::ELEMENT_TABLE,
    element::{ElementFamily, LookupElement},
};
use chemrust_misctools::{write_server_script, ServerScriptType};
use crystal_cif_io::to_cif_document;

use crate::edit::Model;

#[derive(Debug)]
pub struct ExportPackage {
    filename_stem: String,
    full_cell_doc: CellDocument,
    param: CastepParams,
    potentials_loc: PathBuf,
}

#[derive(Debug)]
pub struct SeedFilePackage<P: AsRef<Path>> {
    geom_pack: ExportPackage,
    bs_pack: ExportPackage,
    parent_dir: P,
    output_dir: P,
}

impl<P: AsRef<Path>> SeedFilePackage<P> {
    pub fn new(
        geom_pack: ExportPackage,
        bs_pack: ExportPackage,
        parent_dir: P,
        output_dir: P,
    ) -> Self {
        Self {
            geom_pack,
            bs_pack,
            parent_dir,
            output_dir,
        }
    }

    fn directory_check(&self) -> Result<(), io::Error> {
        if !self.parent_dir.as_ref().exists() {
            create_dir(&self.parent_dir)?;
        }
        if !self.output_dir.as_ref().exists() {
            create_dir(&self.output_dir)
        } else {
            Ok(())
        }
    }

    pub fn write_to_dir(&self) -> Result<(), io::Error> {
        self.directory_check()?;
        self.geom_pack.write_to_path(&self.output_dir)?;
        self.geom_pack.write_cif_to_path(&self.output_dir)?;
        self.geom_pack.write_script(&self.output_dir)?;
        self.bs_pack.write_to_path(&self.output_dir)
    }

    pub fn copy_potentials(&self) -> Result<(), io::Error> {
        self.geom_pack
            .full_cell_doc
            .get_potential_paths(self.geom_pack.potentials_loc())
            .iter()
            .try_for_each(|pot_src| {
                let pot_dest = self
                    .output_dir
                    .as_ref()
                    .join(pot_src.file_name().expect("No filename in path"));
                if pot_dest.exists() {
                    Ok(())
                } else {
                    fs::copy(pot_src, pot_dest)?;
                    Ok(())
                }
            })
    }
}

impl ExportPackage {
    pub fn new(
        model_stem: String,
        full_cell_doc: CellDocument,
        param: CastepParams,
        potentials_loc: PathBuf,
    ) -> Self {
        Self {
            filename_stem: model_stem,
            full_cell_doc,
            param,
            potentials_loc,
        }
    }

    pub fn write_cif_to_path<P: AsRef<Path>>(&self, dest_dir: P) -> Result<(), io::Error> {
        let cif_document = to_cif_document(&self.full_cell_doc, &self.filename_stem);
        let cif_path = Path::new(dest_dir.as_ref())
            .join(&self.filename_stem)
            .with_extension("cif");
        fs::write(cif_path, format!("{}", cif_document))
    }

    pub fn write_to_path<P: AsRef<Path>>(&self, dest_dir: P) -> Result<(), io::Error> {
        let cell_path = Path::new(dest_dir.as_ref())
            .join(&self.filename_stem)
            .with_extension("cell");
        let param_path = Path::new(dest_dir.as_ref())
            .join(&self.filename_stem)
            .with_extension("param");
        fs::write(cell_path, format!("{}", self.full_cell_doc))?;
        fs::write(param_path, format!("{}", self.param))
    }

    pub fn write_script<P: AsRef<Path>>(&self, dest_dir: P) -> Result<(), io::Error> {
        let cell_path = Path::new(dest_dir.as_ref())
            .join(&self.filename_stem)
            .with_extension("cell");
        write_server_script(&cell_path, 8, ServerScriptType::Pbs)
    }

    pub fn potentials_loc(&self) -> &PathBuf {
        &self.potentials_loc
    }
}
pub struct SeedfileGenerator {
    task: CastepTask,
    cell_doc: CellDocument,
    use_edft: Option<bool>,
    kpoint_quality: Option<KpointQuality>,
}

impl SeedfileGenerator {
    pub fn use_edft(&mut self, use_edft: bool) {
        self.use_edft = Some(use_edft);
    }

    pub fn set_kpoint_quality(&mut self, kpoint_quality: KpointQuality) {
        self.kpoint_quality = Some(kpoint_quality);
    }

    pub fn new(task: CastepTask, cell_doc: CellDocument) -> Self {
        let use_edft = cell_doc.get_elements().iter().any(|elm| {
            matches!(elm.family(), ElementFamily::RareEarthLa)
                || matches!(elm.family(), ElementFamily::RareEarthAc)
        });
        Self {
            task,
            cell_doc,
            use_edft: Some(use_edft),
            kpoint_quality: None,
        }
    }

    fn get_total_spin(&self) -> u32 {
        self.cell_doc.total_spin()
    }

    fn get_cutoff_energy<P: AsRef<Path>>(&self, potentials_loc: P) -> Result<f64, io::Error> {
        Ok(self
            .cell_doc
            .get_elements()
            .iter()
            .map(|&elm| -> Result<f64, io::Error> {
                let potential_file = ELEMENT_TABLE.get_by_symbol(elm).potential();
                let potential_path = Path::new(potentials_loc.as_ref()).join(potential_file);
                let file = File::open(potential_path)?;
                let reader = BufReader::new(file);
                let fine_energy = reader
                    .lines()
                    .find(|line| line.as_ref().unwrap().contains("FINE"))
                    .map(|line| {
                        let num_str = line.as_ref().unwrap().split_whitespace().next().unwrap();
                        num_str.parse::<u32>().expect("Can't parse into `u32`")
                    })
                    .expect("Failed to parse fine energy from pseudopotential file.");
                let round_bigger_tenth = |num: u32| -> f64 {
                    match num % 10 {
                        0 => num as f64,
                        _ => ((num / 10 + 1) * 10) as f64,
                    }
                };
                let ultra_fine_energy = round_bigger_tenth((fine_energy as f64 * 1.1) as u32);
                Ok(ultra_fine_energy)
            })
            .filter_map(|res| res.ok())
            .reduce(|prev, next| if next > prev { next } else { prev })
            .expect("Error in comparing the largest cutoff energy"))
    }

    fn geom_opt_cell(&self) -> CellDocument {
        let elements = self.cell_doc.get_elements();
        let entries = vec![
            CellEntries::KpointSettings(KpointSettings::MPSpacing(KpointMPSpacing::default())),
            CellEntries::FixAllCell(FixAllCell::new(true)),
            CellEntries::FixCom(FixCom::new(false)),
            CellEntries::IonicConstraints(IonicConstraintsBlock::default()),
            CellEntries::ExtEfield(ExtEFieldBlock::default()),
            CellEntries::ExtPressure(ExtPressureBlock::default()),
            CellEntries::SpeciesMass(SpeciesMassBlock::from_elements(&elements)),
            CellEntries::SpeciesPot(SpeciesPotBlock::from_elements(&elements)),
            CellEntries::SpeciesLCAOStates(SpeciesLCAOStatesBlock::from_elememts(&elements)),
        ];
        let mut geom_cell = self.cell_doc.clone();
        geom_cell.set_entries(Some(entries));
        geom_cell
    }

    fn bs_cell(&self) -> CellDocument {
        let mut bs_cell = self.cell_doc.clone();
        let elements = self.cell_doc.get_elements();
        let entries = vec![
            CellEntries::KpointSettings(KpointSettings::MPSpacing(KpointMPSpacing::default())),
            CellEntries::NCKpointSettings(NCKpointSettings::PathSpacing(BSKpointPathSpacing::new(
                InvLengthUnit::Ang,
                0.07,
            ))),
            CellEntries::FixAllCell(FixAllCell::new(true)),
            CellEntries::FixCom(FixCom::new(false)),
            CellEntries::IonicConstraints(IonicConstraintsBlock::default()),
            CellEntries::ExtEfield(ExtEFieldBlock::default()),
            CellEntries::ExtPressure(ExtPressureBlock::default()),
            CellEntries::SpeciesMass(SpeciesMassBlock::from_elements(&elements)),
            CellEntries::SpeciesPot(SpeciesPotBlock::from_elements(&elements)),
            CellEntries::SpeciesLCAOStates(SpeciesLCAOStatesBlock::from_elememts(&elements)),
        ];
        bs_cell.set_entries(Some(entries));
        bs_cell
    }

    pub fn generate_cell_file(&self) -> CellDocument {
        match self.task {
            CastepTask::BandStructure => self.bs_cell(),
            CastepTask::GeometryOptimization => self.geom_opt_cell(),
        }
    }

    fn geom_opt_param<P: AsRef<Path>>(&self, potentials_loc: P) -> CastepParams {
        CastepParams::geom_opt(
            self.get_cutoff_energy(potentials_loc)
                .expect("Failed to get cutoff energy"),
            self.get_total_spin(),
            self.use_edft.unwrap_or(false),
        )
    }

    fn bs_param<P: AsRef<Path>>(&self, potentials_loc: P) -> CastepParams {
        CastepParams::band_structure(
            self.get_cutoff_energy(potentials_loc)
                .expect("Failed to get cutoff energy"),
            self.get_total_spin(),
            self.use_edft.unwrap_or(false),
        )
    }

    pub fn generate_castep_param<P: AsRef<Path>>(&self, potentials_loc: P) -> CastepParams {
        match self.task {
            CastepTask::BandStructure => self.bs_param(potentials_loc),
            CastepTask::GeometryOptimization => self.geom_opt_param(potentials_loc),
        }
    }
    // TODO: kptaux
}

impl Model {
    pub fn export_files<P: AsRef<Path>>(
        &self,
        potentials_loc: P,
        task: CastepTask,
    ) -> ExportPackage {
        let seed_gen = SeedfileGenerator::new(task, self.basic_cell().clone());
        let cell = seed_gen.generate_cell_file();
        let param = seed_gen.generate_castep_param(potentials_loc.as_ref());
        let output_name = match task {
            CastepTask::BandStructure => format!("{}_DOS", self.output_name()),
            CastepTask::GeometryOptimization => self.output_name(),
        };
        ExportPackage::new(output_name, cell, param, potentials_loc.as_ref().into())
    }
}
