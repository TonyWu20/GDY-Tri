use std::{
    fs::{self, create_dir},
    io::{self},
    path::Path,
};

use castep_cell_io::{CastepParams, CastepTask, CellDocument, SeedfileGenerator};
use crystal_cif_io::{
    data_dict::{
        core_cif::{
            atom_site::LoopAtomSiteData,
            cell::CellDataSection,
            space_group::{
                CrystalSystem, ITNumber, SpaceGroupItem, SpaceGroupLoopData, SpaceGroupLoopItem,
                SpaceGroupSection,
            },
        },
        CifData, DataBlock, LoopDataEntry,
    },
    CifFile,
};

use crate::edit::Model;

#[derive(Debug)]
pub struct ExportPackage {
    filename_stem: String,
    full_cell_doc: CellDocument,
    param: CastepParams,
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
        self.geom_pack.write_cif_to_path(&self.output_dir)?;
        self.geom_pack.write_to_path(&self.output_dir)?;
        self.bs_pack.write_to_path(&self.output_dir)
    }
}

impl ExportPackage {
    pub fn new(model_stem: String, full_cell_doc: CellDocument, param: CastepParams) -> Self {
        Self {
            filename_stem: model_stem,
            full_cell_doc,
            param,
        }
    }

    pub fn write_cif_to_path<P: AsRef<Path>>(&self, dest_dir: P) -> Result<(), io::Error> {
        let space_group_section = SpaceGroupSection::init_builder()
            .add_entry(SpaceGroupItem::Crystal_system(CrystalSystem::Triclinic))
            .add_entry(SpaceGroupItem::IT_number(ITNumber::new(1)))
            .finish();
        let mut data_block = DataBlock::init(self.filename_stem.to_string());
        data_block
            .add_section(CifData::SpaceGroup(space_group_section))
            .add_section(CifData::SpaceGroupLoop(
                SpaceGroupLoopData::init_builder()
                    .add_entry(
                        LoopDataEntry::init_builder()
                            .add_entry(SpaceGroupLoopItem::Symop_operation_xyz("x,y,z".to_string()))
                            .finish(),
                    )
                    .finish(),
            ))
            .add_section(CifData::CellData(CellDataSection::from(
                self.full_cell_doc.model_description().lattice_block(),
            )))
            .add_section(CifData::AtomSiteLoop(LoopAtomSiteData::from(
                self.full_cell_doc.model_description().ionic_pos_block(),
            )));
        let cif_file = CifFile::new(vec![data_block]);
        let cif_path = Path::new(dest_dir.as_ref())
            .join(&self.filename_stem)
            .with_extension("cif");
        fs::write(cif_path, format!("{}", cif_file))
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
}

impl Model {
    pub fn export_files<P: AsRef<Path>>(
        &self,
        potentials_loc: P,
        task: CastepTask,
    ) -> ExportPackage {
        let seed_gen = SeedfileGenerator::new(task, self.basic_cell().clone());
        let cell = seed_gen.generate_cell_file();
        let param = seed_gen.generate_castep_param(potentials_loc);
        let output_name = match task {
            CastepTask::BandStructure => format!("{}_DOS", self.output_name()),
            CastepTask::GeometryOptimization => self.output_name(),
        };
        ExportPackage::new(output_name, cell, param)
    }
}
