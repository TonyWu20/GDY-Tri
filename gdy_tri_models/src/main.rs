use std::{fs::create_dir, io, path::Path};

use castep_cell_io::CastepTask;
use edit::generate_models;
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;

use crate::export::SeedFilePackage;

mod edit;
mod export;
mod template;

fn execution() -> Result<(), io::Error> {
    let potentials_loc = Path::new(env!("CARGO_MANIFEST_DIR")).join("Potentials");
    let dest_dir = Path::new("GDY_TAC_Models");
    if !dest_dir.exists() {
        create_dir(dest_dir)?;
    }
    generate_models()
        .par_iter()
        .progress()
        .map(|model| {
            let model_family_dir = dest_dir.join(model.dir_name());
            let model_dir = model_family_dir.join(model.output_name());
            SeedFilePackage::new(
                model.export_files(&potentials_loc, CastepTask::GeometryOptimization),
                model.export_files(&potentials_loc, CastepTask::BandStructure),
                model_family_dir,
                model_dir,
            )
        })
        .try_for_each(|seed| seed.write_to_dir())
}

fn main() -> Result<(), io::Error> {
    execution()
}
