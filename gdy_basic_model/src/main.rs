extern crate nalgebra as na;
use std::error::Error;

use gdy_model::{editor, param_writer};
fn main() -> Result<(), Box<dyn Error>> {
    task_gen_all()?;
    Ok(())
}

fn task_gen_all() -> Result<(), Box<dyn Error>> {
    println!("Generate all base models");
    editor::msi_editor::generate_all_base_models("./resources/GDY_tri.msi")?;
    println!("Generate all seed files for base models");
    param_writer::param_writer::generate_all_seed_files("GDY_TAC_models");
    param_writer::param_writer::to_xsd_scripts("GDY_TAC_models");
    Ok(())
}
