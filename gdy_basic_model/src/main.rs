extern crate nalgebra as na;
use std::error::Error;

use gdy_model::{
    editor,
    external_info::project::{load_project_info, ProjectInfo},
    param_writer,
};
fn main() -> Result<(), Box<dyn Error>> {
    let project_info = load_project_info("./resources/project.yaml")?;
    task_gen_all(&project_info)?;
    Ok(())
}

fn task_gen_all(project_info: &ProjectInfo) -> Result<(), Box<dyn Error>> {
    println!("Generate all base models");
    editor::msi_editor::generate_all_base_models(project_info.base_model_loc())?;
    println!("Generate all seed files for base models");
    param_writer::param_writer::generate_all_seed_files(project_info.export_loc())?;
    param_writer::param_writer::to_xsd_scripts(project_info.export_loc());
    Ok(())
}
