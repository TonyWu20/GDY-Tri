use std::error::Error;

use castep_model_generator_backend::external_info::project::{load_project_info, ProjectInfo};
fn main() -> Result<(), Box<dyn Error>> {
    let project_info = load_project_info("./resources/project.yaml")?;
    task_gen_all(&project_info)?;
    Ok(())
}

fn task_gen_all(project_info: &ProjectInfo) -> Result<(), Box<dyn Error>> {
    println!("Generate all base models");
    gdy_tri_basic_models::editor::gdy_tri_editor::generate_all_base_models(
        project_info.base_model_loc(),
        project_info.export_loc(),
    )?;
    Ok(())
}
