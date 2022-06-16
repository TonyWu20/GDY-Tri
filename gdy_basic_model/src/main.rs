extern crate nalgebra as na;
use gdy_model::{editor, param_writer};
fn main() {
    task_gen_all();
}

fn task_gen_all() {
    println!("Generate all base models");
    editor::msi_editor::generate_all_base_models("./resources/GDY_tri.msi");
    println!("Generate all seed files for base models");
    param_writer::param_writer::generate_all_seed_files("GDY_TAC_models");
}
