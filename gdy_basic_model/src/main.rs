extern crate nalgebra as na;
use gdy_model::editor;
fn main() {
    task_gen_all();
}

fn task_gen_all() {
    println!("Generate all base models");
    editor::msi_editor::generate_all_base_models("./resources/GDY_tri.msi");
}
