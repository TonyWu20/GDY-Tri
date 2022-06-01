extern crate gdy_model;
extern crate nalgebra as na;
mod editor;
mod parser;
mod test;
fn main() {
    task_gen_all();
}

fn task_gen_all() {
    println!("Generate all base models");
    editor::msi_editor::generate_all_base_models("GDY_tri.msi");
}
