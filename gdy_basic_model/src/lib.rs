#![allow(dead_code)]
pub mod atom;
pub mod cell;
pub mod editor;
pub mod external_info;
pub mod lattice;
pub mod molecule;
pub mod param_writer;
pub mod parser;
pub mod test;

extern crate nalgebra as na;

// Shared behaviour
pub trait Export {
    fn format_output(&self) -> String;
}
pub trait Transformation {
    fn rotate(&mut self, rotate_quatd: na::UnitQuaternion<f64>);
    fn translate(&mut self, translate_matrix: na::Translation<f64, 3>);
}
// trait ends
