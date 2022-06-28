use crate::Export;
use ::core::fmt;
use std::cmp::Ordering;

use na::Point3;
// Atom
#[derive(Debug, Clone)]
pub struct Atom {
    element_name: String,
    element_id: u8,
    xyz: Point3<f64>,
    atom_id: u8,
    // lcao: u8,
    // mass: f64,
    // pot: String,
    // spin: u8,
}

impl Atom {
    pub fn new(
        element_name: String,
        element_id: u8,
        xyz: Point3<f64>,
        atom_id: u8,
        // lcao: u8,
        // mass: f64,
        // pot: String,
        // spin: u8,
    ) -> Self {
        Self {
            element_name,
            element_id,
            xyz,
            atom_id,
            // lcao,
            // mass,
            // pot,
            // spin,
        }
    }

    pub fn element_name(&self) -> &str {
        &self.element_name
    }
    pub fn set_element_name(&mut self, new_name: &str) {
        self.element_name = new_name.to_string();
    }
    pub fn element_id(&self) -> u8 {
        self.element_id
    }
    pub fn set_element_id(&mut self, new_id: u8) {
        self.element_id = new_id;
    }
    pub fn xyz(&self) -> &Point3<f64> {
        &self.xyz
    }
    pub fn set_xyz(&mut self, new_xyz: Point3<f64>) {
        self.xyz = new_xyz;
    }
    pub fn atom_id(&self) -> u8 {
        self.atom_id
    }
    // pub fn lcao(&self) -> u8 {
    //     self.lcao
    // }
    // pub fn spin(&self) -> u8 {
    //     self.spin
    // }
    // pub fn potential_file(&self) -> String {
    //     self.pot.to_string()
    // }
    // pub fn mass(&self) -> f64 {
    //     self.mass
    // }
}

impl Export for Atom {
    fn format_output(&self) -> String {
        let msi_output: String = format!(
            r#"  ({item_id} Atom
    (A C ACL "{elm_id} {elm}")
    (A C Label "{elm}")
    (A D XYZ ({x:.12} {y:.12} {z:.12}))
    (A I Id {atom_id})
  )
"#,
            item_id = self.atom_id() + 1,
            elm_id = self.element_id(),
            elm = self.element_name(),
            x = self.xyz().x,
            y = self.xyz().y,
            z = self.xyz().z,
            atom_id = self.atom_id(),
        );
        msi_output
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Element: {}\nElement ID: {}\ncoord: {}\nAtom ID: {}",
            self.element_name, self.element_id, self.xyz, self.atom_id
        )
    }
}

impl Ord for Atom {
    fn cmp(&self, other: &Self) -> Ordering {
        self.atom_id.cmp(&other.atom_id)
    }
}

impl PartialOrd for Atom {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Atom {
    fn eq(&self, other: &Self) -> bool {
        self.atom_id == other.atom_id
    }
}

impl Eq for Atom {}
// End Atom
