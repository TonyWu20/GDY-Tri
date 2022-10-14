/**
Module for loading `project.yaml`, which has necessary information to control
the program.
*/
extern crate serde;
use std::{collections::HashMap, error::Error, fs, path::Path};

use serde::{Deserialize, Serialize};

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct ProjectInfo {
    base_model_loc: String,
    element_table_loc: String,
    adsorbate_table_loc: String,
    potentials_loc: String,
    export_loc: String,
    coord_sites: Vec<CoordSite>,
    coord_cases: Vec<CoordCase>,
}

impl ProjectInfo {
    pub fn element_table_loc(&self) -> &str {
        self.element_table_loc.as_ref()
    }

    pub fn adsorbate_table_loc(&self) -> &str {
        self.adsorbate_table_loc.as_ref()
    }

    pub fn potentials_loc(&self) -> &str {
        self.potentials_loc.as_ref()
    }

    pub fn coord_sites(&self) -> &[CoordSite] {
        self.coord_sites.as_ref()
    }
    pub fn hash_coord_site(&self) -> HashMap<u32, String> {
        let mut hash_tab: HashMap<u32, String> = HashMap::new();
        self.coord_sites().iter().for_each(|coord_site| {
            hash_tab.insert(coord_site.atom_id, coord_site.name.to_string());
        });
        hash_tab
    }

    pub fn export_loc(&self) -> &str {
        self.export_loc.as_ref()
    }

    pub fn base_model_loc(&self) -> &str {
        self.base_model_loc.as_ref()
    }

    pub fn coord_cases(&self) -> &[CoordCase] {
        self.coord_cases.as_ref()
    }
}
#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct CoordSite {
    name: String,
    atom_id: u32,
}

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct CoordCase {
    name: String,
    cases: Vec<(u32, Option<u32>)>,
}

impl CoordCase {
    pub fn get_cases(&self, reverse: bool) -> Vec<(u32, Option<u32>)> {
        if reverse == true {
            let reversed_cases: Vec<(u32, Option<u32>)> = self
                .cases
                .iter()
                .map(|(s1, s2)| (s2.unwrap(), Some(*s1)))
                .collect();
            reversed_cases
        } else {
            self.cases.clone()
        }
    }
}

/**
Load project `yaml` from given path
# Arguments
* `filepath` - Type that has trait `AsRef<Path>`
# Returns
* `Result<ProjectInfo, Box<dyn Error>`
*/
pub fn load_project_info<P: AsRef<Path>>(filepath: P) -> Result<ProjectInfo, Box<dyn Error>> {
    let project_yaml = fs::File::open(filepath)?;
    let project_table = serde_yaml::from_reader(project_yaml)?;
    Ok(project_table)
}
