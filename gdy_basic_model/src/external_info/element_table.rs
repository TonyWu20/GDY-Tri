extern crate serde;
use serde::{Deserialize, Serialize};
use std::{collections::HashMap, error::Error, fs, ops::Deref, path::Path};

use super::YamlTable;

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct Element {
    pub element: String,
    #[serde(rename = "atomic_num")]
    pub atomic_number: u8,
    #[serde(rename = "LCAO")]
    pub lcao: u8,
    pub mass: f64,
    pub pot: String,
    pub spin: u8,
}
#[derive(Debug, PartialEq, Serialize, Deserialize)]
pub struct ElmTab {
    #[serde(rename = "Element_info")]
    pub elements: Option<Vec<Element>>,
}

impl YamlTable for ElmTab {
    type Table = ElmTab;
    type TableItem = Element;
    fn load_table<P: AsRef<Path>>(filepath: P) -> Result<Self::Table, Box<dyn std::error::Error>> {
        let yaml_table = fs::File::open(filepath)?;
        let table: ElmTab = serde_yaml::from_reader(yaml_table)?;
        Ok(table)
    }
    fn hash_table<P: AsRef<Path>>(
        filepath: P,
    ) -> Result<HashMap<String, Self::TableItem>, Box<dyn Error>> {
        let table = Self::load_table(filepath)?;
        let mut hash_tab: HashMap<String, Element> = HashMap::new();
        table.elements.unwrap().iter().for_each(|elm: &Element| {
            hash_tab.insert(elm.element.to_string(), elm.deref().clone());
        });
        Ok(hash_tab)
    }
}
