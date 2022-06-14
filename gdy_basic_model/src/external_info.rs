pub mod element_table {
    extern crate serde;
    use serde::{Deserialize, Serialize};
    use std::{collections::HashMap, fs, ops::Deref};

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
    pub struct ElmInfo {
        #[serde(rename = "Element_info")]
        pub elements: Option<Vec<Element>>,
    }

    pub fn load_table() -> Result<ElmInfo, serde_yaml::Error> {
        let yaml_table = fs::File::open("/Users/tonywu/Library/Mobile Documents/com~apple~CloudDocs/Programming/GDY-Tri/gdy_basic_model/resources/element_table.yaml")
            .expect("Something wrong in reading element_table.yaml");
        let table: ElmInfo = serde_yaml::from_reader(yaml_table)?;
        Ok(table)
    }
    pub fn hash_table() -> HashMap<String, Element> {
        let table = load_table().unwrap();
        let mut hash_tab: HashMap<String, Element> = HashMap::new();
        table.elements.unwrap().iter().for_each(|elm: &Element| {
            hash_tab.insert(elm.element.to_string(), elm.deref().clone());
        });
        hash_tab
    }
}

#[test]
fn test_yaml() -> Result<(), serde_yaml::Error> {
    use crate::external_info::element_table::hash_table;
    let table = self::element_table::load_table().unwrap();
    println!("{}", table.elements.is_some());
    table
        .elements
        .unwrap()
        .iter()
        .for_each(|elm| println!("{:#?}", elm));
    let hashtab = hash_table();
    println!("{:#?}", hashtab.get("C").unwrap());
    Ok(())
}
