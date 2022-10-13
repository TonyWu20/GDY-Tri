use std::{collections::HashMap, error::Error, path::Path};

pub mod adsorbate_table;
pub mod element_table;
pub mod project;

pub trait YamlTable {
    type Table;
    type TableItem;
    fn load_table<P: AsRef<Path>>(filepath: P) -> Result<Self::Table, Box<dyn Error>>;
    fn hash_table<P: AsRef<Path>>(
        filepath: P,
    ) -> Result<HashMap<String, Self::TableItem>, Box<dyn Error>>;
}

#[test]
fn test_yaml() -> Result<(), Box<dyn Error>> {
    use crate::external_info::adsorbate_table::AdsTab;
    use crate::external_info::element_table::ElmTab;
    use crate::external_info::project::load_project_info;
    let project_info = load_project_info("resources/project.yaml")?;
    let table = ElmTab::load_table(project_info.element_table_loc()).unwrap();
    println!("{}", table.elements.is_some());
    table
        .elements
        .unwrap()
        .iter()
        .for_each(|elm| println!("{:#?}", elm));
    let hashtab = ElmTab::hash_table(project_info.element_table_loc())?;
    println!("{:#?}", hashtab.get("C").unwrap());
    let ads_table = AdsTab::load_table(project_info.adsorbate_table_loc())?;
    assert!(ads_table.adsorbates().is_some());
    let hash_ads_tab = AdsTab::hash_table(project_info.adsorbate_table_loc())?;
    assert_eq!(
        vec![1, 2, 3],
        hash_ads_tab.get("CH2CHOH").unwrap().plane_atom_ids()
    );
    println!("{:?}", project_info);
    let coord_sites_dict = project_info.hash_coord_site();
    println!("{}", coord_sites_dict.get(&41).unwrap());
    Ok(())
}
