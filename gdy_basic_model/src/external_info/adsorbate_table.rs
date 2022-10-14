/// Parse adsorbate info `yaml` file to get necessary information for loading adsorbate models.
extern crate serde;
use std::{
    collections::HashMap,
    error::Error,
    fs,
    ops::Deref,
    path::{Path, PathBuf},
    str::FromStr,
};

use serde::{Deserialize, Serialize};

use super::YamlTable;

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct AdsInfo {
    name: String,
    #[serde(rename = "coordAtomIds")]
    coord_atom_ids: Vec<u32>,
    #[serde(rename = "stemAtomIds")]
    stem_atom_ids: Vec<u32>,
    #[serde(rename = "planeAtomIds")]
    plane_atom_ids: Vec<u32>,
    vertical: bool,
    #[serde(rename = "bSym")]
    symmetric: bool,
    #[serde(rename = "upperAtomId")]
    upper_atom_id: u32,
    #[serde(rename = "atomNums")]
    atom_nums: u32,
    #[serde(rename = "pathName")]
    path_name: String,
}

impl AdsInfo {
    pub fn name(&self) -> &str {
        self.name.as_ref()
    }

    pub fn coord_atom_ids(&self) -> &[u32] {
        self.coord_atom_ids.as_ref()
    }

    pub fn stem_atom_ids(&self) -> &[u32] {
        self.stem_atom_ids.as_ref()
    }

    pub fn plane_atom_ids(&self) -> &[u32] {
        self.plane_atom_ids.as_ref()
    }

    pub fn vertical(&self) -> bool {
        self.vertical
    }

    pub fn symmetric(&self) -> bool {
        self.symmetric
    }

    pub fn upper_atom_id(&self) -> u32 {
        self.upper_atom_id
    }

    pub fn atom_nums(&self) -> u32 {
        self.atom_nums
    }

    pub fn path_name(&self) -> &str {
        self.path_name.as_ref()
    }
    pub fn file_path(&self, parent_dir: &str) -> Result<PathBuf, std::convert::Infallible> {
        PathBuf::from_str(&format!(
            "{}/{}_path/{}.msi",
            parent_dir,
            self.path_name(),
            self.name()
        ))
    }
}

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct AdsTab {
    directory: String,
    #[serde(rename = "Adsorbates")]
    adsorbates: Option<Vec<AdsInfo>>,
}

impl AdsTab {
    pub fn adsorbates(&self) -> Option<&Vec<AdsInfo>> {
        self.adsorbates.as_ref()
    }

    pub fn directory(&self) -> &str {
        self.directory.as_ref()
    }
}

impl YamlTable for AdsTab {
    type Table = AdsTab;

    type TableItem = AdsInfo;

    fn load_table<P: AsRef<Path>>(filepath: P) -> Result<Self::Table, Box<dyn Error>> {
        let ads_table_src = fs::File::open(filepath)?;
        let ads_table = serde_yaml::from_reader(ads_table_src)?;
        Ok(ads_table)
    }

    fn hash_table(
        &self,
    ) -> Result<std::collections::HashMap<String, Self::TableItem>, Box<dyn Error>> {
        let mut hash_tab: HashMap<String, AdsInfo> = HashMap::new();
        let adsinfo_vec: &Vec<AdsInfo> = self.adsorbates.as_ref().unwrap();
        adsinfo_vec.iter().for_each(|ads: &AdsInfo| {
            hash_tab.insert(ads.name.to_string(), ads.deref().clone());
        });
        Ok(hash_tab)
    }
}

pub fn load_ads_table() -> Result<AdsTab, Box<dyn Error>> {
    let ads_table_src = fs::File::open("./resources/ads_table.yaml")?;
    let ads_table = serde_yaml::from_reader(ads_table_src)?;
    Ok(ads_table)
}
