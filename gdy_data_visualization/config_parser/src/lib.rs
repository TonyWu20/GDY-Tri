#![allow(dead_code)]
use std::{collections::HashMap, fmt::Display};

/**
Define Structs to deserialize config toml files to get necessary parameters for data processing
*/
use serde::Deserialize;

#[derive(Deserialize, Debug)]
#[serde(rename_all = "camelCase")]
pub enum Part {
    Elements { element_symbols: Vec<String> },
    Sites(Sites),
    Adsorbates(AdsorbatesMap),
}

#[derive(Deserialize, Debug)]
pub struct EnergyConfig {
    dir_prefix: String,
    element_symbols: Vec<String>,
    sites: Sites,
    adsorbates: AdsorbatesMap,
    pathways: Pathways,
}

impl EnergyConfig {
    pub fn sites(&self) -> &Sites {
        &self.sites
    }

    pub fn adsorbates(&self) -> &AdsorbatesMap {
        &self.adsorbates
    }

    pub fn pathways(&self) -> &Pathways {
        &self.pathways
    }

    pub fn dir_prefix(&self) -> &str {
        self.dir_prefix.as_ref()
    }

    pub fn element_symbols(&self) -> &[String] {
        self.element_symbols.as_ref()
    }

    pub fn construct_paths(&self) -> Vec<String> {
        todo!();
    }
}

#[derive(Deserialize, Debug)]
pub struct Sites {
    site_names: Vec<Vec<String>>,
    site_series: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct SitesLengthError;

impl Display for SitesLengthError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Inconsistent lengths of site_names and site_series")
    }
}

impl Sites {
    pub fn site_names(&self) -> &[Vec<String>] {
        self.site_names.as_ref()
    }

    pub fn site_series(&self) -> &[String] {
        self.site_series.as_ref()
    }
    /// Output all sites in a flattened vector.
    pub fn all_sites(&self) -> Vec<String> {
        self.site_names.clone().into_iter().flatten().collect()
    }
    pub fn hashmap(&self) -> Result<HashMap<String, Vec<String>>, SitesLengthError> {
        if self.site_names().len() != self.site_series().len() {
            return Err(SitesLengthError);
        }
        let mut hash_map: HashMap<String, Vec<String>> = HashMap::new();
        self.site_series()
            .iter()
            .zip(self.site_names())
            .for_each(|(site_name, sites)| {
                hash_map.insert(site_name.to_string(), sites.clone());
            });
        Ok(hash_map)
    }
}

#[derive(Deserialize, Debug)]
pub struct AdsorbatesMap {
    adsorbates: Vec<Adsorbate>,
}

impl AdsorbatesMap {
    pub fn adsorbates(&self) -> &[Adsorbate] {
        self.adsorbates.as_ref()
    }
    pub fn ads_name_site_hashmap(&self) -> HashMap<String, String> {
        let mut hash_map: HashMap<String, String> = HashMap::new();
        self.adsorbates.iter().for_each(|ads| {
            hash_map.insert(ads.name().to_owned(), ads.sites().to_owned());
        });
        hash_map
    }
}

#[derive(Deserialize, Debug)]
pub struct Adsorbate {
    name: String,
    #[serde(rename = "sites")]
    site_series: String,
}

impl Adsorbate {
    pub fn name(&self) -> &str {
        self.name.as_ref()
    }

    pub fn sites(&self) -> &str {
        self.site_series.as_ref()
    }
}
#[derive(Deserialize, Debug)]
pub struct Pathways {
    item: Vec<PathwayItem>,
}

impl Pathways {
    pub fn item(&self) -> &[PathwayItem] {
        self.item.as_ref()
    }
}
#[derive(Deserialize, Debug)]
pub struct PathwayItem {
    name: String,
    path: Vec<String>,
}

impl PathwayItem {
    pub fn name(&self) -> &str {
        self.name.as_ref()
    }

    pub fn path(&self) -> &[String] {
        self.path.as_ref()
    }
}

#[cfg(test)]
#[test]
fn test_read_config() {
    use std::fs;

    let config_text = fs::read_to_string("../config.toml").unwrap();
    let config: Result<EnergyConfig, toml::de::Error> = toml::from_str(&config_text);
    assert!(config.is_ok());
    let pathways = config.as_ref().unwrap().pathways();
    let ads_hashmap = config
        .as_ref()
        .unwrap()
        .adsorbates()
        .ads_name_site_hashmap();
    pathways.item().iter().for_each(|pathway| {
        pathway.path().iter().for_each(|item| {
            println!(
                "ads: {}, site series: {}",
                item,
                ads_hashmap.get(item).unwrap()
            )
        })
    })
}
