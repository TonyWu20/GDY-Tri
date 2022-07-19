/**
Define Structs to deserialize config toml files to get necessary parameters for data processing
*/
use serde::Deserialize;

#[derive(Deserialize, Debug)]
pub struct EnergyConfig {
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
}

#[derive(Deserialize, Debug)]
pub struct Sites {
    site_names: Vec<Vec<String>>,
    site_series: Vec<String>,
}

impl Sites {
    pub fn site_names(&self) -> &[Vec<String>] {
        self.site_names.as_ref()
    }

    pub fn site_series(&self) -> &[String] {
        self.site_series.as_ref()
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
}

#[derive(Deserialize, Debug)]
pub struct Adsorbate {
    name: String,
    sites: String,
}

impl Adsorbate {
    pub fn name(&self) -> &str {
        self.name.as_ref()
    }

    pub fn sites(&self) -> &str {
        self.sites.as_ref()
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

    let config_text = fs::read_to_string("config.toml").unwrap();
    let config: Result<EnergyConfig, toml::de::Error> = toml::from_str(&config_text);
    assert!(config.is_ok())
}
