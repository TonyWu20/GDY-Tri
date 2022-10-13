use std::{collections::HashMap, error::Error};

use crate::molecule::adsorbate::Adsorbate;

pub trait AdsAddition {
    fn append_mol_name(
        &mut self,
        ads: &Adsorbate,
        site_1: u32,
        site_2: Option<u32>,
        coord_site_dict: &HashMap<u32, String>,
    ) -> Result<(), Box<dyn Error>>;
    fn init_ads_direction(
        &self,
        ads: &mut Adsorbate,
        site_1: u32,
        site_2: Option<u32>,
        flip_upright: bool,
    ) -> Result<(), Box<dyn Error>>;
    fn add_ads(
        &mut self,
        ads: &mut Adsorbate,
        site_1: u32,
        site_2: Option<u32>,
        height: f64,
        flip_upright: bool,
        coord_site_dict: &HashMap<u32, String>,
    ) -> Result<(), Box<dyn Error>>;
}
