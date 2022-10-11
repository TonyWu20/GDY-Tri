use crate::{lattice::Lattice, molecule::Adsorbate};

pub trait AdsAddition {
    fn append_mol_name(&mut self, ads: Adsorbate, site_1: u32, site_2: u32);
    fn init_ads_direction(&self, ads: Adsorbate, site_1: u32, site_2: u32);
    fn add_ads(&self, ads: Adsorbate, site_1: u32, site_2: u32, height: f64);
}
