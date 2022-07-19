use super::config_parser::EnergyConfig;

/// Organizations of all energy data
trait OrganizeEnergy {
    /// Decide type of output. Possibly: Vec<f64> (flatten array), .csv, .dat (polar coord for gnuplot)
    type Output;
    /**
    - basic unit: energy per site
    - theta: chunks of sectors: adsorbate, each unit sector: e_ads, radial: stack of pathway branches of each elements
    */
    fn sites_pathway_element_view(&self) -> Self::Output;
    /**
    - basic unit: energy per site,
    - theta: chunks of sectors: elements, each unit sector: pathway
    - radial: along pathway, greatest common length
    */
    fn sites_element_pathway_view(&self) -> Self::Output;
}
