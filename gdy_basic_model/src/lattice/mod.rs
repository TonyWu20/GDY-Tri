use periodic_table_on_an_enum as pt_enum;
use std::collections::HashMap;
use std::{collections::HashSet, error::Error};

use na::{Matrix3, Translation3, Unit, UnitQuaternion, Vector, Vector3};

use crate::molecule::adsorbate::Adsorbate;
use crate::{
    assemble::AdsAddition,
    atom::{Atom, AtomArrayRef},
    Export, Transformation,
};

#[derive(Clone)]
pub struct Lattice {
    lattice_name: String,
    atoms_vec: Vec<Atom>,
    lattice_vectors: Matrix3<f64>,
    metal_sites: Vec<u32>,
    adsorbate: Option<String>,
    sorted: bool,
    pathway: Option<String>,
}

impl Lattice {
    pub fn new(
        lattice_name: String,
        atoms_vec: Vec<Atom>,
        lattice_vectors: Matrix3<f64>,
        metal_sites: Vec<u32>,
        adsorbate: Option<String>,
        sorted: bool,
    ) -> Self {
        Self {
            lattice_name,
            atoms_vec,
            lattice_vectors,
            metal_sites,
            adsorbate,
            sorted,
            pathway: None,
        }
    }

    pub fn get_lattice_vectors(&self) -> &Matrix3<f64> {
        &self.lattice_vectors
    }
    pub fn set_lattice_vectors(&mut self, new_lattice_vec: Matrix3<f64>) {
        self.lattice_vectors = new_lattice_vec;
    }
    pub fn get_metal_sites(&self) -> &Vec<u32> {
        &self.metal_sites
    }
    // Get element list sorted by atomic number
    pub fn get_element_list(&self) -> Vec<String> {
        let mut elm_list: Vec<String> = vec![];
        elm_list.extend(
            self.atoms_vec
                .iter()
                .map(|atom| atom.element_name().to_string())
                .collect::<Vec<String>>()
                .drain(..)
                .collect::<HashSet<String>>()
                .into_iter(),
        );
        elm_list.sort_unstable_by(|a, b| {
            pt_enum::Element::from_symbol(&a)
                .unwrap()
                .get_atomic_number()
                .cmp(
                    &pt_enum::Element::from_symbol(&b)
                        .unwrap()
                        .get_atomic_number(),
                )
        });
        elm_list
    }
    pub fn get_adsorbate_name(&self) -> Option<&str> {
        match &self.adsorbate {
            Some(x) => Some(x.as_str()),
            None => None,
        }
    }
    pub fn set_adsorbate_name(&mut self, ads_name: String) {
        self.adsorbate = Some(ads_name);
    }
    pub fn rotate_to_standard_orientation(&mut self) {
        let x_axis: Vector3<f64> = Vector::x();
        let a_vec = &self.get_lattice_vectors().column(0);
        let a_to_x_angle: f64 = a_vec.angle(&x_axis);
        if a_to_x_angle == 0.0 {
            return;
        }
        let rot_axis = a_vec.cross(&x_axis).normalize();
        let rot_quatd: UnitQuaternion<f64> = UnitQuaternion::new(rot_axis * a_to_x_angle);
        self.rotate(rot_quatd);
    }
    // pub fn update_base_name(&mut self) {
    //     // Collect all metal's symbols
    //     let metal_names: Vec<String> = self
    //         .metal_sites
    //         .iter()
    //         .map(|metal_id| -> String {
    //             self.atoms_vec
    //                 .get_atom_by_id(*metal_id)
    //                 .unwrap()
    //                 .element_name()
    //                 .to_string()
    //         })
    //         .collect::<Vec<String>>();
    //     // Because we have only 3 metal elements
    //     let new_name = format!(
    //         "GDY_{}_{}_{}",
    //         metal_names[0], metal_names[1], metal_names[2]
    //     );
    //     self.set_lattice_name(new_name);
    // }

    pub fn lattice_name(&self) -> &str {
        self.lattice_name.as_ref()
    }

    pub fn set_lattice_name(&mut self, lattice_name: String) {
        self.lattice_name = lattice_name;
    }

    pub fn atoms_vec(&self) -> &[Atom] {
        self.atoms_vec.as_ref()
    }

    pub fn atoms_vec_mut(&mut self) -> &mut Vec<Atom> {
        &mut self.atoms_vec
    }

    pub fn sorted(&self) -> bool {
        self.sorted
    }

    pub fn set_sorted(&mut self, sorted: bool) {
        self.sorted = sorted;
    }
    pub fn sort_atoms_by_elements(&mut self) {
        self.atoms_vec_mut().sort();
        self.set_sorted(true);
    }
    pub fn lattice_vector_str(&self) -> (String, String) {
        let vectors = self.get_lattice_vectors();
        let mut vectors_string = String::new();
        vectors.column_iter().for_each(|col| {
            vectors_string.push_str(&format!("{:24.18}{:24.18}{:24.18}\n", col.x, col.y, col.z));
        });
        ("LATTICE_CART".to_string(), vectors_string)
    }
    pub fn fractional_coord_matrix(&self) -> Matrix3<f64> {
        let lattice_vectors = self.get_lattice_vectors();
        let vec_a = lattice_vectors.column(0);
        let vec_b = lattice_vectors.column(1);
        let vec_c = lattice_vectors.column(2);
        let len_a: f64 = vec_a.norm();
        let len_b: f64 = vec_b.norm();
        let len_c: f64 = vec_c.norm();
        let (alpha, beta, gamma) = (
            vec_b.angle(&vec_c),
            vec_a.angle(&vec_c),
            vec_a.angle(&vec_b),
        );
        let vol = vec_a.dot(&vec_b.cross(&vec_c));
        let to_cart = Matrix3::new(
            len_a,
            len_b * gamma.cos(),
            len_c * beta.cos(),
            0.0,
            len_b * gamma.sin(),
            len_c * (alpha.cos() - beta.cos() * gamma.cos()) / gamma.sin(),
            0.0,
            0.0,
            vol / (len_a * len_b * gamma.sin()),
        );
        let to_frac = to_cart.try_inverse().unwrap();
        to_frac
    }
    pub fn add_atoms(&mut self, new_atoms: &[Atom]) -> Result<(), Box<dyn Error>> {
        let lat_last_atom_id = self
            .atoms_vec()
            .get_atom_by_id(self.atoms_vec().len() as u32)?
            .atom_id();
        new_atoms.to_vec().into_iter().for_each(|mut atom| {
            atom.set_atom_id(lat_last_atom_id + atom.atom_id());
            self.atoms_vec_mut().push(atom);
        });
        Ok(())
    }

    pub fn pathway(&self) -> Option<&String> {
        self.pathway.as_ref()
    }

    pub fn set_pathway(&mut self, pathway: Option<String>) {
        self.pathway = pathway;
    }
}

impl Export for Lattice {
    fn format_output(&self) -> String {
        let lattice_vectors = self.get_lattice_vectors();
        let vec_a_line = format!(
            "  (A D A3 ({:.12} {:.12} {:.12}))\n",
            lattice_vectors.column(0).x,
            lattice_vectors.column(0).y,
            lattice_vectors.column(0).z
        );
        let vec_b_line = format!(
            "  (A D B3 ({:.12} {:.12} {:.12}))\n",
            lattice_vectors.column(1).x,
            lattice_vectors.column(1).y,
            lattice_vectors.column(1).z
        );
        let vec_c_line = format!(
            "  (A D C3 ({:.12} {:.12} {:.12}))\n",
            lattice_vectors.column(2).x,
            lattice_vectors.column(2).y,
            lattice_vectors.column(2).z
        );
        let headers: String = vec![
            "# MSI CERIUS2 DataModel File Version 4 0\n",
            "(1 Model\n",
            "  (A I CRY/DISPLAY (192 256))\n",
            "  (A I PeriodicType 100)\n",
            "  (A C SpaceGroup \"1 1\")\n",
            &vec_a_line,
            &vec_b_line,
            &vec_c_line,
            "  (A D CRY/TOLERANCE 0.05)\n",
        ]
        .join("");
        let atom_strings: String = self.atoms_vec.format_output();
        let contents: String = format!("{}{})", headers, atom_strings);
        contents
    }
}
impl Transformation for Lattice {
    fn rotate(&mut self, rotate_quatd: na::UnitQuaternion<f64>) {
        self.atoms_vec.rotate(rotate_quatd);
        let rotation_matrix = rotate_quatd.to_rotation_matrix();
        let new_lat_vec: Matrix3<f64> = rotation_matrix * self.get_lattice_vectors();
        self.set_lattice_vectors(new_lat_vec);
    }
    fn translate(&mut self, translate_matrix: na::Translation<f64, 3>) {
        self.atoms_vec.translate(translate_matrix);
    }
}

impl AdsAddition for Lattice {
    fn append_mol_name(
        &mut self,
        ads: &Adsorbate,
        site_1: u32,
        site_2: Option<u32>,
        coord_site_dict: &HashMap<u32, String>,
    ) -> Result<(), Box<dyn Error>> {
        let new_name = match site_2 {
            Some(site_2) => {
                let suf_1 = coord_site_dict.get(&site_1).unwrap().to_owned();
                let suf_2 = coord_site_dict.get(&site_2).unwrap().to_owned();
                format!(
                    "{}_{}_{}_{}",
                    self.lattice_name(),
                    ads.mol_name(),
                    suf_1,
                    suf_2
                )
            }
            None => {
                format!(
                    "{}_{}_{}",
                    self.lattice_name(),
                    ads.mol_name(),
                    coord_site_dict.get(&site_1).unwrap().to_owned()
                )
            }
        };
        self.set_lattice_name(new_name);
        Ok(())
    }

    fn init_ads_direction(
        &self,
        ads: &mut Adsorbate,
        site_1: u32,
        site_2: Option<u32>,
        flip_upright: bool,
    ) -> Result<(), Box<dyn Error>> {
        let ads_stem_vec = ads.get_stem_vector()?;
        let (target_site_1, target_site_2) = match site_2 {
            Some(site_2) => (site_1, site_2),
            None => (41, 42),
        };
        let direction_vec = self
            .atoms_vec()
            .get_vector_ab(target_site_1, target_site_2)?;
        let angle_stem_chain = ads_stem_vec.angle(&direction_vec);
        let rot_axis = Unit::new_normalize(ads_stem_vec.cross(&direction_vec));
        let rot_quatd = UnitQuaternion::from_axis_angle(&rot_axis, angle_stem_chain);
        ads.atoms_vec_mut().rotate(rot_quatd);
        if flip_upright == true {
            ads.make_upright()?;
        }
        Ok(())
    }

    fn add_ads(
        &mut self,
        ads: &Adsorbate,
        site_1: u32,
        site_2: Option<u32>,
        height: f64,
        flip_upright: bool,
        coord_site_dict: &HashMap<u32, String>,
    ) -> Result<(), Box<dyn Error>> {
        let mut ads_clone = ads.clone();
        Self::init_ads_direction(&self, &mut ads_clone, site_1, site_2, flip_upright)?;
        /*
        If site_2 exists or the adsorbate has two coordination atoms,
        the coordinate sites follows the adsorbate info.
        If site_2 does not exists and the adsorbate has one coordinate atom,
        both sites are assigned to the only one coordinate atom.
        */
        let (cd_1, cd_2) = if site_2.is_some() || ads_clone.coord_atom_nums() == 2 {
            (ads_clone.coord_atom_ids()[0], ads_clone.coord_atom_ids()[1])
        } else {
            (ads_clone.coord_atom_ids()[0], ads_clone.coord_atom_ids()[0])
        };
        let (lat_site_1, lat_site_2) = if site_2.is_some() || ads.coord_atom_nums() == 2 {
            (site_1, site_2.unwrap())
        } else {
            (site_1, site_1)
        };
        // Get the center coordinates of the two carbon sites.
        let lat_sites = (
            self.atoms_vec().get_atom_by_id(lat_site_1)?.xyz().clone(),
            self.atoms_vec().get_atom_by_id(lat_site_2)?.xyz().clone(),
        );
        let cd_sites = (
            ads_clone.atoms_vec().get_atom_by_id(cd_1)?.xyz().clone(),
            ads_clone.atoms_vec().get_atom_by_id(cd_2)?.xyz().clone(),
        );
        let lat_sites_centroid = na::center(&lat_sites.0, &lat_sites.1);
        let cd_sites_centroid = na::center(&cd_sites.0, &cd_sites.1);
        let mut trans_matrix = Translation3::from(lat_sites_centroid - cd_sites_centroid);
        trans_matrix.vector.z += height;
        ads_clone.atoms_vec_mut().translate(trans_matrix);
        self.add_atoms(ads_clone.atoms_vec())?;
        self.set_adsorbate_name(ads_clone.mol_name().to_string());
        self.set_pathway(Some(ads_clone.pathway_name().to_string()));
        self.append_mol_name(&ads_clone, site_1, site_2, coord_site_dict)?;
        Ok(())
    }
}
