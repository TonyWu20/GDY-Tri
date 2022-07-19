use std::f64::consts::PI;

use polars::prelude::*;
use polars_lazy::prelude::*;

use super::{export_format::ExportGnuData, plot_data_struct::PolarMeshData};

const PATHWAY_CH3OH_1: [&str; 7] = ["CO2", "COOH", "CO", "CHO", "CH2O", "CH3O", "CH3OH"];
const PATHWAY_CH3OH_2: [&str; 7] = ["CO2", "COOH", "CO", "CHO", "CHOH", "CH2OH", "CH3OH"];
const PATHWAY_CH4_1: [&str; 9] = [
    "CO2", "COOH", "CO", "COH", "CHOH", "CH", "CH2", "CH3", "CH4",
];
const PATHWAY_CH4_2: [&str; 9] = ["CO2", "COOH", "CO", "COH", "C", "CH", "CH2", "CH3", "CH4"];
const PATHWAY_HCOOH: [&str; 3] = ["CO2", "HCOO", "HCOOH"];
const MARK: f64 = 0.0;
fn get_columns_with_prefix(csv_name: &str, column_kw: &str) -> Result<DataFrame> {
    let col_kw = format!("^{}_.*$", column_kw);
    let col = LazyCsvReader::new(csv_name.into())
        .finish()
        .unwrap()
        .select(&[col(&col_kw)])
        .collect();
    col
}
fn get_columns_with_prefix_n_elm(csv_name: &str, column_kw: &str, elm: &str) -> Result<DataFrame> {
    let col_kw = format!("^{}_.*$", column_kw);
    let col = LazyCsvReader::new(csv_name.into())
        .finish()
        .unwrap()
        .filter(col("SAC_GDY_X").eq(lit(elm)))
        .select(&[col(&col_kw)])
        .collect();
    col
}

fn get_cols_by_pathway_n_elm(csv_name: &str, pathway_arr: &[&str], elm: &str) -> Result<DataFrame> {
    let ads_columns = pathway_arr
        .iter()
        .map(|ads| col(&format!("^{}_.*$", ads)))
        .collect::<Vec<Expr>>();
    // let mut selections = vec![col("SAC_GDY_X")];
    // selections.append(&mut ads_columns);
    LazyCsvReader::new(csv_name.to_string())
        .finish()
        .unwrap()
        .filter(col("SAC_GDY_X").eq(lit(elm)))
        .select(ads_columns)
        .collect()
}

fn get_energy_by_pathway_site_elm(
    csv_name: &str,
    pathway_arr: &[&str],
    elm: &str,
    site_name: &str,
) -> Result<DataFrame> {
    let ads_site_columns: Vec<Expr> = pathway_arr
        .iter()
        .map(|ads| col(&format!("{}_{}", ads, site_name)))
        .collect();
    LazyCsvReader::new(csv_name.to_string())
        .finish()
        .unwrap()
        .filter(col("SAC_GDY_X").eq(lit(elm)))
        .select(ads_site_columns)
        .collect()
}

fn compute_pathway_eneregy_for_elm(
    csv_name: &str,
    pathway_arr: &[&str],
    elm: &str,
) -> Result<DataFrame> {
    let ads_columns = pathway_arr
        .iter()
        .map(|ads| col(&format!("^{}_.*$", ads)))
        .collect::<Vec<Expr>>();
    // let mut selections = vec![col("SAC_GDY_X")];
    // selections.append(&mut ads_columns);
    LazyCsvReader::new(csv_name.to_string())
        .finish()
        .unwrap()
        .filter(col("SAC_GDY_X").eq(lit(elm)))
        .select(ads_columns)
        .collect()
}

fn reshape(csv_name: &str) -> Result<DataFrame> {
    let _elements = LazyCsvReader::new(csv_name.to_string())
        .finish()
        .unwrap()
        .select(&[col("SAC_GDY_X")])
        .collect()
        .unwrap();
    todo!();
}
pub fn dataframe_to_vec_f64(df: DataFrame) -> Vec<f64> {
    df.get_row(0)
        .0
        .into_iter()
        .map(|v| {
            v.to_string()
                .replace("\\", "")
                .parse::<f64>()
                .unwrap_or_else(|e| panic!("parsing {} error, {}", v, e))
        })
        .collect::<Vec<f64>>()
}

pub fn c1_routine_pathway_rad() {
    let csv_file = "gdy_c1.csv";
    let elements = LazyCsvReader::new(csv_file.to_string())
        .finish()
        .unwrap()
        .select(&[col("SAC_GDY_X")])
        .collect()
        .unwrap();
    let elements: Vec<String> = elements
        .column("SAC_GDY_X")
        .unwrap()
        .iter()
        .map(|val| val.to_string().replace("\\", "").replace(r#"""#, ""))
        .collect();
    let site_names = ["c1", "c2", "c3", "c4", "c5", "metal"];
    let pathways: Vec<Vec<&str>> = vec![
        PATHWAY_CH3OH_1.to_vec(),
        PATHWAY_CH3OH_2.to_vec(),
        PATHWAY_CH4_1.to_vec(),
        PATHWAY_CH4_2.to_vec(),
        PATHWAY_HCOOH.to_vec(),
    ];
    let pathway_per_site_elm_data: Vec<f64> = elements
        .iter()
        // Three loops, inner to outer: site, pathway, element
        // The vectors are aligned in length
        .flat_map(|elm| -> Vec<f64> {
            pathways
                .iter()
                .enumerate()
                .flat_map(|(i, pathway)| -> Vec<f64> {
                    site_names
                        .iter()
                        .flat_map(|site| -> Vec<f64> {
                            match i {
                                0..=1 => {
                                    let mut data = dataframe_to_vec_f64(
                                        get_energy_by_pathway_site_elm(
                                            csv_file, &pathway, &elm, &site,
                                        )
                                        .unwrap(),
                                    );
                                    let size = data.len();
                                    let (f1, f2) = (data[size - 2], data[size - 1]);
                                    data[size - 1] = f1;
                                    data.push(f2);
                                    data.push(f2);
                                    let base_energy = data[0];
                                    let pathway_energy: Vec<f64> =
                                        data.iter().map(|&v| v - base_energy).collect();
                                    let size = pathway_energy.len();
                                    pathway_energy[1..size].to_vec()
                                }
                                2..=3 => {
                                    let data = dataframe_to_vec_f64(
                                        get_energy_by_pathway_site_elm(
                                            csv_file, &pathway, &elm, &site,
                                        )
                                        .unwrap(),
                                    );
                                    let base_energy = data[0];
                                    let pathway_energy: Vec<f64> =
                                        data.iter().map(|&v| v - base_energy).collect();
                                    let size = pathway_energy.len();
                                    pathway_energy[1..size].to_vec()
                                }
                                4 => {
                                    let data = dataframe_to_vec_f64(
                                        get_energy_by_pathway_site_elm(
                                            csv_file, &pathway, &elm, &site,
                                        )
                                        .unwrap(),
                                    );
                                    (0..8)
                                        .into_iter()
                                        .map(|i| match i {
                                            0..=3 => data[1] - data[0],
                                            4..=8 => data[2] - data[0],
                                            _ => MARK,
                                        })
                                        .collect::<Vec<f64>>()
                                }
                                _ => vec![0.0],
                            }
                        })
                        .collect::<Vec<f64>>()
                })
                .collect()
        })
        .collect();
    let num_sectors = 44 * 5 * 6;
    let polar_mesh_data = PolarMeshData::new(
        (0.0, 7.0),
        8,
        (0.0, PI * 2.0),
        num_sectors as usize,
        pathway_per_site_elm_data,
    );
    polar_mesh_data.to_gnu_data("c1_pr_eth.dat").unwrap();
}

#[cfg(test)]
#[test]
fn test_pathway_rad() {
    use std::f64::consts::PI;

    use crate::data_processing::{export_format::ExportGnuData, plot_data_struct::PolarMeshData};

    let csv_file = "gdy_c1.csv";
    let elements = LazyCsvReader::new(csv_file.to_string())
        .finish()
        .unwrap()
        .select(&[col("SAC_GDY_X")])
        .collect()
        .unwrap();
    let elements: Vec<String> = elements
        .column("SAC_GDY_X")
        .unwrap()
        .iter()
        .map(|val| val.to_string().replace("\\", "").replace(r#"""#, ""))
        .collect();
    let site_names = ["c1", "c2", "c3", "c4", "c5", "metal"];
    let pathways: Vec<Vec<&str>> = vec![
        PATHWAY_CH3OH_1.to_vec(),
        PATHWAY_CH3OH_2.to_vec(),
        PATHWAY_CH4_1.to_vec(),
        PATHWAY_CH4_2.to_vec(),
        PATHWAY_HCOOH.to_vec(),
    ];
    let pathway_per_site_elm_data: Vec<f64> = elements
        .iter()
        // Three loops, inner to outer: site, pathway, element
        // The vectors are aligned in length
        .flat_map(|elm| -> Vec<f64> {
            pathways
                .iter()
                .enumerate()
                .flat_map(|(i, pathway)| -> Vec<f64> {
                    site_names
                        .iter()
                        .flat_map(|site| -> Vec<f64> {
                            match i {
                                0..=1 => {
                                    let mut data = dataframe_to_vec_f64(
                                        get_energy_by_pathway_site_elm(
                                            csv_file, &pathway, &elm, &site,
                                        )
                                        .unwrap(),
                                    );
                                    let size = data.len();
                                    let (f1, f2) = (data[size - 2], data[size - 1]);
                                    data[size - 1] = f1;
                                    data.push(f2);
                                    data.push(f2);
                                    let base_energy = data[0];
                                    let pathway_energy: Vec<f64> =
                                        data.iter().map(|&v| v - base_energy).collect();
                                    let size = pathway_energy.len();
                                    pathway_energy[1..size].to_vec()
                                }
                                2..=3 => {
                                    let data = dataframe_to_vec_f64(
                                        get_energy_by_pathway_site_elm(
                                            csv_file, &pathway, &elm, &site,
                                        )
                                        .unwrap(),
                                    );
                                    let base_energy = data[0];
                                    let pathway_energy: Vec<f64> =
                                        data.iter().map(|&v| v - base_energy).collect();
                                    let size = pathway_energy.len();
                                    pathway_energy[1..size].to_vec()
                                }
                                4 => {
                                    let data = dataframe_to_vec_f64(
                                        get_energy_by_pathway_site_elm(
                                            csv_file, &pathway, &elm, &site,
                                        )
                                        .unwrap(),
                                    );
                                    (0..8)
                                        .into_iter()
                                        .map(|i| match i {
                                            0..=3 => data[1] - data[0],
                                            4..=8 => data[2] - data[0],
                                            _ => MARK,
                                        })
                                        .collect::<Vec<f64>>()
                                }
                                _ => vec![0.0],
                            }
                        })
                        .collect::<Vec<f64>>()
                })
                .collect()
        })
        .collect();
    let num_sectors = 44 * 5 * 6;
    let polar_mesh_data = PolarMeshData::new(
        (0.0, 8.0),
        8,
        (0.0, PI * 2.0),
        num_sectors as usize,
        pathway_per_site_elm_data,
    );
    polar_mesh_data.to_gnu_data("c1_pr_eth.dat").unwrap();
}
#[test]
fn test_use() {
    let csv_file = "gdy_c1.csv";
    let elements = LazyCsvReader::new(csv_file.to_string())
        .finish()
        .unwrap()
        .select(&[col("SAC_GDY_X")])
        .collect()
        .unwrap();
    let elements: Vec<String> = elements
        .column("SAC_GDY_X")
        .unwrap()
        .iter()
        .map(|val| val.to_string().replace("\\", "").replace(r#"""#, ""))
        .collect();
    println!("{:?}", elements);
    let all_data_flatten: Vec<f64> = elements
        .iter()
        .flat_map(|elm| -> Vec<f64> {
            let p1 = get_cols_by_pathway_n_elm(csv_file, &PATHWAY_CH3OH_1, &elm).unwrap();
            let p2 = get_cols_by_pathway_n_elm(csv_file, &PATHWAY_CH3OH_2, &elm).unwrap();
            let p3 = get_cols_by_pathway_n_elm(csv_file, &PATHWAY_CH4_1, &elm).unwrap();
            let p4 = get_cols_by_pathway_n_elm(csv_file, &PATHWAY_CH4_2, &elm).unwrap();
            let p5 = get_cols_by_pathway_n_elm(csv_file, &PATHWAY_HCOOH, &elm).unwrap();
            let to_vec = |df: DataFrame| -> Vec<f64> {
                df.get_row(0)
                    .0
                    .into_iter()
                    .map(|v| {
                        v.to_string()
                            .replace("\\", "")
                            .parse::<f64>()
                            .unwrap_or_else(|e| panic!("parsing {} error, {}", v, e))
                    })
                    .collect()
            };
            let reaction_energy = |ads_e: Vec<f64>| -> Vec<f64> {
                ads_e
                    .chunks(6)
                    .into_iter()
                    .flat_map(|chunk| -> Vec<f64> {
                        chunk
                            .into_iter()
                            .zip(ads_e[0..6].into_iter())
                            .map(|(val, base)| *val - *base)
                            .collect()
                    })
                    .collect()
            };
            let mut p1_vec: Vec<f64> = reaction_energy(to_vec(p1));
            let mut p2_vec: Vec<f64> = reaction_energy(to_vec(p2));
            let p3_vec = reaction_energy(to_vec(p3));
            let p4_vec = reaction_energy(to_vec(p4));
            let mut p5_vec = reaction_energy(to_vec(p5));
            let p1_2_vec_marks = vec![MARK; p3_vec.len() - p1_vec.len()];
            p1_2_vec_marks.into_iter().for_each(|m| {
                p1_vec.push(m);
                p2_vec.push(m);
            });
            let p5_vec_marks = vec![MARK; p3_vec.len() - p5_vec.len()];
            p5_vec_marks.into_iter().for_each(|m| p5_vec.push(m));
            assert_eq!(p3_vec.len(), p1_vec.len());
            let all = vec![p1_vec, p2_vec, p3_vec, p4_vec, p5_vec];
            all.into_iter().flatten().collect()
        })
        .collect();
    println!("{:?}", all_data_flatten);
    // let array: Array3<f64> = Array3::from_shape_vec((220, 9, 6), all_data_flatten).unwrap();
    // println!("{:?}", array);
}
