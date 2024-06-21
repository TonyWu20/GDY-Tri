use castep_cell_io::CellDocument;
use castep_periodic_table::element::ElementSymbol;
use itertools::Itertools;
use rayon::prelude::*;

use crate::template::load_template;

#[derive(Debug, Clone)]
pub struct Model {
    basic_cell: CellDocument,
    new_elements: [ElementSymbol; 3],
}

impl Model {
    pub fn new(basic_cell: CellDocument, new_elements: [ElementSymbol; 3]) -> Self {
        Self {
            basic_cell,
            new_elements,
        }
    }
    pub fn output_name(&self) -> String {
        format!(
            "GDY_{}",
            self.new_elements
                .iter()
                .map(|symbol| format!("{symbol}"))
                .collect::<Vec<String>>()
                .join("_")
        )
    }

    pub fn dir_name(&self) -> String {
        format!(
            "GDY_{}",
            self.new_elements
                .iter()
                .map(|symbol| {
                    match symbol.family() {
                        castep_periodic_table::element::ElementFamily::TransitionMetal3d => {
                            "3d".to_string()
                        }
                        castep_periodic_table::element::ElementFamily::TransitionMetal4d => {
                            "4d".to_string()
                        }
                        castep_periodic_table::element::ElementFamily::TransitionMetal5d => {
                            "5d".to_string()
                        }
                        castep_periodic_table::element::ElementFamily::RareEarthLa => {
                            "RE".to_string()
                        }
                        _ => "else".to_string(),
                    }
                })
                .collect::<Vec<String>>()
                .join("_")
        )
    }

    pub fn basic_cell(&self) -> &CellDocument {
        &self.basic_cell
    }
}

pub fn generate_models() -> Vec<Model> {
    let template = load_template();
    tri_elm_combos()
        .par_iter()
        .map(|combo| {
            let new_doc = edit_template(&template, combo);
            Model::new(new_doc, *combo)
        })
        .collect()
}

fn element_lists() -> Vec<ElementSymbol> {
    [(21..=30), (39..=48), (72..=80), (57..=71)]
        .into_iter()
        .flat_map(|rg| {
            rg.map(|i| ElementSymbol::try_from(i as u8).unwrap())
                .collect::<Vec<ElementSymbol>>()
        })
        .collect()
}

fn tri_elm_combos() -> Vec<[ElementSymbol; 3]> {
    let list = element_lists();
    // let mut used_combos: HashSet<[ElementSymbol; 3]> = HashSet::new();
    list.into_iter()
        .combinations_with_replacement(3)
        .map(|combo| combo.try_into().unwrap())
        .collect()
}

fn edit_template(template: &CellDocument, combo: &[ElementSymbol; 3]) -> CellDocument {
    let mut new_template = template.clone();
    let positions = new_template
        .model_description_mut()
        .ionic_pos_block_mut()
        .positions_mut();
    combo
        .iter()
        .zip(72..75)
        .for_each(|(&symbol, i)| positions[i].set_symbol(symbol));
    new_template
}

#[cfg(test)]
mod test {
    use castep_periodic_table::element::ElementSymbol;

    use crate::template::load_template;

    use super::edit_template;

    #[test]
    fn template_edit() {
        let cell_template = load_template();
        let new = edit_template(
            &cell_template,
            &[ElementSymbol::Sc, ElementSymbol::Ti, ElementSymbol::Au],
        );
        new.model_description()
            .ionic_pos_block()
            .positions()
            .last_chunk::<3>()
            .unwrap()
            .iter()
            .for_each(|pos| println!("{}", pos.symbol()));
    }
}
