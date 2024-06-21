use std::{fs::read_to_string, path::Path};

use castep_cell_io::CellDocument;
use castep_cell_io::CellParser;

pub fn load_template() -> CellDocument {
    let cwd = Path::new(env!("CARGO_MANIFEST_DIR"));
    let template_path = cwd.join("resources").join("GDY_Sc_Sc_Sc.cell");
    let read_to_string =
        read_to_string(template_path).expect("Error when loading the template file.");
    CellParser::from(&read_to_string)
        .parse()
        .expect("Error when parsing the template")
}

#[cfg(test)]
mod test {
    use castep_periodic_table::element::ElementSymbol;
    use itertools::Itertools;

    #[test]
    fn test_lists() {
        // element_lists().iter().for_each(|elm| println!("{elm}"));
        let list = [ElementSymbol::Sc, ElementSymbol::Ti, ElementSymbol::V];
        list.iter()
            .combinations_with_replacement(3)
            .for_each(|combo| println!("{:?}", combo));
        // dbg!(element_lists().len());
        // dbg!(tri_elm_combos().len());
        // dbg!(tri_elm_combos()
        // .iter()
        // .any(|[a, b, c]| { a == b || a == b || b == c }));
    }
    #[test]
    fn rare_earth_regex() {
        let expr = (57..=71)
            .map(|i| ElementSymbol::try_from(i as u8).unwrap())
            .map(|sym| format!("{sym}"))
            .collect::<Vec<String>>()
            .join("|");
        println!("{expr}");
    }
}
