pub mod element_table;

#[test]
fn test_yaml() -> Result<(), serde_yaml::Error> {
    use crate::external_info::element_table::hash_table;
    let table = self::element_table::load_table().unwrap();
    println!("{}", table.elements.is_some());
    table
        .elements
        .unwrap()
        .iter()
        .for_each(|elm| println!("{:#?}", elm));
    let hashtab = hash_table();
    println!("{:#?}", hashtab.get("C").unwrap());
    Ok(())
}
