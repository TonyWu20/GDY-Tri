use std::io::Error;
use std::path::Path;

pub trait MyOutput {
    fn formatted_output(&self) -> String;
}

impl MyOutput for f64 {
    fn formatted_output(&self) -> String {
        format!("{:.8}", *self)
    }
}

impl MyOutput for &str {
    fn formatted_output(&self) -> String {
        self.to_string()
    }
}

pub trait ExportGnuData {
    fn to_gnu_data<P: AsRef<Path>>(&self, filename: P) -> Result<(), Error>;
}
