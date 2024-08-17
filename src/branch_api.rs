/// Lower-level Python API implementation for sourmash_plugin_branchwater
use pyo3::prelude::*;

use crate::utils::build_selection;
use crate::utils::load_collection;
use crate::utils::ReportType;

#[pyclass]
pub struct BranchCollection {
    #[pyo3(get)]
    pub val: i32,
}


#[pyfunction]
pub fn api_load_collection(
    location: String,
    ksize: u8,
    scaled: usize,
    moltype: String,
) -> PyResult<Py<BranchCollection>> {
    let selection = build_selection(ksize, scaled, &moltype);

/*
    match load_collection(&location, &selection, ReportType::Query, true) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
}
     */
    let obj =
        Python::with_gil(|py| Py::new(py, BranchCollection { val: 1001 }).unwrap());
    Ok(obj)
}


