/// Lower-level Python API implementation for sourmash_plugin_branchwater
use pyo3::prelude::*;

use sourmash::collection::Collection;
use crate::utils::build_selection;
use crate::utils::load_collection;
use crate::utils::ReportType;

#[pyclass]
pub struct BranchCollection {
    #[pyo3(get)]
    pub val: i32,
    collection: Collection
}

#[pymethods]
impl BranchCollection {
    pub fn __len__(&self) -> PyResult<usize> {
        Ok(self.collection.len())
    }
}


#[pyfunction]
pub fn api_load_collection(
    location: String,
    ksize: u8,
    scaled: usize,
    moltype: String,
) -> PyResult<Py<BranchCollection>> {
    let selection = build_selection(ksize, scaled, &moltype);

    let collection = load_collection(&location, &selection, ReportType::Query, true).unwrap();
    let obj =
        Python::with_gil(|py| Py::new(py, BranchCollection { val: 1001, collection }).unwrap());
    Ok(obj)
}
