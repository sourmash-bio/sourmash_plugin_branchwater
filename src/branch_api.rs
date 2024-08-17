/// Lower-level Python API implementation for sourmash_plugin_branchwater
use pyo3::prelude::*;

use crate::utils::build_selection;
use crate::utils::load_collection;
use crate::utils::ReportType;
use sourmash::collection::Collection;
use sourmash::manifest::{Manifest, Record};

#[pyclass]
pub struct BranchRecord {
    record: Record,
}

#[pymethods]
impl BranchRecord {
    pub fn get_name(&self) -> PyResult<String> {
        Ok(self.record.name().clone())
    }

    pub fn to_row_dict(&self, py: Python) -> PyResult<String> {
        Ok("hello".to_string())
    }
}

#[pyclass]
pub struct BranchManifest {
    manifest: Manifest,
}

#[pymethods]
impl BranchManifest {
    pub fn __len__(&self) -> PyResult<usize> {
        Ok(self.manifest.len())
    }
    pub fn _check_row_values(&self) -> PyResult<bool> {
        Ok(true)
    }
}

#[pyclass]
pub struct BranchCollection {
    #[pyo3(get)]
    pub location: String,

    #[pyo3(get)]
    pub is_database: bool,

    #[pyo3(get)]
    pub has_manifest: bool,

    collection: Collection,
}

#[pymethods]
impl BranchCollection {
    pub fn __len__(&self) -> PyResult<usize> {
        Ok(self.collection.len())
    }

    #[getter]
    pub fn get_manifest(&self) -> PyResult<Py<BranchManifest>> {
        let manifest: Manifest = self.collection.manifest().clone();
        let obj =
            Python::with_gil(|py| Py::new(py, BranchManifest { manifest: manifest }).unwrap());
        Ok(obj)
    }
    pub fn get_first_record(&self) -> PyResult<Py<BranchRecord>> {
        let records: Vec<_> = self.collection.iter().collect();
        let first_record = records.first().unwrap().1;

        // @CTB: can I turn this into something automatic?
        let obj = Python::with_gil(|py| {
            Py::new(
                py,
                BranchRecord {
                    record: first_record.clone(),
                },
            )
            .unwrap()
        });
        Ok(obj)
    }

    pub fn get_records(&self) -> PyResult<Vec<BranchRecord>> {
        let records: Vec<_> = self.collection.iter().collect();

        let obj = records
                .iter()
                .map(|x| {
                    BranchRecord {
                        record: x.1.clone(),
                    }
                })
            .collect();

        // @CTB: this does the GIL grabbing as needed?
        Ok(obj)
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
    let obj = Python::with_gil(|py| {
        Py::new(
            py,
            BranchCollection {
                location: location,
                collection,
                is_database: false,
                has_manifest: true,
            },
        )
        .unwrap()
    });
    Ok(obj)
}
