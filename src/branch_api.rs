/// Lower-level Python API implementation for sourmash_plugin_branchwater
use pyo3::prelude::*;

use crate::utils::build_selection;
use crate::utils::load_collection;
use crate::utils::ReportType;
use crate::utils::multicollection::MultiCollection;
use sourmash::collection::Collection;
use sourmash::manifest::{Manifest, Record};
use sourmash::index::revindex::{RevIndex, RevIndexOps};
use pyo3::types::{IntoPyDict, PyDict, PyList};
use pyo3::IntoPyObjectExt;
use camino::Utf8PathBuf as PathBuf;

#[pyclass]
pub struct BranchRevIndex {
    db: RevIndex
}

#[pymethods]
impl BranchRevIndex {
    pub fn open(location: PathBuf) -> PyResult<Py<Self>> {
        let db = RevIndex::open(location, true, None).expect("foo");

        let obj = Python::with_gil(|py| {
            Py::new(py, BranchRevIndex { db }).unwrap()
        });
        Ok(obj)
    }
}

#[pyclass]
pub struct BranchRecord {
    record: Record,
}

#[pymethods]
impl BranchRecord {
    pub fn get_name(&self) -> PyResult<String> {
        Ok(self.record.name().clone())
    }

    #[getter]
    pub fn get_as_row<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = {
            let key_vals: Vec<(&str, PyObject)> = vec![
                ("ksize", self.record.ksize().into_bound_py_any(py)?.unbind()),
                ("moltype", self.record.moltype().to_string().into_bound_py_any(py)?.unbind()),
                ("scaled", self.record.scaled().into_bound_py_any(py)?.unbind()),
                ("num", self.record.num().into_bound_py_any(py)?.unbind()),
                ("with_abundance", self.record.with_abundance().into_bound_py_any(py)?.unbind()),
                ("n_hashes", self.record.n_hashes().into_bound_py_any(py)?.unbind()),
            ];
            key_vals.into_py_dict_bound(py)
        };
        Ok(dict)
    }
}

/*
impl<T, I> IntoPyDict for I
where
    T: PyDictItem
    I: IntoIterator<Item = T>
fn into_py_dict(self, py: Python<'_>) -> Bound<'_, PyDict> {
    let dict = PyDict::new(py);
    for item in self {
        dict.set_item(item.key(), item.value())
            .expect("Failed to set_item on dict");
    }
    dict
}
}
*/

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
    #[getter]
    pub fn get_rows<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        let res: Vec<_> = self.manifest.iter().map(|x| { BranchRecord {
            record: x.clone(),
        }.get_as_row(py).unwrap()
        }).collect();

        Ok(res)
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

    collection: MultiCollection,
}

#[pymethods]
impl BranchCollection {
    pub fn __len__(&self) -> PyResult<usize> {
        Ok(self.collection.len())
    }
/*
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

    #[getter]
    pub fn get_rows(&self) -> PyResult<Vec<BranchRecord>> {
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
*/
}

#[pyfunction]
pub fn api_load_collection(
    location: String,
    ksize: u8,
    scaled: u32,
    moltype: String,
) -> PyResult<Py<BranchCollection>> {
    let selection = build_selection(ksize, Some(scaled), &moltype);

    let collection = load_collection(&location, &selection, ReportType::Query, true)?;

    let obj = Python::with_gil(|py| {
        Py::new(
            py,
            BranchCollection {
                location: location,
                collection,
                is_database: false,
                has_manifest: true,
            },
        ).expect("cannot convert collection into py object")
    });
    Ok(obj)
}
