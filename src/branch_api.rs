/// Lower-level Python API implementation for sourmash_plugin_branchwater
use pyo3::prelude::*;
use sourmash::prelude::*;

use crate::fastmultigather_rocksdb::fastmultigather_rocksdb_obj;
use crate::utils::build_selection;
use crate::utils::load_collection;
use crate::utils::multicollection::MultiCollection;
use crate::utils::ReportType;
use pyo3::types::{IntoPyDict, PyDict, PyList};
use pyo3::IntoPyObjectExt;
use sourmash::index::revindex::{RevIndex, RevIndexOps};
use sourmash::manifest::{Manifest, Record};

#[pyclass]
pub struct BranchSelection {
    selection: Selection,
}

impl BranchSelection {
    pub fn new(selection: Selection) -> Self {
        Self { selection }
    }
}

#[pymethods]
impl BranchSelection {
    #[new]
    #[pyo3(signature = (ksize, scaled, moltype))]
    pub fn build(ksize: u8, scaled: u32, moltype: &str) -> BranchSelection {
        let selection = build_selection(ksize, Some(scaled), moltype);
        Self { selection }
    }
}

#[pyclass]
pub struct BranchRevIndex {
    db: RevIndex,
}

#[pymethods]
impl BranchRevIndex {
    #[new]
    pub fn open(location: &str) -> PyResult<Py<Self>> {
        let db = RevIndex::open(location, true, None).expect("foo");

        let obj = Python::with_gil(|py| Py::new(py, BranchRevIndex { db }).unwrap());
        Ok(obj)
    }

    pub fn __len__(&self) -> usize {
        self.db.collection().len()
    }

    pub fn min_max_scaled(&self) -> PyResult<(u32, u32)> {
        let (min_scaled, max_scaled) = self
            .db
            .collection()
            .min_max_scaled()
            .expect("no records in db");
        Ok((*min_scaled, *max_scaled))
    }

    pub fn ksize(&self) -> PyResult<u32> {
        let ksize = self
            .db
            .collection()
            .manifest()
            .first()
            .map(|first| first.ksize())
            .expect("no records in db");
        Ok(ksize)
    }

    pub fn moltype(&self) -> PyResult<String> {
        let moltype = self
            .db
            .collection()
            .manifest()
            .first()
            .map(|first| first.moltype())
            .expect("no records in db");
        Ok(moltype.to_string())
    }

    #[pyo3(signature = (queries_file, ksize, scaled, moltype, threshold_bp, output))]
    pub fn fastmultigather_against(
        &self,
        queries_file: String,
        ksize: u8,
        scaled: u32,
        moltype: String,
        threshold_bp: u32,
        output: String,
    ) -> anyhow::Result<u8> {
        let selection = build_selection(ksize, Some(scaled), &moltype);

        match fastmultigather_rocksdb_obj(
            queries_file,
            &self.db,
            selection,
            threshold_bp,
            Some(output),
            true,
        ) {
            // @CTB
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
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
                (
                    "moltype",
                    self.record
                        .moltype()
                        .to_string()
                        .into_bound_py_any(py)?
                        .unbind(),
                ),
                (
                    "scaled",
                    self.record.scaled().into_bound_py_any(py)?.unbind(),
                ),
                ("num", self.record.num().into_bound_py_any(py)?.unbind()),
                (
                    "with_abundance",
                    self.record.with_abundance().into_bound_py_any(py)?.unbind(),
                ),
                (
                    "n_hashes",
                    self.record.n_hashes().into_bound_py_any(py)?.unbind(),
                ),
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
        let res: Vec<_> = self
            .manifest
            .iter()
            .map(|x| BranchRecord { record: x.clone() }.get_as_row(py).unwrap())
            .collect();

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
        )
        .expect("cannot convert collection into py object")
    });
    Ok(obj)
}
