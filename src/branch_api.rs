/// Lower-level Python API implementation for sourmash_plugin_branchwater
use pyo3::prelude::*;
use sourmash::prelude::*;

use crate::fastmultigather::fastmultigather_obj;
use crate::fastmultigather_rocksdb::fastmultigather_rocksdb_obj2;

use crate::utils::build_selection;
use crate::utils::load_collection;
use crate::utils::multicollection::MultiCollection;
use crate::utils::ReportType;
use pyo3::types::{IntoPyDict, PyDict};
use pyo3::IntoPyObjectExt;
use sourmash::index::revindex::{RevIndex, RevIndexOps};
use sourmash::manifest::{Manifest, Record};

#[pyclass]
pub struct BranchSelection {
    pub selection: Selection,
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

    pub fn ksize(&self) -> PyResult<u32> {
        Ok(self.selection.ksize().expect("ksize not set"))
    }

    pub fn moltype(&self) -> PyResult<String> {
        Ok(self
            .selection
            .moltype()
            .expect("moltype not set")
            .to_string())
    }

    pub fn scaled(&self) -> PyResult<u32> {
        Ok(self.selection.scaled().expect("scaled not set"))
    }
}

#[pyclass]
pub struct BranchRevIndex {
    location: String,
    db: RevIndex,
}

#[pymethods]
impl BranchRevIndex {
    #[new]
    pub fn open(location: &str) -> PyResult<Py<Self>> {
        let db = RevIndex::open(location, true, None).expect("foo");

        let obj = Python::with_gil(|py| Py::new(py, BranchRevIndex {
            location: location.to_string(),
            db
        }).unwrap());
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

    pub fn ksize(&self) -> PyResult<u8> {
        let ksize = self
            .db
            .collection()
            .manifest()
            .first()
            .map(|first| first.ksize())
            .expect("no records in db");
        Ok(ksize.try_into().unwrap())
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

    pub fn selection(&self) -> PyResult<BranchSelection> {
        let ksize = self.ksize().expect("foo 1");
        let moltype = self.moltype().expect("foo 2");
        let (_, max_scaled) = self.min_max_scaled().expect("foo 3");
        let selection = BranchSelection::build(ksize, max_scaled, moltype.as_ref());

        Ok(selection)
    }

    pub fn to_collection(&self) -> Py<BranchCollection> {
        let cs = self.db.collection().clone();
        let mc = MultiCollection::new(vec![cs.into_inner()], true);

        let obj = Python::with_gil(|py| {
            Py::new(
                py,
                BranchCollection {
                    location: self.location.clone(),
                    collection: mc,
                    is_database: true,
                    has_manifest: true,
                },
            )
            .expect("cannot convert collection into py object")
        });
        obj
    }

    #[pyo3(signature = (query_collection, selection, threshold_bp, output))]
    pub fn fastmultigather_against(
        &self,
        query_collection: &BranchCollection,
        selection: &BranchSelection,
        threshold_bp: u32,
        output: String,
    ) -> anyhow::Result<u8> {
        match fastmultigather_rocksdb_obj2(
            &query_collection.collection,
            &self.db,
            &selection.selection,
            threshold_bp,
            Some(output),
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

    pub collection: MultiCollection,
}

/// Wraps MultiCollection

#[pymethods]
impl BranchCollection {
    pub fn __len__(&self) -> PyResult<usize> {
        Ok(self.collection.len())
    }

    pub fn select(&self, selection: &BranchSelection) -> PyResult<Self> {
        let collection = self
            .collection
            .clone()
            .select(&selection.selection)
            .expect("selection failed");

        let obj = BranchCollection {
            location: self.location.clone(),
            collection,
            is_database: false,
            has_manifest: true,
        };
        Ok(obj)
    }

    #[pyo3(signature = (query_collection, threshold_bp, scaled, output))]
    pub fn fastmultigather_against(
        &self,
        query_collection: &BranchCollection,
        threshold_bp: u32,
        scaled: u32,
        output: String,
    ) -> anyhow::Result<u8> {
        let threshold_hashes: u64 = (threshold_bp / scaled).into();
        eprintln!("foo: {} {}", threshold_hashes, scaled);

        match fastmultigather_obj(
            &query_collection.collection,
            &self.collection,
            false,
            Some(output),
            threshold_hashes,
            scaled,
        ) {
            // @CTB
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
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

// @CTB move into struct?
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
