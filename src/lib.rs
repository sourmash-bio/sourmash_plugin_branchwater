/// Python interface Rust code for sourmash_plugin_branchwater.
use pyo3::prelude::*;

#[macro_use]
extern crate simple_error;

mod utils;
use crate::utils::build_template;
use crate::utils::is_revindex_database;
mod check;
mod fastgather;
mod fastmultigather;
mod index;
mod manysearch;
mod manysketch;
mod mastiff_manygather;
mod mastiff_manysearch;
mod multisearch;

#[pyfunction]
fn do_manysearch(
    querylist_path: String,
    siglist_path: String,
    threshold: f64,
    ksize: u8,
    scaled: usize,
    moltype: String,
    output_path: Option<String>,
) -> anyhow::Result<u8> {
    // if siglist_path is revindex, run mastiff_manysearch; otherwise run manysearch
    let template = build_template(ksize, scaled, &moltype);
    if is_revindex_database(siglist_path.as_ref()) {
        match mastiff_manysearch::mastiff_manysearch(
            querylist_path,
            siglist_path,
            template,
            threshold,
            output_path,
        ) {
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
    } else {
        match manysearch::manysearch(
            querylist_path,
            siglist_path,
            template,
            threshold,
            output_path,
        ) {
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
    }
}

#[pyfunction]
fn do_fastgather(
    query_filename: String,
    siglist_path: String,
    threshold_bp: usize,
    ksize: u8,
    scaled: usize,
    moltype: String,
    output_path_prefetch: Option<String>,
    output_path_gather: Option<String>,
) -> anyhow::Result<u8> {
    let template = build_template(ksize, scaled, &moltype);
    match fastgather::fastgather(
        query_filename,
        siglist_path,
        threshold_bp,
        ksize,
        scaled,
        template,
        output_path_prefetch,
        output_path_gather,
    ) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
fn do_fastmultigather(
    query_filenames: String,
    siglist_path: String,
    threshold_bp: usize,
    ksize: u8,
    scaled: usize,
    moltype: String,
    output_path: Option<String>,
) -> anyhow::Result<u8> {
    // if a siglist path is a revindex, run mastiff_manygather. If not, run multigather
    let template = build_template(ksize, scaled, &moltype);
    if is_revindex_database(siglist_path.as_ref()) {
        match mastiff_manygather::mastiff_manygather(
            query_filenames,
            siglist_path,
            template,
            threshold_bp,
            output_path,
        ) {
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
    } else {
        match fastmultigather::fastmultigather(
            query_filenames,
            siglist_path,
            threshold_bp,
            scaled,
            template,
        ) {
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
    }
}

#[pyfunction]
fn set_global_thread_pool(num_threads: usize) -> PyResult<usize> {
    if std::panic::catch_unwind(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
    })
    .is_ok()
    {
        Ok(rayon::current_num_threads())
    } else {
        Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
            "Could not set the number of threads. Global thread pool might already be initialized.",
        ))
    }
}

#[pyfunction]
fn do_index(
    siglist: String,
    ksize: u8,
    scaled: usize,
    moltype: String,
    output: String,
    save_paths: bool,
    colors: bool,
) -> anyhow::Result<u8> {
    // build template from ksize, scaled
    let template = build_template(ksize, scaled, &moltype);
    match index::index(siglist, template, output, save_paths, colors) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
fn do_check(index: String, quick: bool) -> anyhow::Result<u8> {
    match check::check(index, quick) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
fn do_multisearch(
    querylist_path: String,
    siglist_path: String,
    threshold: f64,
    ksize: u8,
    scaled: usize,
    moltype: String,
    output_path: Option<String>,
) -> anyhow::Result<u8> {
    let template = build_template(ksize, scaled, &moltype);
    match multisearch::multisearch(
        querylist_path,
        siglist_path,
        threshold,
        template,
        output_path,
    ) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
fn do_manysketch(filelist: String, param_str: String, output: String) -> anyhow::Result<u8> {
    match manysketch::manysketch(filelist, param_str, output) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pymodule]
fn sourmash_plugin_branchwater(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(do_manysearch, m)?)?;
    m.add_function(wrap_pyfunction!(do_fastgather, m)?)?;
    m.add_function(wrap_pyfunction!(do_fastmultigather, m)?)?;
    m.add_function(wrap_pyfunction!(do_index, m)?)?;
    m.add_function(wrap_pyfunction!(do_check, m)?)?;
    m.add_function(wrap_pyfunction!(do_manysketch, m)?)?;
    m.add_function(wrap_pyfunction!(set_global_thread_pool, m)?)?;
    m.add_function(wrap_pyfunction!(do_multisearch, m)?)?;
    Ok(())
}
