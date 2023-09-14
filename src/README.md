# Source code guide

The pyo3 Rust/Python interface code is in `lib.rs`. The top level Rust functions called by the pyo3 Rust code are in individual files named for the function, and common utility code is in `utils.rs`.

The Python source code is under `python/`, and tests under `python/tests/`.
