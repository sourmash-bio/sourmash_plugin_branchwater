# 2022-pymagsearch

Python wrapper around
[core sra_search code](https://github.com/sourmash-bio/sra_search).

See Rust code in `src/` and Python wrapper in `python/`.

To try out, use a branch of sourmash that contains [sourmash#2438](https://github.com/sourmash-bio/sourmash/pull/2438). Then, `pip install -e .` and
use `sourmash scripts manysearch` or `sourmash scripts manygather`.

Uses [pyo3](https://github.com/PyO3/pyo3) for wrapping.

CTB 02/2023
