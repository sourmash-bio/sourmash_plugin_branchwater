# 2022-pymagsearch

Python wrapper around
[core sra_search code](https://github.com/sourmash-bio/sra_search).

See Rust code in `src/` and Python wrapper in `python/`.

To try out, use a branch of sourmash that contains [sourmash#2438](https://github.com/sourmash-bio/sourmash/pull/2438). Then, `pip install -e .` and
use `sourmash scripts manysearch` or `sourmash scripts fastgather`.

Uses [pyo3](https://github.com/PyO3/pyo3) for wrapping.

## Quickstart for `manysearch`.

This quickstart demonstrates `manysearch` using
[the 64 genomes from Awad et al., 2017](https://osf.io/vk4fa/).

### First, install this code.

Install this repo in developer mode:
```
pip install -e .
```

### Second, download sketches.

The following commands will download sourmash sketches for them and
unpack them into the directory `podar-ref/`:

```
mkdir -p podar-ref
curl -JLO https://osf.io/4t6cq/download
unzip -u podar-reference-genomes-updated-sigs-2017.06.10.zip
```

### Third, create lists of query and subject files.

`manysearch` takes in lists of signatures to search, so we need to
create those files:

```
ls -1 podar-ref/{2,47,63}.* > query-list.txt
ls -1 podar-ref/* > podar-ref-list.txt
```

### Fourth: Execute!

Now run `manysearch`:
```
sourmash scripts manysearch query-list.txt podar-ref-list.txt -o results.csv
```

You will (hopefully ;) see a set of results in `results.csv`.

## Debugging help

If your file lists are not working properly, try running:
```
sourmash sig summarize query-list.txt
sourmash sig summarize podar-ref-list.txt
```


---

CTB Feb 2023
