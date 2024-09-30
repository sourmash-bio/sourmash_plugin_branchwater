# sourmash_plugin_branchwater

[![PyPI](https://img.shields.io/pypi/v/sourmash_plugin_branchwater)](https://pypi.org/project/sourmash_plugin_branchwater/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/sourmash_plugin_branchwater.svg)](https://anaconda.org/conda-forge/sourmash_plugin_branchwater)

tl;dr Do faster and lower-memory sourmash functions via this plugin.

## Details

[sourmash](https://sourmash.readthedocs.io/en/latest/) is a
command-line tool and Python/Rust library for metagenome analysis and
genome comparison using k-mers.  While sourmash is fast and low
memory, sourmash v4 and lower work in single-threaded mode with Python
containers.

The branchwater plugin for sourmash (this plugin!) provides faster and
lower-memory implementations of several important sourmash features -
sketching, searching, and gather (metagenome decomposition). It does
so by implementing higher-level functions in Rust on top of the core
Rust library of sourmash.  As a result it provides some of the same
functionality as sourmash, but 10-100x faster and in 10x lower memory.

This code is still in prototype mode, and does not have all of the
features of sourmash. As we add features we will move it back into the
core sourmash code base; eventually, much of the code in this
repository will be integrated into sourmash directly.

If you're intrigued but not sure where to start with this plugin, we
suggest first identifying what sourmash functionality you need to run
to accomplish your goals. Once you have your sourmash commands
working, revisit these docs and see if there is a faster implementation
available in this plugin!

This repo originated as a [PyO3-based](https://github.com/PyO3/pyo3)
Python wrapper around the
[core branchwater code](https://github.com/sourmash-bio/sra_search).
[Branchwater](https://www.biorxiv.org/content/10.1101/2022.11.02.514947v1)
is a fast, low-memory and multithreaded application for searching very
large collections of
[FracMinHash sketches](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2)
as generated by [sourmash](https://sourmash.readthedocs.io/).

For technical details, see the Rust code in `src/` and Python wrapper
in `src/python/`.

## Documentation

There is a quickstart below, as well as
[more user documentation here](doc/README.md). Nascent
[developer docs](doc/developer.md) are also available!

The
[betterplot plugin](https://github.com/sourmash-bio/sourmash_plugin_betterplot/)
supplies a number of commands that work with branchwater output. In
particular,

* `mds2` and `tsne2` generate MDS and tSNE plots from `pairwise`
  output;
* `clustermap1` generates [seaborn clustermaps](https://seaborn.pydata.org/generated/seaborn.clustermap.html) from `pairwise` and `multisearch` output;
* `clusters_to_categories` uses the output of the `cluster` command
  to generate categories for coloring and labeling plots;
* `pairwise_to_matrix` converts the output of `pairwise` to a sourmash
  comparison matrix;

See
[the betterplot README](https://github.com/sourmash-bio/sourmash_plugin_betterplot/)
for example figures and commands!

## Quickstart demonstrating `multisearch`.

This quickstart demonstrates `multisearch` using
[the 64 genomes from Awad et al., 2017](https://osf.io/vk4fa/).

### 1. Install the branchwater plugin

On Linux and Mac OS X, you can install the latest release of the
branchwater plugin from conda-forge:
```
conda install sourmash_plugin_branchwater
```
Please [see the developer docs](doc/developer.md) for information on
installing the latest development version.

### 2. Download sketches.

The following commands will download sourmash sketches for the podar genomes into the file `podar-ref.zip`:

```
curl -L https://osf.io/4t6cq/download -o podar-ref.zip
```

### 3. Execute!

Now run `multisearch` to search all the sketches against each other:
```
sourmash scripts multisearch podar-ref.zip podar-ref.zip -o results.csv --cores 4
```

You will (hopefully ;)) see a set of results in `results.csv`. These are comparisons of each query against all matching genomes.

## Debugging help

If your collections aren't loading properly, try running `sourmash sig
summarize` on them, like so:

```
sourmash sig summarize podar-ref.zip
```

If this doesn't work, then you're running into problems creating the
collection. Please [ask for help on the sourmash issue tracker!](https://github.com/dib-lab/sourmash/issues)

## Code of Conduct

This project is under the [sourmash Code of Conduct](https://github.com/sourmash-bio/sourmash/blob/latest/CODE_OF_CONDUCT.rst).

## License

This software is under the AGPL license. Please see [LICENSE.txt](LICENSE.txt).

## Authors

* Luiz Irber
* C. Titus Brown
* Mohamed Abuelanin
* N. Tessa Pierce-Ward
* Olga Botvinnik