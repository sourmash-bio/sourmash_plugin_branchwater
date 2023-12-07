# Developer notes

## Developer installation of sourmash_plugin_branchwater

### 1. Install necessary dependencies

You'll need Rust, Python, and maturin to build the branchwater plugin,
and sourmash to run to run it.

There is a list of the necessary conda packages in [`environment.yml`](../environment.yml).  To install them in a new conda environment, run:

```
mamba env create -n branchwater -f environment.yml
```

and then activate the conda environment:
```
conda activate branchwater
```

### 2. Install sourmash_plugin_branchwater.

Install this repo in editable mode:
```
pip install -e .
```

### 3. Fun, profit.

You can now use the plugin in "developer" mode, where any changes you make
to the Python source code will be reflected in the installed version. To
update the installed version with changes to the Rust code, you'll need to
compile the Rust code with:
```
maturin develop
```

## Running the tests locally

Executing:
```
make test
```
will run the Python tests.

## Generating a release

1. Bump version number in `Cargo.toml` and run `make` to update `Cargo.lock`.
Then commit and push to `origin/main`.

2. Make a new release on github with a matching version tag.

3. Then pull, and:

```
make sdist
make upload_sdist
```

to create a new release on PyPI.

## Building wheels

You can build a release wheel for your current platform with:
```
make wheel
```
and it will be placed under `target/wheels/`.

