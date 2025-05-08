.PHONY: all install format test wheel sdist upload_dist

PYTHON ?= python

all:
	maturin develop

install:
	$(PYTHON) -m pip install -e .

format:
	cargo fmt
	python -m black src

test:
	$(PYTHON) -m pytest

wheel:
	$(PYTHON) -m maturin build -r

sdist:
	rm -f target/wheels/sourmash_plugin_branchwater-*.tar.gz
	$(PYTHON) -m maturin sdist

upload_sdist: sdist
	twine upload target/wheels/sourmash_plugin_branchwater-*.tar.gz
