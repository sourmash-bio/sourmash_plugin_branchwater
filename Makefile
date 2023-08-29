PYTHON ?= python

all:
	maturin develop

install:
	$(PYTHON) -m pip install -e .

clean:
	$(PYTHON) -m pip uninstall .

test:
	$(PYTHON) -m pytest

wheel:
	$(PYTHON) -m maturin build -r

sdist:
	rm -f target/wheels/pyo3_branchwater-*.tar.gz
	$(PYTHON) -m maturin sdist

upload_sdist: sdist
	twine upload target/wheels/pyo3_branchwater-*.tar.gz
