PYTHON ?= python

all:
	maturin develop

install:
	$(PYTHON) -m pip install -e .

clean:
	$(PYTHON) -m pip uninstall .

test:
	$(PYTHON) -m pytest
