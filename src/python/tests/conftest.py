import pytest

from .sourmash_tst_utils import TempDirectory, RunnerContext

@pytest.fixture
def runtmp():
    with TempDirectory() as location:
        yield RunnerContext(location)
