import pytest

from .sourmash_tst_utils import TempDirectory, RunnerContext


@pytest.fixture
def runtmp():
    with TempDirectory() as location:
        yield RunnerContext(location)


@pytest.fixture(params=["--internal-storage", "--no-internal-storage"])
def toggle_internal_storage(request):
    return request.param


@pytest.fixture(params=[True, False])
def zip_query(request):
    return request.param


@pytest.fixture(params=[True, False])
def zip_db(request):
    return request.param


@pytest.fixture(params=[True, False])
def zip_against(request):
    return request.param


@pytest.fixture(params=[True, False])
def indexed(request):
    return request.param


@pytest.fixture(params=[True, False])
def indexed_query(request):
    return request.param


@pytest.fixture(params=[True, False])
def indexed_against(request):
    return request.param
