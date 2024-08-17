import sourmash_plugin_branchwater as branch
from . import sourmash_tst_utils as utils
from .sourmash_tst_utils import get_test_data


def test_basic():
    sigfile = get_test_data('SRR606249.sig.gz')
    res = branch.api.api_load_collection(sigfile, 31, 100_000, 'DNA')
    assert res.location == sigfile
    assert len(res) == 1


def test_fail():
    # try to load a (nonexistent) collection
    sigfile = get_test_data('XXX_SRR606249.sig.gz')
    try:
        res = branch.api.api_load_collection(sigfile, 31, 100_000, 'DNA')
    except:
        pass
    # @CTB should do something better here ;)


def test_basic_get_manifest():
    sigfile = get_test_data('SRR606249.sig.gz')
    res = branch.api.api_load_collection(sigfile, 31, 100_000, 'DNA')
    mf = res.manifest
    print(mf, dir(mf))
    assert len(mf) == 1

    rec = res.get_first_record()
    print(rec, dir(rec))

    print(rec.get_name())

    print(res.get_records())
    for rec in res.get_records():
        print(rec.get_name())
    assert 0
