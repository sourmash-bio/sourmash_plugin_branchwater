import os
import pytest
import sourmash

from . import sourmash_tst_utils as utils

def get_test_data(filename):
    thisdir = os.path.dirname(__file__)
    return os.path.join(thisdir, 'test-data', filename)

def make_file_list(filename, paths):
    with open(filename, 'wt') as fp:
        fp.write("\n".join(paths))
        fp.write("\n")

def test_installed(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'sigcat')

    assert 'usage:  sigcat' in runtmp.last_result.err

def zip_siglist(runtmp, siglist, db):
    runtmp.sourmash('sig', 'cat', siglist,
                    '-o', db)
    return db


def test_simple_sigcat(runtmp):
    # test basic execution!
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2])
    make_file_list(against_list, [sig47, sig63])

    output = runtmp.output('out.zip')

    runtmp.sourmash('scripts', 'sigcat', '--signatures', query_list, against_list,
                    '-o', output)
    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    
    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 3

def test_simple_sigcat_with_zip(runtmp):
    # test basic execution with a zip involved
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2])
    make_file_list(against_list, [sig47, sig63])

    output = runtmp.output('out.zip')
    # zip the against_list
    against_list = zip_siglist(runtmp, against_list, runtmp.output('db.zip'))

    runtmp.sourmash('scripts', 'sigcat', '--signatures', query_list, against_list,
                    '-o', output)
    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    
    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 3

