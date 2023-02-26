import os
import pytest
import pandas

import sourmash_tst_utils as utils


def get_test_data(filename):
    thisdir = os.path.dirname(__file__)
    return os.path.join(thisdir, 'test-data', filename)


def make_file_list(filename, paths):
    with open(filename, 'wt') as fp:
        fp.write("\n".join(paths))
        fp.write("\n")


def test_installed(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysearch')

    assert 'usage:  manysearch' in runtmp.last_result.err


def test_simple(runtmp):
    # test basic execution!
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    runtmp.sourmash('scripts', 'manysearch', query_list, against_list,
                    '-o', output)
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 5


def test_simple_threshold(runtmp):
    # test with a simple threshold => only 3 results
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    runtmp.sourmash('scripts', 'manysearch', query_list, against_list,
                    '-o', output, '-t', '0.5')
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 3


def test_bad_query(runtmp):
    # test with a bad query list
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    #make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysearch', query_list, against_list,
                        '-o', output)

    print(runtmp.last_result.err)


def test_bad_against(runtmp):
    # test with a bad (nonexistent) against list
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    #make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysearch', query_list, against_list,
                        '-o', output)

    err = runtmp.last_result.err
    print('XXX', err)

    #assert "Reading search file paths from" in err
    #assert "Error: No such file or directory (os error 2)" in err


def test_empty_query(runtmp):
    # test with an empty query list
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysearch', query_list, against_list,
                        '-o', output)

    print(runtmp.last_result.err)


def test_empty_against(runtmp):
    # test with an empty (nonexistent) against list
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [])

    output = runtmp.output('out.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysearch', query_list, against_list,
                        '-o', output)

    err = runtmp.last_result.err
    print('XXX', err)

    #assert "Reading search file paths from" in err
    #assert "Error: No such file or directory (os error 2)" in err
