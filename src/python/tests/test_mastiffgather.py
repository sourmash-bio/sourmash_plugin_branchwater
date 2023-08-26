"""
Test 'sourmash scripts gather'
"""
import os
import pytest
import pandas

import sourmash
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
        runtmp.sourmash('scripts', 'gather')

    assert 'usage:  gather' in runtmp.last_result.err

def test_gather_simple(runtmp):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')
    
    db = runtmp.output('db')
    output = runtmp.output('out.csv')

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    # build index
    runtmp.sourmash('scripts', 'index', against_list,
                    '-o', db)
    
    # run gather
    runtmp.sourmash('scripts', 'gather', query_list, db,
                    '-o', output, '-t', '0', '--scaled', '100000')

    assert os.path.exists(output)
    
    df = pandas.read_csv(output)
    assert len(df) == 3

    dd = df.to_dict(orient='index')
    print(dd)

    for idx, row in dd.items():
        assert row.keys() == {'query_name', 'query_md5', 'match_name', 'match_md5', 'f_match', 'intersect_bp'}
        assert row['query_name'] == 'SRR606249'
        assert row['query_md5'] == 'dec29ca72e68db0f15de0b1b46f82fc5'