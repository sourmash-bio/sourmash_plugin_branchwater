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
        runtmp.sourmash('scripts', 'search')

    assert 'usage:  search' in runtmp.last_result.err

def test_index(runtmp):
    # test basic index!
    siglist = runtmp.output('db-sigs.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(siglist, [sig2, sig47, sig63])

    output = runtmp.output('db.rocksdb')

    runtmp.sourmash('scripts', 'index', siglist,
                    '-o', output)
    assert os.path.exists(output)
    print(runtmp.last_result.err)

    assert 'index is done' in runtmp.last_result.err

# def test_search(runtmp):

def test_simple_search(runtmp):
    # test basic execution!
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    db = runtmp.output('db.rocksdb')
    output = runtmp.output('out.csv')

    # build index
    runtmp.sourmash('scripts', 'index', against_list,
                    '-o', db)
    # search the index
    runtmp.sourmash('scripts', 'search', query_list, db,
                    '-o', output)
    
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 5

    dd = df.to_dict(orient='index')
    print(dd)

    for idx, row in dd.items():
        # identical?
        # if row['query_md5'] == row['match_md5']:
        if row['match_name'] == row['query_name']:
            assert float(row['containment'] == 1.0)
            # assert float(row['jaccard'] == 1.0)
        else:
            # confirm hand-checked numbers
            q = row['query_name'].split()[0]
            m = row['match_name'].split()[0]
            # jaccard = float(row['jaccard'])
            cont = float(row['containment'])
            # intersect_hashes = int(row['intersect_hashes'])

            # jaccard = round(jaccard, 4)
            cont = round(cont, 4)
            print(q, m, f"{cont:.04}")
            # print(q, m, f"{jaccard:.04}", f"{cont:.04}")

            if q == 'NC_011665.1' and m == 'NC_009661.1':
                # assert jaccard == 0.3207
                assert cont == 0.4828

            if q == 'NC_009661.1' and m == 'NC_011665.1':
                # assert jaccard == 0.3207
                assert cont == 0.4885

def test_simple_containment_threshold(runtmp):
    # test with a simple threshold => only 3 results
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    db = runtmp.output('db.rocksdb')
    output = runtmp.output('out.csv')

    # build index
    runtmp.sourmash('scripts', 'index', against_list,
                    '-o', db)
    runtmp.sourmash('scripts', 'search', query_list, db,
                    '-o', output, '-c', '0.5')
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 3
