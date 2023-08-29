import os
import pytest
import pandas
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


def test_index_bad_siglist(runtmp, capfd):
    # test index with a bad siglist (.sig.gz file instead of pathlist)

    sig2 = get_test_data('2.fa.sig.gz')
    output = runtmp.output('out.db')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'index', sig2,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: invalid line in fromfile" in captured.err
    print(runtmp.last_result.err)


def test_index_missing_siglist(runtmp, capfd):
    # test with a bad against list (a missing file)
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')
    make_file_list(against_list, [sig2, "no-exist"])

    db = runtmp.output('db.rocksdb')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'index', against_list,
                        '-o', db)

    captured = capfd.readouterr()
    print(captured.err)
    assert 'Error processing "no-exist"' in captured.err
    # assert "WARNING: could not load sketches from path 'no-exist'" in captured.err
    # assert "WARNING: 1 signature paths failed to load. See error messages above." in captured.err


def test_index_check(runtmp):
    # test check index
    siglist = runtmp.output('db-sigs.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')

    make_file_list(siglist, [sig2, sig47])

    output = runtmp.output('db.rocksdb')

    runtmp.sourmash('scripts', 'index', siglist,
                    '-o', output)

    runtmp.sourmash('scripts', 'check', output)
    print(runtmp.last_result.err)

    assert 'index is ok' in runtmp.last_result.err


def test_search_simple(runtmp):
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
        if row['match_name'] == row['query_name']:
            assert float(row['containment'] == 1.0)
        else:
            # confirm hand-checked numbers
            q = row['query_name'].split()[0]
            m = row['match_name'].split()[0]
            cont = float(row['containment'])
            intersect_hashes = int(row['intersect_hashes'])

            cont = round(cont, 4)
            print(q, m, f"{cont:.04}")

            if q == 'NC_011665.1' and m == 'NC_009661.1':
                assert cont == 0.4828
                assert intersect_hashes == 2529

            if q == 'NC_009661.1' and m == 'NC_011665.1':
                assert cont == 0.4885
                assert intersect_hashes == 2529


def test_search_containment_threshold(runtmp):
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
                    '-o', output, '--containment-threshold', '0.5')
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 3


# todo: write better error handling for revindex opening errors
def test_search_bad_db_format(runtmp, capfd):
    # index improperly formatted
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    # db = runtmp.output('db.rocksdb')
    output = runtmp.output('out.csv')

    # # build index
    # runtmp.sourmash('scripts', 'index', against_list,
    #                 '-o', db)

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'search', query_list, against_list,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    assert 'Not a directory' in captured.err # this is the error from the index reading...


def test_search_missing_query(runtmp, capfd):
    # test with a missing query list
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    #make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')
    db = runtmp.output('db.rocksdb')

    # build index
    runtmp.sourmash('scripts', 'index', against_list, '-o', db)

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'search', query_list, db,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: No such file or directory ' in captured.err


def test_search_bad_query(runtmp, capfd):
    # test with a bad query (a .sig.gz file instead of a pathlist)
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')
    db = runtmp.output('db.rocksdb')

    # build index
    runtmp.sourmash('scripts', 'index', against_list, '-o', db)

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'search', sig2, db,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: invalid line in fromfile ' in captured.err


# todo: write code in mastiff_manysearch for handling query errors
def test_search_bad_query_2(runtmp, capfd):
    return

    # test with a bad query list (a missing file)
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')
    make_file_list(query_list, [sig2, "no-exist"])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')
    db = runtmp.output('db.rocksdb')

    # build index
    runtmp.sourmash('scripts', 'index', against_list, '-o', db)

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'search', query_list, db,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)

    assert "No such file or directory" in captured.err
    # assert "WARNING: could not load sketches from path 'no-exist'" in captured.err
    # assert "WARNING: 1 signature paths failed to load. See error messages above." in captured.err


# todo: write code in mastiff_manysearch for handling query errors
def test_search_empty_query(runtmp):
    # test with an empty query list
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')
    db = runtmp.output('db.rocksdb')

    # build index
    runtmp.sourmash('scripts', 'index', against_list, '-o', db)

    #with pytest.raises(utils.SourmashCommandFailed):
    runtmp.sourmash('scripts', 'search', query_list, db,
                        '-o', output)

    print(runtmp.last_result.err)
    assert "search is done!" in runtmp.last_result.err


def test_search_nomatch_query(runtmp, capfd):
    # test a non-matching (diff ksize) in query; do we get warning message?
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig1 = get_test_data('1.fa.k21.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63, sig1])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')
    db = runtmp.output('db.rocksdb')

    # build index
    runtmp.sourmash('scripts', 'index', against_list, '-o', db)

    runtmp.sourmash('scripts', 'search', query_list, db,
                    '-o', output)
    assert os.path.exists(output)

    captured = capfd.readouterr()
    print(captured.err)

    assert 'WARNING: skipped 1 paths - no compatible signatures.' in captured.err


def test_search_load_only_one_bug(runtmp, capfd):
    # check that we behave properly when presented with multiple against
    # sketches
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig1_k31 = get_test_data('1.fa.k31.sig.gz')

    # note: this was created as a 3-sketch-in-one-signature directly
    # via sourmash sketch dna -p k=21,k=31,k=51.
    sig1_all = get_test_data('1.combined.sig.gz')

    make_file_list(query_list, [sig1_k31])
    make_file_list(against_list, [sig1_all])

    output = runtmp.output('out.csv')
    db = runtmp.output('db.rocksdb')

    # build index
    runtmp.sourmash('scripts', 'index', against_list, '-o', db)

    runtmp.sourmash('scripts', 'search', query_list, db,
                    '-o', output)
    assert os.path.exists(output)

    captured = capfd.readouterr()
    print(captured.err)

    assert not 'WARNING: skipped 1 paths - no compatible signatures.' in captured.err
    assert not 'WARNING: no compatible sketches in path ' in captured.err


def test_search_load_only_one_bug_as_query(runtmp, capfd):
    # check that we behave properly when presented with multiple query
    # sketches in one file, with only one matching.
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig1_k31 = get_test_data('1.fa.k31.sig.gz')

    # note: this was created as a 3-sketch-in-one-signature directly
    # via sourmash sketch dna -p k=21,k=31,k=51.
    sig1_all = get_test_data('1.combined.sig.gz')

    make_file_list(query_list, [sig1_all])
    make_file_list(against_list, [sig1_k31])

    output = runtmp.output('out.csv')
    db = runtmp.output('db.rocksdb')

    # build index
    runtmp.sourmash('scripts', 'index', against_list, '-o', db)

    runtmp.sourmash('scripts', 'search', query_list, db,
                    '-o', output)
    assert os.path.exists(output)

    captured = capfd.readouterr()
    print(captured.err)
    print(runtmp.last_result.err)

    assert not 'WARNING: skipped 1 paths - no compatible signatures.' in captured.err
    assert not 'WARNING: no compatible sketches in path ' in captured.err


# todo: get match md5s into output so we can check those too
def test_search_md5(runtmp):
    # test that md5s match what was in the original files, not downsampled etc.
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')
    db = runtmp.output('db.rocksdb')

    # build index
    runtmp.sourmash('scripts', 'index', against_list, '-o', db)

    runtmp.sourmash('scripts', 'search', query_list, db,
                    '-o', output)
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 5

    # md5s = list(df['match_md5'])
    # print(md5s)

    # for against_file in (sig2, sig47, sig63):
    #     for ss in sourmash.load_file_as_signatures(against_file, ksize=31):
    #         assert ss.md5sum() in md5s

    md5s = list(df['query_md5'])
    print(md5s)

    for query_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(query_file, ksize=31):
            assert ss.md5sum() in md5s
