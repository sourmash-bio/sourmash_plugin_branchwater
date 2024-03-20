"""
Test 'sourmash scripts fastmultigather'
"""
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


def index_siglist(runtmp, siglist, db, ksize=31, scaled=1000, moltype='DNA'):
    # build index
    runtmp.sourmash('scripts', 'index', siglist,
                    '-o', db, '-k', str(ksize), '--scaled', str(scaled),
                    '--moltype', moltype)
    return db


def test_installed(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastmultigather')

    assert 'usage:  fastmultigather' in runtmp.last_result.err


def zip_siglist(runtmp, siglist, db):
    runtmp.sourmash('sig', 'cat', siglist,
                    '-o', db)
    return db

@pytest.mark.parametrize('zip_against', [False, True])
def test_simple(runtmp, zip_against):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('against.zip'))

    cwd = os.getcwd()
    try:
        os.chdir(runtmp.output(''))
        runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                        '-s', '100000', '-t', '0')
    finally:
        os.chdir(cwd)

    print(os.listdir(runtmp.output('')))

    g_output = runtmp.output('SRR606249.gather.csv')
    p_output = runtmp.output('SRR606249.prefetch.csv')
    assert os.path.exists(p_output)

    # check prefetch output (only non-indexed gather)
    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'intersect_bp'}

    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    print(df)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}


def test_simple_space_in_signame(runtmp):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    renamed_query = runtmp.output('in.zip')
    name = 'my-favorite-signame has spaces'
    # rename signature
    runtmp.sourmash('sig', 'rename', query, name, '-o', renamed_query)

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    against_list = runtmp.output('against.txt')

    make_file_list(against_list, [sig2, sig47, sig63])

    cwd = os.getcwd()
    try:
        os.chdir(runtmp.output(''))
        runtmp.sourmash('scripts', 'fastmultigather', renamed_query, against_list,
                        '-s', '100000', '-t', '0')
    finally:
        os.chdir(cwd)

    print(os.listdir(runtmp.output('')))

    g_output = runtmp.output('my-favorite-signame.gather.csv')
    p_output = runtmp.output('my-favorite-signame.prefetch.csv')
    assert os.path.exists(p_output)
    assert os.path.exists(g_output)


def test_simple_zip_query(runtmp):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    query_list = zip_siglist(runtmp, query_list, runtmp.output('query.zip'))

    cwd = os.getcwd()
    try:
        os.chdir(runtmp.output(''))
        runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                        '-s', '100000', '-t', '0')
    finally:
        os.chdir(cwd)

    print(os.listdir(runtmp.output('')))

    g_output = runtmp.output('SRR606249.gather.csv')
    p_output = runtmp.output('SRR606249.prefetch.csv')

    # check prefetch output (only non-indexed gather)
    assert os.path.exists(p_output)
    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'intersect_bp'}

    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}


def test_simple_read_manifests(runtmp):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    against_list = runtmp.output('against.txt')
    against_mf = runtmp.output("against.csv")
    query_mf = runtmp.output("query.csv")

    make_file_list(against_list, [sig2, sig47, sig63])

    runtmp.sourmash("sig","manifest", query, "-o", query_mf)
    runtmp.sourmash("sig", "manifest", against_list, "-o", against_mf)

    cwd = os.getcwd()
    try:
        os.chdir(runtmp.output(''))
        runtmp.sourmash('scripts', 'fastmultigather', query_mf, against_list,
                        '-s', '100000', '-t', '0')
    finally:
        os.chdir(cwd)

    print(os.listdir(runtmp.output('')))

    g_output = runtmp.output('SRR606249.gather.csv')
    p_output = runtmp.output('SRR606249.prefetch.csv')

    # check prefetch output (only non-indexed gather)
    assert os.path.exists(p_output)
    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'intersect_bp'}

    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}


@pytest.mark.parametrize('zip_query', [False, True])
def test_simple_indexed(runtmp, zip_query):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output('query.zip'))

    g_output = runtmp.output('out.csv')
    against_db = index_siglist(runtmp, against_list, runtmp.output('test.rocksdb'))
    runtmp.sourmash('scripts', 'fastmultigather', query_list,
                        against_db, '-s', '100000', '-t', '0',
                        '-o', g_output)

    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    expected_keys = {'match_name', 'query_filename', 'query_n_hashes', 'match_filename', 'f_match_orig',
            'query_bp', 'query_abundance', 'match_containment_ani', 'intersect_bp', 'total_weighted_hashes',
            'n_unique_weighted_found', 'query_name', 'gather_result_rank', 'moltype',
            'query_containment_ani', 'sum_weighted_found', 'f_orig_query', 'ksize', 'max_containment_ani',
            'std_abund', 'scaled', 'average_containment_ani', 'f_match', 'f_unique_to_query',
            'average_abund', 'unique_intersect_bp', 'median_abund', 'query_md5', 'match_md5', 'remaining_bp',
            'f_unique_weighted'}
    assert  keys == expected_keys


def test_simple_indexed_query_manifest(runtmp):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    query_mf = runtmp.output('query.csv')
    against_list = runtmp.output('against.txt')

    make_file_list(against_list, [sig2, sig47, sig63])
    runtmp.sourmash("sig", "manifest", query, "-o", query_mf)

    g_output = runtmp.output('out.csv')
    against_db = index_siglist(runtmp, against_list, runtmp.output('db'))
    runtmp.sourmash('scripts', 'fastmultigather', query_mf,
                        against_db, '-s', '100000', '-t', '0',
                        '-o', g_output)

    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    expected_keys = {'match_name', 'query_filename', 'query_n_hashes', 'match_filename', 'f_match_orig',
            'query_bp', 'query_abundance', 'match_containment_ani', 'intersect_bp', 'total_weighted_hashes',
            'n_unique_weighted_found', 'query_name', 'gather_result_rank', 'moltype',
            'query_containment_ani', 'sum_weighted_found', 'f_orig_query', 'ksize', 'max_containment_ani',
            'std_abund', 'scaled', 'average_containment_ani', 'f_match', 'f_unique_to_query',
            'average_abund', 'unique_intersect_bp', 'median_abund', 'query_md5', 'match_md5', 'remaining_bp',
            'f_unique_weighted'}
    assert  keys == expected_keys


@pytest.mark.parametrize('zip_query', [False, True])
@pytest.mark.parametrize('indexed', [False, True])
def test_missing_querylist(runtmp, capfd, indexed, zip_query):
    # test missing querylist
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    if zip_query:
        query_list = runtmp.output('query.zip')
    # do not make query_list!
    make_file_list(against_list, [sig2, sig47, sig63])

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output('db'))

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)
    assert 'Error: No such file or directory' in captured.err


@pytest.mark.parametrize('indexed', [False, True])
def test_sig_query(runtmp, capfd, indexed):
    # sig file is now fine as a query
    query = get_test_data('SRR606249.sig.gz')

    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output('db'))
        g_output = runtmp.output('out.csv')
    else:
        g_output = runtmp.output('SRR606249.gather.csv')
        p_output = runtmp.output('SRR606249.prefetch.csv')

    runtmp.sourmash('scripts', 'fastmultigather', query, against_list,
                        '-s', '100000', '-o', g_output)

    captured = capfd.readouterr()
    print(captured.err)
    if not indexed:
        # check prefetch output (only non-indexed gather)
        assert os.path.exists(p_output)
        df = pandas.read_csv(p_output)
        assert len(df) == 3
        keys = set(df.keys())
        assert {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'intersect_bp'}.issubset(keys)

    # check gather output (both)
    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    if indexed:
        assert {'query_name', 'query_md5', 'match_name', 'match_md5', 'f_match', 'intersect_bp'}.issubset(keys)
    else:
        assert {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}.issubset(keys)


@pytest.mark.parametrize('indexed', [False, True])
def test_bad_query(runtmp, capfd, indexed):
    # test with a bad query (a .sig.gz file renamed as zip file)
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    query_zip = runtmp.output('query.zip')
    # cp sig2 into query_zip
    with open(query_zip, 'wb') as fp:
        with open(sig2, 'rb') as fp2:
            fp.write(fp2.read())

    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output('db'))

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastmultigather', query_zip, against_list,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)

    assert "InvalidArchive" in captured.err


@pytest.mark.parametrize('indexed', [False, True])
def test_missing_query(runtmp, capfd, indexed):
    # test missing query
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, 'no-exist'])
    make_file_list(against_list, [sig2, sig47, sig63])

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output('db'))

    runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                    '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)
    assert "WARNING: could not load sketches from path 'no-exist'" in captured.err
    assert "WARNING: 1 query paths failed to load. See error messages above."


@pytest.mark.parametrize('indexed', [False, True])
@pytest.mark.parametrize("zip_query", [False, True])
def test_nomatch_query(runtmp, capfd, indexed, zip_query):
    # test nomatch file in querylist
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')
    badsig1 = get_test_data('1.fa.k21.sig.gz')

    make_file_list(query_list, [sig2, badsig1])
    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output('query.zip'))
    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output('db'))

    runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                    '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)
    assert "WARNING: skipped 1 query paths - no compatible signatures." in captured.err


@pytest.mark.parametrize('zip_against', [False, True])
def test_missing_against(runtmp, capfd, zip_against):
    # test missing against
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])

    if zip_against:
        against_list = runtmp.output('against.zip')
    # do not make against_list

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: No such file or directory' in captured.err


def test_sig_against(runtmp, capfd):
    # against file can be a sig now
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')

    g_output = runtmp.output('SRR606249.gather.csv')
    p_output = runtmp.output('SRR606249.prefetch.csv')
    runtmp.sourmash('scripts', 'fastmultigather', query, sig2,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    # check prefetch output (only non-indexed gather)
    assert os.path.exists(p_output)
    df = pandas.read_csv(p_output)
    assert len(df) == 1
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'intersect_bp'}

    # check gather output
    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 1
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}


def test_bad_against(runtmp, capfd):
    # test bad 'against' file - in this case, one containing a nonexistent file
    query = get_test_data('SRR606249.sig.gz')
    query_list = runtmp.output('query.txt')
    make_file_list(query_list, [query])

    against_list = runtmp.output('against.txt')
    sig2 = get_test_data('2.fa.sig.gz')
    make_file_list(against_list, [sig2, "no exist"])

    # should succeed, but with error output.
    runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                    '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert "WARNING: could not load sketches from path 'no exist'" in captured.err
    assert "WARNING: 1 search paths failed to load. See error messages above." in captured.err


@pytest.mark.parametrize('zip_query', [False, True])
def test_bad_against_2(runtmp, capfd, zip_query):
    # test with a bad against (a .sig.gz file renamed as zip file)
    query = get_test_data('SRR606249.sig.gz')
    query_list = runtmp.output('query.txt')
    make_file_list(query_list, [query])

    sig2 = get_test_data('2.fa.sig.gz')
    against_zip = runtmp.output('against.zip')
    # cp sig2 into query_zip
    with open(against_zip, 'wb') as fp:
        with open(sig2, 'rb') as fp2:
            fp.write(fp2.read())

    output = runtmp.output('out.csv')
    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output('query.zip'))

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastmultigather', query_list, against_zip,
                        '-s', '100000', '-o', output)

    captured = capfd.readouterr()
    print(captured.err)

    assert 'InvalidArchive' in captured.err


def test_empty_against(runtmp, capfd):
    # test bad 'against' file - in this case, an empty one
    query = get_test_data('SRR606249.sig.gz')
    query_list = runtmp.output('query.txt')
    make_file_list(query_list, [query])

    against_list = runtmp.output('against.txt')
    make_file_list(against_list, [])

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert "Sketch loading error: No such file or directory" in captured.err
    assert "No search signatures loaded, exiting." in captured.err


@pytest.mark.parametrize('zip_against', [False, True])
def test_nomatch_in_against(runtmp, capfd, zip_against):
    # test an against file that has a non-matching ksize sig in it
    query = get_test_data('SRR606249.sig.gz')
    query_list = runtmp.output('query.txt')
    make_file_list(query_list, [query])

    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig1 = get_test_data('1.fa.k21.sig.gz')
    make_file_list(against_list, [sig2, sig1])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('against.zip'))

    runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                    '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'WARNING: skipped 1 search paths - no compatible signatures.' in captured.err


@pytest.mark.parametrize('zip_query', [False, True])
def test_md5(runtmp, zip_query):
    # test correct md5s present in output
    query = get_test_data('SRR606249.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output('query.zip'))

    cwd = os.getcwd()
    try:
        os.chdir(runtmp.output(''))
        runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                        '-s', '100000', '-t', '0')
    finally:
        os.chdir(cwd)

    print(os.listdir(runtmp.output('')))

    g_output = runtmp.output('SRR606249.gather.csv')
    p_output = runtmp.output('SRR606249.prefetch.csv')

    # check prefetch output (only non-indexed gather)
    assert os.path.exists(p_output)
    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'intersect_bp'}

    md5s = set(df['match_md5'])
    for against_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(against_file, ksize=31):
            assert ss.md5sum() in md5s

    # check gather output (mostly same for indexed vs non-indexed version)
    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}

    md5s = set(df['match_md5'])
    for against_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(against_file, ksize=31):
            assert ss.md5sum() in md5s


@pytest.mark.parametrize('zip_query', [False, True])
def test_md5_indexed(runtmp, zip_query):
    # test correct md5s present in output
    query = get_test_data('SRR606249.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output('query.zip'))

    g_output = runtmp.output('out.csv')
    against_list = index_siglist(runtmp, against_list, runtmp.output('db'))
    runtmp.sourmash('scripts', 'fastmultigather', query_list,
                    against_list, '-s', '100000', '-t', '0',
                    '-o', g_output)

    # check gather output (mostly same for indexed vs non-indexed version)
    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    expected_keys = {'match_name', 'query_filename', 'query_n_hashes', 'match_filename', 'f_match_orig',
            'query_bp', 'query_abundance', 'match_containment_ani', 'intersect_bp', 'total_weighted_hashes',
            'n_unique_weighted_found', 'query_name', 'gather_result_rank', 'moltype',
            'query_containment_ani', 'sum_weighted_found', 'f_orig_query', 'ksize', 'max_containment_ani',
            'std_abund', 'scaled', 'average_containment_ani', 'f_match', 'f_unique_to_query',
            'average_abund', 'unique_intersect_bp', 'median_abund', 'query_md5', 'match_md5', 'remaining_bp',
            'f_unique_weighted'}
    assert keys == expected_keys

    md5s = set(df['match_md5'])
    for against_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(against_file, ksize=31):
            assert ss.md5sum() in md5s


@pytest.mark.parametrize('zip_query', [False, True])
@pytest.mark.parametrize('zip_against', [False, True])
def test_csv_columns_vs_sourmash_prefetch(runtmp, zip_query, zip_against):
    # the column names should be strict subsets of sourmash prefetch cols
    query = get_test_data('SRR606249.sig.gz')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    query_list = runtmp.output('query.txt')
    make_file_list(query_list, [query])
    against_list = runtmp.output('against.txt')
    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output('query.zip'))
    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('against.zip'))

    cwd = os.getcwd()
    try:
        os.chdir(runtmp.output(''))
        runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                        '-s', '100000', '-t', '0')
    finally:
        os.chdir(cwd)

    g_output = runtmp.output('SRR606249.gather.csv')
    p_output = runtmp.output('SRR606249.prefetch.csv')

    assert os.path.exists(p_output)
    assert os.path.exists(g_output)
    # now run sourmash prefetch
    sp_output = runtmp.output('sourmash-prefetch.csv')
    runtmp.sourmash('prefetch', query, against_list,
                    '-o', sp_output, '--scaled', '100000')

    gather_df = pandas.read_csv(g_output)
    g_keys = set(gather_df.keys())
    assert {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}.issubset(g_keys)
    g_keys.remove('rank')       # 'rank' is not in sourmash prefetch!

    sourmash_prefetch_df = pandas.read_csv(sp_output)
    sp_keys = set(sourmash_prefetch_df.keys())
    print(g_keys - sp_keys)
    assert not g_keys - sp_keys, g_keys - sp_keys


def test_csv_columns_vs_sourmash_gather_indexed(runtmp):
    # the column names should be identical to sourmash gather cols
    query = get_test_data('SRR606249.sig.gz')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    query_list = runtmp.output('query.txt')
    make_file_list(query_list, [query])
    against_list = runtmp.output('against.txt')
    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('out.csv')
    against_db = index_siglist(runtmp, against_list, runtmp.output('db'))
    runtmp.sourmash('scripts', 'fastmultigather', query_list,
                    against_db, '-s', '100000', '-t', '0',
                    '-o', g_output)

    assert os.path.exists(g_output)
    # now run sourmash gather
    sg_output = runtmp.output('sourmash-gather.csv')
    runtmp.sourmash('gather', query, against_list,
                    '-o', sg_output, '--scaled', '100000')

    gather_df = pandas.read_csv(g_output)
    g_keys = set(gather_df.keys())
    expected_keys = {'match_name', 'query_filename', 'query_n_hashes', 'match_filename', 'f_match_orig',
            'query_bp', 'query_abundance', 'match_containment_ani', 'intersect_bp', 'total_weighted_hashes',
            'n_unique_weighted_found', 'query_name', 'gather_result_rank', 'moltype',
            'query_containment_ani', 'sum_weighted_found', 'f_orig_query', 'ksize', 'max_containment_ani',
            'std_abund', 'scaled', 'average_containment_ani', 'f_match', 'f_unique_to_query',
            'average_abund', 'unique_intersect_bp', 'median_abund', 'query_md5', 'match_md5', 'remaining_bp',
            'f_unique_weighted'}
    assert g_keys == expected_keys

    sourmash_gather_df = pandas.read_csv(sg_output)
    sg_keys = set(sourmash_gather_df.keys())
    print(sg_keys)
    modified_keys = ["match_md5", "match_name", "match_filename"]
    sg_keys.update(modified_keys) # fastmultigather is more explicit (match_md5 instead of md5, etc)
    print('g_keys - sg_keys:', g_keys - sg_keys)
    assert not g_keys - sg_keys, g_keys - sg_keys


def test_simple_protein(runtmp):
    # test basic protein execution
    sigs = get_test_data('protein.zip')

    sig_names = ["GCA_001593935", "GCA_001593925"]

    runtmp.sourmash('scripts', 'fastmultigather', sigs, sigs,
                    '-s', '100', '--moltype', 'protein', '-k', '19')

    for qsig in sig_names:
        g_output = runtmp.output(os.path.join(qsig + '.gather.csv'))
        p_output = runtmp.output(os.path.join(qsig + '.prefetch.csv'))
        print(g_output)
        assert os.path.exists(g_output)
        assert os.path.exists(p_output)

        df = pandas.read_csv(g_output)
        assert len(df) == 1
        keys = set(df.keys())
        assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}
        print(df)
        # since we're just matching to identical sigs, the md5s should be the same
        assert df['query_md5'][0] == df['match_md5'][0]


def test_simple_dayhoff(runtmp):
    # test basic protein execution
    sigs = get_test_data('dayhoff.zip')

    sig_names = ["GCA_001593935", "GCA_001593925"]

    runtmp.sourmash('scripts', 'fastmultigather', sigs, sigs,
                    '-s', '100', '--moltype', 'dayhoff', '-k', '19')

    for qsig in sig_names:
        g_output = runtmp.output(os.path.join(qsig + '.gather.csv'))
        p_output = runtmp.output(os.path.join(qsig + '.prefetch.csv'))
        print(g_output)
        assert os.path.exists(g_output)
        assert os.path.exists(p_output)

        df = pandas.read_csv(g_output)
        assert len(df) == 1
        keys = set(df.keys())
        assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}
        print(df)
        # since we're just matching to identical sigs, the md5s should be the same
        assert df['query_md5'][0] == df['match_md5'][0]


def test_simple_hp(runtmp):
    # test basic protein execution
    sigs = get_test_data('hp.zip')

    sig_names = ["GCA_001593935", "GCA_001593925"]

    runtmp.sourmash('scripts', 'fastmultigather', sigs, sigs,
                    '-s', '100', '--moltype', 'hp', '-k', '19')

    for qsig in sig_names:
        g_output = runtmp.output(os.path.join(qsig + '.gather.csv'))
        p_output = runtmp.output(os.path.join(qsig + '.prefetch.csv'))
        print(g_output)
        assert os.path.exists(g_output)
        assert os.path.exists(p_output)

        df = pandas.read_csv(g_output)
        assert len(df) == 1
        keys = set(df.keys())
        assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}
        print(df)
        # since we're just matching to identical sigs, the md5s should be the same
        assert df['query_md5'][0] == df['match_md5'][0]


def test_simple_protein_indexed(runtmp):
    # test basic protein execution
    sigs = get_test_data('protein.zip')

    sigs_db = index_siglist(runtmp, sigs, runtmp.output('db'), ksize=19, moltype='protein', scaled=100)

    out_csv = runtmp.output('out.csv')
    runtmp.sourmash('scripts', 'fastmultigather', sigs, sigs_db,
                    '-s', '100', '--moltype', 'protein', '-k', '19',
                    '-o', out_csv)

    assert os.path.exists(out_csv)

    df = pandas.read_csv(out_csv)
    assert len(df) == 2
    keys = set(df.keys())
    expected_keys = {'match_name', 'query_filename', 'query_n_hashes', 'match_filename', 'f_match_orig',
            'query_bp', 'query_abundance', 'match_containment_ani', 'intersect_bp', 'total_weighted_hashes',
            'n_unique_weighted_found', 'query_name', 'gather_result_rank', 'moltype',
            'query_containment_ani', 'sum_weighted_found', 'f_orig_query', 'ksize', 'max_containment_ani',
            'std_abund', 'scaled', 'average_containment_ani', 'f_match', 'f_unique_to_query',
            'average_abund', 'unique_intersect_bp', 'median_abund', 'query_md5', 'match_md5', 'remaining_bp',
            'f_unique_weighted'}
    assert keys == expected_keys
    print(df)
    # since we're just matching to identical sigs, the md5s should be the same
    assert df['query_md5'][0] == df['match_md5'][0]
    assert df['query_md5'][1] == df['match_md5'][1]


def test_simple_dayhoff_indexed(runtmp):
    # test basic protein execution
    sigs = get_test_data('dayhoff.zip')

    sigs_db = index_siglist(runtmp, sigs, runtmp.output('db'), ksize=19, moltype='dayhoff', scaled=100)

    out_csv = runtmp.output('out.csv')
    runtmp.sourmash('scripts', 'fastmultigather', sigs, sigs_db,
                    '-s', '100', '--moltype', 'dayhoff', '-k', '19',
                    '-o', out_csv)

    assert os.path.exists(out_csv)

    df = pandas.read_csv(out_csv)
    assert len(df) == 2
    keys = set(df.keys())
    expected_keys = {'match_name', 'query_filename', 'query_n_hashes', 'match_filename', 'f_match_orig',
            'query_bp', 'query_abundance', 'match_containment_ani', 'intersect_bp', 'total_weighted_hashes',
            'n_unique_weighted_found', 'query_name', 'gather_result_rank', 'moltype',
            'query_containment_ani', 'sum_weighted_found', 'f_orig_query', 'ksize', 'max_containment_ani',
            'std_abund', 'scaled', 'average_containment_ani', 'f_match', 'f_unique_to_query',
            'average_abund', 'unique_intersect_bp', 'median_abund', 'query_md5', 'match_md5', 'remaining_bp',
            'f_unique_weighted'}
    assert keys == expected_keys
    print(df)
    # since we're just matching to identical sigs, the md5s should be the same
    assert df['query_md5'][0] == df['match_md5'][0]
    assert df['query_md5'][1] == df['match_md5'][1]


def test_simple_hp_indexed(runtmp):
    # test basic protein execution
    sigs = get_test_data('hp.zip')

    sigs_db = index_siglist(runtmp, sigs, runtmp.output('db'), ksize=19, moltype='hp', scaled=100)

    out_csv = runtmp.output('out.csv')
    runtmp.sourmash('scripts', 'fastmultigather', sigs, sigs_db,
                    '-s', '100', '--moltype', 'hp', '-k', '19',
                    '-o', out_csv)

    assert os.path.exists(out_csv)

    df = pandas.read_csv(out_csv)
    assert len(df) == 2
    keys = set(df.keys())
    expected_keys = {'match_name', 'query_filename', 'query_n_hashes', 'match_filename', 'f_match_orig',
            'query_bp', 'query_abundance', 'match_containment_ani', 'intersect_bp', 'total_weighted_hashes',
            'n_unique_weighted_found', 'query_name', 'gather_result_rank', 'moltype',
            'query_containment_ani', 'sum_weighted_found', 'f_orig_query', 'ksize', 'max_containment_ani',
            'std_abund', 'scaled', 'average_containment_ani', 'f_match', 'f_unique_to_query',
            'average_abund', 'unique_intersect_bp', 'median_abund', 'query_md5', 'match_md5', 'remaining_bp',
            'f_unique_weighted'}
    assert keys == expected_keys
    print(df)
    # since we're just matching to identical sigs, the md5s should be the same
    assert df['query_md5'][0] == df['match_md5'][0]
    assert df['query_md5'][1] == df['match_md5'][1]


def test_indexed_full_output(runtmp):
    # test correct md5s present in output
    query = get_test_data('SRR606249.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('out.csv')
    against_db = index_siglist(runtmp, against_list, runtmp.output('rocksdb'))
    runtmp.sourmash('scripts', 'fastmultigather', query_list,
                    against_db, '-s', '100000', '-t', '0',
                    '-o', g_output)

    # check full gather output
    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    expected_keys = {'match_name', 'query_filename', 'query_n_hashes', 'match_filename', 'f_match_orig',
            'query_bp', 'query_abundance', 'match_containment_ani', 'intersect_bp', 'total_weighted_hashes',
            'n_unique_weighted_found', 'query_name', 'gather_result_rank', 'moltype',
            'query_containment_ani', 'sum_weighted_found', 'f_orig_query', 'ksize', 'max_containment_ani',
            'std_abund', 'scaled', 'average_containment_ani', 'f_match', 'f_unique_to_query',
            'average_abund', 'unique_intersect_bp', 'median_abund', 'query_md5', 'match_md5', 'remaining_bp',
            'f_unique_weighted'}
    assert keys == expected_keys
    results = df.values.tolist()

    # check a few columns
    average_ani = set(df['average_containment_ani'])
    avg_ani = set([round(x, 4) for x in average_ani])
    assert avg_ani == {0.8602, 0.8504, 0.8361}

    f_unique_weighted = set(df['f_unique_weighted'])
    f_unique_weighted = set([round(x, 4) for x in f_unique_weighted])
    assert f_unique_weighted == {0.0063, 0.002, 0.0062}

    unique_intersect_bp = set(df['unique_intersect_bp'])
    unique_intersect_bp = set([round(x,4) for x in unique_intersect_bp])
    assert unique_intersect_bp == {44000, 18000, 22000}
