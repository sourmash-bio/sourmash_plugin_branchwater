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

    g_output = runtmp.output('SRR606249.sig.gz.gather.csv')
    p_output = runtmp.output('SRR606249.sig.gz.prefetch.csv')
    assert os.path.exists(p_output)

    # check prefetch output (only non-indexed gather)
    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'intersect_bp'}

    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}


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

    # outputs are based on md5sum, e.g. "{md5}.sig.gz.gather.csv"
    g_output = runtmp.output('dec29ca72e68db0f15de0b1b46f82fc5.sig.gz.gather.csv')
    p_output = runtmp.output('dec29ca72e68db0f15de0b1b46f82fc5.sig.gz.prefetch.csv')

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
    against_db = index_siglist(runtmp, against_list, runtmp.output('db'))
    runtmp.sourmash('scripts', 'fastmultigather', query_list,
                        against_db, '-s', '100000', '-t', '0',
                        '-o', g_output)

    assert os.path.exists(g_output)
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_name', 'query_md5', 'match_name', 'match_md5', 'f_match_query', 'intersect_bp'}


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

    assert 'Error: No such file or directory ' in captured.err


@pytest.mark.parametrize('indexed', [False, True])
def test_bad_query(runtmp, capfd, indexed):
    # test bad querylist (a sig file)
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    if indexed:
        against_list = index_siglist(runtmp, against_list, runtmp.output('db'))

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastmultigather', sig2, against_list,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: invalid line in fromfile ' in captured.err


@pytest.mark.parametrize('indexed', [False, True])
def test_bad_query_2(runtmp, capfd, indexed):
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

    assert 'Error: invalid Zip archive: Could not find central directory end' in captured.err


@pytest.mark.parametrize('indexed', [False, True])
def test_missing_query(runtmp, capfd, indexed):
    # test missingfile in querylist
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

    if zip_query:
        assert "WARNING: no compatible sketches in path " not in captured.err
    else:
        assert "WARNING: no compatible sketches in path " in captured.err
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

    assert 'Error: No such file or directory ' in captured.err


def test_bad_against(runtmp, capfd):
    # test bad 'against' file - in this case, use a .sig.gz file.
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastmultigather', query, sig2,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: invalid line in fromfile ' in captured.err


def test_bad_against_2(runtmp, capfd):
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
def test_bad_against_3(runtmp, capfd, zip_query):
    # test with a bad query (a .sig.gz file renamed as zip file)
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
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: invalid Zip archive: Could not find central directory end' in captured.err


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

    assert "Loaded 0 search signature(s)" in captured.err
    assert "Error: No search signatures loaded, exiting." in captured.err


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

    g_output = runtmp.output('SRR606249.sig.gz.gather.csv')
    p_output = runtmp.output('SRR606249.sig.gz.prefetch.csv')
    if zip_query:
        g_output = runtmp.output('dec29ca72e68db0f15de0b1b46f82fc5.sig.gz.gather.csv')
        p_output = runtmp.output('dec29ca72e68db0f15de0b1b46f82fc5.sig.gz.prefetch.csv')

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
    assert keys == {'query_name', 'query_md5', 'match_name', 'match_md5', 'f_match_query', 'intersect_bp'}

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

    g_output = runtmp.output('SRR606249.sig.gz.gather.csv')
    p_output = runtmp.output('SRR606249.sig.gz.prefetch.csv')
    if zip_query:
        g_output = runtmp.output('dec29ca72e68db0f15de0b1b46f82fc5.sig.gz.gather.csv')
        p_output = runtmp.output('dec29ca72e68db0f15de0b1b46f82fc5.sig.gz.prefetch.csv')

    assert os.path.exists(p_output)
    assert os.path.exists(g_output)
    # now run sourmash prefetch
    sp_output = runtmp.output('sourmash-prefetch.csv')
    runtmp.sourmash('prefetch', query, against_list,
                    '-o', sp_output, '--scaled', '100000')

    gather_df = pandas.read_csv(g_output)
    g_keys = set(gather_df.keys())
    assert g_keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}
    g_keys.remove('rank')       # 'rank' is not in sourmash prefetch!

    sourmash_prefetch_df = pandas.read_csv(sp_output)
    sp_keys = set(sourmash_prefetch_df.keys())
    print(g_keys - sp_keys)
    assert not g_keys - sp_keys, g_keys - sp_keys


@pytest.mark.parametrize('zip_query', [False, True])
def test_csv_columns_vs_sourmash_prefetch_indexed(runtmp, zip_query):
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

    g_output = runtmp.output('out.csv')
    against_db = index_siglist(runtmp, against_list, runtmp.output('db'))
    runtmp.sourmash('scripts', 'fastmultigather', query_list,
                    against_db, '-s', '100000', '-t', '0',
                    '-o', g_output)

    assert os.path.exists(g_output)
    # now run sourmash prefetch
    sp_output = runtmp.output('sourmash-prefetch.csv')
    runtmp.sourmash('prefetch', query, against_list,
                    '-o', sp_output, '--scaled', '100000')

    gather_df = pandas.read_csv(g_output)
    g_keys = set(gather_df.keys())
    assert g_keys == {'query_name', 'query_md5', 'match_name', 'match_md5', 'f_match_query', 'intersect_bp'}

    sourmash_prefetch_df = pandas.read_csv(sp_output)
    sp_keys = set(sourmash_prefetch_df.keys())
    print(g_keys - sp_keys)
    assert not g_keys - sp_keys, g_keys - sp_keys


def test_simple_protein(runtmp):
    # test basic protein execution
    sigs = get_test_data('protein.zip')

    sig_names = ["GCA_001593935.1_ASM159393v1_protein.faa.gz", "GCA_001593925.1_ASM159392v1_protein.faa.gz"]

    runtmp.sourmash('scripts', 'fastmultigather', sigs, sigs,
                    '-s', '100', '--moltype', 'protein', '-k', '19')

    for qsig in sig_names:
        g_output = runtmp.output(os.path.join(qsig + '.sig.gather.csv'))
        p_output = runtmp.output(os.path.join(qsig + '.sig.prefetch.csv'))
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

    sig_names = ["GCA_001593935.1_ASM159393v1_protein.faa.gz", "GCA_001593925.1_ASM159392v1_protein.faa.gz"]

    runtmp.sourmash('scripts', 'fastmultigather', sigs, sigs,
                    '-s', '100', '--moltype', 'dayhoff', '-k', '19')

    for qsig in sig_names:
        g_output = runtmp.output(os.path.join(qsig + '.sig.gather.csv'))
        p_output = runtmp.output(os.path.join(qsig + '.sig.prefetch.csv'))
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

    sig_names = ["GCA_001593935.1_ASM159393v1_protein.faa.gz", "GCA_001593925.1_ASM159392v1_protein.faa.gz"]

    runtmp.sourmash('scripts', 'fastmultigather', sigs, sigs,
                    '-s', '100', '--moltype', 'hp', '-k', '19')

    for qsig in sig_names:
        g_output = runtmp.output(os.path.join(qsig + '.sig.gather.csv'))
        p_output = runtmp.output(os.path.join(qsig + '.sig.prefetch.csv'))
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
    assert keys == {'query_name', 'query_md5', 'match_name', 'match_md5', 'f_match_query', 'intersect_bp'}
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
    assert keys == {'query_name', 'query_md5', 'match_name', 'match_md5', 'f_match_query', 'intersect_bp'}
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
    assert keys == {'query_name', 'query_md5', 'match_name', 'match_md5', 'f_match_query', 'intersect_bp'}
    print(df)
    # since we're just matching to identical sigs, the md5s should be the same
    assert df['query_md5'][0] == df['match_md5'][0]
    assert df['query_md5'][1] == df['match_md5'][1]
