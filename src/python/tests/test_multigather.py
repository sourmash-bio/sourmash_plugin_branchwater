"""
Test 'sourmash scripts fastmultigather'
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
        runtmp.sourmash('scripts', 'fastmultigather')

    assert 'usage:  fastmultigather' in runtmp.last_result.err


def test_simple(runtmp):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    cwd = os.getcwd()
    try:
        os.chdir(runtmp.output(''))
        runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                        '-s', '100000', '-t', '0')
    finally:
        os.chdir(cwd)

    print(os.listdir(runtmp.output('')))

    g_output = runtmp.output('SRR606249.sig.gz.gather.csv')
    assert os.path.exists(g_output)
    p_output = runtmp.output('SRR606249.sig.gz.prefetch.csv')
    assert os.path.exists(p_output)

    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_file', 'match', 'match_md5', 'rank', 'overlap'}

    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_file', 'match', 'match_md5', 'overlap'}


def test_missing_querylist(runtmp, capfd):
    # test missing querylist
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    # do not make query_list!
    make_file_list(against_list, [sig2, sig47, sig63])

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: No such file or directory ' in captured.err


def test_bad_query(runtmp, capfd):
    # test bad querylist (a sig file)
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastmultigather', sig2, against_list,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: invalid line in fromfile ' in captured.err


def test_missing_query(runtmp, capfd):
    # test missingfile in querylist
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, 'no-exist'])
    make_file_list(against_list, [sig2, sig47, sig63])

    runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                    '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert "WARNING: could not load sketches from path 'no-exist'" in captured.err
    assert "WARNING: 1 query paths failed to load. See error messages above."


def test_nomatch_query(runtmp, capfd):
    # test nomatch file in querylist
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')
    badsig1 = get_test_data('1.fa.k21.sig.gz')

    make_file_list(query_list, [sig2, badsig1])
    make_file_list(against_list, [sig2, sig47, sig63])

    runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                    '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert "WARNING: no compatible sketches in path " in captured.err
    assert "WARNING: skipped 1 query paths - no compatible signatures." in captured.err


def test_missing_against(runtmp, capfd):
    # test missing against
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
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

    assert "Loaded 0 sketches to search against." in captured.err
    assert "Error: No sketches loaded to search against!?" in captured.err


def test_nomatch_in_against(runtmp, capfd):
    # test an against file that has a non-matching ksize sig in it
    query = get_test_data('SRR606249.sig.gz')
    query_list = runtmp.output('query.txt')
    make_file_list(query_list, [query])

    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig1 = get_test_data('1.fa.k21.sig.gz')
    make_file_list(against_list, [sig2, sig1])

    runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                    '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'WARNING: skipped 1 search paths - no compatible signatures.' in captured.err


def test_md5(runtmp):
    # test correct md5s present in output
    query = get_test_data('SRR606249.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    make_file_list(query_list, [query])
    make_file_list(against_list, [sig2, sig47, sig63])

    cwd = os.getcwd()
    try:
        os.chdir(runtmp.output(''))
        runtmp.sourmash('scripts', 'fastmultigather', query_list, against_list,
                        '-s', '100000', '-t', '0')
    finally:
        os.chdir(cwd)

    print(os.listdir(runtmp.output('')))

    g_output = runtmp.output('SRR606249.sig.gz.gather.csv')
    assert os.path.exists(g_output)
    p_output = runtmp.output('SRR606249.sig.gz.prefetch.csv')
    assert os.path.exists(p_output)

    # check gather output
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_file', 'match', 'match_md5', 'rank', 'overlap'}

    md5s = set(df['match_md5'])
    for against_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(against_file, ksize=31):
            assert ss.md5sum() in md5s

    # check prefetch output
    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_file', 'match', 'match_md5', 'overlap'}

    md5s = set(df['match_md5'])
    for against_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(against_file, ksize=31):
            assert ss.md5sum() in md5s
