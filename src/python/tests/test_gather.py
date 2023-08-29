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
        runtmp.sourmash('scripts', 'fastgather')

    assert 'usage:  fastgather' in runtmp.last_result.err


def test_simple(runtmp):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '-s', '100000')
    assert os.path.exists(g_output)

    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'match_name', 'match_md5', 'rank', 'intersect_bp'}


def test_simple_with_prefetch(runtmp):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')
    assert os.path.exists(g_output)
    assert os.path.exists(p_output)

    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'match_name', 'match_md5', 'rank', 'intersect_bp'}

    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'match_name', 'match_md5', 'intersect_bp'}


def test_missing_query(runtmp, capfd):
    # test missing query
    query = runtmp.output('no-such-file')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather', query, against_list,
                        '-o', g_output, '--output-prefetch', p_output,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: No such file or directory ' in captured.err


def test_bad_query(runtmp, capfd):
    # test non-sig query
    query = runtmp.output('no-such-file')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    # since 'query' needs to be a sig, this breaks it.
    make_file_list(query, [sig2])

    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather', query, against_list,
                        '-o', g_output, '--output-prefetch', p_output,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: expected value at line 1' in captured.err


def test_missing_against(runtmp, capfd):
    # test missing against
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    #make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather', query, against_list,
                        '-o', g_output, '--output-prefetch', p_output,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: No such file or directory ' in captured.err


def test_bad_against(runtmp, capfd):
    # test bad 'against' file - in this case, use a .sig.gz file.
    query = get_test_data('SRR606249.sig.gz')

    sig2 = get_test_data('2.fa.sig.gz')

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather', query, sig2,
                        '-o', g_output, '--output-prefetch', p_output,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: invalid line in fromfile ' in captured.err


def test_bad_against_2(runtmp, capfd):
    # test bad 'against' file - in this case, one containing a bad filename.
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    make_file_list(against_list, [sig2, 'no-exist'])


    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert "WARNING: could not load sketches from path 'no-exist'" in captured.err
    assert "WARNING: 1 signature paths failed to load. See error messages above." in captured.err


def test_bad_against_3(runtmp, capfd):
    # test bad 'against' file - in this case, one containing an empty file
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    empty_file = runtmp.output('empty.sig')
    with open(empty_file, 'wb') as fp:
        pass
    make_file_list(against_list, [sig2, empty_file])


    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert "Sketch loading error: File is too short, less than five bytes" in captured.err
    assert "WARNING: could not load sketches from path" in captured.err

    assert "WARNING: 1 signature paths failed to load. See error messages above." in captured.err


def test_against_multisigfile(runtmp):
    # test against a sigfile that contains multiple sketches
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    combined = runtmp.output('combined.sig.gz')
    runtmp.sourmash('sig', 'cat', sig2, sig47, sig63, '-o', combined)
    make_file_list(against_list, [combined])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')
    df = pandas.read_csv(g_output)
    assert len(df) == 1
    # @CTB this is a bug :(. It should load multiple sketches properly!


def test_query_multisigfile(runtmp):
    # test with a sigfile that contains multiple sketches
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    combined = runtmp.output('combined.sig.gz')
    runtmp.sourmash('sig', 'cat', sig2, sig47, sig63, '-o', combined)

    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', combined, against_list,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')
    # @CTB this should fail, not succeed :(.
    df = pandas.read_csv(g_output)
    assert len(df) == 1


def test_against_nomatch(runtmp, capfd):
    # test with 'against' file containing a non-matching ksize
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig1 = get_test_data('1.fa.k21.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig1, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'WARNING: skipped 1 paths - no compatible signatures.' in captured.err


def test_md5s(runtmp):
    # check that the correct md5sums (of the original sketches) are in
    # the output files
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')
    assert os.path.exists(g_output)
    assert os.path.exists(p_output)

    # test gather output!
    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'match_name', 'match_md5', 'rank', 'intersect_bp'}

    md5s = list(df['match_md5'])
    print(md5s)

    for against_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(against_file, ksize=31):
            assert ss.md5sum() in md5s

    # test prefetch output!
    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())

    # prefetch output has no rank.
    assert keys == {'query_filename', 'match_name', 'match_md5', 'intersect_bp'}

    md5s = list(df['match_md5'])
    print(md5s)

    for against_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(against_file, ksize=31):
            assert ss.md5sum() in md5s


def test_csv_columns_vs_sourmash_prefetch(runtmp):
    # the column names should be strict subsets of sourmash prefetch cols
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    # first run fastgather
    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')
    assert os.path.exists(g_output)
    assert os.path.exists(p_output)

    # now run sourmash prefetch
    sp_output = runtmp.output('sourmash-prefetch.csv')
    runtmp.sourmash('prefetch', query, against_list,
                    '-o', sp_output, '--scaled', '100000')

    gather_df = pandas.read_csv(g_output)
    g_keys = set(gather_df.keys())
    assert g_keys == {'query_filename', 'match_name', 'match_md5', 'rank', 'intersect_bp'}
    g_keys.remove('rank')       # 'rank' is not in sourmash prefetch!

    sourmash_prefetch_df = pandas.read_csv(sp_output)
    sp_keys = set(sourmash_prefetch_df.keys())
    print(g_keys - sp_keys)
    assert not g_keys - sp_keys, g_keys - sp_keys


def test_fastgather_gatherout_as_picklist(runtmp):
    # should be able to use fastgather gather output as picklist
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    # first run fastgather
    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')
    assert os.path.exists(g_output)
    assert os.path.exists(p_output)

    # now run sourmash gather using as picklist as picklist
    gather_picklist_output = runtmp.output('sourmash-gather+picklist.csv')
    runtmp.sourmash('gather', query, against_list,
                    '-o', gather_picklist_output, '--scaled', '100000',
                    '--picklist', f'{g_output}:match_name:ident')

    # finally, run sourmash gather using fastgather gather output as picklist
    full_gather_output = runtmp.output('sourmash-gather.csv')
    runtmp.sourmash('gather', query, against_list,
                    '-o', full_gather_output, '--scaled', '100000')

    picklist_df = pandas.read_csv(gather_picklist_output)
    full_df = pandas.read_csv(full_gather_output)

    assert picklist_df.equals(full_df)


def test_fastgather_prefetchout_as_picklist(runtmp):
    # should be able to use fastgather prefetch output as picklist
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    # first run fastgather
    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')
    assert os.path.exists(g_output)
    assert os.path.exists(p_output)

    # now run sourmash gather using fastgather prefetch output as picklist
    gather_picklist_output = runtmp.output('sourmash-gather+picklist.csv')
    runtmp.sourmash('gather', query, against_list,
                    '-o', gather_picklist_output, '--scaled', '100000',
                    '--picklist', f'{p_output}:match_name:ident')

    # finally, run sourmash gather using as picklist as picklist
    full_gather_output = runtmp.output('sourmash-gather.csv')
    runtmp.sourmash('gather', query, against_list,
                    '-o', full_gather_output, '--scaled', '100000')

    picklist_df = pandas.read_csv(gather_picklist_output)
    full_df = pandas.read_csv(full_gather_output)

    assert picklist_df.equals(full_df)
