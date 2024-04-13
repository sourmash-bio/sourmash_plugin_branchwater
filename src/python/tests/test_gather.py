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


def zip_siglist(runtmp, siglist, db):
    runtmp.sourmash('sig', 'cat', siglist,
                    '-o', db)
    return db


def test_installed(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather')

    assert 'usage:  fastgather' in runtmp.last_result.err


@pytest.mark.parametrize('zip_against', [False, True])
def test_simple(runtmp, zip_against):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('against.zip'))

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '-s', '100000')
    assert os.path.exists(g_output)

    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}


@pytest.mark.parametrize('zip_against', [False, True])
def test_simple_with_prefetch(runtmp, zip_against):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('against.zip'))

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
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}

    df = pandas.read_csv(p_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'intersect_bp'}


@pytest.mark.parametrize('zip_against', [False, True])
def test_missing_query(runtmp, capfd, zip_against):
    # test missing query
    query = runtmp.output('no-such-file')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('against.zip'))

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather', query, against_list,
                        '-o', g_output, '--output-prefetch', p_output,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: No such file or directory' in captured.err


@pytest.mark.parametrize('zip_against', [False, True])
def test_bad_query(runtmp, capfd, zip_against):
    # test non-sig query
    query = runtmp.output('no-such-file')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    # query doesn't need to be a sig anymore - sig, zip, or pathlist welcome
    # as long as there's only one sketch that matches params
    make_file_list(query, [sig2,sig47])
    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('against.zip'))

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather', query, against_list,
                        '-o', g_output, '--output-prefetch', p_output,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: Fastgather requires a single query sketch. Check input:' in captured.err


@pytest.mark.parametrize('zip_against', [False, True])
def test_missing_against(runtmp, capfd, zip_against):
    # test missing against
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    # don't make against list
    if zip_against:
        against_list = runtmp.output('against.zip')

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather', query, against_list,
                        '-o', g_output, '--output-prefetch', p_output,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: No such file or directory' in captured.err


def test_sig_against(runtmp, capfd):
    # sig file is ok as against file now 
    query = get_test_data('SRR606249.sig.gz')

    sig2 = get_test_data('2.fa.sig.gz')

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, sig2,
                        '-o', g_output, '--output-prefetch', p_output,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert os.path.exists(g_output)

    df = pandas.read_csv(g_output)
    assert len(df) == 1
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}


def test_bad_against(runtmp, capfd):
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
    assert "WARNING: 1 search paths failed to load. See error messages above." in captured.err


def test_bad_against_2(runtmp, capfd):
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

    assert "WARNING: 1 search paths failed to load. See error messages above." in captured.err


def test_bad_against_3(runtmp, capfd):
    # test with a bad against (a .sig.gz file renamed as zip file)
    query = get_test_data('SRR606249.sig.gz')

    sig2 = get_test_data('2.fa.sig.gz')
    against_zip = runtmp.output('against.zip')
    # cp sig2 into against_zip
    with open(against_zip, 'wb') as fp:
        with open(sig2, 'rb') as fp2:
            fp.write(fp2.read())

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather', query, against_zip,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'InvalidArchive' in captured.err


@pytest.mark.parametrize('zip_against', [False, True])
def test_against_multisigfile(runtmp, zip_against):
    # test against a sigfile that contains multiple sketches
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    combined = runtmp.output('combined.sig.gz')
    runtmp.sourmash('sig', 'cat', sig2, sig47, sig63, '-o', combined)
    make_file_list(against_list, [combined])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('against.zip'))

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')
    df = pandas.read_csv(g_output)
    if zip_against:
        assert len(df) == 3
        print(df)
    else:
        print(df)
        assert len(df) == 1
    # @CTB this is a bug :(. It should load multiple sketches properly!


@pytest.mark.parametrize('zip_against', [False, True])
def test_query_multisigfile(runtmp, capfd, zip_against):
    # test with a sigfile that contains multiple sketches
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    combined = runtmp.output('combined.sig.gz')
    runtmp.sourmash('sig', 'cat', sig2, sig47, sig63, '-o', combined)

    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('against.zip'))

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather', combined, against_list,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')
    # this fails now :)
    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: Fastgather requires a single query sketch. Check input:" in captured.err


@pytest.mark.parametrize('zip_against', [False, True])
def test_against_nomatch(runtmp, capfd, zip_against):
    # test with 'against' file containing a non-matching ksize
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig1 = get_test_data('1.fa.k21.sig.gz')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig1, sig47, sig63])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('against.zip'))

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against_list,
                    '-o', g_output, '--output-prefetch', p_output,
                    '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert 'WARNING: skipped 1 search paths - no compatible signatures.' in captured.err


@pytest.mark.parametrize('zip_against', [False, True])
def test_md5s(runtmp, zip_against):
    # check that the correct md5sums (of the original sketches) are in
    # the output files
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('against.zip'))

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
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}

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
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'intersect_bp'}

    md5s = list(df['match_md5'])
    print(md5s)

    for against_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(against_file, ksize=31):
            assert ss.md5sum() in md5s


@pytest.mark.parametrize('zip_against', [False, True])
def test_csv_columns_vs_sourmash_prefetch(runtmp, zip_against):
    # the column names should be strict subsets of sourmash prefetch cols
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('against.zip'))

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
    assert g_keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}
    g_keys.remove('rank')       # 'rank' is not in sourmash prefetch!

    sourmash_prefetch_df = pandas.read_csv(sp_output)
    sp_keys = set(sourmash_prefetch_df.keys())
    print(g_keys - sp_keys)
    assert not g_keys - sp_keys, g_keys - sp_keys


@pytest.mark.parametrize('zip_against', [False, True])
def test_fastgather_gatherout_as_picklist(runtmp, zip_against):
    # should be able to use fastgather gather output as picklist
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('against.zip'))

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


@pytest.mark.parametrize('zip_against', [False, True])
def test_fastgather_prefetchout_as_picklist(runtmp, zip_against):
    # should be able to use fastgather prefetch output as picklist
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    if zip_against:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('against.zip'))

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


def test_simple_protein(runtmp):
    # test basic protein execution
    sigs = get_test_data('protein.zip')

    query = runtmp.output('query.zip')
    against = runtmp.output('against.zip')
    # extract query from zip file
    runtmp.sourmash('sig', 'extract', sigs, '--name', 'GCA_001593935', '-o', query)
    # extract against from zip file
    runtmp.sourmash('sig', 'extract', sigs, '--name', 'GCA_001593925', '-o', against)

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against,
                    '-o', g_output, '-s', '100', '--moltype', 'protein', '-k', '19',
                    '--threshold', '0')
    assert os.path.exists(g_output)

    df = pandas.read_csv(g_output)
    assert len(df) == 1
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}
    print(df)
    assert df['match_md5'][0] == "16869d2c8a1d29d1c8e56f5c561e585e"


def test_simple_dayhoff(runtmp):
    # test basic protein execution
    sigs = get_test_data('dayhoff.zip')

    query = runtmp.output('query.zip')
    against = runtmp.output('against.zip')
    # extract query from zip file
    runtmp.sourmash('sig', 'extract', sigs, '--name', 'GCA_001593935', '-o', query)
    # extract against from zip file
    runtmp.sourmash('sig', 'extract', sigs, '--name', 'GCA_001593925', '-o', against)

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against,
                    '-o', g_output, '-s', '100', '--moltype', 'dayhoff', '-k', '19',
                    '--threshold', '0')
    assert os.path.exists(g_output)

    df = pandas.read_csv(g_output)
    assert len(df) == 1
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}
    print(df)
    assert df['match_md5'][0] == "fbca5e5211e4d58427997fd5c8343e9a"


def test_simple_hp(runtmp):
    # test basic protein execution
    sigs = get_test_data('hp.zip')

    query = runtmp.output('query.zip')
    against = runtmp.output('against.zip')
    # extract query from zip file
    runtmp.sourmash('sig', 'extract', sigs, '--name', 'GCA_001593935', '-o', query)
    # extract against from zip file
    runtmp.sourmash('sig', 'extract', sigs, '--name', 'GCA_001593925', '-o', against)

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query, against,
                    '-o', g_output, '-s', '100', '--moltype', 'hp', '-k', '19',
                    '--threshold', '0')
    assert os.path.exists(g_output)

    df = pandas.read_csv(g_output)
    assert len(df) == 1
    keys = set(df.keys())
    assert keys == {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}
    print(df)
    assert df['match_md5'][0] == "ea2a1ad233c2908529d124a330bcb672"


def test_indexed_against(runtmp, capfd):
    # do not accept rocksdb for now
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')

    make_file_list(against_list, [sig2])
    db_against = runtmp.output('against.rocksdb')

    ## index against
    runtmp.sourmash('scripts', 'index', against_list,
                    '-o', db_against, '-k', str(31), '--scaled', str(1000),
                    '--moltype', "DNA")

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'fastgather', query, db_against,
                        '-o', g_output, '--output-prefetch', p_output,
                        '-s', '100000')

    captured = capfd.readouterr()
    print(captured.err)

    assert "Cannot load search signatures from a 'rocksdb' database. Please use sig, zip, or pathlist." in captured.err


def test_simple_with_manifest_loading(runtmp):
    # test basic execution!
    query = get_test_data('SRR606249.sig.gz')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])
    query_manifest = runtmp.output("query-manifest.csv")
    against_manifest = runtmp.output("against-manifest.csv")

    runtmp.sourmash("sig", "manifest", query, "-o", query_manifest)
    runtmp.sourmash("sig", "manifest", against_list, "-o", against_manifest)

    g_output = runtmp.output('gather.csv')
    p_output = runtmp.output('prefetch.csv')

    runtmp.sourmash('scripts', 'fastgather', query_manifest, against_manifest,
                    '-o', g_output, '-s', '100000')
    assert os.path.exists(g_output)

    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    assert {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'rank', 'intersect_bp'}.issubset(keys)


def test_simple_full_output(runtmp):
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
                    '-o', g_output, '-s', '100000', '--full-results')
    assert os.path.exists(g_output)

    df = pandas.read_csv(g_output)
    assert len(df) == 3
    keys = set(df.keys())
    print(keys)
    print(df)
    assert {'query_filename', 'query_name', 'query_md5', 'match_name', 'match_md5', 'gather_result_rank', 'intersect_bp'}.issubset(keys)
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


    intersect_bp = set(df['intersect_bp'])
    assert intersect_bp == set([4400000, 4100000, 2200000])
    f_unique_to_query = set([round(x,4) for x in df['f_unique_to_query']])
    assert f_unique_to_query == set([0.0053, 0.0105, 0.0044])
    query_containment_ani = set([round(x,4) for x in df['query_containment_ani']])
    assert query_containment_ani == set([0.8632, 0.8444, 0.8391])
    print(query_containment_ani)
    for index, row in df.iterrows():
        print(row.to_dict())


def test_fullres_vs_sourmash_gather(runtmp):
    # the column names should be identical to sourmash gather cols
    query = get_test_data('SRR606249.sig.gz')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    query_list = runtmp.output('query.txt')
    make_file_list(query_list, [query])
    against_list = runtmp.output('against.txt')
    make_file_list(against_list, [sig2, sig47, sig63])

    g_output = runtmp.output('SRR606249.gather.csv')
    runtmp.sourmash('scripts', 'fastgather', query_list,
                    against_list, '-s', '100000', '-t', '0',
                    '--full-results', '-o', g_output)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert os.path.exists(g_output)
    # now run sourmash gather
    sg_output = runtmp.output('.csv')
    runtmp.sourmash('gather', query, against_list,
                    '-o', sg_output, '--scaled', '100000')

    gather_df = pandas.read_csv(g_output)
    g_keys = set(gather_df.keys())

    sourmash_gather_df = pandas.read_csv(sg_output)
    sg_keys = set(sourmash_gather_df.keys())
    print(sg_keys)
    modified_keys = ["match_md5", "match_name", "match_filename"]
    sg_keys.update(modified_keys) # fastgather is more explicit (match_md5 instead of md5, etc)
    print('g_keys - sg_keys:', g_keys - sg_keys)
    assert not g_keys - sg_keys, g_keys - sg_keys

    for index, row in sourmash_gather_df.iterrows():
        print(row.to_dict())

    fmg_intersect_bp = set(gather_df['intersect_bp'])
    g_intersect_bp = set(sourmash_gather_df['intersect_bp'])
    assert fmg_intersect_bp == g_intersect_bp == set([4400000, 4100000, 2200000])

    fmg_f_orig_query =  set([round(x,4) for x in gather_df['f_orig_query']])
    g_f_orig_query =  set([round(x,4) for x in sourmash_gather_df['f_orig_query']])
    assert fmg_f_orig_query == g_f_orig_query == set([0.0098, 0.0105, 0.0052])

    fmg_f_unique_to_query =  set([round(x,3) for x in gather_df['f_unique_to_query']]) # rounding to 4 --> slightly different!
    g_f_unique_to_query =  set([round(x,3) for x in sourmash_gather_df['f_unique_to_query']])
    assert fmg_f_unique_to_query == g_f_unique_to_query == set([0.004, 0.01, 0.005])

    fmg_f_unique_weighted =  set([round(x,4) for x in gather_df['f_unique_weighted']])
    g_f_unique_weighted =  set([round(x,4) for x in sourmash_gather_df['f_unique_weighted']])
    assert fmg_f_unique_weighted== g_f_unique_weighted == set([0.0063, 0.002, 0.0062])

    fmg_average_abund =  set([round(x,4) for x in gather_df['average_abund']])
    g_average_abund =  set([round(x,4) for x in sourmash_gather_df['average_abund']])
    assert fmg_average_abund== g_average_abund == set([8.2222, 10.3864, 21.0455])

    fmg_median_abund =  set([round(x,4) for x in gather_df['median_abund']])
    g_median_abund =  set([round(x,4) for x in sourmash_gather_df['median_abund']])
    assert fmg_median_abund== g_median_abund == set([8.0, 10.5, 21.5])

    fmg_std_abund =  set([round(x,4) for x in gather_df['std_abund']])
    g_std_abund =  set([round(x,4) for x in sourmash_gather_df['std_abund']])
    assert fmg_std_abund== g_std_abund == set([3.172, 5.6446, 6.9322])

    # we can't get match filename in FMG yet.
    # fmg_match_filename =  set(gather_df['match_filename'])
    # assert fmg_match_filename == [] all are nans rn..{nan, nan, nan}
    g_match_filename =  sourmash_gather_df['filename']
    g_match_filename_basename = [os.path.basename(filename) for filename in sourmash_gather_df['filename']]
    assert all([x in g_match_filename_basename for x in ['2.fa.sig.gz', '63.fa.sig.gz', '47.fa.sig.gz']])

    assert list(sourmash_gather_df['name']) == list(gather_df['match_name'])
    assert list(sourmash_gather_df['md5']) == list(gather_df['match_md5'])

    fmg_f_match_orig =  set([round(x,4) for x in gather_df['f_match_orig']])
    g_f_match_orig =  set([round(x,4) for x in sourmash_gather_df['f_match_orig']])
    assert fmg_f_match_orig == g_f_match_orig == set([1.0])

    fmg_f_match =  set([round(x,4) for x in gather_df['f_match']])
    g_f_match =  set([round(x,4) for x in sourmash_gather_df['f_match']])
    assert fmg_f_match == g_f_match == set([0.439, 1.0])

    fmg_unique_intersect_bp = set(gather_df['unique_intersect_bp'])
    g_unique_intersect_bp = set(sourmash_gather_df['unique_intersect_bp'])
    assert fmg_unique_intersect_bp == g_unique_intersect_bp == set([4400000, 1800000, 2200000])

    fmg_total_weighted_hashes= set(gather_df['total_weighted_hashes'])
    g_total_weighted_hashes = set(sourmash_gather_df['total_weighted_hashes'])
    assert fmg_total_weighted_hashes == g_total_weighted_hashes == set([73489])

    fmg_gather_result_rank= set(gather_df['gather_result_rank'])
    g_gather_result_rank = set(sourmash_gather_df['gather_result_rank'])
    assert fmg_gather_result_rank == g_gather_result_rank == set([0,1,2])

    # FIX remaining_bp
    fmg_remaining_bp = set(gather_df['remaining_bp'])
    g_remaining_bp = set(sourmash_gather_df['remaining_bp'])
    print(g_remaining_bp) #{4000000, 0, 1800000}
    print(fmg_remaining_bp) # {415600000, 411600000, 413400000}
    # assert fmg_remaining_bp == g_remaining_bp == set([])

    # f_unique_to_query = set([round(x,4) for x in df['f_unique_to_query']])
    # assert f_unique_to_query == set([0.0053, 0.0105, 0.0044])
    # query_containment_ani = set([round(x,4) for x in df['query_containment_ani']])
    # assert query_containment_ani == set([0.8632, 0.8444, 0.8391])
