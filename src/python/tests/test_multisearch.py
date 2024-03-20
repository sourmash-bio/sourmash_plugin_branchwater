import os
import csv
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
        runtmp.sourmash('scripts', 'multisearch')

    assert 'usage:  multisearch' in runtmp.last_result.err

def zip_siglist(runtmp, siglist, db):
    runtmp.sourmash('sig', 'cat', siglist,
                    '-o', db)
    return db

@pytest.mark.parametrize("zip_query", [False, True])
@pytest.mark.parametrize("zip_db", [False, True])
def test_simple_no_ani(runtmp, zip_query, zip_db):
    # test basic execution!
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    if zip_db:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('db.zip'))
    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output('query.zip'))

    runtmp.sourmash('scripts', 'multisearch', query_list, against_list,
                    '-o', output)
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 5

    dd = df.to_dict(orient='index')
    print(dd)

    for idx, row in dd.items():
        # identical?
        if row['match_name'] == row['query_name']:
            assert row['query_md5'] == row['match_md5'], row
            assert float(row['containment'] == 1.0)
            assert float(row['jaccard'] == 1.0)
            assert float(row['max_containment'] == 1.0)
            assert 'query_containment_ani' not in row
            assert 'match_containment_ani' not in row
            assert 'average_containment_ani' not in row
            assert 'max_containment_ani' not in row

        else:
            # confirm hand-checked numbers
            q = row['query_name'].split()[0]
            m = row['match_name'].split()[0]
            cont = float(row['containment'])
            jaccard = float(row['jaccard'])
            maxcont = float(row['max_containment'])
            intersect_hashes = int(row['intersect_hashes'])

            jaccard = round(jaccard, 4)
            cont = round(cont, 4)
            maxcont = round(maxcont, 4)
            print(q, m, f"{jaccard:.04}", f"{cont:.04}", f"{maxcont:.04}")

            if q == 'NC_011665.1' and m == 'NC_009661.1':
                assert jaccard == 0.3207
                assert cont == 0.4828
                assert maxcont == 0.4885
                assert intersect_hashes == 2529

            if q == 'NC_009661.1' and m == 'NC_011665.1':
                assert jaccard == 0.3207
                assert cont == 0.4885
                assert maxcont == 0.4885
                assert intersect_hashes == 2529


@pytest.mark.parametrize("zip_query", [False, True])
@pytest.mark.parametrize("zip_db", [False, True])
def test_simple_ani(runtmp, zip_query, zip_db):
    # test basic execution!
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    if zip_db:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('db.zip'))
    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output('query.zip'))

    runtmp.sourmash('scripts', 'multisearch', query_list, against_list,
                    '-o', output, '--ani')
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 5

    dd = df.to_dict(orient='index')
    print(dd)

    for idx, row in dd.items():
        # identical?
        if row['match_name'] == row['query_name']:
            assert row['query_md5'] == row['match_md5'], row
            assert float(row['containment'] == 1.0)
            assert float(row['jaccard'] == 1.0)
            assert float(row['max_containment'] == 1.0)
            assert float(row['query_containment_ani'] == 1.0)
            assert float(row['match_containment_ani'] == 1.0)
            assert float(row['average_containment_ani'] == 1.0)
            assert float(row['max_containment_ani'] == 1.0)

        else:
            # confirm hand-checked numbers
            q = row['query_name'].split()[0]
            m = row['match_name'].split()[0]
            cont = float(row['containment'])
            jaccard = float(row['jaccard'])
            maxcont = float(row['max_containment'])
            intersect_hashes = int(row['intersect_hashes'])
            q1_ani = float(row['query_containment_ani'])
            q2_ani = float(row['match_containment_ani'])
            avg_ani = float(row['average_containment_ani'])
            max_ani = float(row['max_containment_ani'])


            jaccard = round(jaccard, 4)
            cont = round(cont, 4)
            maxcont = round(maxcont, 4)
            q1_ani = round(q1_ani, 4)
            q2_ani = round(q2_ani, 4)
            avg_ani = round(avg_ani, 4)
            max_ani = round(max_ani, 4)
            print(q, m, f"{jaccard:.04}", f"{cont:.04}", f"{maxcont:.04}", f"{q1_ani:.04}", f"{q2_ani:.04}", f"{avg_ani:.04}", f"{max_ani:.04}")

            if q == 'NC_011665.1' and m == 'NC_009661.1':
                assert jaccard == 0.3207
                assert cont == 0.4828
                assert maxcont == 0.4885
                assert intersect_hashes == 2529
                assert q1_ani == 0.9768
                assert q2_ani == 0.9772
                assert avg_ani == 0.977
                assert max_ani == 0.9772

            if q == 'NC_009661.1' and m == 'NC_011665.1':
                assert jaccard == 0.3207
                assert cont == 0.4885
                assert maxcont == 0.4885
                assert intersect_hashes == 2529
                assert q1_ani == 0.9772
                assert q2_ani == 0.9768
                assert avg_ani == 0.977
                assert max_ani == 0.9772


@pytest.mark.parametrize("zip_query", [False, True])
@pytest.mark.parametrize("zip_db", [False, True])
def test_simple_threshold(runtmp, zip_query, zip_db):
    # test with a simple threshold => only 3 results
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    if zip_db:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('db.zip'))
    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output('query.zip'))

    runtmp.sourmash('scripts', 'multisearch', query_list, against_list,
                    '-o', output, '-t', '0.5')
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 3


def test_simple_manifest(runtmp):
    # test with a simple threshold => only 3 results
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    query_mf = runtmp.output('qmf.csv')
    against_mf = runtmp.output('amf.csv')

    runtmp.sourmash("sig", "manifest", query_list, "-o", query_mf)
    runtmp.sourmash("sig", "manifest", against_list, "-o", against_mf)

    output = runtmp.output('out.csv')

    runtmp.sourmash('scripts', 'multisearch', query_mf, against_mf,
                    '-o', output, '-t', '0.5')
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 3


@pytest.mark.parametrize("zip_query", [False, True])
def test_missing_query(runtmp, capfd, zip_query):
    # test with a missing query list
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    if zip_query:
        query_list = runtmp.output('query.zip')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'multisearch', query_list, against_list,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: No such file or directory' in captured.err


def test_sig_query(runtmp, capfd):
    # sig is ok as query now
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    runtmp.sourmash('scripts', 'multisearch', sig2, against_list,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)

    assert os.path.exists(output)
    df = pandas.read_csv(output)
    assert len(df) == 1


def test_bad_query(runtmp, capfd):
    # test with a bad query list (a missing file)
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')
    make_file_list(query_list, [sig2, "no-exist"])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    runtmp.sourmash('scripts', 'multisearch', query_list, against_list,
                    '-o', output)

    captured = capfd.readouterr()
    print(captured.err)

    assert "WARNING: could not load sketches from path 'no-exist'" in captured.err
    assert "WARNING: 1 query paths failed to load. See error messages above." in captured.err


def test_bad_query_3(runtmp, capfd):
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

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'multisearch', query_zip, against_list,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)

    assert 'InvalidArchive' in captured.err


@pytest.mark.parametrize("zip_db", [False, True])
def test_missing_against(runtmp, capfd, zip_db):
    # test with a missing against list
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    # do not create against_list

    if zip_db:
        #specify .zip but don't create the file
        against_list = runtmp.output('db.zip')

    output = runtmp.output('out.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'multisearch', query_list, against_list,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)

    assert 'Error: No such file or directory' in captured.err


def test_sig_against(runtmp, capfd):
    # against can be sig now
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    runtmp.sourmash('scripts', 'multisearch', query_list, sig2,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)

    assert os.path.exists(output)
    df = pandas.read_csv(output)
    assert len(df) == 1


def test_bad_against(runtmp, capfd):
    # test with a bad against list (a missing file)
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')
    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, "no-exist"])

    output = runtmp.output('out.csv')

    runtmp.sourmash('scripts', 'multisearch', query_list, against_list,
                    '-o', output)

    captured = capfd.readouterr()
    print(captured.err)

    assert "WARNING: could not load sketches from path 'no-exist'" in captured.err
    assert "WARNING: 1 search paths failed to load. See error messages above." in captured.err


def test_empty_query(runtmp, capfd):
    # test with an empty query list - fail with error
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'multisearch', query_list, against_list,
                        '-o', output)

    print(runtmp.last_result.err)
    captured = capfd.readouterr()
    print(captured.err)
    assert "No query signatures loaded, exiting." in captured.err
    # @CTB


@pytest.mark.parametrize("zip_query", [False, True])
def test_nomatch_query(runtmp, capfd, zip_query):
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

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output('query.zip'))

    runtmp.sourmash('scripts', 'multisearch', query_list, against_list,
                    '-o', output)
    assert os.path.exists(output)

    captured = capfd.readouterr()
    print(captured.err)

    assert 'WARNING: skipped 1 query paths - no compatible signatures' in captured.err


@pytest.mark.parametrize("zip_db", [False, True])
def test_load_only_one_bug(runtmp, capfd, zip_db):
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

    if zip_db:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('db.zip'))

    output = runtmp.output('out.csv')

    runtmp.sourmash('scripts', 'multisearch', query_list, against_list,
                    '-o', output)
    assert os.path.exists(output)

    captured = capfd.readouterr()
    print(captured.err)

    assert not 'WARNING: skipped 1 paths - no compatible signatures.' in captured.err
    assert not 'WARNING: no compatible sketches in path' in captured.err


@pytest.mark.parametrize("zip_query", [False, True])
def test_load_only_one_bug_as_query(runtmp, capfd, zip_query):
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

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output('query.zip'))

    output = runtmp.output('out.csv')

    runtmp.sourmash('scripts', 'multisearch', query_list, against_list,
                    '-o', output)
    assert os.path.exists(output)

    captured = capfd.readouterr()
    print(captured.err)

    assert not 'WARNING: skipped 1 paths - no compatible signatures.' in captured.err
    assert not 'WARNING: no compatible sketches in path ' in captured.err


@pytest.mark.parametrize("zip_query", [False, True])
@pytest.mark.parametrize("zip_db", [False, True])
def test_md5(runtmp, zip_query, zip_db):
    # test that md5s match what was in the original files, not downsampled etc.
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output('query.zip'))
    if zip_db:
        against_list = zip_siglist(runtmp, against_list, runtmp.output('db.zip'))

    runtmp.sourmash('scripts', 'multisearch', query_list, against_list,
                    '-o', output)
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 5

    md5s = list(df['query_md5'])
    print(md5s)

    for query_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(query_file, ksize=31):
            assert ss.md5sum() in md5s

    md5s = list(df['match_md5'])
    print(md5s)

    for against_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(against_file, ksize=31):
            assert ss.md5sum() in md5s


def test_simple_prot(runtmp):
    # test basic execution with protein sigs
    sigs = get_test_data('protein.zip')

    output = runtmp.output('out.csv')

    runtmp.sourmash('scripts', 'multisearch', sigs, sigs,
                    '-o', output, '--moltype', 'protein',
                    '-k', '19', '--scaled', '100', '--ani')
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 4

    dd = df.to_dict(orient='index')
    print(dd)

    for idx, row in dd.items():
        # identical?
        if row['match_name'] == row['query_name']:
            assert row['query_md5'] == row['match_md5'], row
            assert float(row['containment'] == 1.0)
            assert float(row['jaccard'] == 1.0)
            assert float(row['max_containment'] == 1.0)
            assert float(row['query_containment_ani'] == 1.0)
            assert float(row['match_containment_ani'] == 1.0)
            assert float(row['average_containment_ani'] == 1.0)
            assert float(row['max_containment_ani'] == 1.0)

        else:
            # confirm hand-checked numbers
            q = row['query_name'].split()[0]
            m = row['match_name'].split()[0]
            cont = float(row['containment'])
            jaccard = float(row['jaccard'])
            maxcont = float(row['max_containment'])
            intersect_hashes = int(row['intersect_hashes'])
            q1_ani = float(row['query_containment_ani'])
            q2_ani = float(row['match_containment_ani'])
            avg_ani = float(row['average_containment_ani'])
            max_ani = float(row['max_containment_ani'])

            jaccard = round(jaccard, 4)
            cont = round(cont, 4)
            maxcont = round(maxcont, 4)
            q1_ani = round(q1_ani, 4)
            q2_ani = round(q2_ani, 4)
            avg_ani = round(avg_ani, 4)
            max_ani = round(max_ani, 4)
            print(q, m, f"{jaccard:.04}", f"{cont:.04}", f"{maxcont:.04}", intersect_hashes, f"{q1_ani:.04}", f"{q2_ani:.04}", f"{avg_ani:.04}", f"{max_ani:.04}")

            if q == 'GCA_001593925' and m == 'GCA_001593935':
                assert jaccard == 0.0434
                assert cont == 0.1003
                assert maxcont == 0.1003
                assert intersect_hashes == 342
                assert q1_ani == 0.886
                assert q2_ani == 0.8702
                assert avg_ani == 0.8781
                assert max_ani == 0.886

            if q == 'GCA_001593935' and m == 'GCA_001593925':
                assert jaccard == 0.0434
                assert cont == 0.0712
                assert maxcont == 0.1003
                assert intersect_hashes == 342
                assert q1_ani == 0.8702
                assert q2_ani == 0.886
                assert avg_ani == 0.8781
                assert max_ani == 0.886


def test_simple_dayhoff(runtmp):
    # test basic execution with dayhoff sigs
    sigs = get_test_data('dayhoff.zip')

    output = runtmp.output('out.csv')

    runtmp.sourmash('scripts', 'multisearch', sigs, sigs,
                    '-o', output, '--moltype', 'dayhoff',
                    '-k', '19', '--scaled', '100', '--ani')
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 4

    dd = df.to_dict(orient='index')
    print(dd)

    for idx, row in dd.items():
        # identical?
        if row['match_name'] == row['query_name']:
            assert row['query_md5'] == row['match_md5'], row
            assert float(row['containment'] == 1.0)
            assert float(row['jaccard'] == 1.0)
            assert float(row['max_containment'] == 1.0)
            assert float(row['query_containment_ani'] == 1.0)
            assert float(row['match_containment_ani'] == 1.0)
            assert float(row['average_containment_ani'] == 1.0)
            assert float(row['max_containment_ani'] == 1.0)

        else:
            # confirm hand-checked numbers
            q = row['query_name'].split()[0]
            m = row['match_name'].split()[0]
            cont = float(row['containment'])
            jaccard = float(row['jaccard'])
            maxcont = float(row['max_containment'])
            intersect_hashes = int(row['intersect_hashes'])
            q1_ani = float(row['query_containment_ani'])
            q2_ani = float(row['match_containment_ani'])
            avg_ani = float(row['average_containment_ani'])
            max_ani = float(row['max_containment_ani'])

            jaccard = round(jaccard, 4)
            cont = round(cont, 4)
            maxcont = round(maxcont, 4)
            q1_ani = round(q1_ani, 4)
            q2_ani = round(q2_ani, 4)
            avg_ani = round(avg_ani, 4)
            max_ani = round(max_ani, 4)
            print(q, m, f"{jaccard:.04}", f"{cont:.04}", f"{maxcont:.04}", intersect_hashes, f"{q1_ani:.04}", f"{q2_ani:.04}", f"{avg_ani:.04}", f"{max_ani:.04}")

            if q == 'GCA_001593925' and m == 'GCA_001593935':
                assert jaccard == 0.1326
                assert cont == 0.2815
                assert maxcont == 0.2815
                assert intersect_hashes == 930
                assert q1_ani == 0.9355
                assert q2_ani == 0.9189
                assert avg_ani == 0.9272
                assert max_ani == 0.9355

            if q == 'GCA_001593935' and m == 'GCA_001593925':
                assert jaccard == 0.1326
                assert cont == 0.2004
                assert maxcont == 0.2815
                assert intersect_hashes == 930
                assert q1_ani == 0.9189
                assert q2_ani == 0.9355
                assert avg_ani == 0.9272
                assert max_ani == 0.9355


def test_simple_hp(runtmp):
    # test basic execution with hp sigs
    sigs = get_test_data('hp.zip')

    output = runtmp.output('out.csv')

    runtmp.sourmash('scripts', 'multisearch', sigs, sigs,
                    '-o', output, '--moltype', 'hp',
                    '-k', '19', '--scaled', '100', '--ani')
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 4

    dd = df.to_dict(orient='index')
    print(dd)

    for idx, row in dd.items():
        # identical?
        if row['match_name'] == row['query_name']:
            assert row['query_md5'] == row['match_md5'], row
            assert float(row['containment'] == 1.0)
            assert float(row['jaccard'] == 1.0)
            assert float(row['max_containment'] == 1.0)
            assert float(row['query_containment_ani'] == 1.0)
            assert float(row['match_containment_ani'] == 1.0)
            assert float(row['average_containment_ani'] == 1.0)
            assert float(row['max_containment_ani'] == 1.0)

        else:
            # confirm hand-checked numbers
            q = row['query_name'].split()[0]
            m = row['match_name'].split()[0]
            cont = float(row['containment'])
            jaccard = float(row['jaccard'])
            maxcont = float(row['max_containment'])
            intersect_hashes = int(row['intersect_hashes'])
            q1_ani = float(row['query_containment_ani'])
            q2_ani = float(row['match_containment_ani'])
            avg_ani = float(row['average_containment_ani'])
            max_ani = float(row['max_containment_ani'])

            jaccard = round(jaccard, 4)
            cont = round(cont, 4)
            maxcont = round(maxcont, 4)
            q1_ani = round(q1_ani, 4)
            q2_ani = round(q2_ani, 4)
            avg_ani = round(avg_ani, 4)
            max_ani = round(max_ani, 4)
            print(q, m, f"{jaccard:.04}", f"{cont:.04}", f"{maxcont:.04}", intersect_hashes, f"{q1_ani:.04}", f"{q2_ani:.04}", f"{avg_ani:.04}", f"{max_ani:.04}")

            if q == 'GCA_001593925' and m == 'GCA_001593935':
                assert jaccard == 0.4983
                assert cont == 0.747
                assert maxcont == 0.747
                assert intersect_hashes == 1724
                assert q1_ani == 0.9848
                assert q2_ani == 0.9734
                assert avg_ani == 0.9791
                assert max_ani == 0.9848

            if q == 'GCA_001593935' and m == 'GCA_001593925':
                assert jaccard == 0.4983
                assert cont == 0.5994
                assert maxcont == 0.747
                assert intersect_hashes == 1724
                assert q1_ani == 0.9734
                assert q2_ani == 0.9848
                assert avg_ani == 0.9791
                assert max_ani == 0.9848


def test_simple_below_threshold(runtmp):
    # test basic execution!
    query_list = runtmp.output('query.txt')
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])
    make_file_list(against_list, [sig2, sig47, sig63])

    output = runtmp.output('out.csv')

    runtmp.sourmash('scripts', 'multisearch', query_list, against_list,
                    '-o', output, '--ani', '--threshold', '0.5')
    assert os.path.exists(output)

    with open(output, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)
        assert len(rows) == 3
        for row in rows:
            # only identical reported
            print(row)
            assert row['query_md5'] == row['match_md5']
            assert float(row['containment']) == 1.0
            assert float(row['jaccard']) == 1.0
            assert float(row['max_containment']) == 1.0
            assert float(row['query_containment_ani']) == 1.0
            assert float(row['match_containment_ani']) == 1.0
            assert float(row['average_containment_ani']) == 1.0
            assert float(row['max_containment_ani']) == 1.0
