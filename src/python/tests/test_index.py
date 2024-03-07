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
        runtmp.sourmash('scripts', 'index')

    assert 'usage:  index' in runtmp.last_result.err


def test_index(runtmp):
    # test basic index!
    siglist = runtmp.output('db-sigs.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(siglist, [sig2, sig47, sig63])

    output = runtmp.output('db.rdb')

    runtmp.sourmash('scripts', 'index', siglist,
                    '-o', output)
    assert os.path.exists(output)
    print(runtmp.last_result.err)

    assert 'index is done' in runtmp.last_result.err


def test_index_protein(runtmp):
    sigs = get_test_data('protein.zip')
    output = runtmp.output('db.rocksdb')

    runtmp.sourmash('scripts', 'index', sigs, '-k', '19', '-s', '100',
                    '--moltype', 'protein', '-o', output)
    assert os.path.exists(output)
    print(runtmp.last_result.err)
    assert 'index is done' in runtmp.last_result.err


def test_index_dayhoff(runtmp):
    sigs = get_test_data('dayhoff.zip')
    output = runtmp.output('db.rocksdb')

    runtmp.sourmash('scripts', 'index', sigs, '-k', '19', '-s', '100',
                    '--moltype', 'dayhoff', '-o', output)
    assert os.path.exists(output)
    print(runtmp.last_result.err)
    assert 'index is done' in runtmp.last_result.err


def test_index_protein(runtmp):
    sigs = get_test_data('hp.zip')
    output = runtmp.output('db.rocksdb')

    runtmp.sourmash('scripts', 'index', sigs, '-k', '19', '-s', '100',
                    '--moltype', 'hp', '-o', output)
    assert os.path.exists(output)
    print(runtmp.last_result.err)
    assert 'index is done' in runtmp.last_result.err


def test_index_missing_siglist(runtmp, capfd):
    # test missing siglist file
    siglist = runtmp.output('db-sigs.txt')
    output = runtmp.output('out.db')
    # make_file_list(siglist, []) # don't make siglist file

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'index', siglist,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    assert 'Error: No such file or directory' in captured.err


def test_index_sig(runtmp, capfd):
    # test index with a .sig.gz file instead of pathlist
    # (should work now)
    sig2 = get_test_data('2.fa.sig.gz')
    output = runtmp.output('out.db')

    runtmp.sourmash('scripts', 'index', sig2,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    print(runtmp.last_result.err)
    assert 'index is done' in runtmp.last_result.err


def test_index_manifest(runtmp, capfd):
    # test index with a manifest file
    sig2 = get_test_data('2.fa.sig.gz')
    output = runtmp.output('out.db')
    sig_mf = runtmp.output('mf.csv')
    runtmp.sourmash("sig", "manifest", sig2, "-o", sig_mf)

    runtmp.sourmash('scripts', 'index', sig_mf,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    print(runtmp.last_result.err)
    assert 'index is done' in runtmp.last_result.err


def test_index_bad_siglist_2(runtmp, capfd):
    # test with a bad siglist (containing a missing file)
    against_list = runtmp.output('against.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')
    make_file_list(against_list, [sig2, "no-exist"])

    db = runtmp.output('db.rdb')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'index', against_list,
                        '-o', db)

    captured = capfd.readouterr()
    print(captured.err)
    assert "WARNING: could not load sketches from path 'no-exist'" in captured.err


def test_index_empty_siglist(runtmp, capfd):
    # test empty siglist file
    siglist = runtmp.output('db-sigs.txt')
    output = runtmp.output('out.db')
    make_file_list(siglist, []) # empty

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'index', siglist,
                        '-o', output)

    captured = capfd.readouterr()
    assert not os.path.exists(output) # do we want an empty file, or no file?
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    print(captured.err)
    assert "Error: Signatures failed to load. Exiting." in captured.err


def test_index_nomatch(runtmp, capfd):
    # test index with a siglist file that has (only) a non-matching ksize sig
    siglist = runtmp.output('against.txt')
    db = runtmp.output('db.rdb')

    sig1 = get_test_data('1.fa.k21.sig.gz')
    make_file_list(siglist, [sig1])

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'index', siglist,
                        '-o', db)

    captured = capfd.readouterr()
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    print(captured.out)
    print(captured.err)
    assert not os.path.exists(db)


def test_index_nomatch_sig_in_siglist(runtmp, capfd):
    # test index with a siglist file that has both matching and non-matching sigs
    siglist = runtmp.output('against.txt')
    db = runtmp.output('db.rdb')

    sig2 = get_test_data('2.fa.sig.gz')
    sig1 = get_test_data('1.fa.k21.sig.gz')
    make_file_list(siglist, [sig2, sig1])

    runtmp.sourmash('scripts', 'index', siglist,
                        '-o', db)

    captured = capfd.readouterr()
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    print(captured.out)
    print(captured.err)
    assert os.path.exists(db)


def test_index_zipfile(runtmp, capfd):
    # test basic index from sourmash zipfile
    siglist = runtmp.output('db-sigs.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(siglist, [sig2, sig47, sig63])

    zipf = runtmp.output('sigs.zip')

    runtmp.sourmash('sig', 'cat', siglist, '-o', zipf)

    output = runtmp.output('db.rdb')

    runtmp.sourmash('scripts', 'index', zipf,
                    '-o', output)
    assert os.path.exists(output)
    print(runtmp.last_result.err)

    assert 'index is done' in runtmp.last_result.err
    captured = capfd.readouterr()
    print(captured.err)


def test_index_zipfile_repeated_md5sums(runtmp, capfd):
    # test that we're reading all files, including repeated md5sums
    siglist = runtmp.output('db-sigs.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig2a = runtmp.output('sig2a.sig.gz')
    sig2b = runtmp.output('sig2b.sig.gz')
    runtmp.sourmash('sig', 'rename', sig2, 'name2', '-o', sig2a)
    runtmp.sourmash('sig', 'rename', sig2, 'name3', '-o', sig2b)

    make_file_list(siglist, [sig2, sig2a, sig2b])

    zipf = runtmp.output('sigs.zip')
    runtmp.sourmash('sig', 'cat', siglist, '-o', zipf)

    output = runtmp.output('db.rdb')

    runtmp.sourmash('scripts', 'index', zipf,
                    '-o', output)
    assert os.path.exists(output)
    print(runtmp.last_result.err)

    captured = capfd.readouterr()
    print(captured.err)

    assert 'index is done' in runtmp.last_result.err


def test_index_zipfile_multiparam(runtmp, capfd):
    # test index from sourmash zipfile with multiple ksizes / scaled /moltype
    siglist = runtmp.output('db-sigs.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')
    sig1 = get_test_data('1.combined.sig.gz')
    srr = get_test_data('SRR606249.sig.gz')
    prot = get_test_data('protein.zip')

    make_file_list(siglist, [sig2, sig47, sig63, sig1, srr, prot])

    zipf = runtmp.output('sigs.zip')

    runtmp.sourmash('sig', 'cat', siglist, '-o', zipf)

    output = runtmp.output('db.rdb')

    runtmp.sourmash('scripts', 'index', zipf,
                    '-o', output)
    assert os.path.exists(output)
    print(runtmp.last_result.err)

    assert 'index is done' in runtmp.last_result.err
    captured = capfd.readouterr()
    print(captured.err)


def test_index_zipfile_bad(runtmp, capfd):
    # test with a bad input zipfile (a .sig.gz file renamed as zip file)
    sig2 = get_test_data('2.fa.sig.gz')

    query_zip = runtmp.output('query.zip')
    # cp sig2 into query_zip
    with open(query_zip, 'wb') as fp:
        with open(sig2, 'rb') as fp2:
            fp.write(fp2.read())

    output = runtmp.output('out.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'index', query_zip,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)

    assert "Couldn't find End Of Central Directory Record" in captured.err


def test_index_check(runtmp):
    # test check index
    siglist = runtmp.output('db-sigs.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')

    make_file_list(siglist, [sig2, sig47])

    output = runtmp.output('db.rdb')

    runtmp.sourmash('scripts', 'index', siglist,
                    '-o', output)

    runtmp.sourmash('scripts', 'check', output)
    print(runtmp.last_result.err)

    assert 'index is ok' in runtmp.last_result.err


def test_index_check_quick(runtmp):
    # test check index
    siglist = runtmp.output('db-sigs.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')

    make_file_list(siglist, [sig2, sig47])

    output = runtmp.output('db.rdb')

    runtmp.sourmash('scripts', 'index', siglist,
                    '-o', output)

    runtmp.sourmash('scripts', 'check', '--quick', output)
    print(runtmp.last_result.err)

    assert 'index is ok' in runtmp.last_result.err
