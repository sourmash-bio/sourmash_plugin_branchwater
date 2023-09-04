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
        runtmp.sourmash('scripts', 'manysketch')

    assert 'usage:  manysketch' in runtmp.last_result.err


def test_manysketch(runtmp):
    falist = runtmp.output('db-fa.txt')

    fa1 = get_test_data('short.fa')
    fa2 = get_test_data('short2.fa')
    fa3 = get_test_data('short3.fa')

    make_file_list(falist, [fa1, fa2, fa3])

    output = runtmp.output('db.zip')

    runtmp.sourmash('scripts', 'manysketch', falist, '-o', output)

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    sigs = list(sourmash.load_file_as_signatures(output))
    assert len(sigs) == 3


def test_sketch_missing_falist(runtmp, capfd):
    # test missing falist file
    falist = runtmp.output('falist.txt')
    output = runtmp.output('out.zip')
    # make_file_list(falist, []) # don't make falist file

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch', falist,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    assert 'Error: No such file or directory ' in captured.err


def test_manysketch_bad_falist(runtmp, capfd):
    # test sketch with fasta provided instead of falist
    fa1 = get_test_data('short.fa')
    output = runtmp.output('db.zip')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch', fa1, '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "ERROR: Could not load fasta files: no signatures created." in captured.err


def test_manysketch_bad_falist_2(runtmp, capfd):
    # test index with a bad siglist (.sig.gz file instead of fasta)
    siglist = runtmp.output('db-sigs.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(siglist, [sig2, sig47, sig63])

    output = runtmp.output('db.zip')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch', siglist, '-o', output) 

    captured = capfd.readouterr()
    print(captured.err)
    assert "ERROR: Could not load fasta files: no signatures created." in captured.err


def test_manysketch_empty_falist(runtmp, capfd):
    # test empty falist file
    falist = runtmp.output('fa.txt')
    output = runtmp.output('out.zip')
    make_file_list(falist, []) # empty

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch', falist,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "ERROR: Could not load fasta files: no signatures created." in captured.err

