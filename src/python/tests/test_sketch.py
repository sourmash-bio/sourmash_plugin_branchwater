import os
import pytest
import pandas
import sourmash
from sourmash import index

from . import sourmash_tst_utils as utils


def get_test_data(filename):
    thisdir = os.path.dirname(__file__)
    return os.path.join(thisdir, 'test-data', filename)


def make_file_csv(filename, genome_paths, protein_paths = []):
    names = [x.split('.fa')[0] for x in genome_paths]
    if not protein_paths:
        protein_paths = ["" for x in genome_paths]
    with open(filename, 'wt') as fp:
        fp.write("name,genome_filename,protein_filename\n")
        for name, genome_path, protein_path in zip(names, genome_paths, protein_paths):

            fp.write("{},{},{}\n".format(name, genome_path, protein_path))


def test_installed(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch')

    assert 'usage:  manysketch' in runtmp.last_result.err


def test_manysketch(runtmp):
    falist = runtmp.output('db-fa.txt')

    fa1 = get_test_data('short.fa')
    fa2 = get_test_data('short2.fa')
    fa3 = get_test_data('short3.fa')

    make_file_csv(falist, [fa1, fa2, fa3])

    output = runtmp.output('db.zip')

    runtmp.sourmash('scripts', 'manysketch', falist, '-o', output,
                    '--param-str', "dna,k=31,scaled=1")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 3


def test_manysketch_mult_k(runtmp):
    falist = runtmp.output('db-fa.txt')

    fa1 = get_test_data('short.fa')
    fa2 = get_test_data('short2.fa')
    fa3 = get_test_data('short3.fa')

    make_file_csv(falist, [fa1, fa2, fa3])

    output = runtmp.output('db.zip')

    runtmp.sourmash('scripts', 'manysketch', falist, '-o', output,
                    '--param-str', "dna,k=21,k=31,scaled=1")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 6


def test_manysketch_mult_k_2(runtmp):
    falist = runtmp.output('db-fa.txt')

    fa1 = get_test_data('short.fa')
    fa2 = get_test_data('short2.fa')
    fa3 = get_test_data('short3.fa')

    make_file_csv(falist, [fa1, fa2, fa3])

    output = runtmp.output('db.zip')

    runtmp.sourmash('scripts', 'manysketch', falist, '-o', output,
                    '--param-str', "dna,k=21,scaled=1",
                    '--param-str', "dna,k=31,scaled=1",
                    '--param-str', "dna,k=21,scaled=1")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 6


def test_manysketch_missing_falist(runtmp, capfd):
    # test missing falist file
    falist = runtmp.output('falist.txt')
    output = runtmp.output('out.zip')
    # make_file_list(falist, []) # don't make falist file

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch', falist,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "Could not load fromfile csv" in captured.err


def test_manysketch_bad_falist(runtmp, capfd):
    # siglist instead of fastalist
    siglist = runtmp.output('db-sigs.txt')

    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_csv(siglist, [sig2, sig47, sig63])

    output = runtmp.output('db.zip')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch', siglist, '-o', output) 

    captured = capfd.readouterr()
    print(captured.err)
    assert "Could not load fasta files: no signatures created." in captured.err


def test_manysketch_bad_falist_2(runtmp, capfd):
    # test sketch with fasta provided instead of falist
    output = runtmp.output('out.zip')
    fa1 = get_test_data('short.fa')
    print(fa1)

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch', fa1,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "Could not load fromfile csv" in captured.err


def test_manysketch_empty_falist(runtmp, capfd):
    # test empty falist file
    falist = runtmp.output('fa.txt')
    output = runtmp.output('out.zip')
    make_file_csv(falist, []) # empty

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch', falist,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: No files to load, exiting." in captured.err


def test_zip_manifest(runtmp, capfd):
    # test basic manifest-generating functionality.
    falist = runtmp.output('db-fa.txt')

    fa1 = get_test_data('short.fa')
    fa2 = get_test_data('short2.fa')
    fa3 = get_test_data('short3.fa')

    make_file_csv(falist, [fa1, fa2, fa3])
    output = runtmp.output('db.zip')

    runtmp.sourmash('scripts', 'manysketch', falist, '-o', output,
                    '--param-str', "dna,k=31,scaled=1")

    loader = sourmash.load_file_as_index(output)

    rows = []
    siglist = []
    for (sig, loc) in loader._signatures_with_internal():
        row = index.CollectionManifest.make_manifest_row(sig, loc)
        rows.append(row)
        siglist.append(sig)

    manifest = index.CollectionManifest(rows)

    assert len(manifest) == len(rows)
    assert len(manifest) == 3

    md5_list = [ row['md5'] for row in manifest.rows ]
    assert '9191284a3a23a913d8d410f3d53ce8f0' in md5_list
    assert 'd663bb55b2a0f8782c53c8af89f20fff' in md5_list
    assert 'bf752903d635b1eb83c53fe4aae951db' in md5_list

    for sig in siglist:
        assert sig in manifest
        assert sig.minhash.ksize == 31
        assert sig.minhash.moltype == 'DNA'
        assert sig.minhash.scaled == 1
