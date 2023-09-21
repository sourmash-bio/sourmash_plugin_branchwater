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
    names = [os.path.basename(x).split('.fa')[0] for x in genome_paths]
    # Check if the number of protein paths is less than genome paths
    # and fill in the missing paths with "".
    if len(protein_paths) < len(genome_paths):
        protein_paths.extend(["" for _ in range(len(genome_paths) - len(protein_paths))])
    with open(filename, 'wt') as fp:
        fp.write("name,genome_filename,protein_filename\n")
        for name, genome_path, protein_path in zip(names, genome_paths, protein_paths):
            fp.write("{},{},{}\n".format(name, genome_path, protein_path))


def test_installed(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch')

    assert 'usage:  manysketch' in runtmp.last_result.err


def test_manysketch_simple(runtmp):
    fa_csv = runtmp.output('db-fa.txt')

    fa1 = get_test_data('short.fa')
    fa2 = get_test_data('short2.fa')
    fa3 = get_test_data('short3.fa')

    make_file_csv(fa_csv, [fa1, fa2, fa3])

    output = runtmp.output('db.zip')

    runtmp.sourmash('scripts', 'manysketch', fa_csv, '-o', output,
                    '--param-str', "dna,k=31,scaled=1")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 3


def test_manysketch_mult_k(runtmp):
    fa_csv = runtmp.output('db-fa.txt')

    fa1 = get_test_data('short.fa')
    fa2 = get_test_data('short2.fa')
    fa3 = get_test_data('short3.fa')

    make_file_csv(fa_csv, [fa1, fa2, fa3])

    output = runtmp.output('db.zip')

    runtmp.sourmash('scripts', 'manysketch', fa_csv, '-o', output,
                    '--param-str', "dna,k=21,k=31,scaled=1")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 6


def test_manysketch_mult_k_2(runtmp):
    fa_csv = runtmp.output('db-fa.txt')

    fa1 = get_test_data('short.fa')
    fa2 = get_test_data('short2.fa')
    fa3 = get_test_data('short3.fa')

    make_file_csv(fa_csv, [fa1, fa2, fa3])

    output = runtmp.output('db.zip')

    runtmp.sourmash('scripts', 'manysketch', fa_csv, '-o', output,
                    '--param-str', "dna,k=21,scaled=1",
                    '--param-str', "dna,k=31,scaled=1",
                    '--param-str', "dna,k=21,scaled=1")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 6


def test_manysketch_mult_moltype(runtmp):
    fa_csv = runtmp.output('db-fa.csv')

    fa1 = get_test_data('short.fa')
    fa2 = get_test_data('short2.fa')
    fa3 = get_test_data('short3.fa')
    protfa1 = get_test_data('short-protein.fa')

    make_file_csv(fa_csv, [fa1, fa2, fa3], [protfa1])

    output = runtmp.output('db.zip')

    runtmp.sourmash('scripts', 'manysketch', fa_csv, '-o', output,
                    '--param-str', "dna,k=21,scaled=1",
                    '--param-str', "protein,k=10,scaled=1")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 4
    # check moltypes, etc!
    for sig in sigs:
        if sig.name == 'short':
            if sig.minhash.is_dna:
                assert sig.minhash.ksize == 21
                assert sig.minhash.scaled == 1
                assert sig.md5sum() == "1474578c5c46dd09da4c2df29cf86621"
            else:
                assert sig.minhash.ksize == 10
                assert sig.minhash.scaled == 1
                assert sig.md5sum() == "eb4467d11e0ecd2dbde4193bfc255310"
        else:
            assert sig.minhash.ksize == 21
            assert sig.minhash.scaled == 1
            assert sig.minhash.is_dna
            assert sig.md5sum() in ["4efeebd26644278e36b9553e018a851a","f85747ac4f473c4a71c1740d009f512b"]


def test_manysketch_only_incompatible_fastas(runtmp, capfd):
    # provide dna, protein fastas, but only sketch protein (skip protein fastas!)
    fa_csv = runtmp.output('db-fa.csv')

    fa1 = get_test_data('short.fa')
    fa2 = get_test_data('short2.fa')
    fa3 = get_test_data('short3.fa')

    make_file_csv(fa_csv, [fa1, fa2, fa3])

    output = runtmp.output('db.zip')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch', fa_csv, '-o', output,
                        '--param-str', "protein,k=10,scaled=1")

    assert os.path.exists(output) # output will still exist - is this desired?
    assert not runtmp.last_result.out # stdout should be empty

    captured = capfd.readouterr()
    print(captured.err)

    assert 'DONE. Processed 3 fasta files' in captured.err
    assert 'Error: No fasta files compatible with provided sketch parameters: no signatures created.' in captured.err


def test_manysketch_skip_incompatible_fastas(runtmp, capfd):
    # provide dna, protein fastas, but only sketch protein (skip protein fastas!)
    fa_csv = runtmp.output('db-fa.csv')

    fa1 = get_test_data('short.fa')
    fa2 = get_test_data('short2.fa')
    fa3 = get_test_data('short3.fa')
    protfa1 = get_test_data('short-protein.fa')

    make_file_csv(fa_csv, [fa1, fa2, fa3], [protfa1])

    output = runtmp.output('db.zip')

    runtmp.sourmash('scripts', 'manysketch', fa_csv, '-o', output,
                    '--param-str', "protein,k=10,scaled=1")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    captured = capfd.readouterr()
    print(captured.err)

    assert len(sigs) == 1
    # check moltypes, etc!
    for sig in sigs:
        assert sig.minhash.ksize == 10
        assert sig.minhash.scaled == 1
        assert sig.md5sum() == "eb4467d11e0ecd2dbde4193bfc255310"
    assert 'Starting file 2/4 (50%)' in captured.err
    assert 'Starting file 3/4 (75%)' in captured.err
    assert 'Starting file 4/4 (100%)' in captured.err
    assert 'DONE. Processed 4 fasta files' in captured.err
    assert 'WARNING: 3 fasta files skipped - no compatible signatures.' in captured.err


def test_manysketch_missing_fa_csv(runtmp, capfd):
    # test missing fa_csv file
    fa_csv = runtmp.output('fa_csv.txt')
    output = runtmp.output('out.zip')
    # make_file_list(fa_csv, []) # don't make fa_csv file

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch', fa_csv,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "Could not load fromfile csv" in captured.err


def test_manysketch_bad_fa_csv(runtmp, capfd):
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


def test_manysketch_bad_fa_csv_2(runtmp, capfd):
    # test sketch with fasta provided instead of fa_csv
    output = runtmp.output('out.zip')
    fa1 = get_test_data('short.fa')
    print(fa1)

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch', fa1,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "Invalid header" in captured.err
    assert "Could not load fromfile csv" in captured.err


def test_manysketch_bad_fa_csv_3(runtmp, capfd):
    # test sketch with improperly formatted fa_csv
    fa_csv = runtmp.output('db-fa.csv')

    fa1 = get_test_data('short.fa')
    fa2 = get_test_data('short2.fa')
    fa3 = get_test_data('short3.fa')
    protfa1 = get_test_data('short-protein.fa')

    # make file csv but don't fill empty protein rows with ,""
    make_file_csv(fa_csv, [fa1, fa2, fa3], [protfa1])
    g_fa = [fa1, fa2, fa3]
    p_fa = [protfa1]
    with open(fa_csv, 'wt') as fp:
        fp.write("name,genome_filename,protein_filename\n")
        for i, g in enumerate(g_fa):
            name = os.path.basename(g).split('.fa')[0]
            if i < len(p_fa):
                p = p_fa[i]
                fp.write("{},{},{}\n".format(name, g, p))
            else:
                fp.write("{},{}\n".format(name, g)) # missing prot path, no trailing comma

    output = runtmp.output('db.zip')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch', fa_csv,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    assert 'found record with 2 fields' in captured.err
    assert "Could not load fromfile csv" in captured.err
 

def test_manysketch_empty_fa_csv(runtmp, capfd):
    # test empty fa_csv file
    fa_csv = runtmp.output('fa.txt')
    output = runtmp.output('out.zip')
    make_file_csv(fa_csv, []) # empty

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysketch', fa_csv,
                        '-o', output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: No files to load, exiting." in captured.err


def test_manysketch_duplicated_rows(runtmp, capfd):
    fa_csv = runtmp.output('db-fa.csv')

    fa1 = get_test_data('short.fa')
    fa2 = get_test_data('short2.fa')
    fa3 = get_test_data('short3.fa')
    protfa1 = get_test_data('short-protein.fa')

    make_file_csv(fa_csv, [fa1, fa1, fa1, fa3])

    output = runtmp.output('db.zip')

    runtmp.sourmash('scripts', 'manysketch', fa_csv, '-o', output,
                    '--param-str', "dna,k=21,scaled=1",
                    '--param-str', "protein,k=10,scaled=1")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    captured = capfd.readouterr()
    print(captured.err)
    assert len(sigs) == 2
    assert "DONE. Processed 2 fasta files" in captured.err


def test_manysketch_N_in_dna(runtmp):
    # make sure we can handle Ns in DNA sequences
    fa_csv = runtmp.output('db-fa.txt')
    fa1 = runtmp.output('bad.fa')
    with open (fa1, 'wt') as fp:
        fp.write(">bad\n")
        fp.write("ACAGTN\n")

    make_file_csv(fa_csv, [fa1])

    output = runtmp.output('db.zip')

    runtmp.sourmash('scripts', 'manysketch', fa_csv, '-o', output,
                    '--param-str', "dna,k=4,scaled=1")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 1


def test_zip_manifest(runtmp, capfd):
    # test basic manifest-generating functionality.
    fa_csv = runtmp.output('db-fa.txt')

    fa1 = get_test_data('short.fa')
    fa2 = get_test_data('short2.fa')
    fa3 = get_test_data('short3.fa')

    make_file_csv(fa_csv, [fa1, fa2, fa3])
    output = runtmp.output('db.zip')

    runtmp.sourmash('scripts', 'manysketch', fa_csv, '-o', output,
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
