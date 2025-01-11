import os
import pytest
import csv
import pandas
import sourmash
import subprocess
from sourmash import index, sourmash_args
import io
from . import sourmash_tst_utils as utils


def get_test_data(filename):
    thisdir = os.path.dirname(__file__)
    return os.path.join(thisdir, "test-data", filename)


def make_assembly_csv(filename, genome_paths, protein_paths=[]):
    # equalize path lengths by adding "".
    names = [os.path.basename(x).split(".fa")[0] for x in genome_paths]
    if len(protein_paths) < len(genome_paths):
        protein_paths.extend(
            ["" for _ in range(len(genome_paths) - len(protein_paths))]
        )
    elif len(genome_paths) < len(protein_paths):
        genome_paths.extend(["" for _ in range(len(protein_paths) - len(genome_paths))])
        names = [os.path.basename(x).split(".fa")[0] for x in protein_paths]

    with open(filename, "wt") as fp:
        fp.write("name,genome_filename,protein_filename\n")
        for name, genome_path, protein_path in zip(names, genome_paths, protein_paths):
            fp.write("{},{},{}\n".format(name, genome_path, protein_path))


def make_reads_csv(filename, reads_tuples=[]):
    # reads tuples should be (name,read1,read2)
    with open(filename, "wt") as fp:
        fp.write("name,read1,read2\n")
        for name, read1, read2 in reads_tuples:
            print(f"{name},{read1},{read2}")
            fp.write("{},{},{}\n".format(name, read1, read2))


def test_installed(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "manysketch")

    assert "usage:  manysketch" in runtmp.last_result.err


def test_manysketch_simple(runtmp):
    fa_csv = runtmp.output("db-fa.txt")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")

    make_assembly_csv(fa_csv, [fa1, fa2, fa3])

    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "dna,k=31,scaled=1",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 3


def test_manysketch_mult_k(runtmp):
    fa_csv = runtmp.output("db-fa.txt")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")

    make_assembly_csv(fa_csv, [fa1, fa2, fa3])

    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "dna,k=21,k=31,scaled=1",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 6

    names = [sig.name for sig in sigs]
    print(names)
    assert names.count("short") == 2
    assert names.count("short2") == 2
    assert names.count("short3") == 2


def test_manysketch_mult_k_2(runtmp):
    fa_csv = runtmp.output("db-fa.txt")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")

    make_assembly_csv(fa_csv, [fa1, fa2, fa3])

    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "dna,k=21,scaled=1",
        "--param-str",
        "dna,k=31,scaled=1",
        "--param-str",
        "dna,k=21,scaled=1",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 6

    names = [sig.name for sig in sigs]
    print(names)
    assert names.count("short") == 2
    assert names.count("short2") == 2
    assert names.count("short3") == 2


def test_manysketch_mult_moltype(runtmp):
    fa_csv = runtmp.output("db-fa.csv")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")
    protfa1 = get_test_data("short-protein.fa")

    make_assembly_csv(fa_csv, [fa1, fa2, fa3], [protfa1])

    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "dna,k=21,scaled=1",
        "--param-str",
        "protein,k=10,scaled=1",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 4
    # check moltypes, etc!
    for sig in sigs:
        if sig.name == "short":
            if sig.minhash.is_dna:
                assert sig.minhash.ksize == 21
                assert sig.minhash.scaled == 1
                assert sig.md5sum() == "1474578c5c46dd09da4c2df29cf86621"
            else:
                assert sig.name == "short"
                assert sig.minhash.ksize == 10
                assert sig.minhash.scaled == 1
                assert sig.md5sum() == "eb4467d11e0ecd2dbde4193bfc255310"
        else:
            assert sig.name in ["short", "short2", "short3"]
            assert sig.minhash.ksize == 21
            assert sig.minhash.scaled == 1
            assert sig.minhash.is_dna
            assert sig.md5sum() in [
                "4efeebd26644278e36b9553e018a851a",
                "f85747ac4f473c4a71c1740d009f512b",
            ]


def test_manysketch_mult_moltype_protein(runtmp):
    fa_csv = runtmp.output("db-fa.csv")

    protfa1 = get_test_data("short-protein.fa")

    make_assembly_csv(fa_csv, [], [protfa1])

    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "dayhoff,k=10,scaled=1",
        "--param-str",
        "hp,k=24,scaled=1",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 2
    # check moltypes, etc!
    total_checked = 0
    for sig in sigs:
        print(sig.name)
        assert sig.name == "short-protein"
        if sig.minhash.dayhoff:
            assert sig.md5sum() == "320464775fe704d9f938a8c63d8dd722"
            total_checked += 1
        elif sig.minhash.hp:
            assert sig.md5sum() == "e8ccc6ca7ad560072f51be631d1c39c0"
            total_checked += 1
    assert total_checked == 2


def test_manysketch_only_incompatible_fastas(runtmp, capfd):
    # provide dna, protein fastas, but only sketch protein (skip protein fastas!)
    fa_csv = runtmp.output("db-fa.csv")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")

    make_assembly_csv(fa_csv, [fa1, fa2, fa3])

    output = runtmp.output("db.zip")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts",
            "manysketch",
            fa_csv,
            "-o",
            output,
            "--param-str",
            "protein,k=10,scaled=1",
        )

    assert os.path.exists(output)  # output will still exist - is this desired?
    assert not runtmp.last_result.out  # stdout should be empty

    captured = capfd.readouterr()
    print(captured.err)

    assert "DONE. Processed 3 fasta files" in captured.err
    assert (
        "Error: No fasta files compatible with provided sketch parameters: no signatures created."
        in captured.err
    )


def test_manysketch_skip_incompatible_fastas(runtmp, capfd):
    # provide dna, protein fastas, but only sketch protein (skip protein fastas!)
    fa_csv = runtmp.output("db-fa.csv")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")
    protfa1 = get_test_data("short-protein.fa")

    make_assembly_csv(fa_csv, [fa1, fa2, fa3], [protfa1])

    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "protein,k=10,scaled=1",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty

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
    assert "DONE. Processed 4 fasta files" in captured.err
    assert "WARNING: 3 fasta files skipped - no compatible signatures." in captured.err


def test_manysketch_missing_fa_csv(runtmp, capfd):
    # test missing fa_csv file
    fa_csv = runtmp.output("fa_csv.txt")
    output = runtmp.output("out.zip")
    # make_file_list(fa_csv, []) # don't make fa_csv file

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "manysketch", fa_csv, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "Could not load fromfile csv" in captured.err


def test_manysketch_bad_fa_csv(runtmp, capfd):
    # siglist instead of fastalist
    siglist = runtmp.output("db-sigs.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_assembly_csv(siglist, [sig2, sig47, sig63])

    output = runtmp.output("db.zip")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "manysketch", siglist, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "Could not load fasta files: no signatures created." in captured.err


def test_manysketch_bad_fa_csv_2(runtmp, capfd):
    # bad file within filelist
    siglist = runtmp.output("bad.txt")

    # fa_file = runtmp.output("bad.fa")
    make_assembly_csv(siglist, ["bad2.fa"])

    output = runtmp.output("db.zip")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "manysketch", siglist, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "Could not load fasta files: no signatures created." in captured.err
    assert "Error building signatures from file: bad2.fa" in captured.err


def test_manysketch_bad_fa_csv_3(runtmp, capfd):
    # test sketch with fasta provided instead of fa_csv
    output = runtmp.output("out.zip")
    fa1 = get_test_data("short.fa")
    print(fa1)

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "manysketch", fa1, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "Invalid header" in captured.err
    assert "Could not load fromfile csv" in captured.err


def test_manysketch_bad_fa_csv_4(runtmp, capfd):
    # test sketch with improperly formatted fa_csv
    fa_csv = runtmp.output("db-fa.csv")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")
    protfa1 = get_test_data("short-protein.fa")

    # make file csv but don't fill empty protein rows with ,""
    make_assembly_csv(fa_csv, [fa1, fa2, fa3], [protfa1])
    g_fa = [fa1, fa2, fa3]
    p_fa = [protfa1]
    with open(fa_csv, "wt") as fp:
        fp.write("name,genome_filename,protein_filename\n")
        for i, g in enumerate(g_fa):
            name = os.path.basename(g).split(".fa")[0]
            if i < len(p_fa):
                p = p_fa[i]
                fp.write("{},{},{}\n".format(name, g, p))
            else:
                fp.write(
                    "{},{}\n".format(name, g)
                )  # missing prot path, no trailing comma

    output = runtmp.output("db.zip")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "manysketch", fa_csv, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "found record with 2 fields" in captured.err
    assert "Could not load fromfile csv" in captured.err


def test_manysketch_bad_param_str_moltype(runtmp, capfd):
    # no moltype provided in param str
    fa_csv = runtmp.output("db-fa.txt")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")

    make_assembly_csv(fa_csv, [fa1, fa2, fa3])
    output = runtmp.output("out.zip")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts", "manysketch", fa_csv, "-o", output, "-p", "k=31,scaled=100"
        )

    captured = capfd.readouterr()
    print(captured.err)
    assert (
        "Error parsing params string 'k=31,scaled=100': No moltype provided"
        in captured.err
    )
    assert "Failed to parse params string" in captured.err


def test_manysketch_empty_fa_csv(runtmp, capfd):
    # test empty fa_csv file
    fa_csv = runtmp.output("fa.txt")
    output = runtmp.output("out.zip")
    make_assembly_csv(fa_csv, [])  # empty

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "manysketch", fa_csv, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: No files to load, exiting." in captured.err


def test_manysketch_duplicated_rows(runtmp, capfd):
    fa_csv = runtmp.output("db-fa.csv")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")
    protfa1 = get_test_data("short-protein.fa")

    make_assembly_csv(fa_csv, [fa1, fa1, fa1, fa3])

    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "dna,k=21,scaled=1",
        "--param-str",
        "protein,k=10,scaled=1",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    captured = capfd.readouterr()
    print(captured.err)
    assert len(sigs) == 2
    assert "DONE. Processed 2 fasta files" in captured.err


def test_manysketch_N_in_dna(runtmp):
    # make sure we can handle Ns in DNA sequences
    fa_csv = runtmp.output("db-fa.txt")
    fa1 = runtmp.output("bad.fa")
    with open(fa1, "wt") as fp:
        fp.write(">bad\n")
        fp.write("ACAGTN\n")

    make_assembly_csv(fa_csv, [fa1])

    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts", "manysketch", fa_csv, "-o", output, "--param-str", "dna,k=4,scaled=1"
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 1


def test_zip_manifest(runtmp, capfd):
    # test basic manifest-generating functionality.
    fa_csv = runtmp.output("db-fa.txt")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")

    make_assembly_csv(fa_csv, [fa1, fa2, fa3])
    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "dna,k=31,scaled=1",
    )

    idx = sourmash.load_file_as_index(output)
    manifest = sourmash_args.get_manifest(idx)
    assert len(manifest) == 3
    md5_nhashes = [(row["md5"], row["n_hashes"]) for row in manifest.rows]
    assert ("9191284a3a23a913d8d410f3d53ce8f0", 970) in md5_nhashes
    assert ("d663bb55b2a0f8782c53c8af89f20fff", 925) in md5_nhashes
    assert ("bf752903d635b1eb83c53fe4aae951db", 955) in md5_nhashes


def test_protein_zip_manifest(runtmp, capfd):
    # test basic manifest-generating functionality.
    fa_csv = runtmp.output("db-fa.csv")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short-protein.fa")

    make_assembly_csv(fa_csv, [fa1], [fa2])
    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "protein,k=10,scaled=1",
    )

    idx = sourmash.load_file_as_index(output)
    manifest = sourmash_args.get_manifest(idx)
    assert len(manifest) == 1
    # for row in manifest.rows:
    row = manifest.rows[0]
    assert row["md5"] == "eb4467d11e0ecd2dbde4193bfc255310"
    assert row["ksize"] == 10
    assert row["scaled"] == 1
    assert row["moltype"] == "protein"
    assert row["n_hashes"] == 902


def test_manysketch_singleton(runtmp):
    fa_csv = runtmp.output("db-fa.txt")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")

    make_assembly_csv(fa_csv, [fa1, fa2, fa3])

    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "dna,k=31,scaled=1",
        "--singleton",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 4
    singleton_sketch = runtmp.output("short3.sig")
    runtmp.sourmash(
        "sketch",
        "dna",
        fa3,
        "-o",
        singleton_sketch,
        "--param-str",
        "dna,k=31,scaled=1",
        "--singleton",
    )
    ss_sketch = sourmash.load_signatures(singleton_sketch)
    ss_sketch1 = next(ss_sketch)
    ss_sketch2 = next(ss_sketch)

    expected_signames = ["shortName", "tr1 4", "firstname", "other"]
    for sig in sigs:
        assert sig.name in expected_signames
        if sig.name == "firstname":
            assert sig == ss_sketch1
        if sig.name == "other":
            assert sig == ss_sketch2


def test_manysketch_reads(runtmp, capfd):
    fa_csv = runtmp.output("db-fa.csv")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")

    make_reads_csv(
        fa_csv, [("short", fa1, fa2), ("short3", fa3, "")]
    )  # make sure we can just do read1 alone

    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "dna,k=31,scaled=1",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)
    assert "Found 'reads' CSV, assuming all files are DNA." in captured.out
    assert "Starting file 3/3 (100%)" in captured.err
    assert "DONE. Processed 3 fasta files" in captured.err

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 2
    s1 = runtmp.output("short.sig")
    runtmp.sourmash(
        "sketch",
        "dna",
        fa1,
        fa2,
        "-o",
        s1,
        "--param-str",
        "k=31,scaled=1",
        "--name",
        "short",
    )
    sig1 = sourmash.load_one_signature(s1)
    s3 = runtmp.output("short3.sig")
    runtmp.sourmash(
        "sketch",
        "dna",
        fa3,
        "-o",
        s3,
        "--param-str",
        "k=31,scaled=1",
        "--name",
        "short3",
    )
    sig2 = sourmash.load_one_signature(s3)

    expected_signames = ["short", "short3"]
    for sig in sigs:
        assert sig.name in expected_signames
        if sig.name == "short":
            assert sig == sig1
        if sig.name == "short3":
            assert sig == sig2


def test_manysketch_reads_singleton(runtmp, capfd):
    fa_csv = runtmp.output("db-fa.csv")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")

    make_reads_csv(
        fa_csv,
        [
            ("short", fa2, fa3),
        ],
    )

    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "dna,k=31,scaled=1",
        "--singleton",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)
    assert "Found 'reads' CSV, assuming all files are DNA." in captured.out
    assert "Starting file 2/2 (100%)" in captured.err
    assert "DONE. Processed 2 fasta files" in captured.err

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 3
    s1 = runtmp.output("singleton.sig")
    runtmp.sourmash(
        "sketch",
        "dna",
        fa2,
        fa3,
        "-o",
        s1,
        "--param-str",
        "k=31,scaled=1",
        "--singleton",
    )
    ss = sourmash.load_signatures(s1)

    ss_sketch1 = next(ss)
    ss_sketch2 = next(ss)
    ss_sketch3 = next(ss)

    expected_signames = ["tr1 4", "firstname", "other"]
    for sig in sigs:
        assert sig.name in expected_signames
        if sig.name == "tr1 4":
            assert sig == ss_sketch1
        elif sig.name == "firstname":
            assert sig == ss_sketch2
        elif sig.name == "other":
            assert sig == ss_sketch3


def test_manysketch_prefix(runtmp, capfd):
    fa_csv = runtmp.output("db-fa.csv")

    fa1 = get_test_data("short.fa")

    fa_path = os.path.dirname(fa1)
    dna_prefix = os.path.join(
        fa_path, "short*fa"
    )  # need to avoid matching short-protein.fa
    prot_prefix = os.path.join(fa_path, "*protein.fa")

    # make prefix input file
    with open(fa_csv, "wt") as fp:
        fp.write("name,input_moltype,prefix,exclude\n")
        fp.write(
            f"short,DNA,{dna_prefix},{prot_prefix}\n"
        )  # short.fa, short2.fa, short3.fa, short-protein.fa
        fp.write(f"short_protein,protein,{prot_prefix},\n")  # short-protein.fa only

    output = runtmp.output("prefix.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "dna,k=31,scaled=1",
        "-p",
        "protein,k=10,scaled=1",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)
    assert (
        "Found 'prefix' CSV. Using 'glob' to find files based on 'prefix' column."
        in captured.out
    )
    assert "DONE. Processed 4 fasta files" in captured.err

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 2

    # make same sigs with sourmash
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")
    fa4 = get_test_data("short-protein.fa")
    s1 = runtmp.output("short.sig")
    runtmp.sourmash(
        "sketch",
        "dna",
        fa1,
        fa2,
        fa3,
        "-o",
        s1,
        "--param-str",
        "dna,k=31,scaled=1",
        "--name",
        "short",
    )
    sig1 = sourmash.load_one_signature(s1)
    s2 = runtmp.output("short-protein.sig")
    runtmp.sourmash(
        "sketch",
        "protein",
        fa4,
        "-o",
        s2,
        "--param-str",
        "protein,k=10,scaled=1",
        "--name",
        "short_protein",
    )
    sig2 = sourmash.load_one_signature(s2)

    expected_signames = ["short", "short_protein"]
    for sig in sigs:
        assert sig.name in expected_signames
        if sig.name == "short":
            assert sig, minhash.hashes == sig1.minhash.hashes
        if sig.name == "short_protein":
            assert sig == sig2


def test_manysketch_prefix2(runtmp, capfd):
    fa_csv = runtmp.output("db-fa.csv")

    fa1 = get_test_data("short.fa")

    fa_path = os.path.dirname(fa1)
    # test without '*'
    dna_prefix = os.path.join(
        fa_path, "short"
    )  # need to avoid matching short-protein.fa
    prot_prefix = os.path.join(fa_path, "*protein.fa")
    zip_exclude = os.path.join(fa_path, "*zip")

    # make prefix input file
    with open(fa_csv, "wt") as fp:
        fp.write("name,input_moltype,prefix,exclude\n")
        fp.write(
            f"short,DNA,{dna_prefix},{prot_prefix}\n"
        )  # short.fa, short2.fa, short3.fa, short-protein.fa
        fp.write(
            f"short_protein,protein,{prot_prefix},{zip_exclude}\n"
        )  # short-protein.fa only

    output = runtmp.output("prefix.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "dna,k=31,scaled=1",
        "-p",
        "protein,k=10,scaled=1",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)
    assert (
        "Found 'prefix' CSV. Using 'glob' to find files based on 'prefix' column."
        in captured.out
    )
    assert "DONE. Processed 4 fasta files" in captured.err

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 2

    # make same sigs with sourmash
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")
    fa4 = get_test_data("short-protein.fa")
    s1 = runtmp.output("short.sig")
    runtmp.sourmash(
        "sketch",
        "dna",
        fa1,
        fa2,
        fa3,
        "-o",
        s1,
        "--param-str",
        "dna,k=31,scaled=1",
        "--name",
        "short",
    )
    sig1 = sourmash.load_one_signature(s1)
    s2 = runtmp.output("short-protein.sig")
    runtmp.sourmash(
        "sketch",
        "protein",
        fa4,
        "-o",
        s2,
        "--param-str",
        "protein,k=10,scaled=1",
        "--name",
        "short_protein",
    )
    sig2 = sourmash.load_one_signature(s2)

    expected_signames = ["short", "short_protein"]
    for sig in sigs:
        assert sig.name in expected_signames
        if sig.name == "short":
            # minhash is not defined? How does this test work? - @olgabot
            assert sig, minhash.hashes == sig1.minhash.hashes
        if sig.name == "short_protein":
            assert sig == sig2


def test_manysketch_prefix_duplicated_fail(runtmp, capfd):
    fa_csv = runtmp.output("db-fa.csv")

    fa1 = get_test_data("short.fa")

    fa_path = os.path.dirname(fa1)
    # test without '*'
    dna_prefix = os.path.join(
        fa_path, "short"
    )  # need to avoid matching short-protein.fa
    prot_prefix = os.path.join(fa_path, "*protein")
    zip_exclude = os.path.join(fa_path, "*zip")

    # make prefix input file
    with open(fa_csv, "wt") as fp:
        fp.write("name,input_moltype,prefix,exclude\n")
        fp.write(
            f"short,DNA,{dna_prefix},{prot_prefix}\n"
        )  # short.fa, short2.fa, short3.fa, short-protein.fa
        fp.write(
            f"short,DNA,{dna_prefix},{prot_prefix}\n"
        )  # duplicate of row one -- this should just be skipped
        fp.write(
            f"short_protein,protein,{prot_prefix},{zip_exclude}\n"
        )  # short-protein.fa only
        # ALSO short-protein.fa, but different name. should raise err without force
        fp.write(f"second_protein,protein,{prot_prefix},{zip_exclude}\n")

    output = runtmp.output("prefix.zip")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts",
            "manysketch",
            fa_csv,
            "-o",
            output,
            "--param-str",
            "dna,k=31,scaled=1",
            "-p",
            "protein,k=10,scaled=1",
        )

    assert not os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)
    assert (
        "Found 'prefix' CSV. Using 'glob' to find files based on 'prefix' column."
        in captured.out
    )
    assert "Found identical FASTA paths in more than one row!" in captured.err
    assert "Duplicated paths:" in captured.err
    assert "short-protein.fa" in captured.err
    assert (
        "Duplicated FASTA files found. Please use --force to bypass this check"
        in captured.err
    )


def test_manysketch_prefix_duplicated_force(runtmp, capfd):
    fa_csv = runtmp.output("db-fa.csv")

    fa1 = get_test_data("short.fa")

    fa_path = os.path.dirname(fa1)
    # test without '*'
    dna_prefix = os.path.join(
        fa_path, "short"
    )  # need to avoid matching short-protein.fa
    prot_prefix = os.path.join(fa_path, "*protein*fa")
    zip_exclude = os.path.join(fa_path, "*zip")

    # make prefix input file
    with open(fa_csv, "wt") as fp:
        fp.write("name,input_moltype,prefix,exclude\n")
        fp.write(
            f"short,DNA,{dna_prefix},{prot_prefix}\n"
        )  # short.fa, short2.fa, short3.fa, short-protein.fa
        fp.write(
            f"short,DNA,{dna_prefix},{prot_prefix}\n"
        )  # duplicate of row one -- this should just be skipped
        fp.write(
            f"short_protein,protein,{prot_prefix},{zip_exclude}\n"
        )  # short-protein.fa only
        # ALSO short-protein.fa, but different name. should raise err without force
        fp.write(f"second_protein,protein,{prot_prefix},{zip_exclude}\n")

    output = runtmp.output("prefix.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "dna,k=31,scaled=1",
        "-p",
        "protein,k=10,scaled=1",
        "--force",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)
    assert (
        "Found 'prefix' CSV. Using 'glob' to find files based on 'prefix' column."
        in captured.out
    )
    assert (
        "Loaded 3 rows in total (3 DNA FASTA and 2 protein FASTA), 1 duplicate rows skipped."
        in captured.out
    )
    assert "Found identical FASTA paths in more than one row!" in captured.err
    assert "Duplicated paths:" in captured.err
    assert "short-protein.fa" in captured.err
    assert "--force is set. Continuing..." in captured.err

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    assert len(sigs) == 3


def test_singlesketch_simple(runtmp):
    """Test basic single sketching with default parameters."""
    fa1 = get_test_data("short.fa")
    output = runtmp.output("short.sig")

    # Run the singlesketch command
    runtmp.sourmash("scripts", "singlesketch", fa1, "-o", output, "-p", "scaled=10")

    # Check if the output exists and contains the expected data
    assert os.path.exists(output)
    sig = sourmash.load_one_signature(output)

    assert sig.name == "short.fa"
    assert sig.minhash.ksize == 31
    assert sig.minhash.is_dna
    assert sig.minhash.scaled == 10
    print("HASHES", sig.minhash.hashes)

    # validate against sourmash sketch
    output2 = runtmp.output("short2.sig")
    runtmp.sourmash("sketch", "dna", fa1, "-o", output2, "-p", "scaled=10")
    sig2 = sourmash.load_one_signature(output2)
    assert sig.minhash.hashes == sig2.minhash.hashes


def test_singlesketch_with_name(runtmp):
    """Test single sketching with a custom name."""
    fa1 = get_test_data("short.fa")
    output = runtmp.output("short_named.sig")

    # Run the singlesketch command with the --name option
    runtmp.sourmash("scripts", "singlesketch", fa1, "-o", output, "-n", "custom_name")

    # Check if the output exists and contains the expected data
    assert os.path.exists(output)
    sig = sourmash.load_one_signature(output)

    assert sig.name == "custom_name"
    assert sig.minhash.ksize == 31
    assert sig.minhash.is_dna
    assert sig.minhash.scaled == 1000


def test_singlesketch_mult_k(runtmp):
    """Test single sketching with multiple k-mer sizes."""
    fa1 = get_test_data("short.fa")
    output = runtmp.output("short_mult_k.sig")

    # Run the singlesketch command with multiple k sizes
    runtmp.sourmash(
        "scripts",
        "singlesketch",
        fa1,
        "-o",
        output,
        "-p",
        "k=21,scaled=100",
        "-p",
        "k=31,scaled=100",
    )

    # Check if the output exists and contains the expected data
    assert os.path.exists(output)
    sigs = list(sourmash.load_signatures(output))

    # Verify that two signatures with different k-mer sizes exist
    assert len(sigs) == 2
    assert any(sig.minhash.ksize == 21 for sig in sigs)
    assert any(sig.minhash.ksize == 31 for sig in sigs)


def test_singlesketch_mult_k_2(runtmp):
    """Test single sketching with multiple k-mer sizes in one param string"""
    fa1 = get_test_data("short.fa")
    output = runtmp.output("short_mult_k.sig")

    # Run the singlesketch command with multiple k sizes
    runtmp.sourmash(
        "scripts",
        "singlesketch",
        fa1,
        "-o",
        output,
        "-p",
        "k=21,k=31,scaled=100",
    )

    # Check if the output exists and contains the expected data
    assert os.path.exists(output)
    sigs = list(sourmash.load_signatures(output))

    # Verify that two signatures with different k-mer sizes exist
    assert len(sigs) == 2
    assert any(sig.minhash.ksize == 21 for sig in sigs)
    assert any(sig.minhash.ksize == 31 for sig in sigs)


def test_singlesketch_explicit_dna(runtmp):
    """Test single sketching with explicit DNA in name"""
    fa1 = get_test_data("short.fa")
    output = runtmp.output("short_dna.sig")

    # Run the singlesketch command with multiple k sizes
    runtmp.sourmash(
        "scripts",
        "singlesketch",
        fa1,
        "-o",
        output,
        "-p",
        "k=21,k=31,scaled=100,dna",
    )

    # Check if the output exists and contains the expected data
    assert os.path.exists(output)
    sigs = list(sourmash.load_signatures(output))

    # Verify that two signatures with different k-mer sizes exist
    assert len(sigs) == 2
    assert any(sig.minhash.ksize == 21 for sig in sigs)
    assert any(sig.minhash.ksize == 31 for sig in sigs)


def test_singlesketch_protein_moltype(runtmp):
    """Test single sketching with different molecule types."""
    fa1 = get_test_data("short-protein.fa")
    output = runtmp.output("short_mult_moltype.sig")

    # Run the singlesketch command with prot molecule types
    runtmp.sourmash(
        "scripts",
        "singlesketch",
        fa1,
        "-o",
        output,
        "-p",
        "protein,k=10,scaled=100",
        "--input-moltype",
        "protein",
    )

    # Check if the output exists and contains the expected data
    assert os.path.exists(output)
    sig = sourmash.load_one_signature(output)

    # Verify that the signature has the correct molecule type and parameters
    assert sig.minhash.ksize == 10
    assert sig.minhash.is_protein
    assert sig.minhash.scaled == 100
    print("HASHES:", sig.minhash.hashes)

    # validate against sourmash sketch
    output2 = runtmp.output("short2.sig")
    runtmp.sourmash("sketch", "protein", fa1, "-p", "k=10,scaled=100", "-o", output2)
    sig2 = sourmash.load_one_signature(output2)
    assert sig.minhash.hashes == sig2.minhash.hashes


def test_singlesketch_invalid_params(runtmp, capfd):
    """Test singlesketch command with invalid parameters."""
    fa1 = get_test_data("short.fa")
    output = runtmp.output("short_invalid.sig")

    # Run the singlesketch command with an invalid parameter string
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts", "singlesketch", fa1, "-o", output, "-p", "invalid_param"
        )

    # Check that the error message is correct
    captured = capfd.readouterr()
    assert "Failed to parse params string" in captured.err


@pytest.mark.xfail(reason="needs to be implemented")
def test_singlesketch_translate(runtmp):
    """Test basic single sketching with input = DNA, output = protein"""
    fa1 = get_test_data("short.fa")
    output = runtmp.output("short.sig")

    # Run the singlesketch command
    runtmp.sourmash(
        "scripts",
        "singlesketch",
        fa1,
        "-o",
        output,
        "--input-moltype",
        "dna",
        "-p",
        "protein,k=7",
    )

    # Check if the output exists and contains the expected data
    assert os.path.exists(output)
    sig = sourmash.load_one_signature(output)

    assert sig.name == "short.fa"
    assert sig.minhash.ksize == 31
    assert sig.minhash.is_dna
    assert sig.minhash.scaled == 1000


@pytest.mark.xfail(reason="needs to be implemented")
def test_singlesketch_multimoltype_fail(runtmp):
    """Test failure with multiple moltype"""
    fa1 = get_test_data("short.fa")
    output = runtmp.output("short.sig")

    # Run the singlesketch command
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash(
            "scripts",
            "singlesketch",
            fa1,
            "-o",
            output,
            "--input-moltype",
            "dna",
            "-p",
            "protein,dna,k=7",
        )


def test_singlesketch_gzipped_output(runtmp):
    """Test singlesketch with gzipped output."""
    fa1 = get_test_data("short.fa")
    output = runtmp.output("short.sig.gz")

    # Run the singlesketch command
    runtmp.sourmash("scripts", "singlesketch", fa1, "-o", output, "-p", "scaled=10")

    # Check if the output exists and contains the expected data
    assert os.path.exists(output)

    # Verify the file is gzipped
    import gzip

    try:
        with gzip.open(output, "rt") as f:
            f.read(1)  # Try to read a single character to ensure it's valid gzip
    except gzip.BadGzipFile:
        assert False, f"Output file {output} is not a valid gzipped file."

    # check the signatures
    sig = sourmash.load_one_signature(output)
    print("HASHES:", sig.minhash.hashes)

    assert sig.name == "short.fa"
    assert sig.minhash.ksize == 31
    assert sig.minhash.is_dna
    assert sig.minhash.scaled == 10

    # validate against sourmash sketch
    output2 = runtmp.output("short2.sig")
    runtmp.sourmash("sketch", "dna", fa1, "-o", output2, "-p", "scaled=10")
    sig2 = sourmash.load_one_signature(output2)
    assert sig.minhash.hashes == sig2.minhash.hashes


def test_singlesketch_zip_output(runtmp):
    """Test singlesketch with zip output."""
    fa1 = get_test_data("short.fa")
    output = runtmp.output("short.zip")

    # Run the singlesketch command
    runtmp.sourmash("scripts", "singlesketch", fa1, "-o", output, "-p", "scaled=10")

    # Check if the output exists and contains the expected data
    assert os.path.exists(output)
    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    assert len(sigs) == 1
    print(sigs)
    sig = sigs[0]
    print("HASHES:", sig.minhash.hashes)

    assert sig.name == "short.fa"
    assert sig.minhash.ksize == 31
    assert sig.minhash.is_dna
    assert sig.minhash.scaled == 10

    # validate against sourmash sketch
    output2 = runtmp.output("short2.sig")
    runtmp.sourmash("sketch", "dna", fa1, "-o", output2, "-p", "scaled=10")
    sig2 = sourmash.load_one_signature(output2)
    assert sig.minhash.hashes == sig2.minhash.hashes


def test_manysketch_skipm2n3(runtmp, capfd):
    fa_csv = runtmp.output("db-fa.csv")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")
    fa3 = get_test_data("short3.fa")
    protfa1 = get_test_data("short-protein.fa")

    make_assembly_csv(fa_csv, [fa1, fa2, fa3], [protfa1])

    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts",
        "manysketch",
        fa_csv,
        "-o",
        output,
        "--param-str",
        "dna,k=21,scaled=1",
        "--param-str",
        "skipm2n3,k=31,scaled=30",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty
    captured = capfd.readouterr()
    print(captured.err)

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)

    # note: requires sourmash v4.8.13 or later.
    assert len(sigs) == 6  # 3 dna, 3 skipmer.

    # check moltypes, etc!
    dna_md5sums = {
        "short": "1474578c5c46dd09da4c2df29cf86621",
        "short2": "4efeebd26644278e36b9553e018a851a",
        "short3": "f85747ac4f473c4a71c1740d009f512b",
    }
    skip_md5sums = {
        "short2": "ec6305f5d82e51659f3914d47fcc32ee",
        "short": "0486fcae73545363da9cd5bfcf18d322",
        "short3": "890557b39ae66d3177035296818de7c6",
    }
    for sig in sigs:
        if sig.minhash.is_dna:
            assert sig.minhash.ksize == 21
            assert sig.minhash.scaled == 1
            print("DNA: ", sig.name, sig.md5sum())
            assert sig.md5sum() == dna_md5sums[sig.name]
        elif sig.minhash.moltype == "skipm2n3":
            print(sig.minhash.ksize, sig.minhash.scaled, sig.name, sig.md5sum())
            assert sig.minhash.ksize == 31
            assert sig.minhash.scaled == 30
            assert sig.md5sum() == skip_md5sums[sig.name]

    # read the file with python and check sigs
    import zipfile, gzip, json

    with zipfile.ZipFile(output, "r") as zf:
        # Check the manifest exists
        assert "SOURMASH-MANIFEST.csv" in zf.namelist()

        expected_signatures = [
            {
                "name": "short",
                "ksize": 31,
                "scaled": 30,
                "moltype": "skipm2n3",
                "md5sum": "0486fcae73545363da9cd5bfcf18d322",
            },
            {
                "name": "short3",
                "ksize": 31,
                "scaled": 30,
                "moltype": "skipm2n3",
                "md5sum": "890557b39ae66d3177035296818de7c6",
            },
            {
                "name": "short2",
                "ksize": 31,
                "scaled": 30,
                "moltype": "skipm2n3",
                "md5sum": "ec6305f5d82e51659f3914d47fcc32ee",
            },
        ]
        expected_signatures_dict = {exp["md5sum"]: exp for exp in expected_signatures}

        # Read and parse the manifest
        with zf.open("SOURMASH-MANIFEST.csv") as manifest_file:
            manifest_data = manifest_file.read().decode("utf-8").splitlines()
            manifest_data = [line for line in manifest_data if not line.startswith("#")]
            manifest_reader = csv.DictReader(manifest_data)

            for row in manifest_reader:
                if row["moltype"] == "skipm2n3":
                    print("Manifest Row:", row)

                    # Validate row fields
                    md5sum = row["md5"]
                    assert (
                        md5sum in expected_signatures_dict
                    ), f"Unexpected md5sum: {md5sum}"
                    expected = expected_signatures_dict[md5sum]
                    assert (
                        row["name"] == expected["name"]
                    ), f"Name mismatch: {row['name']}"
                    assert (
                        int(row["ksize"]) == expected["ksize"]
                    ), f"Ksize mismatch: {row['ksize']}"
                    assert (
                        row["moltype"] == expected["moltype"]
                    ), f"Moltype mismatch: {row['moltype']}"

                    sig_path = row["internal_location"]
                    assert sig_path.startswith("signatures/")

                    # Extract and read the signature file
                    with zf.open(sig_path) as sig_gz:
                        with gzip.open(sig_gz, "rt") as sig_file:
                            sig_contents = json.load(sig_file)
                            print("Signature Contents:", sig_contents)

                            # Validate signature contents
                            sig_data = sig_contents[0]
                            print(sig_data)
                            siginfo = sig_data["signatures"][0]
                            assert (
                                siginfo["md5sum"] == md5sum
                            ), f"MD5 mismatch: {siginfo['md5sum']}"
                            assert (
                                siginfo["ksize"] == expected["ksize"]
                            ), f"Ksize mismatch: {siginfo['ksize']}"
                            assert (
                                siginfo["molecule"] == expected["moltype"]
                            ), f"Moltype mismatch: {siginfo['molecule']}"


def test_singlesketch_skipm2n3(runtmp):
    """Test singlesketch with skipm2n3."""
    fa1 = get_test_data("short.fa")
    output = runtmp.output("short.sig")

    # Run the singlesketch command
    runtmp.sourmash(
        "scripts", "singlesketch", fa1, "-p", "skipm2n3,k=31,scaled=100", "-o", output
    )

    # Check if the output exists and contains the expected data
    assert os.path.exists(output)
    # Load the output signature file
    import json

    with open(output, "r") as f:
        data = json.load(f)

    # Extract the signature part from the JSON
    signatures = data[0]["signatures"]

    # Expected signature fields
    expected_signatures = [
        {
            "name": "short.fa",
            "ksize": 31,
            "moltype": "skipm2n3",
            "md5sum": "387d25c8e4b4878c78872efd13621491",
        }
    ]

    # Check if the signatures match the expected
    assert len(signatures) == len(
        expected_signatures
    ), "Number of signatures does not match."

    for sig, expected in zip(signatures, expected_signatures):
        assert sig["ksize"] == expected["ksize"], f"Unexpected ksize: {sig['ksize']}"
        assert (
            sig["molecule"] == expected["moltype"]
        ), f"Unexpected moltype: {sig['molecule']}"
        assert (
            sig["md5sum"] == expected["md5sum"]
        ), f"Unexpected md5sum: {sig['md5sum']}"
        assert (
            data[0]["name"] == expected["name"]
        ), f"Unexpected name: {data[0]['name']}"


def test_singlesketch_stdin(runtmp):
    """Test basic single sketching with default parameters."""
    fa1 = get_test_data("short.fa")
    output = runtmp.output("short.sig")

    # Run the singlesketch command using subprocess
    cmd = f"cat {fa1} | sourmash scripts singlesketch - --name short -o {output} -p dna,scaled=10"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Check if the command succeeded
    assert result.returncode == 0, f"Command failed: {result.stderr}"

    # Check if the output exists and contains the expected data
    assert os.path.exists(output)
    sig = sourmash.load_one_signature(output)

    assert sig.name == "short"
    assert sig.minhash.ksize == 31
    assert sig.minhash.is_dna
    assert sig.minhash.scaled == 10
    print("HASHES:", sig.minhash.hashes)

    # validate against sourmash sketch
    output2 = runtmp.output("short2.sig")
    runtmp.sourmash("sketch", "dna", fa1, "-o", output2, "-p", "dna,scaled=10")
    sig2 = sourmash.load_one_signature(output2)
    assert sig.minhash.hashes == sig2.minhash.hashes


def test_singlesketch_multifiles(runtmp, capfd):
    # multiple input files to singlesketch
    fa_csv = runtmp.output("db-fa.csv")

    fa1 = get_test_data("short.fa")
    fa2 = get_test_data("short2.fa")

    output = runtmp.output("db.zip")

    runtmp.sourmash(
        "scripts",
        "singlesketch",
        fa1,
        fa2,
        "-o",
        output,
        "--param-str",
        "dna,k=31,scaled=1",
    )

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)
    assert "calculated 1 signatures for 2 sequences in 2 files" in captured.err

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    print(sigs)
    assert len(sigs) == 1
    made_sig = sigs[0]
    assert made_sig.name == "short.fa"

    s1 = runtmp.output("short.sig")
    runtmp.sourmash(
        "sketch",
        "dna",
        fa1,
        fa2,
        "-o",
        s1,
        "--param-str",
        "k=31,scaled=1",
        "--name",
        "short.fa",
    )
    sig1 = sourmash.load_one_signature(s1)

    assert made_sig == sig1
