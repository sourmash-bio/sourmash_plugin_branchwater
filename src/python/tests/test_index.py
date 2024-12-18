import os
import pytest
import pandas
import sourmash
import shutil

from . import sourmash_tst_utils as utils
from .sourmash_tst_utils import get_test_data, make_file_list, zip_siglist


def test_installed(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "index")

    assert "usage:  index" in runtmp.last_result.err


def test_index(runtmp, toggle_internal_storage):
    # test basic index!
    siglist = runtmp.output("db-sigs.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(siglist, [sig2, sig47, sig63])

    output = runtmp.output("db.rocksdb")

    runtmp.sourmash("scripts", "index", siglist, "-o", output, toggle_internal_storage)
    assert os.path.exists(output)
    print(runtmp.last_result.err)

    assert "index is done" in runtmp.last_result.err


def test_index_warning_message(runtmp, capfd):
    # test basic index when it has to load things into memory - see #451.
    siglist = runtmp.output("db-sigs.txt")

    # note: can't use zip w/o breaking index. See sourmash-bio/sourmash#3321.
    sig2 = get_test_data("2.sig.zip")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(siglist, [sig2, sig47, sig63])

    output = runtmp.output("db.rocksdb")

    runtmp.sourmash("scripts", "index", siglist, "-o", output)
    assert os.path.exists(output)
    print(runtmp.last_result.err)

    assert "index is done" in runtmp.last_result.err
    captured = capfd.readouterr()
    print(captured.err)
    assert (
        "WARNING: loading all sketches into memory in order to index." in captured.err
    )


def test_index_error_message(runtmp, capfd):
    # test basic index when it errors out b/c can't load
    siglist = runtmp.output("db-sigs.txt")

    # note: can't use zip w/o breaking index. See sourmash-bio/sourmash#3321.
    sig2 = get_test_data("2.sig.zip")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(siglist, [sig2, sig47, sig63])

    output = runtmp.output("db.rocksdb")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts", "index", siglist, "-o", output, "--no-internal-storage"
        )

    captured = capfd.readouterr()
    print(captured.err)
    assert "cannot index this type of collection with external storage" in captured.err


def test_index_recursive(runtmp, capfd):
    # test index of pathlist containing standalone manifest containing zip.
    # a little ridiculous, but should hit the various branches in
    # MultiCollection::load
    siglist = runtmp.output("db-sigs.txt")

    # our basic list of sketches...
    sig2_zip = get_test_data("2.sig.zip")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    # generate a standalone mf containing a sip
    standalone_mf = runtmp.output("stand-mf.csv")
    runtmp.sourmash("sig", "collect", "-F", "csv", "-o", standalone_mf, sig2_zip)

    # now make a file list containing that mf
    make_file_list(siglist, [standalone_mf, sig47, sig63])

    output = runtmp.output("db.rocksdb")

    runtmp.sourmash("scripts", "index", siglist, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)
    assert (
        "WARNING: loading all sketches into memory in order to index." in captured.err
    )
    assert "index is done" in runtmp.last_result.err
    assert "Indexing 3 sketches." in captured.err


def test_index_protein(runtmp, toggle_internal_storage):
    sigs = get_test_data("protein.zip")
    output = runtmp.output("db.rocksdb")

    runtmp.sourmash(
        "scripts",
        "index",
        sigs,
        "-k",
        "19",
        "-s",
        "100",
        "--moltype",
        "protein",
        "-o",
        output,
        toggle_internal_storage,
    )
    assert os.path.exists(output)
    print(runtmp.last_result.err)
    assert "index is done" in runtmp.last_result.err


def test_index_dayhoff(runtmp, toggle_internal_storage):
    sigs = get_test_data("dayhoff.zip")
    output = runtmp.output("db.rocksdb")

    runtmp.sourmash(
        "scripts",
        "index",
        sigs,
        "-k",
        "19",
        "-s",
        "100",
        "--moltype",
        "dayhoff",
        "-o",
        output,
        toggle_internal_storage,
    )
    assert os.path.exists(output)
    print(runtmp.last_result.err)
    assert "index is done" in runtmp.last_result.err


def test_index_protein(runtmp, toggle_internal_storage):
    sigs = get_test_data("hp.zip")
    output = runtmp.output("db.rocksdb")

    runtmp.sourmash(
        "scripts",
        "index",
        sigs,
        "-k",
        "19",
        "-s",
        "100",
        "--moltype",
        "hp",
        "-o",
        output,
        toggle_internal_storage,
    )
    assert os.path.exists(output)
    print(runtmp.last_result.err)
    assert "index is done" in runtmp.last_result.err


def test_index_missing_siglist(runtmp, capfd, toggle_internal_storage):
    # test missing siglist file
    siglist = runtmp.output("db-sigs.txt")
    output = runtmp.output("out.db")
    # make_file_list(siglist, []) # don't make siglist file

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts", "index", siglist, "-o", output, toggle_internal_storage
        )

    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: No such file or directory: " in captured.err


def test_index_sig(runtmp, capfd, toggle_internal_storage):
    # test index with a .sig.gz file instead of pathlist
    # (should work now)
    sig2 = get_test_data("2.fa.sig.gz")
    output = runtmp.output("out.db")

    runtmp.sourmash("scripts", "index", sig2, "-o", output, toggle_internal_storage)

    captured = capfd.readouterr()
    print(captured.err)
    print(runtmp.last_result.err)
    assert "index is done" in runtmp.last_result.err


def test_index_manifest(runtmp, capfd, toggle_internal_storage):
    # test index with a manifest file
    sig2 = get_test_data("2.fa.sig.gz")
    output = runtmp.output("out.db")
    sig_mf = runtmp.output("mf.csv")
    runtmp.sourmash("sig", "manifest", sig2, "-o", sig_mf)

    runtmp.sourmash("scripts", "index", sig_mf, "-o", output, toggle_internal_storage)

    captured = capfd.readouterr()
    print(captured.err)
    print(runtmp.last_result.err)
    assert "index is done" in runtmp.last_result.err


def test_index_bad_siglist_2(runtmp, capfd):
    # test with a bad siglist (containing a missing file)
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")
    make_file_list(against_list, [sig2, "no-exist"])

    db = runtmp.output("db.rocksdb")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "index", against_list, "-o", db)

    captured = capfd.readouterr()
    print(captured.err)
    assert "WARNING: could not load sketches from path 'no-exist'" in captured.err


def test_index_empty_siglist(runtmp, capfd):
    # test empty siglist file
    siglist = runtmp.output("db-sigs.txt")
    output = runtmp.output("out.db")
    make_file_list(siglist, [])  # empty

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "index", siglist, "-o", output)

    captured = capfd.readouterr()
    assert not os.path.exists(output)  # do we want an empty file, or no file?
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    print(captured.err)
    assert "Error: Signatures failed to load. Exiting." in captured.err


def test_index_nomatch(runtmp, capfd):
    # test index with a siglist file that has (only) a non-matching ksize sig
    siglist = runtmp.output("against.txt")
    db = runtmp.output("db.rocksdb")

    sig1 = get_test_data("1.fa.k21.sig.gz")
    make_file_list(siglist, [sig1])

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "index", siglist, "-o", db)

    captured = capfd.readouterr()
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    print(captured.out)
    print(captured.err)
    assert not os.path.exists(db)


def test_index_nomatch_sig_in_siglist(runtmp, capfd):
    # test index with a siglist file that has both matching and non-matching sigs
    siglist = runtmp.output("against.txt")
    db = runtmp.output("db.rocksdb")

    sig2 = get_test_data("2.fa.sig.gz")
    sig1 = get_test_data("1.fa.k21.sig.gz")
    make_file_list(siglist, [sig2, sig1])

    runtmp.sourmash("scripts", "index", siglist, "-o", db)

    captured = capfd.readouterr()
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    print(captured.out)
    print(captured.err)
    assert os.path.exists(db)


def test_index_zipfile(runtmp, capfd, toggle_internal_storage):
    # test basic index from sourmash zipfile
    siglist = runtmp.output("db-sigs.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(siglist, [sig2, sig47, sig63])

    zipf = runtmp.output("sigs.zip")

    runtmp.sourmash("sig", "cat", siglist, "-o", zipf)

    output = runtmp.output("db.rocksdb")

    runtmp.sourmash("scripts", "index", zipf, "-o", output, toggle_internal_storage)
    assert os.path.exists(output)
    print(runtmp.last_result.err)

    assert "index is done" in runtmp.last_result.err
    captured = capfd.readouterr()
    print(captured.err)


def test_index_zipfile_subdir(runtmp, capfd, toggle_internal_storage):
    # test index from sourmash zipfile in different directory.

    # this was a tough test to get to fail!! have to:
    # * use non-abspath for zip file creation
    # * use non-abspath to zip file for indexing
    # so that the relative path gets things wrong.

    siglist = runtmp.output("db-sigs.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    shutil.copyfile(sig2, runtmp.output("2.fa.sig.gz"))
    shutil.copyfile(sig47, runtmp.output("47.fa.sig.gz"))
    shutil.copyfile(sig63, runtmp.output("63.fa.sig.gz"))

    os.mkdir(runtmp.output("subdir"))

    zipf = "sigs.zip"

    runtmp.sourmash(
        "sig", "cat", "2.fa.sig.gz", "47.fa.sig.gz", "63.fa.sig.gz", "-o", zipf
    )

    output = runtmp.output("subdir/db.rocksdb")

    runtmp.sourmash(
        "scripts",
        "index",
        zipf,
        "-o",
        output,
        in_directory=runtmp.output(""),
        toggle_internal_storage=toggle_internal_storage,
    )
    assert os.path.exists(output)
    print(runtmp.last_result.err)

    assert "index is done" in runtmp.last_result.err
    captured = capfd.readouterr()
    print(captured.err)

    runtmp.sourmash(
        "scripts", "check", "db.rocksdb", in_directory=runtmp.output("subdir")
    )
    runtmp.sourmash(
        "scripts", "check", "subdir/db.rocksdb", in_directory=runtmp.output("")
    )


def test_index_zipfile_repeated_md5sums(runtmp, capfd, toggle_internal_storage):
    # test that we're reading all files, including repeated md5sums
    siglist = runtmp.output("db-sigs.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig2a = runtmp.output("sig2a.sig.gz")
    sig2b = runtmp.output("sig2b.sig.gz")
    runtmp.sourmash("sig", "rename", sig2, "name2", "-o", sig2a)
    runtmp.sourmash("sig", "rename", sig2, "name3", "-o", sig2b)

    make_file_list(siglist, [sig2, sig2a, sig2b])

    zipf = runtmp.output("sigs.zip")
    runtmp.sourmash("sig", "cat", siglist, "-o", zipf)

    output = runtmp.output("db.rocksdb")

    runtmp.sourmash("scripts", "index", zipf, "-o", output, toggle_internal_storage)
    assert os.path.exists(output)
    print(runtmp.last_result.err)

    captured = capfd.readouterr()
    print(captured.err)

    assert "index is done" in runtmp.last_result.err


def test_index_zipfile_multiparam(runtmp, capfd, toggle_internal_storage):
    # test index from sourmash zipfile with multiple ksizes / scaled /moltype
    # SHOULD FAIL.
    siglist = runtmp.output("db-sigs.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")
    sig1 = get_test_data("1.combined.sig.gz")
    srr = get_test_data("SRR606249.sig.gz")
    prot = get_test_data("protein.zip")

    make_file_list(siglist, [sig2, sig47, sig63, sig1, srr, prot])

    zipf = runtmp.output("sigs.zip")

    runtmp.sourmash("sig", "cat", siglist, "-o", zipf)

    output = runtmp.output("db.rocksdb")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "index", zipf, "-o", output, toggle_internal_storage)


def test_index_check(runtmp, toggle_internal_storage):
    # test check index
    siglist = runtmp.output("db-sigs.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")

    make_file_list(siglist, [sig2, sig47])

    output = runtmp.output("db.rocksdb")

    runtmp.sourmash("scripts", "index", siglist, "-o", output, toggle_internal_storage)

    runtmp.sourmash("scripts", "check", output)
    print(runtmp.last_result.err)

    assert "index is ok" in runtmp.last_result.err


def test_index_check_quick(runtmp, toggle_internal_storage):
    # test check index
    siglist = runtmp.output("db-sigs.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")

    make_file_list(siglist, [sig2, sig47])

    output = runtmp.output("db.rocksdb")

    runtmp.sourmash("scripts", "index", siglist, "-o", output, toggle_internal_storage)

    runtmp.sourmash("scripts", "check", "--quick", output)
    print(runtmp.last_result.err)

    assert "index is ok" in runtmp.last_result.err


def test_index_subdir(runtmp, toggle_internal_storage):
    # test basic index & output to subdir
    siglist = runtmp.output("db-sigs.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(siglist, [sig2, sig47, sig63])

    os.mkdir(runtmp.output("subdir"))
    output = runtmp.output("subdir/db.rocksdb")

    runtmp.sourmash("scripts", "index", siglist, "-o", output, toggle_internal_storage)
    assert os.path.exists(output)
    print(runtmp.last_result.err)

    runtmp.sourmash("scripts", "check", output)
