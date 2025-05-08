import os
import csv
import pytest
import pandas
import sourmash

from . import sourmash_tst_utils as utils
from .sourmash_tst_utils import (
    get_test_data,
    make_file_list,
    zip_siglist,
    index_siglist,
)


def test_installed(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "pairwise")

    assert "usage:  pairwise" in runtmp.last_result.err


def test_on_dir(runtmp, capfd):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash(
            "scripts", "pairwise", runtmp.output(""), "-o", runtmp.output("xxx.csv")
        )

    captured = capfd.readouterr()
    print(captured.err)
    assert "arbitrary directories are not supported" in captured.err


def test_simple_no_ani(runtmp, capfd, zip_query, indexed):
    # test basic execution!
    query_list = runtmp.output("query.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    if indexed:
        query_list = index_siglist(runtmp, query_list, runtmp.output("db"))

    runtmp.sourmash("scripts", "pairwise", query_list, "-o", output, "-t", "-1")
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 3

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        # confirm hand-checked numbers
        q = row["query_name"].split()[0]
        m = row["match_name"].split()[0]
        cont = float(row["containment"])
        jaccard = float(row["jaccard"])
        maxcont = float(row["max_containment"])
        intersect_hashes = int(row["intersect_hashes"])
        assert "query_containment_ani" not in row
        assert "match_containment_ani" not in row
        assert "average_containment_ani" not in row
        assert "max_containment_ani" not in row

        jaccard = round(jaccard, 4)
        cont = round(cont, 4)
        maxcont = round(maxcont, 4)
        print(q, m, f"{jaccard:.04}", f"{cont:.04}", f"{maxcont:.04}")

        if q == "NC_011665.1" and m == "NC_009661.1":
            assert jaccard == 0.3207
            assert cont == 0.4828
            assert maxcont == 0.4885
            assert intersect_hashes == 2529

    captured = capfd.readouterr()
    print(captured.err)

    if indexed:
        assert (
            "WARNING: loading all sketches from a RocksDB into memory!" in captured.err
        )


def test_simple_no_ani_output_all(runtmp, capfd, zip_query, indexed):
    # test basic execution!
    query_list = runtmp.output("query.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    if indexed:
        query_list = index_siglist(runtmp, query_list, runtmp.output("db"))

    runtmp.sourmash("scripts", "pairwise", query_list, "-o", output, "-t", "-1", "-A")
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 6


def test_simple_ani(runtmp, zip_query):
    # test basic execution!
    query_list = runtmp.output("query.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    runtmp.sourmash(
        "scripts", "pairwise", query_list, "-o", output, "-t", "-1", "--ani"
    )
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 3

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        # confirm hand-checked numbers
        q = row["query_name"].split()[0]
        m = row["match_name"].split()[0]
        cont = float(row["containment"])
        jaccard = float(row["jaccard"])
        maxcont = float(row["max_containment"])
        intersect_hashes = int(row["intersect_hashes"])
        q1_ani = float(row["query_containment_ani"])
        q2_ani = float(row["match_containment_ani"])
        avg_ani = float(row["average_containment_ani"])
        max_ani = float(row["max_containment_ani"])

        jaccard = round(jaccard, 4)
        cont = round(cont, 4)
        maxcont = round(maxcont, 4)
        q1_ani = round(q1_ani, 4)
        q2_ani = round(q2_ani, 4)
        avg_ani = round(avg_ani, 4)
        max_ani = round(max_ani, 4)
        print(
            q,
            m,
            f"{jaccard:.04}",
            f"{cont:.04}",
            f"{maxcont:.04}",
            f"{q1_ani:.04}",
            f"{q2_ani:.04}",
            f"{avg_ani:.04}",
            f"{max_ani:.04}",
        )

        if q == "NC_011665.1" and m == "NC_009661.1":
            assert jaccard == 0.3207
            assert cont == 0.4828
            assert maxcont == 0.4885
            assert intersect_hashes == 2529
            assert q1_ani == 0.9768
            assert q2_ani == 0.9772
            assert avg_ani == 0.977
            assert max_ani == 0.9772
        elif m == "NC_011665.1" and q == "NC_009661.1":
            assert jaccard == 0.3207
            assert cont == 0.4885
            assert maxcont == 0.4885
            assert intersect_hashes == 2529
            assert q2_ani == 0.9768
            assert q1_ani == 0.9772
            assert avg_ani == 0.977
            assert max_ani == 0.9772
        elif q == m:
            assert jaccard == 1
        else:
            assert jaccard == 0
            assert cont == 0
            assert maxcont == 0
            assert intersect_hashes == 0
            assert q1_ani == 0
            assert q2_ani == 0
            assert avg_ani == 0
            assert max_ani == 0


def test_simple_threshold(runtmp, zip_query):
    # test with a simple threshold => only 3 results
    query_list = runtmp.output("query.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    runtmp.sourmash("scripts", "pairwise", query_list, "-o", output, "-t", "0.1")
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 1


def test_simple_manifest(runtmp):
    # test with a simple threshold => only 3 results
    query_list = runtmp.output("query.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")

    query_mf = runtmp.output("qmf.csv")

    runtmp.sourmash("sig", "manifest", query_list, "-o", query_mf)

    runtmp.sourmash("scripts", "pairwise", query_mf, "-o", output, "-t", "0.1")
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 1


def test_sig_query(runtmp, capfd):
    # sig query is ok now, but fails bc only one sig
    sig2 = get_test_data("2.fa.sig.gz")

    output = runtmp.output("out.csv")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "pairwise", sig2, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: Pairwise requires two or more sketches. Check input" in captured.err


def test_bad_query(runtmp, capfd):
    # test with a bad query list (a missing file)
    query_list = runtmp.output("query.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    make_file_list(query_list, [sig2, sig47, "no-exist"])

    output = runtmp.output("out.csv")

    runtmp.sourmash("scripts", "pairwise", query_list, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)

    assert "WARNING: could not load sketches from path 'no-exist'" in captured.err
    assert (
        "WARNING: 1 analysis paths failed to load. See error messages above."
        in captured.err
    )


def test_missing_query(runtmp, capfd, zip_db):
    # test with a missing query list
    query_list = runtmp.output("query.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    output = runtmp.output("out.csv")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "pairwise", query_list, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)

    assert "Error: No such file or directory" in captured.err


def test_empty_query(runtmp, capfd):
    # test with an empty query list
    query_list = runtmp.output("query.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [])

    output = runtmp.output("out.csv")

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "pairwise", query_list, "-o", output)

    captured = capfd.readouterr()
    assert "Error: No analysis signatures loaded, exiting." in captured.err


def test_nomatch_query_warn(runtmp, capfd, zip_query):
    # test a non-matching (diff ksize) in query; do we get warning message?
    query_list = runtmp.output("query.txt")

    sig1 = get_test_data("1.fa.k21.sig.gz")
    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63, sig1])

    output = runtmp.output("out.csv")

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    runtmp.sourmash("scripts", "pairwise", query_list, "-o", output)
    assert os.path.exists(output)

    captured = capfd.readouterr()
    print(captured.err)

    assert (
        "WARNING: skipped 1 analysis paths - no compatible signatures" in captured.err
    )


def test_nomatch_query_exit(runtmp, capfd, zip_query):
    # test a non-matching (diff ksize) in query; do we get warning message?
    query_list = runtmp.output("query.txt")

    sig1 = get_test_data("1.fa.k21.sig.gz")
    sig2 = get_test_data("2.fa.k21.sig.gz")

    make_file_list(query_list, [sig1, sig2])

    output = runtmp.output("out.csv")

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash("scripts", "pairwise", query_list, "-o", output)

    captured = capfd.readouterr()
    print(captured.err)

    assert (
        "WARNING: skipped 2 analysis paths - no compatible signatures" in captured.err
    )
    assert "Error: No analysis signatures loaded, exiting." in captured.err


def test_load_only_one_bug(runtmp, capfd, zip_db):
    # check that we behave properly when presented with multiple query
    # sketches
    query_list = runtmp.output("query.txt")

    sig1_k31 = get_test_data("1.fa.k31.sig.gz")

    # note: this was created as a 3-sketch-in-one-signature directly
    # via sourmash sketch dna -p k=21,k=31,k=51.
    sig1_all = get_test_data("1.combined.sig.gz")

    make_file_list(query_list, [sig1_all, sig1_k31])

    if zip_db:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("db.zip"))

    output = runtmp.output("out.csv")

    runtmp.sourmash("scripts", "pairwise", query_list, "-o", output)
    assert os.path.exists(output)

    captured = capfd.readouterr()
    print(captured.err)

    assert not "WARNING: skipped 1 paths - no compatible signatures." in captured.err
    assert not "WARNING: no compatible sketches in path " in captured.err


def test_md5(runtmp, zip_query):
    # test that md5s match what was in the original files, not downsampled etc.
    query_list = runtmp.output("query.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")

    if zip_query:
        query_list = zip_siglist(runtmp, query_list, runtmp.output("query.zip"))

    runtmp.sourmash("scripts", "pairwise", query_list, "-o", output, "-t", "-0.1")
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 3

    md5s = list(df["query_md5"]) + list(df["match_md5"])
    print(f"md5s: {md5s}")

    for query_file in (sig2, sig47, sig63):
        for ss in sourmash.load_file_as_signatures(query_file, ksize=31):
            assert ss.md5sum() in md5s

    md5s = list(df["match_md5"])
    print(md5s)


def test_simple_prot_ani(runtmp):
    # test basic execution with protein sigs
    sigs = get_test_data("protein.zip")

    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts",
        "pairwise",
        sigs,
        "-o",
        output,
        "--moltype",
        "protein",
        "-k",
        "19",
        "--scaled",
        "100",
        "--ani",
    )
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 1

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        # confirm hand-checked numbers
        q = row["query_name"].split()[0]
        m = row["match_name"].split()[0]
        cont = float(row["containment"])
        jaccard = float(row["jaccard"])
        maxcont = float(row["max_containment"])
        intersect_hashes = int(row["intersect_hashes"])
        q1_ani = float(row["query_containment_ani"])
        q2_ani = float(row["match_containment_ani"])
        avg_ani = float(row["average_containment_ani"])
        max_ani = float(row["max_containment_ani"])

        jaccard = round(jaccard, 4)
        cont = round(cont, 4)
        maxcont = round(maxcont, 4)
        q1_ani = round(q1_ani, 4)
        q2_ani = round(q2_ani, 4)
        avg_ani = round(avg_ani, 4)
        max_ani = round(max_ani, 4)
        print(
            q,
            m,
            f"{jaccard:.04}",
            f"{cont:.04}",
            f"{maxcont:.04}",
            intersect_hashes,
            f"{q1_ani:.04}",
            f"{q2_ani:.04}",
            f"{avg_ani:.04}",
            f"{max_ani:.04}",
        )

        if q == "GCA_001593925" and m == "GCA_001593935":
            assert jaccard == 0.0434
            assert cont == 0.1003
            assert maxcont == 0.1003
            assert intersect_hashes == 342
            assert q1_ani == 0.8860
            assert q2_ani == 0.8702
            assert avg_ani == 0.8781
            assert max_ani == 0.886


def test_simple_dayhoff_ani(runtmp):
    # test basic execution with dayhoff sigs
    sigs = get_test_data("dayhoff.zip")

    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts",
        "pairwise",
        sigs,
        "-o",
        output,
        "--moltype",
        "dayhoff",
        "-k",
        "19",
        "--scaled",
        "100",
        "--ani",
    )
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 1

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        # confirm hand-checked numbers
        q = row["query_name"].split()[0]
        m = row["match_name"].split()[0]
        cont = float(row["containment"])
        jaccard = float(row["jaccard"])
        maxcont = float(row["max_containment"])
        intersect_hashes = int(row["intersect_hashes"])
        q1_ani = float(row["query_containment_ani"])
        q2_ani = float(row["match_containment_ani"])
        avg_ani = float(row["average_containment_ani"])
        max_ani = float(row["max_containment_ani"])

        jaccard = round(jaccard, 4)
        cont = round(cont, 4)
        maxcont = round(maxcont, 4)
        q1_ani = round(q1_ani, 4)
        q2_ani = round(q2_ani, 4)
        avg_ani = round(avg_ani, 4)
        max_ani = round(max_ani, 4)
        print(
            q,
            m,
            f"{jaccard:.04}",
            f"{cont:.04}",
            f"{maxcont:.04}",
            intersect_hashes,
            f"{q1_ani:.04}",
            f"{q2_ani:.04}",
            f"{avg_ani:.04}",
            f"{max_ani:.04}",
        )

        if q == "GCA_001593925" and m == "GCA_001593935":
            assert jaccard == 0.1326
            assert cont == 0.2815
            assert maxcont == 0.2815
            assert intersect_hashes == 930
            assert q1_ani == 0.9355
            assert q2_ani == 0.9189
            assert avg_ani == 0.9272
            assert max_ani == 0.9355


def test_simple_hp_ani(runtmp):
    # test basic execution with hp sigs
    sigs = get_test_data("hp.zip")

    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts",
        "pairwise",
        sigs,
        "-o",
        output,
        "--moltype",
        "hp",
        "-k",
        "19",
        "--scaled",
        "100",
        "--ani",
    )
    assert os.path.exists(output)

    df = pandas.read_csv(output)
    assert len(df) == 1

    dd = df.to_dict(orient="index")
    print(dd)

    for idx, row in dd.items():
        # confirm hand-checked numbers
        q = row["query_name"].split()[0]
        m = row["match_name"].split()[0]
        cont = float(row["containment"])
        jaccard = float(row["jaccard"])
        maxcont = float(row["max_containment"])
        intersect_hashes = int(row["intersect_hashes"])
        q1_ani = float(row["query_containment_ani"])
        q2_ani = float(row["match_containment_ani"])
        avg_ani = float(row["average_containment_ani"])
        max_ani = float(row["max_containment_ani"])

        jaccard = round(jaccard, 4)
        cont = round(cont, 4)
        maxcont = round(maxcont, 4)
        q1_ani = round(q1_ani, 4)
        q2_ani = round(q2_ani, 4)
        avg_ani = round(avg_ani, 4)
        max_ani = round(max_ani, 4)
        print(
            q,
            m,
            f"{jaccard:.04}",
            f"{cont:.04}",
            f"{maxcont:.04}",
            intersect_hashes,
            f"{q1_ani:.04}",
            f"{q2_ani:.04}",
            f"{avg_ani:.04}",
            f"{max_ani:.04}",
        )

        if q == "GCA_001593925" and m == "GCA_001593935":
            assert jaccard == 0.4983
            assert cont == 0.747
            assert maxcont == 0.747
            assert intersect_hashes == 1724
            assert q1_ani == 0.9848
            assert q2_ani == 0.9734
            assert avg_ani == 0.9791
            assert max_ani == 0.9848


def test_simple_below_threshold(runtmp):
    # test basic execution!
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts", "pairwise", query_list, "-o", output, "--ani", "--threshold", "0.5"
    )
    assert os.path.exists(output)

    with open(output, "r") as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
        print(rows)
        assert len(rows) == 0


def test_simple_below_threshold_write_all(runtmp):
    # test basic execution!
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts",
        "pairwise",
        query_list,
        "-o",
        output,
        "--ani",
        "--threshold",
        "0.5",
        "--write-all",
    )
    assert os.path.exists(output)

    with open(output, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)
        print(rows)
        assert len(rows) == 3
        for row in rows:
            assert float(row["query_containment_ani"]) == 1.0
            assert float(row["match_containment_ani"]) == 1.0
            assert float(row["average_containment_ani"]) == 1.0
            assert float(row["max_containment_ani"]) == 1.0
            assert float(row["containment"]) == 1.0
            assert float(row["max_containment"]) == 1.0
            assert float(row["jaccard"]) == 1.0
            assert row["query_name"] == row["match_name"]
            assert row["query_md5"] == row["match_md5"]


def test_simple_below_threshold_write_all_no_ani(runtmp):
    # test basic execution!
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")

    runtmp.sourmash(
        "scripts",
        "pairwise",
        query_list,
        "-o",
        output,
        "--threshold",
        "0.5",
        "--write-all",
    )
    assert os.path.exists(output)

    with open(output, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)
        print(rows)
        assert len(rows) == 3
        for row in rows:
            assert "query_containment_ani" not in row.keys()
            assert "match_containment_ani" not in row.keys()
            assert "average_containment_ani" not in row.keys()
            assert "max_containment_ani" not in row.keys()
            assert float(row["containment"]) == 1.0
            assert float(row["max_containment"]) == 1.0
            assert float(row["jaccard"]) == 1.0
            assert row["query_name"] == row["match_name"]
            assert row["query_md5"] == row["match_md5"]


def test_simple_scaled(runtmp):
    # test basic execution w/scaled!
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")

    make_file_list(query_list, [sig2, sig47, sig63])

    output = runtmp.output("out.csv")

    runtmp.sourmash("scripts", "pairwise", query_list, "-o", output, "-s", "10_000")
    assert os.path.exists(output)
    df = pandas.read_csv(output)
    assert len(df) == 1
    assert set(list(df["scaled"])) == {10_000}


def test_simple_scaled_heterogenous(runtmp):
    # test basic execution w/heterogeneous scaled
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")
    sig47_ds = runtmp.output("47-10k.sig.zip")

    runtmp.sourmash("sig", "downsample", sig47, "-o", sig47_ds, "--scaled", "10_000")

    make_file_list(query_list, [sig2, sig47_ds, sig63])

    output = runtmp.output("out.csv")

    runtmp.sourmash("scripts", "pairwise", query_list, "-o", output)
    assert os.path.exists(output)
    df = pandas.read_csv(output)
    assert len(df) == 1
    assert set(list(df["scaled"])) == {10_000}


def test_simple_scaled_heterogenous(runtmp):
    # test basic execution w/heterogeneous scaled - specified on command line
    query_list = runtmp.output("query.txt")
    against_list = runtmp.output("against.txt")

    sig2 = get_test_data("2.fa.sig.gz")
    sig47 = get_test_data("47.fa.sig.gz")
    sig63 = get_test_data("63.fa.sig.gz")
    sig47_ds = runtmp.output("47-10k.sig.zip")

    runtmp.sourmash("sig", "downsample", sig47, "-o", sig47_ds, "--scaled", "10_000")

    make_file_list(query_list, [sig2, sig47_ds, sig63])

    output = runtmp.output("out.csv")

    runtmp.sourmash("scripts", "pairwise", query_list, "-o", output, "--scaled=15_000")
    assert os.path.exists(output)
    df = pandas.read_csv(output)
    assert len(df) == 1
    assert set(list(df["scaled"])) == {15_000}
