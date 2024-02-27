import os, csv
import pytest

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
        runtmp.sourmash('scripts', 'cluster')

    assert 'usage:  cluster' in runtmp.last_result.err


def test_cluster_help(runtmp):
    # test sourmash scripts cluster --help /-h
    runtmp.sourmash('scripts', 'cluster', '-h')

    print(runtmp.last_result.err)
    out = runtmp.last_result.out
    print(out)

    assert "usage:  cluster" in out
    assert "positional arguments:" in out
    assert "options:" in out


def test_cluster_containment(runtmp): 
    pairwise_csv = get_test_data('cluster.pairwise.csv')
    output = runtmp.output('clusters.csv')
    sizes = runtmp.output('sizes.csv')
    threshold = '0.5'

    runtmp.sourmash('scripts', 'cluster', pairwise_csv, '-o', output,
                    '--similarity-column', "containment", "--cluster-sizes",
                    sizes, '--threshold', threshold)

    assert os.path.exists(output)

    # check cluster output
    with open(output, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster','nodes']
    assert len(rows) == 1, f"Expected 1 data row but found {len(rows)}"
    assert rows[0]['cluster'] == 'Component_1'
    expected = set("n2;n3;n7;n1;n6;n5;n4".split(';'))
    assert set(rows[0]['nodes'].split(';')) == expected

    # check cluster size histogram
    with open(sizes, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster_size','count']
    assert len(rows) == 1, f"Expected 1 data row but found {len(rows)}"
    assert rows[0]['cluster_size'] == '7'
    assert rows[0]['count'] == '1'


def test_cluster_max_containment_1(runtmp):
    pairwise_csv = get_test_data('cluster.pairwise.csv')
    output = runtmp.output('clusters.csv')
    sizes = runtmp.output('sizes.csv')
    threshold = '0.7'

    runtmp.sourmash('scripts', 'cluster', pairwise_csv, '-o', output,
                    '--similarity-column', "max_containment", "--cluster-sizes",
                    sizes, '--threshold', threshold)

    assert os.path.exists(output)

    # check cluster output
    with open(output, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster','nodes']
    assert len(rows) == 1, f"Expected 1 data row but found {len(rows)}"
    assert rows[0]['cluster'] == 'Component_1'
    expected = set("n2;n3;n7;n1;n6;n5;n4".split(';'))
    assert set(rows[0]['nodes'].split(';')) == expected

    # check cluster size histogram
    with open(sizes, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster_size','count']
    assert len(rows) == 1, f"Expected 1 data row but found {len(rows)}"
    assert rows[0]['cluster_size'] == '7'
    assert rows[0]['count'] == '1'


def test_cluster_max_containment_2(runtmp):
    pairwise_csv = get_test_data('cluster.pairwise.csv')
    output = runtmp.output('clusters.csv')
    sizes = runtmp.output('sizes.csv')
    threshold = '0.9'

    runtmp.sourmash('scripts', 'cluster', pairwise_csv, '-o', output,
                    '--similarity-column', "max_containment", "--cluster-sizes",
                    sizes, '--threshold', threshold)

    assert os.path.exists(output)

    # check cluster output
    with open(output, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster','nodes']
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    assert rows[0]['cluster'] == 'Component_1'
    expected_node_sets = [
    set("n1;n2;n3;n4;n5".split(';')),
    set("n6;n7".split(';')),
    ]
    for row in rows:
        assert set(row['nodes'].split(';')) in expected_node_sets

    # check cluster size histogram
    with open(sizes, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster_size','count']
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    rows_as_tuples = {tuple(row.values()) for row in rows}
    expected = {('5', '1'), ('2', '1')}
    assert rows_as_tuples == expected


def test_cluster_jaccard(runtmp): 
    pairwise_csv = get_test_data('cluster.pairwise.csv')
    output = runtmp.output('clusters.csv')
    sizes = runtmp.output('sizes.csv')
    threshold = '0.6'

    runtmp.sourmash('scripts', 'cluster', pairwise_csv, '-o', output,
                    '--similarity-column', "jaccard", "--cluster-sizes",
                    sizes, '--threshold', threshold)

    assert os.path.exists(output)

    # check cluster output
    with open(output, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster','nodes']
    assert len(rows) == 4, f"Expected 4 data rows but found {len(rows)}"
    assert rows[0]['cluster'] == 'Component_1'
    expected_node_sets = [
    set("n3;n4;n5;n6".split(';')),
    set("n1".split(';')),
    set("n2".split(';')),
    set("n7".split(';'))
    ]
    for row in rows:
        assert set(row['nodes'].split(';')) in expected_node_sets

    # check cluster size histogram
    with open(sizes, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster_size','count']
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    rows_as_tuples = {tuple(row.values()) for row in rows}
    expected = {('1', '3'), ('4', '1')}
    assert rows_as_tuples == expected


def test_cluster_default_similarity(runtmp):
    pairwise_csv = get_test_data('cluster.pairwise.csv')
    output = runtmp.output('clusters.csv')
    sizes = runtmp.output('sizes.csv')
    threshold = '0.9'

    runtmp.sourmash('scripts', 'cluster', pairwise_csv, '-o', output,
                    '--threshold', threshold)

    assert os.path.exists(output)

    # check cluster output
    with open(output, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster','nodes']
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    assert rows[0]['cluster'] == 'Component_1'
    expected_node_sets = [
    set("n1;n2;n3;n4;n5".split(';')),
    set("n6;n7".split(';'))
    ]
    for row in rows:
        assert set(row['nodes'].split(';')) in expected_node_sets

    # check cluster size histogram
    assert not os.path.exists(sizes)


def test_cluster_default_threshold(runtmp):
    # test default threshold (0.95)
    pairwise_csv = get_test_data('cluster.pairwise.csv')
    output = runtmp.output('clusters.csv')
    sizes = runtmp.output('sizes.csv')

    runtmp.sourmash('scripts', 'cluster', pairwise_csv, '-o', output)

    assert os.path.exists(output)

    # check cluster output
    with open(output, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster','nodes']
    assert len(rows) == 5, f"Expected 5 data rows but found {len(rows)}"
    assert rows[0]['cluster'] == 'Component_1'
    expected_node_sets = [
    set("n1".split(';')),
    set("n2;n3;n4".split(';')),
    set("n5".split(';')),
    set("n6".split(';')),
    set("n7".split(';'))
    ]
    for row in rows:
        assert set(row['nodes'].split(';')) in expected_node_sets

    # check cluster size histogram
    assert not os.path.exists(sizes)


def test_cluster_ani(runtmp):
    pairwise_csv = get_test_data('cluster.pairwise.csv')
    output = runtmp.output('clusters.csv')
    sizes = runtmp.output('sizes.csv')
    threshold = '0.9'

    runtmp.sourmash('scripts', 'cluster', pairwise_csv, '-o', output,
                    '--similarity-column', "average_containment_ani", "--cluster-sizes",
                    sizes, '--threshold', threshold)

    assert os.path.exists(output)

    # check cluster output
    with open(output, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster','nodes']
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    assert rows[0]['cluster'] == 'Component_1'
    expected_node_sets = [
    set("n1;n2;n3;n4;n5".split(';')),
    set("n6;n7".split(';'))
    ]
    for row in rows:
        assert set(row['nodes'].split(';')) in expected_node_sets

    # check cluster size histogram
    with open(sizes, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster_size','count']
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    rows_as_tuples = {tuple(row.values()) for row in rows}
    expected = {('5', '1'), ('2', '1')}
    assert rows_as_tuples == expected


def test_cluster_max_ani(runtmp):
    pairwise_csv = get_test_data('cluster.pairwise.csv')
    output = runtmp.output('clusters.csv')
    sizes = runtmp.output('sizes.csv')
    threshold = '0.9'

    runtmp.sourmash('scripts', 'cluster', pairwise_csv, '-o', output,
                    '--similarity-column', "max_containment_ani", "--cluster-sizes",
                    sizes, '--threshold', threshold)

    assert os.path.exists(output) 

    # check cluster output
    with open(output, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster','nodes']
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    assert rows[0]['cluster'] == 'Component_1'
    expected_node_sets = [set("n1;n2;n3;n4;n5".split(';')), set("n6;n7".split(';'))]
    for row in rows:
        assert set(row['nodes'].split(';')) in expected_node_sets

    # check cluster size histogram
    with open(sizes, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster_size','count']
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    rows_as_tuples = {tuple(row.values()) for row in rows}
    expected = {('5', '1'), ('2', '1')}
    assert rows_as_tuples == expected


def test_cluster_ani_pairwise(runtmp):
    pairwise_csv = runtmp.output('pairwise.csv')
    output = runtmp.output('clusters.csv')
    sizes = runtmp.output('sizes.csv')
    cluster_threshold = '0.90'

    query_list = runtmp.output('query.txt')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])

    runtmp.sourmash('scripts', 'pairwise', query_list,
                    '-o', pairwise_csv, "-t", "-0.1", "--ani")

    assert os.path.exists(pairwise_csv)

    runtmp.sourmash('scripts', 'cluster', pairwise_csv, '-o', output,
                    '--similarity-column', "average_containment_ani", "--cluster-sizes",
                    sizes, '--threshold', cluster_threshold)

    assert os.path.exists(output)

    # check cluster output
    with open(output, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster','nodes']
    print(rows)
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    assert rows[0]['cluster'] == 'Component_1'
    expected_node_sets = [set("NC_009661.1;NC_011665.1".split(';')), set("CP001071.1".split(';'))]
    for row in rows:
        assert set(row['nodes'].split(';')) in expected_node_sets

    # check cluster size histogram
    with open(sizes, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster_size','count']
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    rows_as_tuples = {tuple(row.values()) for row in rows}
    expected = {('1', '1'), ('2', '1')}
    assert rows_as_tuples == expected


def test_cluster_avg_ani_no_ani(runtmp, capfd):
    pairwise_csv = runtmp.output('pairwise.csv')
    output = runtmp.output('clusters.csv')
    sizes = runtmp.output('sizes.csv')
    cluster_threshold = '0.9'

    query_list = runtmp.output('query.txt')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])

    runtmp.sourmash('scripts', 'pairwise', query_list,
                    '-o', pairwise_csv, "-t", "-0.1") # do not pass `--ani`

    assert os.path.exists(pairwise_csv)

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'cluster', pairwise_csv, '-o', output,
                    '--similarity-column', "average_containment_ani", "--cluster-sizes",
                    sizes, '--threshold', cluster_threshold)

    print(runtmp.last_result.err)
    captured = capfd.readouterr()
    print(captured.err)
    assert 'average_containment_ani is None. Did you estimate ANI?' in captured.err


def test_cluster_max_ani_no_ani(runtmp, capfd):
    pairwise_csv = runtmp.output('pairwise.csv')
    output = runtmp.output('clusters.csv')
    sizes = runtmp.output('sizes.csv')
    cluster_threshold = '0.9'

    query_list = runtmp.output('query.txt')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])

    runtmp.sourmash('scripts', 'pairwise', query_list,
                    '-o', pairwise_csv, "-t", "-0.1") # do not pass `--ani`

    assert os.path.exists(pairwise_csv)

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'cluster', pairwise_csv, '-o', output,
                    '--similarity-column', "max_containment_ani", "--cluster-sizes",
                    sizes, '--threshold', cluster_threshold)

    print(runtmp.last_result.err)
    captured = capfd.readouterr()
    print(captured.err)
    assert 'max_containment_ani is None. Did you estimate ANI?' in captured.err


def test_cluster_ani_multisearch(runtmp):
    multisearch_csv = runtmp.output('multisearch.csv')
    output = runtmp.output('clusters.csv')
    sizes = runtmp.output('sizes.csv')
    cluster_threshold = '0.90'

    query_list = runtmp.output('query.txt')
    sig2 = get_test_data('2.fa.sig.gz')
    sig47 = get_test_data('47.fa.sig.gz')
    sig63 = get_test_data('63.fa.sig.gz')

    make_file_list(query_list, [sig2, sig47, sig63])

    runtmp.sourmash('scripts', 'multisearch', query_list, query_list,
                    '-o', multisearch_csv, "-t", "-0.1", "--ani")

    assert os.path.exists(multisearch_csv)

    runtmp.sourmash('scripts', 'cluster', multisearch_csv, '-o', output,
                    '--similarity-column', "average_containment_ani", "--cluster-sizes",
                    sizes, '--threshold', cluster_threshold)

    assert os.path.exists(output)

    # check cluster output
    with open(output, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster','nodes']
    print(rows)
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    assert rows[0]['cluster'] == 'Component_1'
    expected_node_sets = [set("NC_009661.1;NC_011665.1".split(';')), set("CP001071.1".split(';'))]
    for row in rows:
        assert set(row['nodes'].split(';')) in expected_node_sets

    # check cluster size histogram
    with open(sizes, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster_size','count']
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    rows_as_tuples = {tuple(row.values()) for row in rows}
    expected = {('1', '1'), ('2', '1')}
    assert rows_as_tuples == expected


def test_empty_file(runtmp, capfd):
    # test with an empty query list
    csv = runtmp.output('empty.csv')

    make_file_list(csv, [])

    output = runtmp.output('out.csv')
    out2 = runtmp.output('counts.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'cluster', csv,
                        '-o', output, '--cluster-sizes', out2)

    print(runtmp.last_result.err)
    captured = capfd.readouterr()
    print(captured.err)

    assert "Error: Failed to build graph" in captured.err


def test_bad_file(runtmp, capfd):
    # test with an empty query list
    csv = runtmp.output('bad.csv')
    with open(csv, 'w') as out:
        out.write('column1,column2')

    make_file_list(csv, [])

    output = runtmp.output('out.csv')
    out2 = runtmp.output('counts.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'cluster', csv,
                        '-o', output, '--cluster-sizes', out2)

    print(runtmp.last_result.err)
    captured = capfd.readouterr()
    print(captured.err)

    assert "Error: Failed to build graph" in captured.err
