import os, csv

from . import sourmash_tst_utils as utils

def get_test_data(filename):
    thisdir = os.path.dirname(__file__)
    return os.path.join(thisdir, 'test-data', filename)

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
    expected = set("n1;n2;n3;n4;n5".split(';'))
    assert set(rows[0]['nodes'].split(';')) == expected
    expected = set("n6;n7".split(';'))
    assert set(rows[1]['nodes'].split(';')) == expected

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
    assert len(rows) == 1, f"Expected 1 data row but found {len(rows)}"
    assert rows[0]['cluster'] == 'Component_1'
    expected = set("n3;n4;n5;n6".split(';'))
    assert set(rows[0]['nodes'].split(';')) == expected

    # check cluster size histogram
    with open(sizes, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster_size','count']
    assert len(rows) == 1, f"Expected 1 data row but found {len(rows)}"
    assert rows[0]['cluster_size'] == '4'
    assert rows[0]['count'] == '1'


def test_cluster_ani(runtmp):
    pairwise_csv = get_test_data('cluster.pairwise.csv')
    output = runtmp.output('clusters.csv')
    sizes = runtmp.output('sizes.csv')
    threshold = '0.9'

    runtmp.sourmash('scripts', 'cluster', pairwise_csv, '-o', output,
                    '--similarity-column', "ani", "--cluster-sizes",
                    sizes, '--threshold', threshold)

    assert os.path.exists(output)

    # check cluster output
    with open(output, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster','nodes']
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    assert rows[0]['cluster'] == 'Component_1'
    expected = set("n1;n2;n3;n4;n5".split(';'))
    assert set(rows[0]['nodes'].split(';')) == expected
    expected = set("n6;n7".split(';'))
    assert set(rows[1]['nodes'].split(';')) == expected

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
                    '--similarity-column', "max_ani", "--cluster-sizes",
                    sizes, '--threshold', threshold)

    assert os.path.exists(output)

    # check cluster output
    with open(output, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster','nodes']
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    assert rows[0]['cluster'] == 'Component_1'
    expected = set("n1;n2;n3;n4;n5".split(';'))
    assert set(rows[0]['nodes'].split(';')) == expected
    expected = set("n6;n7".split(';'))
    assert set(rows[1]['nodes'].split(';')) == expected

    # check cluster size histogram
    with open(sizes, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [row for row in reader]
        assert reader.fieldnames == ['cluster_size','count']
    assert len(rows) == 2, f"Expected 2 data rows but found {len(rows)}"
    rows_as_tuples = {tuple(row.values()) for row in rows}
    expected = {('5', '1'), ('2', '1')}
    assert rows_as_tuples == expected
