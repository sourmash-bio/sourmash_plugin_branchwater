import os, csv

from . import sourmash_tst_utils as utils

def get_test_data(filename):
    thisdir = os.path.dirname(__file__)
    return os.path.join(thisdir, 'test-data', filename)

def test_cluster_containment(runtmp): 
    pairwise_csv = get_test_data('cluster.pairwise.csv')
    output = runtmp.output('clusters.csv')


    runtmp.sourmash('scripts', 'cluster', pairwise_csv, '-o', output,
                    '--similarity-column', "containment")

    assert os.path.exists(output)

    cluster_count = 0
    with open(output, 'r') as file:
        for line in file:
            if line.startswith("Component"):
                cluster_count += 1

    # to do: check to see if this is actually right
    assert cluster_count == 6


def test_cluster_jaccard(runtmp): 
    pairwise_csv = get_test_data('cluster.pairwise.csv')
    output = runtmp.output('clusters.csv')


    runtmp.sourmash('scripts', 'cluster', pairwise_csv, '-o', output,
                    '--similarity-column', "jaccard")

    assert os.path.exists(output)

    cluster_count = 0
    with open(output, 'r') as file:
        for line in file:
            if line.startswith("Component"):
                cluster_count += 1

    # to do: check to see if this is actually right
    assert cluster_count == 8