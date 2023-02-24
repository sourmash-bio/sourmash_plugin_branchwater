import os
import pytest

import sourmash_tst_utils as utils


def get_test_data(filename):
    thisdir = os.path.dirname(__file__)
    return os.path.join(thisdir, filename)


def test_install(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'manysearch')

    assert 'usage:  manysearch' in runtmp.last_result.err
