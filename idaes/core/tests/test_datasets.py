"""
Tests for idaes.core.datasets module
"""
import pytest
from idaes.core import datasets


@pytest.mark.unit
def test_publication_unknown():
    # random, unknown publication
    with pytest.raises(KeyError):
        pub = datasets._Publication("test")


@pytest.mark.unit
def test_publication_known():
    # known publication
    pub = datasets._Publication("Pitzer:1984")
    assert pub
    assert pub.list_tables()
    assert pub.get_table(pub.list_tables()[0])


@pytest.mark.unit
def test_pitzer():
    p = datasets.Pitzer()
    assert p
    assert p.list_tables()
    assert p.get_table(p.list_tables()[0])
