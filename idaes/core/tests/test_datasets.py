"""
Tests for idaes.core.datasets module
"""
import pytest
from idaes.core import datasets


@pytest.mark.unit
def test_pitzer():
    p = datasets.Pitzer()
    assert p
    assert p.list_tables()
    assert p.get_table(p.list_tables()[0])
