##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Tests for idaes.vis.plotbase
"""
# third-party
import pytest
import pandas as pd
# package
from idaes.vis import plotbase


def test_plotbase_construct_abstract():
    assert pytest.raises(TypeError, plotbase.PlotBase, 'foo')


def test_plotregistry_construct():
    reg = plotbase.PlotRegistry()


def test_plotregistry_register():
    reg = plotbase.PlotRegistry()
    reg.remove_all()
    obj = "hello"
    plot_class = plotbase.PlotBase
    reg.register(obj, "test1", plot_class)
    reg.remove_all()


def setup_plot(plot, obj):
    return plot, obj


def test_plotregistry_register_setup():
    reg = plotbase.PlotRegistry()
    reg.remove_all()
    obj = "hello"
    plot_class = plotbase.PlotBase
    reg.register(obj, "test1", plot_class, setup_fn=setup_plot)


def test_plotregistry_register_overwrite():
    reg = plotbase.PlotRegistry()
    reg.remove_all()
    obj = "hello"
    plot_class = plotbase.PlotBase
    reg.register(obj, "test1", plot_class)
    # without overwrite: KeyError
    assert pytest.raises(
        KeyError, reg.register, obj, "test1", plot_class, dict(setup_fn=setup_plot)
    )
    # with overwrite: ok
    reg.register(obj, "test1", plot_class, overwrite=True)


def test_plotregistry_get():
    reg = plotbase.PlotRegistry()
    reg.remove_all()
    obj = "hello"
    plot_class = plotbase.PlotBase
    reg.register(obj, "test1", plot_class, setup_fn=setup_plot)
    result = reg.get(obj, "test1")
    assert result == (plot_class, obj)


class DummyPlot(plotbase.PlotBase):
    def __init__(self):
        super().__init__(None)

    def annotate(self, x, y, label: str):
        return

    def resize(self, height: int = -1, width: int = -1):
        return

    def save(self, destination: str):
        return 'filename'

    def show(self, in_notebook=True):
        return


@pytest.fixture
def plotdf():
    return pd.DataFrame({'x': [1, 2], 'y1': [1, 2], 'y2': [1, 2]})


def test_validate(plotdf):
    # reject
    for bad_args, bad_kwargs in [
        ((None, None, None), {}),
        ((None, None, None), {'legend': 'i_am'}),
        ((None, 'x', []), {}),
        ((plotdf, 'x', ['z']), {}),
        ((plotdf, 'foo', ['y1']), {}),
    ]:
        plot = DummyPlot()
        result, msg = plot.validate(*bad_args, **bad_kwargs)
        assert result is False
    # pass
    for good_args, good_kwargs in [
        ((plotdf, 'x', ['y1']), {}),
        ((plotdf, 'x', ['y1', 'y2']), {}),
        ((plotdf, 'x', ['y1']), {'legend': 'foo'}),
    ]:
        plot = DummyPlot()
        result, msg = plot.validate(*good_args, **good_kwargs)
        assert result is True
