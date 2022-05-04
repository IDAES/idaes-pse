#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Tests the pyomo model created by the NetworkData parser
"""

import pytest
import os
import unittest
from pyomo.common.fileutils import this_file_dir
from urllib.request import urlopen
from zipfile import ZipFile
from io import BytesIO
from socket import timeout
from idaes.models_extra.gas_distribution.parsers.gaslib_parser import parse_gaslib_network

pytest.importorskip("plotly", reason="plotly not available")

try:
    data = urlopen("https://gaslib.zib.de/download/testData/GasLib-11.zip", timeout=3)
    zipfile = ZipFile(BytesIO(data.read()))
    files = [zipfile.open(file) for file in zipfile.namelist()]
    netfile = files[1]
    scnfile = files[2]
    net = parse_gaslib_network(netfile, scnfile)
    m = net.build_model(with_bounds=False)

except timeout:
    raise unittest.SkipTest('An internet connection is required to test gaslib parsing tools')

@pytest.mark.unit
def test_nd_model():
    nodes = 11
    pipes = 8
    compressors = 2

    # Check number of data objects are correct
    assert (len(list(m.fs.nodes)) == nodes)
    assert (len(list(m.fs.pipelines)) == pipes)
    assert (len(list(m.fs.compressors)) == compressors)

    # Get graph nodes to check connections
    e1 = m.fs.nodes['entry01']
    e2 = m.fs.nodes['entry02']
    e3 = m.fs.nodes['entry03']
    n1 = m.fs.nodes['N01']
    n2 = m.fs.nodes['N02']
    n3 = m.fs.nodes['N03']
    n4 = m.fs.nodes['N04']
    n5 = m.fs.nodes['N05']
    x1 = m.fs.nodes['exit01']
    x2 = m.fs.nodes['exit02']
    x3 = m.fs.nodes['exit03']

    # Check correct number of connections
    assert(len(e1.inlets) == 0 and len(e1.outlets) == 1)
    assert(len(e2.inlets) == 0 and len(e2.outlets) == 1)
    assert(len(e3.inlets) == 1 and len(e3.outlets) == 1)
    assert(len(n1.inlets) == 1 and len(n1.outlets) == 2)
    assert(len(n2.inlets) == 1 and len(n2.outlets) == 2)
    assert(len(n3.inlets) == 2 and len(n3.outlets) == 1)
    assert(len(n4.inlets) == 2 and len(n4.outlets) == 1)
    assert(len(n5.inlets) == 1 and len(n5.outlets) == 2)
    assert(len(x1.inlets) == 1 and len(x1.outlets) == 0)
    assert(len(x2.inlets) == 1 and len(x2.outlets) == 0)
    assert(len(x3.inlets) == 1 and len(x3.outlets) == 0)

    # e1_e3 connection
    e1_e3 = m.fs.pipelines['pipe01_entry01_entry03']
    assert(e1.outlets[0].arc.ports[0] == e1_e3.inlet_port)
    assert(e3.inlets[0].arc.ports[0] == e1_e3.outlet_port)
    assert (e1.outlets[0].arc.ports[1] == e1.outlets[0].port and e3.inlets[0].arc.ports[1] == e3.inlets[0].port)

    # e3_n1 connection - compressor
    e3_n1 = m.fs.compressors['CS01_entry03_N01']
    assert(e3.outlets[0].arc.ports[0] == e3_n1.inlet_port)
    assert(n1.inlets[0].arc.ports[0] == e3_n1.outlet_port)
    assert(e3.outlets[0].arc.ports[1] == e3.outlets[0].port and n1.inlets[0].arc.ports[1] == n1.inlets[0].port)


    # n1_n2 connection - pipe
    n1_n2 = m.fs.pipelines['pipe02_N01_N02']
    assert (n1.outlets[0].arc.ports[0] == n1_n2.inlet_port)
    assert (n2.inlets[0].arc.ports[0] == n1_n2.outlet_port)

    # n1_n3 connection - direct connection
    assert(n1.outlets[1].arc.ports[0] == n1.outlets[1].port)
    assert(n1.outlets[1].arc.ports[1] == n3.inlets[1].port)

    # e2_n3 connection - pipe
    e2_n3 = m.fs.pipelines['pipe03_entry02_N03']
    assert(e2.outlets[0].arc.ports[0] == e2_n3.inlet_port)
    assert(n3.inlets[0].arc.ports[0] == e2_n3.outlet_port)

    # n2_x1 connection - pipe
    n2_x1 = m.fs.pipelines['pipe04_N02_exit01']
    assert(n2.outlets[0].arc.ports[0] == n2_x1.inlet_port)
    assert(x1.inlets[0].arc.ports[0] == n2_x1.outlet_port)

    # n2_n4 connection - pipe
    n2_n4 = m.fs.pipelines['pipe05_N02_N04']
    assert(n2.outlets[1].arc.ports[0] == n2_n4.inlet_port)
    assert(n4.inlets[0].arc.ports[0] == n2_n4.outlet_port)

    # n3_n4 connection - pipe
    n3_n4 = m.fs.pipelines['pipe06_N03_N04']
    assert(n3.outlets[0].arc.ports[0] == n3_n4.inlet_port)
    assert(n4.inlets[1].arc.ports[0] == n3_n4.outlet_port)

    # n4_n5 connection - compressor
    n4_n5 = m.fs.compressors['CS02_N04_N05']
    assert(n4.outlets[0].arc.ports[0] == n4_n5.inlet_port)
    assert(n5.inlets[0].arc.ports[0] == n4_n5.outlet_port)

    # n5_x2 connection - pipe
    n5_x2 = m.fs.pipelines['pipe07_N05_exit02']
    assert(n5.outlets[0].arc.ports[0] == n5_x2.inlet_port)
    assert(x2.inlets[0].arc.ports[0] == n5_x2.outlet_port)

    # n5_x3 connection - pipe
    n5_x3 = m.fs.pipelines['pipe08_N05_exit03']
    assert(n5.outlets[1].arc.ports[0] == n5_x3.inlet_port)
    assert(x3.inlets[0].arc.ports[0] == n5_x3.outlet_port)

@pytest.mark.unit
def test_nd_plot():
    net.plot_results_at_time(m, m.fs.time.first(), show_plot=False)



