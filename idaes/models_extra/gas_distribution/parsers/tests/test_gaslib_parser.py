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
Tests the parsing of gaslib XML files
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

try:
    data = urlopen("https://gaslib.zib.de/download/testData/GasLib-11.zip", timeout=3)
    zipfile = ZipFile(BytesIO(data.read()))
    files = [zipfile.open(file) for file in zipfile.namelist()]
    netfile = files[1]
    scnfile = files[2]

except timeout:
    raise unittest.SkipTest('An internet connection is required to test gaslib parsing tools')

@pytest.mark.unit
def test_network_generation():
    net = parse_gaslib_network(netfile, scnfile)
    nodes = 11
    sources = 3
    sinks = 3
    pipes = 8
    compressors = 2
    shortpipes = 1

    # Check number of data objects are correct
    assert(len(list(net.nodes())) == nodes)
    assert(len(list(net.sources())) == sources)
    assert(len(list(net.sinks())) == sinks)
    assert(len(list(net.pipes())) == pipes)
    assert(len(list(net.compressors())) == compressors)
    assert(len(list(net.short_pipes())) == shortpipes)

    # Get connections of the graph
    e1_e3 = net.get_link('pipe01_entry01_entry03')
    e2_n3 = net.get_link('pipe03_entry02_N03')
    e3_n1 = net.get_link('CS01_entry03_N01')

    n1_n2 = net.get_link('pipe02_N01_N02')
    n1_n3 = net.get_link('V01_N01_N03')
    n2_n4 = net.get_link('pipe05_N02_N04')
    n3_n4 = net.get_link('pipe06_N03_N04')
    n4_n5 = net.get_link('CS02_N04_N05')

    n2_x1 = net.get_link('pipe04_N02_exit01')
    n5_x2 = net.get_link('pipe07_N05_exit02')
    n5_x3 = net.get_link('pipe08_N05_exit03')

    # Check entrance connections and flow direction
    assert(e1_e3.from_node_name == 'entry01' and e1_e3.to_node_name == 'entry03')
    assert(e2_n3.from_node_name == 'entry02' and e2_n3.to_node_name == 'N03')
    assert(e3_n1.from_node_name == 'entry03' and e3_n1.to_node_name == 'N01')

    # Check node to node connections and flow direction
    assert(n1_n2.from_node_name == 'N01' and n1_n2.to_node_name == 'N02')
    assert(n1_n3.from_node_name == 'N01' and n1_n3.to_node_name == 'N03')
    assert(n2_n4.from_node_name == 'N02' and n2_n4.to_node_name == 'N04')
    assert(n3_n4.from_node_name == 'N03' and n3_n4.to_node_name == 'N04')
    assert(n4_n5.from_node_name == 'N04' and n4_n5.to_node_name == 'N05')

    # Check exit connections and flow direction
    assert(n2_x1.from_node_name == 'N02' and n2_x1.to_node_name == 'exit01')
    assert(n5_x2.from_node_name == 'N05' and n5_x2.to_node_name == 'exit02')
    assert(n5_x3.from_node_name == 'N05' and n5_x3.to_node_name == 'exit03')






