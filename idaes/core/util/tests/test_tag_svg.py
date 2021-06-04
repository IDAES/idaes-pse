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
This module contains miscalaneous utility functions for use in IDAES models.
"""
__author__ = "John Eslick"
import pytest

from pyomo.environ import ConcreteModel, Set, Block, Var
from pyomo.network import Port

from idaes.core.util.misc import TagReference, svg_tag

# simple svg from drawn in incscape, with text ids:
#   TAGME_4.x, TAGME_4.y, and TAGME_4.f
svg_test_str = \
"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- Created with Inkscape (http://www.inkscape.org/) -->

<svg
   xmlns:dc="http://purl.org/dc/elements/1.1/"
   xmlns:cc="http://creativecommons.org/ns#"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
   width="210mm"
   height="297mm"
   viewBox="0 0 210 297"
   version="1.1"
   id="svg2045"
   inkscape:version="0.92.4 (5da689c313, 2019-01-14)"
   sodipodi:docname="test.svg">
  <defs
     id="defs2039" />
  <sodipodi:namedview
     id="base"
     pagecolor="#ffffff"
     bordercolor="#666666"
     borderopacity="1.0"
     inkscape:pageopacity="0.0"
     inkscape:pageshadow="2"
     inkscape:zoom="0.98994949"
     inkscape:cx="174.80933"
     inkscape:cy="316.53691"
     inkscape:document-units="mm"
     inkscape:current-layer="layer1"
     showgrid="false"
     inkscape:window-width="3840"
     inkscape:window-height="2035"
     inkscape:window-x="-13"
     inkscape:window-y="-13"
     inkscape:window-maximized="1" />
  <metadata
     id="metadata2042">
    <rdf:RDF>
      <cc:Work
         rdf:about="">
        <dc:format>image/svg+xml</dc:format>
        <dc:type
           rdf:resource="http://purl.org/dc/dcmitype/StillImage" />
        <dc:title></dc:title>
      </cc:Work>
    </rdf:RDF>
  </metadata>
  <g
     inkscape:label="Layer 1"
     inkscape:groupmode="layer"
     id="layer1">
    <ellipse
       style="stroke-width:0.26499999;stroke-miterlimit:4;stroke-dasharray:none"
       id="path2590"
       cx="113.18864"
       cy="130.22382"
       rx="30.869631"
       ry="31.003265" />
    <path
       style="fill:none;stroke:#000000;stroke-width:0.26458332px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       d="m 45.70309,75.433559 37.685,33.943231"
       id="path2592"
       inkscape:connector-curvature="0" />
    <path
       style="fill:none;stroke:#000000;stroke-width:0.26458332px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       d="m 136.04019,155.61442 49.17759,45.16855"
       id="path2594"
       inkscape:connector-curvature="0" />
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:normal;font-size:10.58333302px;line-height:1.25;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.26458332"
       x="63.610146"
       y="68.751816"
       id="text2598"><tspan
         sodipodi:role="line"
         id="tspan2596"
         x="63.610146"
         y="68.751816"
         style="stroke-width:0.26458332">x:</tspan></text>
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:normal;font-size:10.58333302px;line-height:1.25;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.26458332"
       x="66.282845"
       y="86.391609"
       id="text2602"><tspan
         sodipodi:role="line"
         x="66.282845"
         y="86.391609"
         style="stroke-width:0.26458332"
         id="tspan2604">y:</tspan></text>
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:normal;font-size:10.58333302px;line-height:1.25;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.26458332"
       x="166.77618"
       y="171.11606"
       id="text2610"><tspan
         sodipodi:role="line"
         id="tspan2608"
         x="166.77618"
         y="171.11606"
         style="stroke-width:0.26458332">f:</tspan></text>
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:normal;font-size:10.58333302px;line-height:1.25;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.26458332"
       x="77.723328"
       y="68.57209"
       id="TAGME_4.x"
       inkscape:label="#text2598-9"><tspan
         sodipodi:role="line"
         id="tspan2596-1"
         x="77.723328"
         y="68.57209"
         style="stroke-width:0.26458332">1</tspan></text>
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:normal;font-size:10.58333302px;line-height:1.25;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.26458332"
       x="80.396027"
       y="86.211876"
       id="TAGME_4.y"
       inkscape:label="#text2602-0"><tspan
         sodipodi:role="line"
         x="80.396027"
         y="86.211876"
         style="stroke-width:0.26458332"
         id="tspan2604-1">1</tspan></text>
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:normal;font-size:10.58333302px;line-height:1.25;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.26458332"
       x="176.65611"
       y="171.99467"
       id="TAGME_4.f"
       inkscape:label="#text2610-0"><tspan
         sodipodi:role="line"
         id="tspan2608-1"
         x="176.65611"
         y="171.99467"
         style="stroke-width:0.26458332"></tspan></text>
    <flowRoot
       xml:space="preserve"
       id="flowRoot2642"
       style="fill:black;fill-opacity:1;stroke:none;font-family:sans-serif;font-style:normal;font-weight:normal;font-size:40px;line-height:1.25;letter-spacing:0px;word-spacing:0px"><flowRegion
         id="flowRegion2644"><rect
           id="rect2646"
           width="533.36053"
           height="577.80725"
           x="714.17786"
           y="232.57529" /></flowRegion><flowPara
         id="flowPara2648"></flowPara></flowRoot>  </g>
</svg>"""


@pytest.mark.unit
def test_tag_reference():
    m = ConcreteModel()
    m.x = Var([0,1], initialize={0:2.22,1:3.33})
    m.y = Var([0,1], initialize={0:4.44,1:5.55})
    m.f = Var([0,1], initialize={0:6.66,1:7.77})
    test_tag = {}
    test_tag["TAGME@4.x"] = TagReference(m.x[:], description="x tag")
    test_tag["TAGME@4.y"] = TagReference(m.y[:], description="y tag")
    test_tag["TAGME@4.f"] = TagReference(m.f[:], description="z tag")
    m.tag = test_tag

    xml_str = svg_tag(m.tag, svg_test_str, idx=0)
    # lazy testing
    assert("2.2200e" in xml_str)
    assert("4.4400e" in xml_str)
    assert("6.6600e" in xml_str)

    xml_str = svg_tag(m.tag, svg_test_str, idx=1)
    # lazy testing
    assert("3.3300e" in xml_str)
    assert("5.5500e" in xml_str)
    assert("7.7700e" in xml_str)

    xml_str = svg_tag(m.tag, svg_test_str, show_tags=True)
    # lazy testing
    assert("TAGME@4.x" in xml_str)
    assert("TAGME@4.y" in xml_str)
    assert("TAGME@4.f" in xml_str)

    tag_data_like = {}
    tag_data_like["TAGME@4.x"] = 1.1212
    tag_data_like["TAGME@4.y"] = 2.1212
    tag_data_like["TAGME@4.f"] = "3.1212 Hello"
    xml_str = svg_tag(tag_data_like, svg_test_str, idx=None)
    assert("1.1212" in xml_str)
    assert("2.1212" in xml_str)
    assert("3.1212 Hello" in xml_str)


@pytest.mark.unit
def test_tag_reference_default_format():
    m = ConcreteModel()
    m.x = Var([0,1], initialize={0:2.22,1:3.33})
    m.y = Var([0,1], initialize={0:4.44,1:5.55})
    m.f = Var([0,1], initialize={0:6.66,1:7.77})
    test_tag = {}
    test_tag["TAGME@4.x"] = TagReference(m.x[:], description="x tag")
    test_tag["TAGME@4.y"] = TagReference(m.y[:], description="y tag")
    test_tag["TAGME@4.f"] = TagReference(m.f[:], description="z tag")
    m.tag = test_tag

    xml_str = svg_tag(m.tag, svg_test_str, idx=0, tag_format_default="{:.2f}X")
    # lazy testing
    assert("2.22X" in xml_str)
    assert("4.44X" in xml_str)
    assert("6.66X" in xml_str)


@pytest.mark.unit
def test_tag_reference_tag_format():
    m = ConcreteModel()
    m.x = Var([0,1], initialize={0:2.22,1:3.33})
    m.y = Var([0,1], initialize={0:4.44,1:5.55})
    m.f = Var([0,1], initialize={0:6.66,1:7.77})
    test_tag = {}
    test_tag["TAGME@4.x"] = TagReference(m.x[:], description="x tag")
    test_tag["TAGME@4.y"] = TagReference(m.y[:], description="y tag")
    test_tag["TAGME@4.f"] = TagReference(m.f[:], description="z tag")
    m.tag = test_tag

    xml_str = svg_tag(m.tag, svg_test_str, idx=0, tag_format={"TAGME@4.y":"{:.5e}"})
    # lazy testing
    assert("2.2200e" in xml_str)
    assert("4.44000e" in xml_str)
    assert("6.6600e" in xml_str)


@pytest.mark.unit
def test_tag_reference_tag_format_conditional():
    m = ConcreteModel()
    m.x = Var([0,1], initialize={0:10000,1:3.33})
    m.y = Var([0,1], initialize={0:1.1,1:5.55})
    m.f = Var([0,1], initialize={0:4.5,1:7.77})
    test_tag = {}
    test_tag["TAGME@4.x"] = TagReference(m.x[:], description="x tag")
    test_tag["TAGME@4.y"] = TagReference(m.y[:], description="y tag")
    test_tag["TAGME@4.f"] = TagReference(m.f[:], description="z tag")
    m.tag = test_tag


    xml_str = svg_tag(
        m.tag,
        svg_test_str,
        idx=0,
        tag_format_default=lambda x: "{:,.0f} kPa" if x >= 100 else "{:.2f} kPa")
    # lazy testing
    assert("10,000 kPa" in xml_str)
    assert("1.10 kPa" in xml_str)
    assert("4.50 kPa" in xml_str)
