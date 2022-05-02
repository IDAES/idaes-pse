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

import pyomo.environ as pe
from idaes.models_extra.gas_distribution.data import network_data as nd
from xml.etree import ElementTree
from typing import MutableMapping, Optional
import logging
from pyomo.core.expr.numeric_expr import NPV_ProductExpression

logger = logging.getLogger(__name__)

def _extract_expression_with_units(attrib: MutableMapping):
    val = float(attrib["value"])
    unit_str = attrib["unit"]
    try:
        unit = getattr(pe.units, unit_str)
    except (AttributeError, ValueError) as e:
        if unit_str == "1000m_cube_per_hour":
            val *= 1000
            unit = pe.units.m**3 / pe.units.hr
        elif unit_str == "Celsius":
            val += 273.15
            unit = pe.units.K
        else:
            raise e

    return NPV_ProductExpression((val, unit))


def process_source(xml_node: ElementTree.Element, net: nd.NetworkData):
    node_attrib = xml_node.attrib
    node_name = node_attrib["id"]
    x = float(node_attrib["x"])
    y = float(node_attrib["y"])
    min_pressure = xml_node.find("{http://gaslib.zib.de/Gas}pressureMin").attrib
    max_pressure = xml_node.find("{http://gaslib.zib.de/Gas}pressureMax").attrib
    min_flow = xml_node.find("{http://gaslib.zib.de/Gas}flowMin").attrib
    max_flow = xml_node.find("{http://gaslib.zib.de/Gas}flowMax").attrib
    temp = xml_node.find("{http://gaslib.zib.de/Gas}gasTemperature").attrib

    source = nd.Source(
        name=node_name,
        x=x,
        y=y,
        min_pressure=_extract_expression_with_units(min_pressure),
        max_pressure=_extract_expression_with_units(max_pressure),
        min_flow=_extract_expression_with_units(min_flow),
        max_flow=_extract_expression_with_units(max_flow),
        temperature=_extract_expression_with_units(temp),
    )

    net.add_node(node_name, source)


def process_sink(xml_node: ElementTree.Element, net: nd.NetworkData):
    node_attrib = xml_node.attrib
    node_name = node_attrib["id"]
    x = float(node_attrib["x"])
    y = float(node_attrib["y"])
    min_pressure = xml_node.find("{http://gaslib.zib.de/Gas}pressureMin").attrib
    max_pressure = xml_node.find("{http://gaslib.zib.de/Gas}pressureMax").attrib
    min_flow = xml_node.find("{http://gaslib.zib.de/Gas}flowMin").attrib
    max_flow = xml_node.find("{http://gaslib.zib.de/Gas}flowMax").attrib

    sink = nd.Sink(
        name=node_name,
        x=x,
        y=y,
        min_pressure=_extract_expression_with_units(min_pressure),
        max_pressure=_extract_expression_with_units(max_pressure),
        target_demand=None,
    )

    net.add_node(node_name, sink)


def process_innode(xml_node: ElementTree.Element, net: nd.NetworkData):
    node_attrib = xml_node.attrib
    node_name = node_attrib["id"]
    x = float(node_attrib["x"])
    y = float(node_attrib["y"])
    min_pressure = xml_node.find("{http://gaslib.zib.de/Gas}pressureMin").attrib
    max_pressure = xml_node.find("{http://gaslib.zib.de/Gas}pressureMax").attrib

    node = nd.Node(
        name=node_name,
        x=x,
        y=y,
        min_pressure=_extract_expression_with_units(min_pressure),
        max_pressure=_extract_expression_with_units(max_pressure),
    )

    net.add_node(node_name, node)


def process_pipe(
    xml_node: ElementTree.Element,
    net: nd.NetworkData,
    use_pipe_roughness_as_rugosity: bool,
):
    node_attrib = xml_node.attrib
    link_name = node_attrib["id"]
    from_node_name = node_attrib["from"]
    to_node_name = node_attrib["to"]
    length = xml_node.find("{http://gaslib.zib.de/Gas}length").attrib
    diameter = xml_node.find("{http://gaslib.zib.de/Gas}diameter").attrib
    if use_pipe_roughness_as_rugosity:
        rugosity = _extract_expression_with_units(
            xml_node.find("{http://gaslib.zib.de/Gas}roughness").attrib
        )
    else:
        rugosity = 0.025 * pe.units.mm

    link = nd.Pipe(
        name=link_name,
        from_node_name=from_node_name,
        to_node_name=to_node_name,
        length=_extract_expression_with_units(length),
        diameter=_extract_expression_with_units(diameter),
        rugosity=rugosity,
    )

    net.add_link(link_name, link)


def process_short_pipe(xml_node: ElementTree.Element, net: nd.NetworkData):
    node_attrib = xml_node.attrib
    link_name = node_attrib["id"]
    from_node_name = node_attrib["from"]
    to_node_name = node_attrib["to"]

    link = nd.ShortPipe(
        name=link_name, from_node_name=from_node_name, to_node_name=to_node_name
    )

    net.add_link(link_name, link)


def process_valve(xml_node: ElementTree.Element, net: nd.NetworkData):
    node_attrib = xml_node.attrib
    link_name = node_attrib["id"]
    from_node_name = node_attrib["from"]
    to_node_name = node_attrib["to"]

    link = nd.ShortPipe(
        name=link_name, from_node_name=from_node_name, to_node_name=to_node_name
    )

    net.add_link(link_name, link)


def process_compressor(xml_node: ElementTree.Element, net: nd.NetworkData):
    node_attrib = xml_node.attrib
    link_name = node_attrib["id"]
    from_node_name = node_attrib["from"]
    to_node_name = node_attrib["to"]
    min_flow = xml_node.find("{http://gaslib.zib.de/Gas}flowMin").attrib
    max_flow = xml_node.find("{http://gaslib.zib.de/Gas}flowMax").attrib
    min_inlet_pressure = xml_node.find("{http://gaslib.zib.de/Gas}pressureInMin").attrib
    max_outlet_pressure = xml_node.find(
        "{http://gaslib.zib.de/Gas}pressureOutMax"
    ).attrib

    link = nd.Compressor(
        name=link_name,
        from_node_name=from_node_name,
        to_node_name=to_node_name,
        min_flow=_extract_expression_with_units(min_flow),
        max_flow=_extract_expression_with_units(max_flow),
        min_inlet_pressure=_extract_expression_with_units(min_inlet_pressure),
        max_outlet_pressure=_extract_expression_with_units(max_outlet_pressure),
    )

    net.add_link(link_name, link)


def process_scenario_demand(xml_node: ElementTree.Element, net: nd.NetworkData):
    assert xml_node.attrib["type"] == "exit"
    node_name = xml_node.attrib["id"]
    node: nd.Sink = net.get_node(node_name)
    assert node.is_sink()
    if len(xml_node) == 1:
        xml_flow = xml_node[0]
        assert xml_flow.attrib["bound"].lower() == "both"
    else:
        assert len(xml_node) == 2
        xml_flow = xml_node[0]
        other_xml_flow = xml_node[1]
        assert (
            xml_flow.attrib["bound"].lower() == "lower"
            and other_xml_flow.attrib["bound"].lower() == "upper"
        ) or (
            xml_flow.attrib["bound"].lower() == "upper"
            and other_xml_flow.attrib["bound"].lower() == "lower"
        )
        if xml_flow.attrib["value"] != other_xml_flow.attrib["value"]:
            raise RuntimeError("scenario demand must be fixed")
        if xml_flow.attrib["unit"] != other_xml_flow.attrib["unit"]:
            raise RuntimeError("scenario demand must be fixed")
    target_demand = _extract_expression_with_units(xml_flow.attrib)
    node.target_demand = target_demand


def parse_gaslib_network(
    net_filename: str,
    demand_scenario_filename: str,
    use_pipe_roughness_as_rugosity: bool = True,
    net: Optional[nd.NetworkData] = None,
) -> nd.NetworkData:
    tree = ElementTree.parse(net_filename)
    root = tree.getroot()
    if net is None:
        net = nd.NetworkData()

    nodes = root.find("{http://gaslib.zib.de/Framework}nodes")
    assert nodes is not None

    for n in nodes:
        n_tag = n.tag.split("}")[1]
        if n_tag == "source":
            process_source(n, net)
        elif n_tag == "sink":
            process_sink(n, net)
        elif n_tag == "innode":
            process_innode(n, net)
        else:
            raise ValueError(f"Unrecognized type of node: {n_tag}")

    links = root.find("{http://gaslib.zib.de/Framework}connections")
    assert links is not None

    for l in links:
        l_tag = l.tag.split("}")[1]
        if l_tag == "pipe":
            process_pipe(
                l, net, use_pipe_roughness_as_rugosity=use_pipe_roughness_as_rugosity
            )
        elif l_tag == "compressorStation":
            process_compressor(l, net)
        elif l_tag == "shortPipe":
            process_short_pipe(l, net)
        elif l_tag == "valve":
            logger.warning("valves are currently fixed to always be open.")
            process_valve(l, net)
        else:
            raise ValueError(f"Unrecognized type of link: {l_tag}")

    tree = ElementTree.parse(demand_scenario_filename)
    root = tree.getroot()
    scenario = root.find("{http://gaslib.zib.de/Gas}scenario")
    for child in scenario:
        child_tag = child.tag.split("}")[1]
        assert child_tag == "node"
        if child.attrib["type"] == "entry":
            continue
        process_scenario_demand(child, net)

    return net
