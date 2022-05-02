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
from typing import Optional, MutableMapping, MutableSet, Generator
from pyomo.core.expr.numeric_expr import ExpressionBase
import logging
from pyomo.common.collections.orderedset import OrderedSet
from pyomo.core.base.block import _BlockData, ScalarBlock
import idaes.core

from idaes.models_extra.gas_distribution.properties.natural_gas import NaturalGasParameterBlock
from idaes.models_extra.gas_distribution.unit_models.pipeline import GasPipeline
from idaes.models_extra.gas_distribution.unit_models.compressor import IsothermalCompressor
from idaes.models_extra.gas_distribution.unit_models.node import PipelineNode

from pyomo.network.arc import Arc
import math
from pyomo.util.check_units import assert_units_consistent
import plotly.graph_objects as go


logger = logging.getLogger(__name__)


class Node(object):
    def __init__(
        self,
        name,
        x=None,
        y=None,
        min_pressure=1 * pe.units.bar,
        max_pressure=None,
    ):
        self.name: str = name
        self.x: Optional[float] = x
        self.y: Optional[float] = y
        self.min_pressure: Optional[ExpressionBase] = min_pressure
        self.max_pressure: Optional[ExpressionBase] = max_pressure

    def is_source(self):
        return False

    def is_sink(self):
        return False


class Source(Node):
    def __init__(
        self,
        name,
        x=None,
        y=None,
        min_pressure=1 * pe.units.bar,
        max_pressure=None,
        min_flow=0 * pe.units.m ** 3 / pe.units.hr,
        max_flow=None,
        temperature=298.15 * pe.units.K,
    ):
        super().__init__(
            name,
            x=x,
            y=y,
            min_pressure=min_pressure,
            max_pressure=max_pressure,
        )
        self.min_flow: Optional[ExpressionBase] = min_flow
        self.max_flow: Optional[ExpressionBase] = max_flow
        self.temperature: ExpressionBase = temperature

    def is_source(self):
        return True


class Sink(Node):
    def __init__(
        self,
        name,
        x=None,
        y=None,
        min_pressure=1 * pe.units.bar,
        max_pressure=None,
        target_demand=None
    ):
        super().__init__(
            name,
            x=x,
            y=y,
            min_pressure=min_pressure,
            max_pressure=max_pressure,
        )
        self.target_demand: Optional[ExpressionBase] = target_demand

    def is_sink(self):
        return True


class Link(object):
    def __init__(self, name, from_node_name, to_node_name):
        self.name: str = name
        self.from_node_name: str = from_node_name
        self.to_node_name: str = to_node_name

    def is_pipe(self):
        return False

    def is_compressor(self):
        return False

    def is_short_pipe(self):
        return False


class Pipe(Link):
    def __init__(
        self,
        name,
        from_node_name,
        to_node_name,
        length=100 * pe.units.km,
        diameter=1 * pe.units.m,
        rugosity=0.025 * pe.units.mm,
    ):
        super().__init__(name, from_node_name, to_node_name)
        self.length: ExpressionBase = length
        self.diameter: ExpressionBase = diameter
        self.rugosity: ExpressionBase = rugosity

    def is_pipe(self):
        return True


class ShortPipe(Link):
    def is_short_pipe(self):
        return True


class Compressor(Link):
    def __init__(
        self,
        name,
        from_node_name,
        to_node_name,
        min_flow=None,
        max_flow=None,
        min_inlet_pressure=None,
        max_outlet_pressure=None,
    ):
        super().__init__(name, from_node_name, to_node_name)
        self.min_flow: Optional[ExpressionBase] = min_flow
        self.max_flow: Optional[ExpressionBase] = max_flow
        self.min_inlet_pressure: Optional[ExpressionBase] = min_inlet_pressure
        self.max_outlet_pressure: Optional[ExpressionBase] = max_outlet_pressure

    def is_compressor(self):
        return True


def _get_link_txt(l: Link, m: _BlockData, t: float):
    txt = f'{l.name}<br>'
    if l.is_pipe():
        txt += f'Type: Pipe<br>'
    elif l.is_compressor():
        txt += f'Type: Compressor<br>'
    elif l.is_short_pipe():
        txt += f'Type: Short Pipe<br>'
    else:
        raise ValueError(f'Unrecognized link type: {str(type(l))}')
    return txt


def _get_pipe_txt(l: Pipe, m: _BlockData, t: float):
    cv = m.fs.pipelines[l.name].control_volume
    ld = cv.length_domain
    f_in = cv.properties[t,ld.first()].flow_mol
    f_out = cv.properties[t, ld.last()].flow_mol
    f_units = pe.units.get_units(f_in)
    txt = f'Flow In: {f_in.value:<.3f} {f_units}<br>'
    txt += f'Flow Out: {f_out.value:<.3f} {f_units}<br>'
    dens_in = cv.properties[t, ld.first()].dens_mol
    txt += f'Inlet Density: {dens_in.value:<.3f} {pe.units.get_units(dens_in)}<br>'
    dens_out = cv.properties[t, ld.last()].dens_mol
    txt += f'Outlet Density: {dens_out.value:<.3f} {pe.units.get_units(dens_out)}<br>'
    return txt


def _get_short_pipe_txt(l: ShortPipe, m: _BlockData, t: float):
    from_node = m.fs.nodes[l.from_node_name]
    outlet_idx = None
    for tmp_outlet_idx, tmp_outlet_link in from_node._outlet_pipelines.items():
        if tmp_outlet_link == l.name:
            outlet_idx = tmp_outlet_idx
    assert outlet_idx is not None
    f = from_node.outlets[outlet_idx].state[t].flow_mol
    f_units = pe.units.get_units(f)
    txt = f'Flow: {f.value:<.3f} {f_units}<br>'
    return txt


def _get_compressor_txt(l: Compressor, m: _BlockData, t: float):
    c = m.fs.compressors[l.name]
    f_in = c.inlet_state[t].flow_mol
    f_out = c.outlet_state[t].flow_mol
    f_units = pe.units.get_units(f_in)
    min_flow = l.min_flow * m.fs.properties.dens_nominal / m.fs.properties.natural_gas.mw
    max_flow = l.max_flow * m.fs.properties.dens_nominal / m.fs.properties.natural_gas.mw
    min_flow = pe.units.convert(min_flow, f_units)
    max_flow = pe.units.convert(max_flow, f_units)
    txt = f'Flow: {f_in.value:<.3f} {f_units}<br>'
    txt += f'Min Flow: {pe.value(min_flow):<.3f} {pe.units.get_units(min_flow)}<br>'
    txt += f'Max Flow: {pe.value(max_flow):<.3f} {pe.units.get_units(max_flow)}<br>'
    txt += f'Inlet Density: {pe.value(c.inlet_state[t].dens_mol):<.3f} ' \
           f'{pe.units.get_units(c.inlet_state[t].dens_mol)}<br>'

    bp = m.fs.compressors[l.name].boost_pressure[t]
    bp_units = pe.units.get_units(bp)
    if bp.has_lb():
        bp_lb = bp.lb
    else:
        bp_lb = -math.inf
    if bp.has_ub():
        bp_ub = bp.ub
    else:
        bp_ub = math.inf
    txt += f'Boost Pressure: {bp.value:<.3f} {bp_units}'
    if bp.fixed:
        txt += ', fixed'
    txt += '<br>'
    txt += f'Min Boost Pressure: {bp_lb:<.3f} {bp_units}<br>'
    txt += f'Max Boost Pressure: {bp_ub:<.3f} {bp_units}<br>'

    p_in = m.fs.compressors[l.name].inlet_state[t].pressure
    p_out = m.fs.compressors[l.name].outlet_state[t].pressure
    p_units = pe.units.get_units(p_in)
    if p_in.has_lb():
        p_lb = p_in.lb
    else:
        p_lb = -math.inf
    if p_out.has_ub():
        p_ub = p_out.ub
    else:
        p_ub = math.inf
    txt += f'Min Pressure In: {p_lb:<.3f} {p_units}<br>'
    txt += f'Max Pressure Out: {p_ub:<.3f} {p_units}<br>'
    return txt


def _get_node_txt(n: Node, m: _BlockData, t: float):
    txt = f'{n.name}<br>'
    if n.is_source():
        txt += f'Type: Supply<br>'
    elif n.is_sink():
        txt += f'Type: Demand<br>'
    else:
        txt += f'Type: Node<br>'

    p = m.fs.nodes[n.name].state[t].pressure
    p_units = str(pe.units.get_units(p))
    txt += f'Pressure: {p.value:<.3f} {p_units}'
    if p.fixed or (n.is_source() and m.fs.nodes[n.name].supplies[0].state[t].pressure.fixed):
        txt += ', fixed'
    txt += '<br>'
    if p.has_lb():
        p_lb = p.lb
    else:
        p_lb = -math.inf
    if p.has_ub():
        p_ub = p.ub
    else:
        p_ub = math.inf
    if n.is_source():
        sp = m.fs.nodes[n.name].supplies[0].state[t].pressure
        if sp.has_lb() and sp.lb > p_lb:
            p_lb = sp.lb
        if sp.has_ub() and sp.ub < p_ub:
            p_ub = sp.ub
    txt += f'Pmin: {p_lb} {p_units}<br>'
    txt += f'Pmax: {p_ub} {p_units}<br>'
    return txt


def _get_supply_txt(n: Source, m: _BlockData, t: float):
    txt = ''
    f = m.fs.nodes[n.name].supplies[0].flow_mol[t]
    f_units = pe.units.get_units(f)
    txt += f'Supply Flow: {f.value:<.3f} {str(f_units)}'
    if f.fixed:
        txt += ', fixed'
    txt += '<br>'
    min_flow = n.min_flow * m.fs.properties.dens_nominal / m.fs.properties.natural_gas.mw
    max_flow = n.max_flow * m.fs.properties.dens_nominal / m.fs.properties.natural_gas.mw
    min_flow = pe.units.convert(min_flow, f_units)
    max_flow = pe.units.convert(max_flow, f_units)
    txt += f'Min Supply: {pe.value(min_flow):<.3f} {pe.units.get_units(min_flow)}<br>'
    txt += f'Max Supply: {pe.value(max_flow):<.3f} {pe.units.get_units(max_flow)}<br>'
    return txt


def _get_demand_txt(n: Sink, m: _BlockData, t: float):
    txt = ''
    f = m.fs.nodes[n.name].demands[0].flow_mol[t]
    f_units = pe.units.get_units(f)
    txt += f'Demand Flow: {f.value:<.3f} {str(f_units)}'
    if f.fixed:
        txt += ', fixed'
    txt += '<br>'
    td = n.target_demand
    td *= m.fs.properties.dens_nominal / m.fs.properties.natural_gas.mw
    td = pe.units.convert(td, f_units)
    txt += f'Target Demand: {pe.value(td):<.3f}' \
           f' {pe.units.get_units(td)}<br>'
    return txt


class NetworkData(object):
    def __init__(self):
        self._nodes: MutableMapping[str, Node] = dict()
        self._links: MutableMapping[str, Link] = dict()
        self._node_inlets: MutableMapping[str, MutableSet[Link]] = dict()
        self._node_outlets: MutableMapping[str, MutableSet[Link]] = dict()

    def build_model(
        self,
        nfex: int = 5,
        xscheme: str = "FORWARD",
        dynamic: bool = False,
        nfet: Optional[int] = None,
        tscheme: str = "BACKWARD",
        horizon: Optional[float] = None,
        with_bounds: bool = True,
        fix_demands: bool = True,
    ) -> _BlockData:
        if dynamic:
            if horizon is None:
                raise ValueError('horizon must be specified in order to build a dynamic model.')
            if nfet is None:
                raise ValueError('nfet must be specified in order to build a dynamic model.')

        default = dict(dynamic=dynamic)
        if dynamic:
            default['time_set'] = [0, horizon]
            default['time_units'] = pe.units.hr

        m = ScalarBlock(concrete=True)
        m.fs = idaes.core.FlowsheetBlock(default=default)
        m.fs.properties = NaturalGasParameterBlock()

        self._build_node_models(m)
        self._build_pipeline_models(m, nfex, xscheme)
        self._connect_pipes_to_nodes(m)
        self._build_short_pipe_models(m)
        self._build_compressor_models(m)

        pe.TransformationFactory('network.expand_arcs').apply_to(m)
        assert_units_consistent(m)

        if dynamic:
            disc = pe.TransformationFactory('dae.finite_difference')
            disc.apply_to(m, nfe=nfet, wrt=m.fs.time, scheme=tscheme)

        self._fix_source_mole_frac_and_temp(m)

        if with_bounds:
            self.apply_compressor_flow_bounds(m)
            self.apply_compressor_pressure_bounds(m)
            self.apply_node_pressure_bounds(m)
            self.apply_source_flow_bounds(m)

        if fix_demands:
            self.fix_demand_flows(m)

        return m

    def _build_short_pipe_models(self, m: _BlockData):
        for sp in self.short_pipes():
            from_node = m.fs.nodes[sp.from_node_name]
            to_node = m.fs.nodes[sp.to_node_name]
            outlet_idx = from_node.n_outlet_pipelines
            inlet_idx = to_node.n_inlet_pipelines
            from_node.outlets[outlet_idx].arc = Arc(ports=(from_node.outlets[outlet_idx].port, to_node.inlets[inlet_idx].port))
            from_node.n_outlet_pipelines += 1
            to_node.n_inlet_pipelines += 1
            from_node.outlets[outlet_idx].has_pipeline = True
            to_node.inlets[inlet_idx].has_pipeline = True
            from_node._outlet_pipelines[outlet_idx] = sp.name
            to_node._inlet_pipelines[inlet_idx] = sp.name

    def _build_compressor_models(self, m: _BlockData):
        compressor_names = [c.name for c in self.compressors()]
        compressor_config = dict(property_package=m.fs.properties)
        m.fs.compressor_set = pe.Set(initialize=compressor_names)
        m.fs.compressors = IsothermalCompressor(m.fs.compressor_set, default=compressor_config)

        for c in self.compressors():
            m.fs.nodes[c.from_node_name].add_pipeline_to_outlet(m.fs.compressors[c.name])
            m.fs.nodes[c.to_node_name].add_pipeline_to_inlet(m.fs.compressors[c.name])

    def apply_compressor_flow_bounds(self, m: _BlockData):
        for c in self.compressors():
            mc = m.fs.compressors[c.name]
            for t in m.fs.time:
                flow_units = pe.units.get_units(mc.inlet_state[t].flow_mol)
                min_flow_mol = c.min_flow * m.fs.properties.dens_nominal / m.fs.properties.natural_gas.mw
                max_flow_mol = c.max_flow * m.fs.properties.dens_nominal / m.fs.properties.natural_gas.mw
                min_flow_mol = pe.units.convert(min_flow_mol, flow_units)
                max_flow_mol = pe.units.convert(max_flow_mol, flow_units)
                mc.inlet_state[t].flow_mol.setlb(pe.value(min_flow_mol))
                mc.inlet_state[t].flow_mol.setub(pe.value(max_flow_mol))

    def apply_compressor_pressure_bounds(self, m: _BlockData):
        t0 = m.fs.time.first()
        for c in self.compressors():
            pressure_units = pe.units.get_units(m.fs.compressors[c.name].inlet_state[t0].pressure)
            pmin = pe.units.convert(c.min_inlet_pressure, pressure_units)
            m.fs.compressors[c.name].inlet_state[:].pressure.setlb(pe.value(pmin))
            pmax = pe.units.convert(c.max_outlet_pressure, pressure_units)
            m.fs.compressors[c.name].outlet_state[:].pressure.setub(pe.value(pmax))

    def _connect_pipes_to_nodes(self, m: _BlockData):
        for p in self.pipes():
            m.fs.nodes[p.from_node_name].add_pipeline_to_outlet(m.fs.pipelines[p.name])
            m.fs.nodes[p.to_node_name].add_pipeline_to_inlet(m.fs.pipelines[p.name])

    def _build_pipeline_models(self, m: _BlockData, nfex: int, xscheme: str):
        pipe_names = [p.name for p in self.pipes()]
        pipeline_config = dict(property_package=m.fs.properties,
                               finite_elements=nfex,
                               transformation_method='dae.finite_difference',
                               transformation_scheme=xscheme,
                               has_holdup=True)
        m.fs.pipeline_set = pe.Set(initialize=pipe_names)
        m.fs.pipelines = GasPipeline(m.fs.pipeline_set, default=pipeline_config)

        for p in self.pipes():
            m.fs.pipelines[p.name].diameter.fix(p.diameter)
            m.fs.pipelines[p.name].control_volume.area.fix(math.pi/4*p.diameter**2)
            d_eq = m.fs.pipelines[p.name].diameter_eqn
            assert d_eq.lb == d_eq.ub
            assert pe.value(d_eq.body) >= d_eq.lb - 1e-6
            d_eq.deactivate()
            m.fs.pipelines[p.name].control_volume.length.fix(p.length)
            m.fs.pipelines[p.name].rugosity.value = p.rugosity

    def _build_node_models(self, m: _BlockData):
        node_configs = dict()
        for node in self.nodes():
            ncd = dict()
            ncd['property_package'] = m.fs.properties
            ncd['n_inlet_pipelines'] = len(self.get_inlet_links_for_node(node.name))
            ncd['n_outlet_pipelines'] = len(self.get_outlet_links_for_node(node.name))
            if node.is_source():
                ncd['n_supplies'] = 1
                ncd['n_demands'] = 0
            elif node.is_sink():
                ncd['n_supplies'] = 0
                ncd['n_demands'] = 1
            node_configs[node.name] = ncd

        m.fs.node_set = pe.Set(initialize=list(node_configs.keys()))
        m.fs.nodes = PipelineNode(m.fs.node_set, initialize=node_configs)

    def _fix_source_mole_frac_and_temp(self, m: _BlockData):
        for n in self.sources():
            m.fs.nodes[n.name].supplies[:].state[:].mole_frac_comp['natural_gas'].fix(1)
            m.fs.nodes[n.name].supplies[:].state[:].temperature.fix(n.temperature)

    def apply_node_pressure_bounds(self, m: _BlockData):
        t0 = m.fs.time.first()
        for node in self.nodes():
            pressure_units = pe.units.get_units(m.fs.nodes[node.name].state[t0].pressure)
            pmin = pe.units.convert(node.min_pressure, pressure_units)
            pmax = pe.units.convert(node.max_pressure, pressure_units)
            m.fs.nodes[node.name].state[:].pressure.setlb(pe.value(pmin))
            m.fs.nodes[node.name].state[:].pressure.setub(pe.value(pmax))

    def apply_source_flow_bounds(self, m: _BlockData):
        for node in self.sources():
            s = m.fs.nodes[node.name].supplies[0]
            for t in m.fs.time:
                flow_units = pe.units.get_units(s.state[t].flow_mol)
                min_flow_mol = node.min_flow * m.fs.properties.dens_nominal / m.fs.properties.natural_gas.mw
                max_flow_mol = node.max_flow * m.fs.properties.dens_nominal / m.fs.properties.natural_gas.mw
                min_flow_mol = pe.units.convert(min_flow_mol, flow_units)
                max_flow_mol = pe.units.convert(max_flow_mol, flow_units)
                s.state[t].flow_mol.setlb(pe.value(min_flow_mol))
                s.state[t].flow_mol.setub(pe.value(max_flow_mol))

    def fix_demand_flows(self, m: _BlockData):
        for node in self.sinks():
            n = m.fs.nodes[node.name]
            d = n.demands[0]
            for t in m.fs.time:
                flow_units = pe.units.get_units(d.flow_mol[t])
                target_flow_mol = node.target_demand * m.fs.properties.dens_nominal / m.fs.properties.natural_gas.mw
                target_flow_mol = pe.units.convert(target_flow_mol, flow_units)
                d.flow_mol[t].fix(pe.value(target_flow_mol))

    def add_node(self, name: str, node: Node):
        if name in self._nodes:
            raise ValueError(f"Network already contains a node named {name}")
        self._nodes[name] = node
        self._node_inlets[name] = OrderedSet()
        self._node_outlets[name] = OrderedSet()

    def add_link(self, name: str, link: Link):
        if name in self._links:
            raise ValueError(f"Network already contains a link named {name}")
        self._links[name] = link
        self._node_inlets[link.to_node_name].add(link)
        self._node_outlets[link.from_node_name].add(link)

    def remove_node(self, name: str):
        if len(self._node_inlets[name]) != 0:
            raise ValueError("Cannot remove node that has inlet links")
        if len(self._node_outlets[name]) != 0:
            raise ValueError("Cannot remove node that has outlet links")
        self._nodes.pop(name)
        self._node_inlets.pop(name)
        self._node_outlets.pop(name)

    def remove_link(self, name: str):
        link = self._links[name]
        self._node_inlets[link.to_node_name].remove(link)
        self._node_outlets[link.from_node_name].remove(link)
        self._links.pop(name)

    def nodes(self) -> Generator[Node, None, None]:
        for n in self._nodes.values():
            yield n

    def links(self) -> Generator[Link, None, None]:
        for l in self._links.values():
            yield l

    def sources(self) -> Generator[Source, None, None]:
        for n in self.nodes():
            if n.is_source():
                yield n

    def sinks(self) -> Generator[Sink, None, None]:
        for n in self.nodes():
            if n.is_sink():
                yield n

    def pipes(self) -> Generator[Pipe, None, None]:
        for n in self.links():
            if n.is_pipe():
                yield n

    def compressors(self) -> Generator[Compressor, None, None]:
        for n in self.links():
            if n.is_compressor():
                yield n

    def short_pipes(self) -> Generator[ShortPipe, None, None]:
        for n in self.links():
            if n.is_short_pipe():
                yield n

    def get_inlet_links_for_node(self, name: str):
        return list(self._node_inlets[name])

    def get_outlet_links_for_node(self, name: str):
        return list(self._node_outlets[name])

    def get_node(self, node_name: str):
        return self._nodes[node_name]

    def get_link(self, link_name: str):
        return self._links[link_name]

    def plot_results_at_time(self, m: _BlockData, t: float):
        # plot the nodes
        node_x = list()
        node_y = list()
        node_color = list()
        node_text = list()
        for n in self.nodes():
            node_x.append(n.x)
            node_y.append(n.y)
            txt = _get_node_txt(n, m, t)
            if n.is_source():
                node_color.append('blue')
                txt += _get_supply_txt(n, m, t)
            elif n.is_sink():
                node_color.append('red')
                txt += _get_demand_txt(n, m, t)
            else:
                node_color.append('black')
            node_text.append(txt)

        node_trace = go.Scatter(
            x=node_x, y=node_y, mode='markers', hoverinfo='text', text=node_text,
            marker=dict(color=node_color, size=15)
        )

        # plot the edges
        edge_x = list()
        edge_y = list()
        for l in self.links():
            from_node = self.get_node(l.from_node_name)
            to_node = self.get_node(l.to_node_name)
            edge_x.extend([from_node.x, to_node.x, None])
            edge_y.extend([from_node.y, to_node.y, None])

        edge_trace = go.Scatter(
            x=edge_x, y=edge_y, line=dict(width=1.5, color='black'), hoverinfo='none',
            mode='lines'
        )

        # now edge text
        edge_x = list()
        edge_y = list()
        edge_text = list()
        for l in self.links():
            from_node = self.get_node(l.from_node_name)
            to_node = self.get_node(l.to_node_name)
            x = 0.5 * (from_node.x + to_node.x)
            y = 0.5 * (from_node.y + to_node.y)
            edge_x.append(x)
            edge_y.append(y)
            txt = _get_link_txt(l, m, t)
            if l.is_pipe():
                txt += _get_pipe_txt(l, m, t)
            elif l.is_compressor():
                txt += _get_compressor_txt(l, m, t)
            elif l.is_short_pipe():
                txt += _get_short_pipe_txt(l, m, t)
            edge_text.append(txt)

        edge_text_trace = go.Scatter(
            x=edge_x, y=edge_y, mode='markers', hoverinfo='text', text=edge_text,
            marker=dict(size=1)
        )

        fig = go.Figure(
            data=[node_trace, edge_trace, edge_text_trace],
            layout=go.Layout(hovermode='closest')
        )
        fig.show()
