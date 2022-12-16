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
import pytest

import pyomo.common.unittest as unittest
import pyomo.environ as pyo
from pyomo.network.arc import Arc

from pyomo.contrib.incidence_analysis import (
    IncidenceGraphInterface,
)
from pyomo.util.check_units import assert_units_consistent

import idaes.core as idaes
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models_extra.gas_distribution.properties.natural_gas import (
    NaturalGasParameterBlock,
)
from idaes.models_extra.gas_distribution.unit_models.pipeline import GasPipeline
from idaes.models_extra.gas_distribution.unit_models.node import PipelineNode
from idaes.models_extra.gas_distribution.unit_models.compressor import (
    IsothermalCompressor,
)

"""
These tests test construction of some slightly less trivial flowsheets
involving nodes, compressors, and pipelines.

These are construction tests only. We construct the flowsheet model,
fix degrees of freedom, and assert that the model is structurally
nonsingular. Tests that actually solve flowsheet models are included
in the gas_distribution/flowsheets/tests/ directory.

"""


@pytest.mark.unit
class TestConstructFlowsheets(unittest.TestCase):
    """ """

    def test_four_nodes_linear(self):
        """
        The network looks something like this:
        s            s
        |            |
        0 (c)-> 1 -> 2 (c)-> 3
                |            |
                d            d

        """
        m = pyo.ConcreteModel()
        fs_config = {
            "dynamic": True,
            "time_set": [0.0, 20.0],
            "time_units": pyo.units.hr,
        }
        m.fs = idaes.FlowsheetBlock(**fs_config)
        m.fs.properties = NaturalGasParameterBlock()

        pipeline_config = {
            "property_package": m.fs.properties,
            "finite_elements": 2,
        }
        m.fs.pipeline_set = pyo.Set(initialize=range(3))
        m.fs.pipeline = GasPipeline(m.fs.pipeline_set, **pipeline_config)

        node_config = [
            {
                "property_package": m.fs.properties,
                "n_inlet_pipelines": 0,
                "n_outlet_pipelines": 1,
                "n_supplies": 1,
                "n_demands": 0,
            },
            {
                "property_package": m.fs.properties,
                "n_inlet_pipelines": 1,
                "n_outlet_pipelines": 1,
                "n_supplies": 0,
                "n_demands": 1,
            },
            {
                "property_package": m.fs.properties,
                "n_inlet_pipelines": 1,
                "n_outlet_pipelines": 1,
                "n_supplies": 1,
                "n_demands": 0,
            },
            {
                "property_package": m.fs.properties,
                "n_inlet_pipelines": 1,
                "n_outlet_pipelines": 0,
                "n_supplies": 0,
                "n_demands": 1,
            },
        ]
        node_config = {i: config for i, config in enumerate(node_config)}
        m.fs.node_set = pyo.Set(initialize=range(4))
        m.fs.node = PipelineNode(m.fs.node_set, initialize=node_config)

        compressor_config = {"property_package": m.fs.properties}
        m.fs.compressor_set = pyo.Set(initialize=range(2))
        m.fs.compressor = IsothermalCompressor(m.fs.compressor_set, **compressor_config)

        # Connect compressors to pipelines
        # Should/could I make this easier?
        pipeline_idx_map = {0: 0, 1: 2}

        def compressor_to_pipeline_rule(fs, i):
            return (
                m.fs.compressor[i].outlet_port,
                m.fs.pipeline[pipeline_idx_map[i]].inlet_port,
            )

        m.fs.compressor_to_pipeline = Arc(
            m.fs.compressor_set, rule=compressor_to_pipeline_rule
        )

        # Note that we are adding a compressor, not a pipeline, to the
        # outlet here.
        m.fs.node[0].add_pipeline_to_outlet(m.fs.compressor[0])
        m.fs.node[1].add_pipeline_to_inlet(m.fs.pipeline[0])
        m.fs.node[1].add_pipeline_to_outlet(m.fs.pipeline[1])
        m.fs.node[2].add_pipeline_to_inlet(m.fs.pipeline[1])
        m.fs.node[2].add_pipeline_to_outlet(m.fs.compressor[1])
        m.fs.node[3].add_pipeline_to_inlet(m.fs.pipeline[2])

        expand_arcs = pyo.TransformationFactory("network.expand_arcs")
        expand_arcs.apply_to(m)

        assert_units_consistent(m)

        ntfe = 2
        disc = pyo.TransformationFactory("dae.finite_difference")
        disc.apply_to(m, wrt=m.fs.time, nfe=ntfe)

        # Fix "design" variables
        for pipeline in m.fs.pipeline.values():
            pipeline.diameter.fix()
            pipeline.control_volume.length.fix()

        # Dynamic inputs:
        pred_dof = 10 * len(m.fs.time)
        # Initial conditions:
        pred_dof += 2 * (len(m.fs.pipeline[0].control_volume.length_domain) - 1)
        pred_dof += 2 * (len(m.fs.pipeline[1].control_volume.length_domain) - 1)
        pred_dof += 2 * (len(m.fs.pipeline[2].control_volume.length_domain) - 1)

        # Make sure we have the expected number of degrees of freedom
        self.assertEqual(degrees_of_freedom(m), pred_dof)

        # Fix predicted degrees of freedom:
        for t in m.fs.time:
            m.fs.node[0].supplies[0].state[t].mole_frac_comp[:].fix()
            m.fs.node[0].state[t].temperature.fix()
            m.fs.node[0].state[t].pressure.fix()
            m.fs.node[2].supplies[0].state[t].mole_frac_comp[:].fix()
            m.fs.node[2].supplies[0].state[t].flow_mol.fix()
            m.fs.node[2].supplies[0].state[t].temperature.fix()
            m.fs.node[1].demands[0].flow_mol[t].fix()
            m.fs.node[3].demands[0].flow_mol[t].fix()
            m.fs.compressor[0].boost_pressure[t].fix()
            m.fs.compressor[1].boost_pressure[t].fix()

        t0 = m.fs.time.first()
        x0 = m.fs.pipeline[0].control_volume.length_domain.first()
        xf = m.fs.pipeline[0].control_volume.length_domain.last()
        for x in m.fs.pipeline[0].control_volume.length_domain:
            # Here I assume that all three pipelines have the same
            # length domain.
            if x != x0:
                m.fs.pipeline[0].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[1].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[2].control_volume.pressure[t0, x].fix()
            if x != xf:
                m.fs.pipeline[0].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[1].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[2].control_volume.flow_mass[t0, x].fix()

        # Make sure we have zero degrees of freedom and a perfect matching
        # between constraints and variables.
        self.assertEqual(degrees_of_freedom(m), 0)
        igraph = IncidenceGraphInterface(m)
        N, M = len(igraph.constraints), len(igraph.variables)
        matching = igraph.maximum_matching()
        self.assertEqual(len(matching), N)

    def test_two_supply_one_demand(self):
        """
        The network looks something like this:
        s-0 
           \
            --> 2 (c)-> 3-d
           /
        s-1
        """
        m = pyo.ConcreteModel()
        fs_config = {
            "dynamic": True,
            "time_set": [0.0, 20.0],
            "time_units": pyo.units.hr,
        }
        m.fs = idaes.FlowsheetBlock(**fs_config)
        m.fs.properties = NaturalGasParameterBlock()

        pipeline_config = {
            "property_package": m.fs.properties,
            "finite_elements": 2,
        }
        m.fs.pipeline_set = pyo.Set(initialize=range(3))
        m.fs.pipeline = GasPipeline(m.fs.pipeline_set, **pipeline_config)

        node_config = [
            {
                "property_package": m.fs.properties,
                "n_inlet_pipelines": 0,
                "n_outlet_pipelines": 1,
                "n_supplies": 1,
                "n_demands": 0,
            },
            {
                "property_package": m.fs.properties,
                "n_inlet_pipelines": 0,
                "n_outlet_pipelines": 1,
                "n_supplies": 1,
                "n_demands": 0,
            },
            {
                "property_package": m.fs.properties,
                "n_inlet_pipelines": 2,
                "n_outlet_pipelines": 1,
                "n_supplies": 0,
                "n_demands": 0,
            },
            {
                "property_package": m.fs.properties,
                "n_inlet_pipelines": 1,
                "n_outlet_pipelines": 0,
                "n_supplies": 0,
                "n_demands": 1,
            },
        ]
        node_config = {i: config for i, config in enumerate(node_config)}
        m.fs.node_set = pyo.Set(initialize=range(4))
        m.fs.node = PipelineNode(m.fs.node_set, initialize=node_config)

        compressor_config = {"property_package": m.fs.properties}
        m.fs.compressor_set = pyo.Set(initialize=range(1))
        m.fs.compressor = IsothermalCompressor(m.fs.compressor_set, **compressor_config)

        # Connect compressors to pipelines
        # Should/could I make this easier?
        pipeline_idx_map = {0: 2}

        def compressor_to_pipeline_rule(fs, i):
            return (
                m.fs.compressor[i].outlet_port,
                m.fs.pipeline[pipeline_idx_map[i]].inlet_port,
            )

        m.fs.compressor_to_pipeline = Arc(
            m.fs.compressor_set, rule=compressor_to_pipeline_rule
        )

        m.fs.node[0].add_pipeline_to_outlet(m.fs.pipeline[0])
        m.fs.node[1].add_pipeline_to_outlet(m.fs.pipeline[1])
        m.fs.node[2].add_pipeline_to_inlet(m.fs.pipeline[0])
        m.fs.node[2].add_pipeline_to_inlet(m.fs.pipeline[1])
        m.fs.node[2].add_pipeline_to_outlet(m.fs.compressor[0])
        m.fs.node[3].add_pipeline_to_inlet(m.fs.pipeline[2])

        expand_arcs = pyo.TransformationFactory("network.expand_arcs")
        expand_arcs.apply_to(m)

        assert_units_consistent(m)

        ntfe = 2
        disc = pyo.TransformationFactory("dae.finite_difference")
        disc.apply_to(m, wrt=m.fs.time, nfe=ntfe)

        # Fix "design" variables
        for pipeline in m.fs.pipeline.values():
            pipeline.diameter.fix()
            pipeline.control_volume.length.fix()

        # Predicted degrees of freedom:
        pred_dof = 8 * len(m.fs.time)
        pred_dof += 2 * (len(m.fs.pipeline[0].control_volume.length_domain) - 1)
        pred_dof += 2 * (len(m.fs.pipeline[1].control_volume.length_domain) - 1)
        pred_dof += 2 * (len(m.fs.pipeline[2].control_volume.length_domain) - 1)

        self.assertEqual(degrees_of_freedom(m), pred_dof)

        # Fix predicted degrees of freedom:
        # Dynamic inputs:
        for t in m.fs.time:
            m.fs.node[0].supplies[0].state[t].mole_frac_comp[:].fix()
            m.fs.node[0].state[t].temperature.fix()
            m.fs.node[0].state[t].pressure.fix()
            m.fs.node[1].supplies[0].state[t].mole_frac_comp[:].fix()
            m.fs.node[1].supplies[0].state[t].flow_mol.fix()
            m.fs.node[1].supplies[0].state[t].temperature.fix()
            m.fs.node[3].demands[0].flow_mol[t].fix()
            m.fs.compressor[0].boost_pressure[t].fix()

        # Initial conditions:
        t0 = m.fs.time.first()
        x0 = m.fs.pipeline[0].control_volume.length_domain.first()
        xf = m.fs.pipeline[0].control_volume.length_domain.last()
        for x in m.fs.pipeline[0].control_volume.length_domain:
            # Here I assume that all three pipelines have the same
            # length domain.
            if x != x0:
                m.fs.pipeline[0].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[1].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[2].control_volume.pressure[t0, x].fix()
            if x != xf:
                m.fs.pipeline[0].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[1].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[2].control_volume.flow_mass[t0, x].fix()

        # Make sure we have zero degrees of freedom and a perfect matching
        # between constraints and variables.
        self.assertEqual(degrees_of_freedom(m), 0)
        igraph = IncidenceGraphInterface(m)
        N, M = len(igraph.constraints), len(igraph.variables)
        matching = igraph.maximum_matching()
        self.assertEqual(len(matching), N)

    def test_one_supply_two_demand(self):
        """
        The network looks something like this:
                   (c)-> 2-d
                  /
        s-0 -> 1 - 
                  \
                   (c)-> 3-d
        """
        m = pyo.ConcreteModel()
        fs_config = {
            "dynamic": True,
            "time_set": [0.0, 20.0],
            "time_units": pyo.units.hr,
        }
        m.fs = idaes.FlowsheetBlock(**fs_config)
        m.fs.properties = NaturalGasParameterBlock()

        pipeline_config = {
            "property_package": m.fs.properties,
            "finite_elements": 2,
        }
        m.fs.pipeline_set = pyo.Set(initialize=range(3))
        m.fs.pipeline = GasPipeline(m.fs.pipeline_set, **pipeline_config)

        node_config = [
            {
                "property_package": m.fs.properties,
                "n_inlet_pipelines": 0,
                "n_outlet_pipelines": 1,
                "n_supplies": 1,
                "n_demands": 0,
            },
            {
                "property_package": m.fs.properties,
                "n_inlet_pipelines": 1,
                "n_outlet_pipelines": 2,
                "n_supplies": 0,
                "n_demands": 0,
            },
            {
                "property_package": m.fs.properties,
                "n_inlet_pipelines": 1,
                "n_outlet_pipelines": 0,
                "n_supplies": 0,
                "n_demands": 1,
            },
            {
                "property_package": m.fs.properties,
                "n_inlet_pipelines": 1,
                "n_outlet_pipelines": 0,
                "n_supplies": 0,
                "n_demands": 1,
            },
        ]
        node_config = {i: config for i, config in enumerate(node_config)}
        m.fs.node_set = pyo.Set(initialize=range(4))
        m.fs.node = PipelineNode(m.fs.node_set, initialize=node_config)

        compressor_config = {"property_package": m.fs.properties}
        m.fs.compressor_set = pyo.Set(initialize=range(2))
        m.fs.compressor = IsothermalCompressor(m.fs.compressor_set, **compressor_config)

        # Connect compressors to pipelines
        # Should/could I make this easier?
        pipeline_idx_map = {0: 1, 1: 2}

        def compressor_to_pipeline_rule(fs, i):
            return (
                m.fs.compressor[i].outlet_port,
                m.fs.pipeline[pipeline_idx_map[i]].inlet_port,
            )

        m.fs.compressor_to_pipeline = Arc(
            m.fs.compressor_set, rule=compressor_to_pipeline_rule
        )

        m.fs.node[0].add_pipeline_to_outlet(m.fs.pipeline[0])
        m.fs.node[1].add_pipeline_to_inlet(m.fs.pipeline[0])
        m.fs.node[1].add_pipeline_to_outlet(m.fs.compressor[0])
        m.fs.node[1].add_pipeline_to_outlet(m.fs.compressor[1])
        m.fs.node[2].add_pipeline_to_inlet(m.fs.pipeline[1])
        m.fs.node[3].add_pipeline_to_inlet(m.fs.pipeline[2])

        expand_arcs = pyo.TransformationFactory("network.expand_arcs")
        expand_arcs.apply_to(m)

        assert_units_consistent(m)

        ntfe = 2
        disc = pyo.TransformationFactory("dae.finite_difference")
        disc.apply_to(m, wrt=m.fs.time, nfe=ntfe)

        # Fix "design" variables
        for pipeline in m.fs.pipeline.values():
            pipeline.diameter.fix()
            pipeline.control_volume.length.fix()

        # Predicted degrees of freedom (assuming all pipelines have same length
        # discretization)
        pred_dof = 7 * len(m.fs.time)
        pred_dof += 6 * (len(m.fs.pipeline[0].control_volume.length_domain) - 1)
        self.assertEqual(degrees_of_freedom(m), pred_dof)

        # Fix predicted degrees of freedom:
        # Dynamic inputs:
        for t in m.fs.time:
            m.fs.node[0].supplies[0].state[t].mole_frac_comp[:].fix()
            m.fs.node[0].state[t].temperature.fix()
            m.fs.node[0].state[t].pressure.fix()
            m.fs.compressor[0].boost_pressure[t].fix()
            m.fs.compressor[1].boost_pressure[t].fix()
            m.fs.node[2].demands[0].flow_mol[t].fix()
            m.fs.node[3].demands[0].flow_mol[t].fix()

        # Initial conditions:
        t0 = m.fs.time.first()
        x0 = m.fs.pipeline[0].control_volume.length_domain.first()
        xf = m.fs.pipeline[0].control_volume.length_domain.last()
        for x in m.fs.pipeline[0].control_volume.length_domain:
            # Here I assume that all three pipelines have the same
            # length domain.
            if x != x0:
                m.fs.pipeline[0].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[1].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[2].control_volume.pressure[t0, x].fix()
            if x != xf:
                m.fs.pipeline[0].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[1].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[2].control_volume.flow_mass[t0, x].fix()

        # Make sure we have zero degrees of freedom and a perfect matching
        # between constraints and variables.
        self.assertEqual(degrees_of_freedom(m), 0)
        igraph = IncidenceGraphInterface(m)
        N, M = len(igraph.constraints), len(igraph.variables)
        matching = igraph.maximum_matching()
        self.assertEqual(len(matching), N)


if __name__ == "__main__":
    unittest.main()
