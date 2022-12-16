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
from pyomo.common.collections import ComponentSet
import pyomo.environ as pyo
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

"""
Tests for the PipelineNode unit model
"""


@pytest.mark.unit
class TestConstructNode(unittest.TestCase):
    """
    Tests for the construction of nodes with various attached pipeline
    configurations.
    """

    def test_construct_node_with_pipelines(self):
        """
        A node with one inlet, one outlet, no supplies, and
        no demands.
        """
        m = pyo.ConcreteModel()
        default = {
            "dynamic": False,
        }
        m.fs = idaes.FlowsheetBlock(**default)
        m.fs.properties = NaturalGasParameterBlock()
        pipeline_config = {
            "property_package": m.fs.properties,
            "finite_elements": 2,
        }
        m.fs.pipeline1 = GasPipeline(**pipeline_config)
        pipeline1 = m.fs.pipeline1
        m.fs.pipeline2 = GasPipeline(**pipeline_config)
        pipeline2 = m.fs.pipeline2

        node_config = {"property_package": m.fs.properties}
        m.fs.node = PipelineNode(**node_config)
        node = m.fs.node

        var_set = ComponentSet(node.component_data_objects(pyo.Var))
        # Three state blocks, four variables each.
        self.assertEqual(len(var_set), 12)

        igraph = IncidenceGraphInterface(node)
        # This test asserts that each of these variables participates
        # in an active constraint (is recognized by the incidence graph)
        self.assertEqual(len(igraph.variables), len(var_set))
        for var in igraph.variables:
            self.assertIn(var, var_set)

        time = m.fs.time
        t0 = time.first()
        j = next(iter(m.fs.properties.component_list))

        states = [node.state, node.inlets[0].state, node.outlets[0].state]
        for state in states:
            # Make sure the inlet and outlet blocks have state blocks
            # with their state variables.
            self.assertTrue(isinstance(state[t0].temperature, pyo.Var))
            self.assertIn(state[t0].temperature, var_set)
            self.assertTrue(isinstance(state[t0].pressure, pyo.Var))
            self.assertIn(state[t0].pressure, var_set)
            self.assertTrue(isinstance(state[t0].mole_frac_comp, pyo.Var))
            self.assertIn(state[t0].mole_frac_comp[j], var_set)
            self.assertTrue(isinstance(state[t0].flow_mol, pyo.Var))
            self.assertIn(state[t0].flow_mol, var_set)

        con_set = ComponentSet(node.component_data_objects(pyo.Constraint))
        # Pressure equality between inlet and node, temperature, pressure,
        # and mole fraction equalities between the outlet and node,
        # plus flow balance, component mixing, and total flow equations.
        self.assertEqual(len(con_set), 8)

        blocks = [node.inlets[0], node.outlets[0]]
        for block in blocks:
            self.assertTrue(isinstance(block.pressure_eq, pyo.Constraint))
            self.assertIn(block.pressure_eq[t0], con_set)

        # Inlet mole fractions and temperatures are already specified,
        # so we don't have this equation.
        block = node.outlets[0]
        self.assertTrue(isinstance(block.mole_frac_comp_eq, pyo.Constraint))
        self.assertIn(block.mole_frac_comp_eq[t0, j], con_set)
        self.assertTrue(isinstance(block.temperature_eq, pyo.Constraint))
        self.assertIn(block.temperature_eq[t0], con_set)

        self.assertTrue(isinstance(node.flow_balance, pyo.Constraint))
        self.assertTrue(isinstance(node.component_mixing_eq, pyo.Constraint))
        self.assertTrue(isinstance(node.total_flow_eq, pyo.Constraint))
        self.assertTrue(isinstance(node.enthalpy_mixing_eq, pyo.Constraint))
        self.assertIn(node.flow_balance[t0], con_set)
        self.assertIn(node.component_mixing_eq[t0, j], con_set)
        self.assertIn(node.total_flow_eq[t0], con_set)
        self.assertIn(node.enthalpy_mixing_eq[t0], con_set)

        node.add_pipeline_to_inlet(pipeline1)
        node.add_pipeline_to_outlet(pipeline2)

        p1 = node.get_inlet_pipeline(0)
        self.assertIs(p1, pipeline1)
        p2 = node.get_outlet_pipeline(0)
        self.assertIs(p2, pipeline2)

        expand_arcs = pyo.TransformationFactory("network.expand_arcs")
        expand_arcs.apply_to(m)

        pipeline1.diameter.fix()
        pipeline1.control_volume.length.fix()
        pipeline2.diameter.fix()
        pipeline2.control_volume.length.fix()

        cv1 = pipeline1.control_volume
        x0 = cv1.length_domain.first()
        cv1.properties[t0, x0].temperature.fix()
        cv1.properties[t0, x0].pressure.fix()
        cv1.properties[t0, x0].mole_frac_comp[j].fix()
        cv1.properties[t0, x0].flow_mol.fix()

        assert_units_consistent(m)
        igraph = IncidenceGraphInterface(m)
        matching = igraph.maximum_matching()
        N, M = len(igraph.constraints), len(igraph.variables)
        self.assertEqual(degrees_of_freedom(m), 0)
        self.assertEqual(N, M)  # Sanity check
        self.assertEqual(len(matching), N)

    def test_pipeline_with_two_nodes(self):
        """
        One pipeline connecting two nodes.
        """
        m = pyo.ConcreteModel()
        default = {
            "dynamic": False,
        }
        m.fs = idaes.FlowsheetBlock(**default)
        m.fs.properties = NaturalGasParameterBlock()
        pipeline_config = {
            "property_package": m.fs.properties,
            "finite_elements": 2,
        }
        m.fs.pipeline = GasPipeline(**pipeline_config)
        pipeline = m.fs.pipeline

        node_configs = [
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
                "n_outlet_pipelines": 0,
                "n_supplies": 0,
                "n_demands": 1,
            },
        ]
        m.fs.node_set = pyo.Set(initialize=list(range(len(node_configs))))
        node_configs = {i: config for i, config in enumerate(node_configs)}
        m.fs.nodes = PipelineNode(m.fs.node_set, initialize=node_configs)

        time = m.fs.time
        t0 = m.fs.time.first()

        m.fs.nodes[0].add_pipeline_to_outlet(pipeline)
        m.fs.nodes[1].add_pipeline_to_inlet(pipeline)
        expand_arcs = pyo.TransformationFactory("network.expand_arcs")
        expand_arcs.apply_to(m)

        # Fix degrees of freedom
        pipeline.control_volume.length.fix()
        pipeline.diameter.fix()

        # Fix state variables at inlet node.
        # Fixing the flow rate here, which represents the total flow rate
        # through the node, implicitly fixes the supply flow.
        # Same goes for the mole fraction.
        #
        # This loop also specifies pressure and temperature of the supply
        for var in m.fs.nodes[0].state[t0].component_data_objects(pyo.Var):
            var.fix()

        assert_units_consistent(m)
        self.assertEqual(degrees_of_freedom(m), 0)
        igraph = IncidenceGraphInterface(m)
        N, M = len(igraph.constraints), len(igraph.variables)
        matching = igraph.maximum_matching()
        self.assertEqual(N, M)
        self.assertEqual(N, len(matching))

        var_set = ComponentSet(igraph.variables)
        self.assertIn(m.fs.nodes[0].supplies[0].flow_mol[t0], var_set)
        self.assertIn(m.fs.nodes[1].demands[0].flow_mol[t0], var_set)

        # Temperature and pressure equality constraints make sure
        # supply variables are the same as node variables.
        self.assertIn(m.fs.nodes[0].supplies[0].state[t0].temperature, var_set)
        self.assertIn(m.fs.nodes[0].supplies[0].state[t0].pressure, var_set)

    def test_multiple_inlets_outlets(self):
        """
        A test with multiple inlets and outlets
        """
        m = pyo.ConcreteModel()
        default = {
            "dynamic": False,
        }
        m.fs = idaes.FlowsheetBlock(**default)
        m.fs.properties = NaturalGasParameterBlock()
        pipeline_config = {
            "property_package": m.fs.properties,
            "finite_elements": 2,
        }
        m.fs.pipeline_set = pyo.Set(initialize=list(range(7)))
        m.fs.pipeline = GasPipeline(m.fs.pipeline_set, **pipeline_config)

        node_config = {
            "property_package": m.fs.properties,
            "n_supplies": 1,
            "n_demands": 2,
            "n_inlet_pipelines": 3,
            "n_outlet_pipelines": 4,
        }
        m.fs.node = PipelineNode(**node_config)

        m.fs.node.add_pipeline_to_inlet(m.fs.pipeline[0])
        m.fs.node.add_pipeline_to_inlet(m.fs.pipeline[1])
        m.fs.node.add_pipeline_to_inlet(m.fs.pipeline[2])
        m.fs.node.add_pipeline_to_outlet(m.fs.pipeline[3])
        m.fs.node.add_pipeline_to_outlet(m.fs.pipeline[4])
        m.fs.node.add_pipeline_to_outlet(m.fs.pipeline[5])
        m.fs.node.add_pipeline_to_outlet(m.fs.pipeline[6])

        expand_arcs = pyo.TransformationFactory("network.expand_arcs")
        expand_arcs.apply_to(m)

        j = next(iter(m.fs.properties.component_list))

        for pipeline in m.fs.pipeline.values():
            pipeline.control_volume.length.fix()
            pipeline.diameter.fix()

        t0 = m.fs.time.first()
        x0 = m.fs.pipeline[0].control_volume.length_domain.first()
        state0 = m.fs.pipeline[0].control_volume.properties[t0, x0]

        # Fix pressure of one inlet
        state0.pressure.fix()

        # Fix flow rates and mole fractions of all inlets
        for pipeline in m.fs.node.inlet_pipelines().values():
            pipeline.control_volume.flow_mass[t0, x0].fix()
            state = pipeline.control_volume.properties[t0, x0]
            state.mole_frac_comp[j].fix()
            state.temperature.fix()

        # Fix supply flows and mole fractions
        m.fs.node.supplies[:].flow_mol[:].fix()
        m.fs.node.supplies[:].state[:].mole_frac_comp[:].fix()
        m.fs.node.supplies[:].state[:].temperature.fix()

        # Fix demand flows, and all outlet flow rates except one
        m.fs.node.demands[:].flow_mol[:].fix()
        m.fs.pipeline[4].control_volume.flow_mass[t0, x0].fix()
        m.fs.pipeline[5].control_volume.flow_mass[t0, x0].fix()
        m.fs.pipeline[6].control_volume.flow_mass[t0, x0].fix()

        assert_units_consistent(m)
        self.assertEqual(degrees_of_freedom(m), 0)
        igraph = IncidenceGraphInterface(m)
        N, M = len(igraph.constraints), len(igraph.variables)
        matching = igraph.maximum_matching()
        self.assertEqual(N, M)
        self.assertEqual(N, len(matching))

    def test_dynamic_two_nodes(self):
        """
        A dynamic instance of the test with two nodes
        """
        m = pyo.ConcreteModel()
        default = {
            "dynamic": True,
            "time_set": [0.0, 20.0],
            "time_units": pyo.units.hr,
        }
        m.fs = idaes.FlowsheetBlock(**default)
        m.fs.properties = NaturalGasParameterBlock()
        pipeline_config = {
            "property_package": m.fs.properties,
            "finite_elements": 2,
        }
        m.fs.pipeline = GasPipeline(**pipeline_config)
        pipeline = m.fs.pipeline

        node_configs = [
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
                "n_outlet_pipelines": 0,
                "n_supplies": 0,
                "n_demands": 1,
            },
        ]
        m.fs.node_set = pyo.Set(initialize=list(range(len(node_configs))))
        node_configs = {i: config for i, config in enumerate(node_configs)}
        m.fs.nodes = PipelineNode(m.fs.node_set, initialize=node_configs)

        time = m.fs.time
        t0 = m.fs.time.first()
        ntfe = 4
        disc = pyo.TransformationFactory("dae.finite_difference")
        disc.apply_to(m, wrt=time, nfe=ntfe, scheme="BACKWARD")

        m.fs.nodes[0].add_pipeline_to_outlet(pipeline)
        m.fs.nodes[1].add_pipeline_to_inlet(pipeline)
        expand_arcs = pyo.TransformationFactory("network.expand_arcs")
        expand_arcs.apply_to(m)

        cv = pipeline.control_volume
        x0 = cv.length_domain.first()
        xf = cv.length_domain.last()
        j = next(iter(m.fs.properties.component_list))

        # Deactivate momentum balance where neither derivative is defined
        cv.momentum_balance[t0, xf].deactivate()

        # Fix geometry variables
        pipeline.control_volume.length.fix()
        pipeline.diameter.fix()

        # Fix dynamic inputs: supply pressures and demand flows
        m.fs.nodes[0].supplies[:].state[:].mole_frac_comp[j].fix()
        m.fs.nodes[0].state[:].pressure.fix()
        m.fs.nodes[0].state[:].temperature.fix()
        m.fs.nodes[1].state[:].flow_mol.fix()

        # Fix initial conditions: flow and pressure except where already
        # defined.
        for x in cv.length_domain:
            if x != x0:
                cv.pressure[t0, x].fix()
            if x != xf:
                cv.flow_mass[t0, x].fix()

        self.assertEqual(degrees_of_freedom(m), 0)
        igraph = IncidenceGraphInterface(m)
        N, M = len(igraph.constraints), len(igraph.variables)
        matching = igraph.maximum_matching()
        self.assertEqual(N, M)
        self.assertEqual(N, len(matching))

        var_set = ComponentSet(igraph.variables)
        for t in time:
            self.assertIn(m.fs.nodes[0].supplies[0].flow_mol[t], var_set)
            self.assertIn(m.fs.nodes[0].supplies[0].state[t].temperature, var_set)
            self.assertIn(m.fs.nodes[0].supplies[0].state[t].pressure, var_set)
            self.assertIn(m.fs.nodes[1].demands[0].flow_mol[t], var_set)

    def test_dynamic_multiple_inlet_outlet(self):
        """
        A dynamic instance of the multiple inlet, multiple outlet test
        """
        m = pyo.ConcreteModel()
        default = {
            "dynamic": True,
            "time_set": [0.0, 20.0],
            "time_units": pyo.units.hr,
        }
        m.fs = idaes.FlowsheetBlock(**default)
        m.fs.properties = NaturalGasParameterBlock()
        pipeline_config = {
            "property_package": m.fs.properties,
            "finite_elements": 2,
        }
        m.fs.pipeline_set = pyo.Set(initialize=list(range(7)))
        m.fs.pipeline = GasPipeline(m.fs.pipeline_set, **pipeline_config)

        node_config = {
            "property_package": m.fs.properties,
            "n_supplies": 1,
            "n_demands": 2,
            "n_inlet_pipelines": 3,
            "n_outlet_pipelines": 4,
        }
        m.fs.node = PipelineNode(**node_config)

        m.fs.node.add_pipeline_to_inlet(m.fs.pipeline[0])
        m.fs.node.add_pipeline_to_inlet(m.fs.pipeline[1])
        m.fs.node.add_pipeline_to_inlet(m.fs.pipeline[2])
        m.fs.node.add_pipeline_to_outlet(m.fs.pipeline[3])
        m.fs.node.add_pipeline_to_outlet(m.fs.pipeline[4])
        m.fs.node.add_pipeline_to_outlet(m.fs.pipeline[5])
        m.fs.node.add_pipeline_to_outlet(m.fs.pipeline[6])

        expand_arcs = pyo.TransformationFactory("network.expand_arcs")
        expand_arcs.apply_to(m)

        time = m.fs.time
        t0 = time.first()
        ntfe = 4
        disc = pyo.TransformationFactory("dae.finite_difference")
        disc.apply_to(m, wrt=time, nfe=ntfe, scheme="BACKWARD")

        x0 = m.fs.pipeline[0].control_volume.length_domain.first()
        xf = m.fs.pipeline[0].control_volume.length_domain.last()
        # Deactivate momentum balance at point where both derivatives
        # are not defined.
        mom_bal = m.fs.pipeline[:].control_volume.momentum_balance[t0, xf]
        mom_bal.deactivate()

        j = next(iter(m.fs.properties.component_list))

        # Fix "static inputs," i.e. design variables
        for pipeline in m.fs.pipeline.values():
            pipeline.control_volume.length.fix()
            pipeline.diameter.fix()

        # Fix dynamic inputs

        state0 = m.fs.pipeline[0].control_volume.properties

        # Pressure of one inlet
        state0[:, x0].pressure.fix()

        # Flow rates, mole fractions, and temperatures of all inlets
        for pipeline in m.fs.node.inlet_pipelines().values():
            pipeline.control_volume.flow_mass[:, x0].fix()
            state = pipeline.control_volume.properties
            state[:, x0].mole_frac_comp[j].fix()
            state[:, x0].temperature.fix()

        # Supply flows, mole fractions, and temperatures
        m.fs.node.supplies[:].flow_mol[:].fix()
        m.fs.node.supplies[:].state[:].mole_frac_comp[:].fix()
        m.fs.node.supplies[:].state[:].temperature.fix()

        # Demand flows, and all outlet flow rates except one
        m.fs.node.demands[:].flow_mol[:].fix()
        m.fs.pipeline[4].control_volume.flow_mass[:, x0].fix()
        m.fs.pipeline[5].control_volume.flow_mass[:, x0].fix()
        m.fs.pipeline[6].control_volume.flow_mass[:, x0].fix()

        # Initial conditions. These are the differential variables of the
        # pipelines, pressure and flow_mass, except where they overlap
        # with "dynamic inputs."
        for pipeline in m.fs.node.inlet_pipelines().values():
            cv = pipeline.control_volume
            state = cv.properties
            for x in cv.length_domain:
                if x != cv.length_domain.first():
                    cv.flow_mass[t0, x].fix()

                if state is state0:
                    # Here, inlet pressure has already been specified
                    if x != cv.length_domain.first():
                        cv.pressure[t0, x].fix()
                else:
                    # Outlet pressure is specified at the initial condition
                    # by node pressure
                    if x != cv.length_domain.last():
                        cv.pressure[t0, x].fix()

        for pipeline in m.fs.node.outlet_pipelines().values():
            cv = pipeline.control_volume
            state = cv.properties
            for x in cv.length_domain:
                if x != cv.length_domain.first():
                    # Only values at the node have already been specified.
                    cv.flow_mass[t0, x].fix()
                    cv.pressure[t0, x].fix()

        igraph = IncidenceGraphInterface(m)
        N, M = len(igraph.constraints), len(igraph.variables)
        matching = igraph.maximum_matching()
        self.assertEqual(degrees_of_freedom(m), 0)
        self.assertEqual(N, M)
        self.assertEqual(N, len(matching))


if __name__ == "__main__":
    unittest.main()
