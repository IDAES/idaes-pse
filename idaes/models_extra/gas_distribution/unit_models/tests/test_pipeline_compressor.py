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
from pyomo.dae.flatten import flatten_dae_components
from pyomo.network.arc import Arc

from pyomo.contrib.incidence_analysis import (
    IncidenceGraphInterface,
)
from pyomo.util.check_units import assert_units_consistent

import idaes.core as idaes
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    large_residuals_set,
)
from idaes.models_extra.gas_distribution.properties.natural_gas import (
    NaturalGasParameterBlock,
)
from idaes.models_extra.gas_distribution.unit_models.pipeline import GasPipeline
from idaes.models_extra.gas_distribution.unit_models.compressor import (
    IsothermalCompressor as Compressor,
)

from idaes.apps.nmpc.dynamic_data import (
    load_inputs_into_model,
    interval_data_from_time_series,
)


@pytest.mark.component
class TestSolveDynamicPipelineCompressor(unittest.TestCase):
    def make_steady_model(
        self,
        nfe=2,
        scheme="FORWARD",
    ):
        m = pyo.ConcreteModel()
        default = {"dynamic": False}
        m.fs = idaes.FlowsheetBlock(**default)
        m.fs.properties = NaturalGasParameterBlock()
        pipeline_config = {
            "property_package": m.fs.properties,
            "finite_elements": nfe,
            "transformation_scheme": scheme,
            "has_holdup": True,
        }
        m.fs.pipeline = GasPipeline(**pipeline_config)
        pipeline = m.fs.pipeline
        compressor_config = {"property_package": m.fs.properties}
        m.fs.compressor = Compressor(**compressor_config)
        compressor = m.fs.compressor
        cv = m.fs.pipeline.control_volume
        # Fix geometry variables
        m.fs.pipeline.diameter.fix(0.92 * pyo.units.m)
        cv.length.fix(300.0 * pyo.units.km)

        m._compressor_to_pipeline = Arc(
            ports=(compressor.outlet_port, pipeline.inlet_port),
        )
        expand_arcs = pyo.TransformationFactory("network.expand_arcs")
        expand_arcs.apply_to(m)
        return m

    def fix_model_inlets(
        self,
        model,
        inlet_flow_mass=3.0e5 * pyo.units.kg / pyo.units.hr,
        inlet_pressure=57.0 * pyo.units.bar,
        inlet_temperature=293.15 * pyo.units.K,
    ):
        cv = model.fs.pipeline.control_volume
        j = next(iter(model.fs.properties.component_list))
        x0 = cv.length_domain.first()
        state = model.fs.compressor.inlet_state
        state[:].mole_frac_comp[j].fix()
        state[:].temperature.fix(inlet_temperature)
        state[:].pressure.fix(inlet_pressure)
        cv.flow_mass[:, x0].fix(inlet_flow_mass)

    def get_scalar_data_from_model(
        self,
        model,
        time,
        scalar_vars=None,
        dae_vars=None,
    ):
        if scalar_vars is None or dae_vars is None:
            scalar_vars, dae_vars = flatten_dae_components(model, time, pyo.Var)
        return {str(pyo.ComponentUID(var)): var.value for var in scalar_vars}

    def get_data_from_model_at_time(
        self,
        model,
        time,
        scalar_vars=None,
        dae_vars=None,
        t0=0,
    ):
        if scalar_vars is None or dae_vars is None:
            scalar_vars, dae_vars = flatten_dae_components(model, time, pyo.Var)
        return {str(pyo.ComponentUID(var.referent)): var[t0].value for var in dae_vars}

    def test_sim(self):
        """
        Inlet pressure and outlet flow rate will be fixed.
        """
        nxfe = 4
        ipopt = pyo.SolverFactory("ipopt")

        m_steady = self.make_steady_model(nfe=nxfe)
        self.fix_model_inlets(m_steady, inlet_pressure=50.0 * pyo.units.bar)
        m_steady.fs.compressor.boost_pressure[:].fix(7.0 * pyo.units.bar)
        ipopt.solve(m_steady, tee=True)
        time_steady = m_steady.fs.time
        scalar_data = self.get_scalar_data_from_model(m_steady, time_steady)
        initial_data = self.get_data_from_model_at_time(m_steady, time_steady)

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
            "finite_elements": nxfe,
        }
        m.fs.pipeline = GasPipeline(**pipeline_config)
        pipeline = m.fs.pipeline
        compressor_config = {"property_package": m.fs.properties}
        m.fs.compressor = Compressor(**compressor_config)
        compressor = m.fs.compressor
        m._compressor_to_pipeline = Arc(
            ports=(compressor.outlet_port, pipeline.inlet_port),
        )
        expand_arcs = pyo.TransformationFactory("network.expand_arcs")
        expand_arcs.apply_to(m)

        cv = m.fs.pipeline.control_volume
        assert_units_consistent(m)

        disc = pyo.TransformationFactory("dae.finite_difference")
        ntfe = 20
        disc.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme="BACKWARD")

        time = m.fs.time
        t0 = m.fs.time.first()
        x0 = cv.length_domain.first()
        xf = cv.length_domain.last()
        j = next(iter(m.fs.properties.component_list))

        # Fix geometry variables
        m.fs.pipeline.diameter.fix(0.92 * pyo.units.m)
        cv.length.fix(300.0 * pyo.units.km)

        # Fix boost pressure
        compressor.boost_pressure[:].fix()

        # Inlets to the compressor are fixed, except for flow, where
        # the outlet is fixed.
        state = compressor.inlet_state
        state[:].pressure.fix()
        state[:].mole_frac_comp[j].fix()
        state[:].temperature.fix()
        cv.flow_mass[:, xf].fix()

        # Fix initial conditions. Here, pressure and volume for all
        # non-specified points.
        for x in cv.length_domain:
            if x != x0:
                cv.pressure[t0, x].fix()
            if x != xf:
                cv.flow_mass[t0, x].fix()

        # I want to deactivate differential equations at (t0, xf)
        # Material balance already doesn't exist here.
        cv.momentum_balance[t0, xf].deactivate()

        self.assertEqual(degrees_of_freedom(m), 0)

        # Load initial steady state into model at all time points.
        for name, val in initial_data.items():
            var = m.find_component(name)
            for t in time:
                var[t].set_value(val)
        # Load scalar data from initial steady state
        # (initialize area, basically)
        for name, val in scalar_data.items():
            var = m.find_component(name)
            var.set_value(val)

        cv.material_accumulation[...].set_value(0.0)
        cv.flow_mass_dt[...].set_value(0.0)

        for con in large_residuals_set(m):
            resid = pyo.value(con.body - con.upper)
            print(resid, con.name)
        ipopt.solve(m, tee=True)

        # Load input sequence into model
        sample_points = [4.0, 20.0]
        input_name = "fs.pipeline.control_volume.flow_mass[*,1.0]"
        nominal_density = 0.72
        val = 12.0 * 1e6 / 24 * nominal_density  # 12 (1e6 SCM)/day
        input_series_data = (
            sample_points,
            {input_name: [val, val]},
        )
        input_interval_data = interval_data_from_time_series(input_series_data)
        load_inputs_into_model(m, time, input_interval_data)
        # Solve with loaded inputs
        res = ipopt.solve(m, tee=True)
        self.assertIs(
            res.solver.termination_condition,
            pyo.TerminationCondition.optimal,
        )

        # These predicted values come from a simulation of a single pipeline
        # model from the Pyomo DAE example. flow_mass has been converted
        # to kg/hr from (1e4 SCM/hr) by a factor of 0.72*1e4, where
        # 0.72 kg/m**3 is the gas density at standard conditions.
        pred_values = (
            list(time),
            {
                "fs.pipeline.control_volume.flow_mass[*,%s]"
                % x0: [
                    3.000e5,
                    2.999e5,
                    2.999e5,
                    2.999e5,
                    3.000e5,
                    3.174e5,
                    3.301e5,
                    3.389e5,
                    3.449e5,
                    3.492e5,
                    3.523e5,
                    3.544e5,
                    3.560e5,
                    3.571e5,
                    3.579e5,
                    3.585e5,
                    3.589e5,
                    3.592e5,
                    3.594e5,
                    3.595e5,
                    3.597e5,
                ],
                "fs.pipeline.control_volume.pressure[*,%s]"
                % xf: [
                    50.90,
                    50.90,
                    50.90,
                    50.90,
                    50.90,
                    49.83,
                    49.31,
                    48.95,
                    48.69,
                    48.51,
                    48.38,
                    48.29,
                    48.22,
                    48.17,
                    48.14,
                    48.11,
                    48.10,
                    48.08,
                    48.07,
                    48.07,
                    48.06,
                ],
                "fs.compressor.power[*]": [
                    1.590e3,
                    1.590e3,
                    1.590e3,
                    1.590e3,
                    1.590e3,
                    1.682e3,
                    1.750e3,
                    1.796e3,
                    1.828e3,
                    1.851e3,
                    1.867e3,
                    1.878e3,
                    1.887e3,
                    1.892e3,
                    1.897e3,
                    1.900e3,
                    1.902e3,
                    1.904e3,
                    1.905e3,
                    1.906e3,
                    1.906e3,
                ],
            },
        )
        output_names = [
            "fs.pipeline.control_volume.flow_mass[*,%s]" % x0,
            "fs.pipeline.control_volume.pressure[*,%s]" % xf,
            "fs.compressor.power[*]",
        ]
        actual_values = (
            list(time),
            {
                name: [var.value for var in m.find_component(name).values()]
                for name in output_names
            },
        )
        # Note: We fail with a reltol of 0.01, due to flow rate discrepancies
        # in positions 6, 7, 8, and 9. A reltol of 0.02 seems reasonable to me.
        self.assertStructuredAlmostEqual(pred_values, actual_values, reltol=0.02)


@pytest.mark.unit
class TestConstructPipelineCompressorFlowsheet(unittest.TestCase):
    """
    Test for construction of a flowsheet consisting of one compressor
    feeding into one pipeline. This is equivalent to an "active link"
    in V. Zavala's 2014 paper.
    """

    def test_steady(self):
        """ """
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
        compressor_default = {"property_package": m.fs.properties}
        m.fs.compressor = Compressor(**compressor_default)
        compressor = m.fs.compressor

        m._compressor_to_pipeline = Arc(
            ports=(compressor.outlet_port, pipeline.inlet_port),
        )
        expand_arcs = pyo.TransformationFactory("network.expand_arcs")
        expand_arcs.apply_to(m)

        cv = m.fs.pipeline.control_volume
        p = "Vap"
        j = next(iter(m.fs.properties.component_list))
        t = m.fs.time.first()
        x0 = cv.length_domain.first()

        # Fix pipeline degrees of freedom:
        # Design variables:
        m.fs.pipeline.diameter.fix(0.92 * pyo.units.m)
        m.fs.pipeline.control_volume.length.fix(300.0 * pyo.units.m)

        # Fix compressor degrees of freedom:
        compressor.boost_pressure[:].fix()
        # Inlet variables:
        state = m.fs.compressor.inlet_state
        state[:].mole_frac_comp[j].fix(1.0)
        state[:].temperature.fix(300.0 * pyo.units.K)
        state[:].pressure.fix(57.0 * pyo.units.bar)

        # With just the pipeline unit model, we've been fixing to
        # 3e5 kg/hr. Need to convert this into kmol/hr
        ng_comp = m.fs.properties.get_component(j)
        # Here I'm assuming that we have a single component
        state[:].flow_mol.fix(3.0e5 * pyo.units.kg / pyo.units.hr / ng_comp.mw)

        # This test asserts:
        # (a) consistent units
        # (b) zero degrees of freedom
        # (c) a perfect matching between constraints and variables
        assert_units_consistent(m)
        self.assertEqual(degrees_of_freedom(m), 0)
        igraph = IncidenceGraphInterface(m)
        matching = igraph.maximum_matching()
        N, M = igraph.incidence_matrix.shape
        self.assertEqual(N, M)
        self.assertEqual(len(matching), N)

    def test_dynamic(self):
        """ """
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
            "finite_elements": 4,
        }
        m.fs.pipeline = GasPipeline(**pipeline_config)
        pipeline = m.fs.pipeline

        compressor_config = {"property_package": m.fs.properties}
        m.fs.compressor = Compressor(**compressor_config)
        compressor = m.fs.compressor

        m._pipeline_to_compressor = Arc(
            ports=(compressor.outlet_port, pipeline.inlet_port),
        )
        expand_arcs = pyo.TransformationFactory("network.expand_arcs")
        expand_arcs.apply_to(m)

        p = "Vap"
        j = next(iter(m.fs.properties.component_list))

        # Assert units consistent before discretization transformation
        assert_units_consistent(m)
        disc = pyo.TransformationFactory("dae.finite_difference")
        ntfe = 5
        disc.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme="BACKWARD")

        cv = m.fs.pipeline.control_volume
        t0 = m.fs.time.first()
        x0 = cv.length_domain.first()
        xf = cv.length_domain.last()

        # Fix geometry variables
        m.fs.pipeline.diameter.fix(0.92 * pyo.units.m)
        cv.length.fix(300.0 * pyo.units.km)

        # Fix boost pressure
        compressor.boost_pressure[:].fix()

        # Fix inlet and outlet degrees of freedom
        state = compressor.inlet_state
        state[:].pressure.fix()
        state[:].mole_frac_comp[j].fix()
        state[:].temperature.fix()
        cv.flow_mass[:, xf].fix()

        # Fix initial conditions. Here, pressure and volume for all
        # non-inlet points.
        for x in cv.length_domain:
            if x != x0:
                cv.pressure[t0, x].fix()
            if x != xf:
                cv.flow_mass[t0, x].fix()

        # I want to deactivate differential equations at (t0, xf)
        # Material balance already doesn't exist here.
        cv.momentum_balance[t0, xf].deactivate()

        igraph = IncidenceGraphInterface(m)
        N, M = igraph.incidence_matrix.shape
        matching = igraph.maximum_matching()
        self.assertEqual(degrees_of_freedom(m), 0)
        self.assertEqual(N, M)  # Sanity check
        self.assertEqual(len(matching), N)


if __name__ == "__main__":
    unittest.main()
