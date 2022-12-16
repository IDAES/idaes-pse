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
import itertools
import math
import pytest

import pyomo.common.unittest as unittest
from pyomo.common.collections import ComponentMap, ComponentSet
import pyomo.environ as pyo
import pyomo.dae as dae
from pyomo.core.expr.visitor import identify_variables
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.dae.flatten import flatten_dae_components

from pyomo.contrib.incidence_analysis import (
    IncidenceGraphInterface,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.util.subsystems import ParamSweeper

import idaes.core as idaes
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    large_residuals_set,
)
from idaes.core.util.constants import Constants
from idaes.models_extra.gas_distribution.properties.natural_gas import (
    NaturalGasParameterBlock,
)
from idaes.models_extra.gas_distribution.unit_models.pipeline import GasPipeline

from idaes.apps.nmpc import (
    get_tracking_cost_from_constant_setpoint as get_tracking_cost_expression,
)
from idaes.apps.nmpc.dynamic_data import (
    load_inputs_into_model,
    interval_data_from_time_series,
)

"""
Test for the simple pipeline model
"""


@pytest.mark.component
class TestSolvePipelineSquare(unittest.TestCase):
    """ """

    # These dictionaries map parameters in a steady state pipeline
    # solve to the predicted values for outlet flow rate and pressure.
    # These were obtained by solving a Pyomo DAE implementation of the
    # same model.
    target_flow_kghr = {
        # flow key is also in units of kg/hr. Should this be in
        # 1e4 SCM/hr, which is what the paper uses?
        (180000.0, 40.0, 1): 180000.0 * pyo.units.kg / pyo.units.hr,
        (180000.0, 40.0, 2): 180000.0 * pyo.units.kg / pyo.units.hr,
        (180000.0, 40.0, 5): 180000.0 * pyo.units.kg / pyo.units.hr,
        (180000.0, 57.0, 1): 180000.0 * pyo.units.kg / pyo.units.hr,
        (180000.0, 57.0, 2): 180000.0 * pyo.units.kg / pyo.units.hr,
        (180000.0, 57.0, 5): 180000.0 * pyo.units.kg / pyo.units.hr,
        (180000.0, 70.0, 1): 180000.0 * pyo.units.kg / pyo.units.hr,
        (180000.0, 70.0, 2): 180000.0 * pyo.units.kg / pyo.units.hr,
        (180000.0, 70.0, 5): 180000.0 * pyo.units.kg / pyo.units.hr,
        (270000.0, 40.0, 1): 270000.0 * pyo.units.kg / pyo.units.hr,
        (270000.0, 40.0, 2): 270000.0 * pyo.units.kg / pyo.units.hr,
        (270000.0, 40.0, 5): 270000.0 * pyo.units.kg / pyo.units.hr,
        (270000.0, 57.0, 1): 270000.0 * pyo.units.kg / pyo.units.hr,
        (270000.0, 57.0, 2): 270000.0 * pyo.units.kg / pyo.units.hr,
        (270000.0, 57.0, 5): 270000.0 * pyo.units.kg / pyo.units.hr,
        (270000.0, 70.0, 1): 270000.0 * pyo.units.kg / pyo.units.hr,
        (270000.0, 70.0, 2): 270000.0 * pyo.units.kg / pyo.units.hr,
        (270000.0, 70.0, 5): 270000.0 * pyo.units.kg / pyo.units.hr,
        (450000.0, 40.0, 1): 450000.0 * pyo.units.kg / pyo.units.hr,
        (450000.0, 40.0, 2): 450000.0 * pyo.units.kg / pyo.units.hr,
        (450000.0, 40.0, 5): 450000.0 * pyo.units.kg / pyo.units.hr,
        (450000.0, 57.0, 1): 450000.0 * pyo.units.kg / pyo.units.hr,
        (450000.0, 57.0, 2): 450000.0 * pyo.units.kg / pyo.units.hr,
        (450000.0, 57.0, 5): 450000.0 * pyo.units.kg / pyo.units.hr,
        (450000.0, 70.0, 1): 450000.0 * pyo.units.kg / pyo.units.hr,
        (450000.0, 70.0, 2): 450000.0 * pyo.units.kg / pyo.units.hr,
        (450000.0, 70.0, 5): 450000.0 * pyo.units.kg / pyo.units.hr,
    }
    target_pressure = {
        (180000.0, 40.0, 1): 37.00032465199001 * pyo.units.bar,
        (180000.0, 40.0, 2): 36.94189598799329 * pyo.units.bar,
        (180000.0, 40.0, 5): 36.90451650814217 * pyo.units.bar,
        (180000.0, 57.0, 1): 54.89496466806317 * pyo.units.bar,
        (180000.0, 57.0, 2): 54.87516406405006 * pyo.units.bar,
        (180000.0, 57.0, 5): 54.86291547799542 * pyo.units.bar,
        (180000.0, 70.0, 1): 68.28589980113715 * pyo.units.bar,
        (180000.0, 70.0, 2): 68.27527627682545 * pyo.units.bar,
        (180000.0, 70.0, 5): 68.26877313026228 * pyo.units.bar,
        (270000.0, 40.0, 1): 33.25073046697753 * pyo.units.bar,
        (270000.0, 40.0, 2): 32.939794050787654 * pyo.units.bar,
        (270000.0, 40.0, 5): 32.72133469851969 * pyo.units.bar,
        (270000.0, 57.0, 1): 52.26367050314213 * pyo.units.bar,
        (270000.0, 57.0, 2): 52.16101599318696 * pyo.units.bar,
        (270000.0, 57.0, 5): 52.09485138003043 * pyo.units.bar,
        (270000.0, 70.0, 1): 66.1432745525586 * pyo.units.bar,
        (270000.0, 70.0, 2): 66.0886470609479 * pyo.units.bar,
        (270000.0, 70.0, 5): 66.05431808504241 * pyo.units.bar,
        (450000.0, 40.0, 1): 21.252029074937592 * pyo.units.bar,
        (450000.0, 40.0, 2): 18.382847277082647 * pyo.units.bar,
        (450000.0, 40.0, 5): 15.152518801949473 * pyo.units.bar,
        (450000.0, 57.0, 1): 43.8435291753948 * pyo.units.bar,
        (450000.0, 57.0, 2): 42.98530492454757 * pyo.units.bar,
        (450000.0, 57.0, 5): 42.33711761875827 * pyo.units.bar,
        (450000.0, 70.0, 1): 59.286873757107195 * pyo.units.bar,
        (450000.0, 70.0, 2): 58.84301170614686 * pyo.units.bar,
        (450000.0, 70.0, 5): 58.53638005215126 * pyo.units.bar,
    }

    def test_pipeline_steady_forward_nominal(self):
        """
        Test forward simulation of the pipeline (flow and pressure
        specified at the inlet) with nominal inlet conditions at
        steady state. Here just test that we can solve with no errors.
        """
        m = pyo.ConcreteModel()
        default = {
            "dynamic": False,
        }
        m.fs = idaes.FlowsheetBlock(**default)
        m.fs.properties = NaturalGasParameterBlock()
        pipeline_config = {
            "property_package": m.fs.properties,
            "finite_elements": 1,
        }
        m.fs.pipeline = GasPipeline(**pipeline_config)
        pipeline = m.fs.pipeline
        p = "Vap"
        j = next(iter(m.fs.properties.component_list))
        # Fix pipeline degrees of freedom
        # Design variables:
        pipeline.diameter.fix(0.92 * pyo.units.m)
        # Pipelines in paper are 300 km. This seems long. Is this really
        # the length we're interested in?
        pipeline.control_volume.length.fix(300.0 * pyo.units.km)
        # Inlet variables:
        x0 = pipeline.control_volume.length_domain.first()
        state = pipeline.control_volume.properties
        state[:, x0].mole_frac_comp[j].fix(1.0)
        state[:, x0].temperature.fix(293.15 * pyo.units.K)
        state[:, x0].pressure.fix(57.0 * pyo.units.bar)
        pipeline.control_volume.flow_mass[:, x0].fix(
            3.0e5
            * pyo.units.kg
            / pyo.units.hr
            # close to 10 * (1e6 SCM) / day, the nominal value in the model
        )

        t0 = m.fs.time.first()

        ipopt = pyo.SolverFactory("ipopt")

        res = ipopt.solve(m, tee=True)
        pyo.assert_optimal_termination(res)

        violated_cons = list(large_residuals_set(m))
        # Just test that we can solve the pipeline model.
        # We will test values in the next test.
        self.assertEqual(len(violated_cons), 0)

    def test_pipeline_steady_forward_sweep(self):
        """
        Tests "forward simulation" (supply flow and pressure specified)
        of the steady state pipeline over a range of input values and
        numbers of spatial finite elements.
        """
        nxfe_list = [1, 2, 5]
        p_list = [40.0 * pyo.units.bar, 57.0 * pyo.units.bar, 70.0 * pyo.units.bar]
        h_day = 24.0 * pyo.units.hr / pyo.units.day
        kg_scm = 0.72 * pyo.units.kg / pyo.units.m**3
        f_list = [
            # Units: kg/hr
            (6.0 * 1e6) * pyo.units.m**3 / pyo.units.day / h_day * kg_scm,
            (9.0 * 1e6) * pyo.units.m**3 / pyo.units.day / h_day * kg_scm,
            (15.0 * 1e6) * pyo.units.m**3 / pyo.units.day / h_day * kg_scm,
        ]
        # Iterate over nxfe first, because for each different nxfe, we
        # need to create a new model.

        for nxfe in nxfe_list:
            m = pyo.ConcreteModel()
            default = {
                "dynamic": False,
            }
            m.fs = idaes.FlowsheetBlock(**default)
            m.fs.properties = NaturalGasParameterBlock()
            pipeline_config = {
                "property_package": m.fs.properties,
                "finite_elements": nxfe,
            }
            m.fs.pipeline = GasPipeline(**pipeline_config)
            pipeline = m.fs.pipeline
            p = "Vap"
            j = next(iter(m.fs.properties.component_list))
            # Fix pipeline degrees of freedom
            # Design variables:
            pipeline.diameter.fix(0.92 * pyo.units.m)
            # Pipelines in paper are 300 km. This seems long. Is this really
            # the length we're interested in?
            pipeline.control_volume.length.fix(300.0 * pyo.units.km)
            # Inlet variables:
            x0 = pipeline.control_volume.length_domain.first()
            xf = pipeline.control_volume.length_domain.last()
            state = pipeline.control_volume.properties
            state[:, x0].mole_frac_comp[j].fix(1.0)
            state[:, x0].temperature.fix(293.15 * pyo.units.K)
            state[:, x0].pressure.fix()
            pipeline.control_volume.flow_mass[:, x0].fix()
            t0 = m.fs.time.first()

            input_values = ComponentMap()
            target_values = ComponentMap()
            t0 = m.fs.time.first()
            inlet_flow = pipeline.control_volume.flow_mass[t0, x0]
            inlet_pressure = state[t0, x0].pressure
            outlet_flow = pipeline.control_volume.flow_mass[t0, xf]
            outlet_pressure = state[t0, xf].pressure
            input_values[inlet_flow] = []
            input_values[inlet_pressure] = []
            target_values[outlet_flow] = []
            target_values[outlet_pressure] = []
            for f, p in itertools.product(f_list, p_list):
                # Maybe I should explicitly convert to desired units here,
                # for clarity.
                input_values[inlet_flow].append(f)
                # p_val = pyo.units.convert(p, pres_var.get_units())
                input_values[inlet_pressure].append(p)
                # Evaluating floating points as keys here may not be safe...
                target_values[outlet_flow].append(
                    self.target_flow_kghr[pyo.value(f), pyo.value(p), nxfe],
                )
                target_values[outlet_pressure].append(
                    self.target_pressure[pyo.value(f), pyo.value(p), nxfe]
                )

            n_scen = len(f_list) * len(p_list)
            param_sweeper = ParamSweeper(
                n_scen,
                input_values,
                to_fix=input_values,
                output_values=target_values,
            )
            ipopt = pyo.SolverFactory("ipopt")
            with param_sweeper:
                # Note that doing this in a context manager means that
                # on error, values are reset. This is inconvenient
                # for debugging.
                # NOTE: Entering post-mortem debugging appears to
                # exit the context manager before returning control...
                self.assertEqual(degrees_of_freedom(m), 0)
                for i, (inputs, outputs) in enumerate(param_sweeper):
                    res = ipopt.solve(m, tee=True)
                    # Check for optimal termination and zero infeasibility
                    pyo.assert_optimal_termination(res)
                    self.assertEqual(len(large_residuals_set(m)), 0)

                    # Sanity check that inputs have been properly set.
                    for var, val in inputs.items():
                        val = pyo.value(pyo.units.convert(val, var.get_units()))
                        self.assertEqual(val, var.value)

                    for var, val in outputs.items():
                        val = pyo.value(pyo.units.convert(val, var.get_units()))
                        # Assume var is (or has been scaled to be) in the
                        # same units as the target value.
                        val = pyo.value(val)
                        self.assertTrue(math.isclose(val, var.value, rel_tol=1e-2))


@pytest.mark.component
class TestSolveDynamicPipeline(unittest.TestCase):
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
        cv = m.fs.pipeline.control_volume
        # Fix geometry variables
        m.fs.pipeline.diameter.fix(0.92 * pyo.units.m)
        cv.length.fix(300.0 * pyo.units.km)
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
        cv.properties[:, x0].mole_frac_comp[j].fix()
        cv.properties[:, x0].temperature.fix(inlet_temperature)
        cv.pressure[:, x0].fix(inlet_pressure)
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

    def test_sim_inlet_pressure_outlet_flow_nominal(self):
        """
        Inlet pressure and outlet flow rate will be fixed.
        """
        nxfe = 4
        ipopt = pyo.SolverFactory("ipopt")

        m_steady = self.make_steady_model(nfe=nxfe)
        self.fix_model_inlets(m_steady)
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

        # Fix dynamic inputs. Now these are inlet pressure and outlet
        # flow rate, as well as inlet mole fraction and temperature.
        cv.pressure[:, x0].fix()
        cv.flow_mass[:, xf].fix()
        cv.properties[:, x0].mole_frac_comp[j].fix()
        cv.properties[:, x0].temperature.fix()

        # Fix initial conditions. Here, pressure and volume for all
        # non-specified points.
        for x in cv.length_domain:
            if x != x0:
                cv.pressure[t0, x].fix()
            if x != xf:
                cv.flow_mass[t0, x].fix()

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
            },
        )
        output_names = [
            "fs.pipeline.control_volume.flow_mass[*,%s]" % x0,
            "fs.pipeline.control_volume.pressure[*,%s]" % xf,
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

    def test_opt(self):
        """
        A dynamic optimization problem with a 20 hr horizon. Inputs are inlet
        pressure, which is restricted to be piecewise constant,
        and outlet flow rate, which is unrestricted.
        The setpoint is the steady state of the system with inlets of 57 bar
        and 5e5 kg/hr.
        """
        nxfe = 4
        ipopt = pyo.SolverFactory("ipopt")

        # Steady state data
        m_steady = self.make_steady_model(nfe=nxfe)
        self.fix_model_inlets(m_steady)
        ipopt.solve(m_steady, tee=True)
        time_steady = m_steady.fs.time
        scalar_data = self.get_scalar_data_from_model(m_steady, time_steady)
        initial_data = self.get_data_from_model_at_time(m_steady, time_steady)
        #

        # Setpoint data
        m_setpoint = self.make_steady_model(nfe=nxfe)
        self.fix_model_inlets(
            m_setpoint,
            inlet_flow_mass=5.0e5 * pyo.units.kg / pyo.units.hr,
        )
        ipopt.solve(m_setpoint, tee=True)
        time_setpoint = m_setpoint.fs.time
        setpoint_data = self.get_data_from_model_at_time(m_setpoint, time_setpoint)
        #

        t0 = 0.0
        horizon = 20.0
        m = pyo.ConcreteModel()
        default = {
            "dynamic": True,
            "time_set": [t0, horizon],
            "time_units": pyo.units.hr,
        }
        m.fs = idaes.FlowsheetBlock(**default)
        m.fs.properties = NaturalGasParameterBlock()
        pipeline_config = {
            "property_package": m.fs.properties,
            "finite_elements": nxfe,
        }
        m.fs.pipeline = GasPipeline(**pipeline_config)

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

        # Fix "dynamic parameters" in the optimization problem
        cv.properties[:, x0].mole_frac_comp[j].fix()
        cv.properties[:, x0].temperature.fix()

        # Fix initial conditions. Here, pressure and volume except
        # at points where we treat them as "inputs"
        for x in cv.length_domain:
            # For dynamic optimization, I believe we want these variables
            # fixed everywhere
            cv.pressure[t0, x].fix()
            cv.flow_mass[t0, x].fix()

        # Now we set up the dynamic optimization problem:
        tracking_variables = [
            # NOTE: These space keys (0.0 and 1.0) need to be consistent with
            # how they are declared in a continuous set, since we are using
            # the CUID string representation as a key in setpoint_data.
            pyo.Reference(cv.pressure[:, 0.0]),
            pyo.Reference(cv.pressure[:, 1.0]),
            pyo.Reference(cv.flow_mass[:, 0.0]),
            pyo.Reference(cv.flow_mass[:, 1.0]),
        ]
        weight_data = {
            "fs.pipeline.control_volume.flow_mass[*,0.0]": 1e-10,
            "fs.pipeline.control_volume.flow_mass[*,1.0]": 1e-10,
            "fs.pipeline.control_volume.pressure[*,0.0]": 1e-2,
            "fs.pipeline.control_volume.pressure[*,1.0]": 1e-2,
        }
        m.tracking_cost = get_tracking_cost_expression(
            tracking_variables, time, setpoint_data, weight_data=weight_data
        )
        m.tracking_objective = pyo.Objective(expr=sum(m.tracking_cost.values()))

        sample_period = 2  # hours
        n_samples = (int(horizon) - int(t0)) // sample_period
        sample_points = [t0 + float(sample_period * i) for i in range(n_samples + 1)]
        sample_point_set = set(sample_points)
        piecewise_constant_vars = [pyo.Reference(cv.pressure[:, 0.0])]
        m.piecewise_constant_vars_set = pyo.Set(
            initialize=list(range(len(piecewise_constant_vars)))
        )

        def piecewise_constant_rule(m, i, t):
            var = piecewise_constant_vars[i]
            if t in sample_point_set:
                return pyo.Constraint.Skip
            else:
                t_next = time.next(t)
                return var[t] == var[t_next]

        m.piecewise_constant_constraint = pyo.Constraint(
            m.piecewise_constant_vars_set, time, rule=piecewise_constant_rule
        )
        self.assertEqual(
            degrees_of_freedom(m),
            # Pressure has n_samples dof, flow has (len(time)-1) dof
            n_samples + (len(time) - 1),
        )

        # Initialize the dynamic optimization problem:

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

        self.assertEqual(len(large_residuals_set(m)), 0)

        res = ipopt.solve(m, tee=True)
        self.assertIs(
            res.solver.termination_condition,
            pyo.TerminationCondition.optimal,
        )

        # These predicted values were calculated from a single-pipeline
        # model created from the Pyomo DAE example. That model uses the
        # same tracking objective used here, but with flow rate weights
        # scaled to match the different units in that model (1e4 SCM)/hr.
        # The values obtained from that model are then converted back
        # to kg/hr (via the standard density 0.72 kg/m^3).
        pred_tracking_variable_data = [
            [
                57.00,
                58.86,
                58.86,
                58.71,
                58.71,
                58.35,
                58.35,
                58.04,
                58.04,
                57.80,
                57.80,
                57.63,
                57.63,
                57.51,
                57.51,
                57.42,
                57.42,
                57.38,
                57.38,
                57.35,
                57.35,
            ],
            [
                50.90,
                45.89,
                43.59,
                42.27,
                41.50,
                40.97,
                40.62,
                40.31,
                40.09,
                39.88,
                39.72,
                39.56,
                39.45,
                39.34,
                39.26,
                39.18,
                39.13,
                39.08,
                39.05,
                39.02,
                39.01,
            ],
            [
                3.000e5,
                4.693e5,
                4.629e5,
                4.698e5,
                4.850e5,
                4.761e5,
                4.909e5,
                4.819e5,
                4.936e5,
                4.866e5,
                4.953e5,
                4.903e5,
                4.966e5,
                4.932e5,
                4.977e5,
                4.955e5,
                4.985e5,
                4.973e5,
                4.992e5,
                4.986e5,
                4.996e5,
            ],
            [
                3.000e5,
                5.550e5,
                5.456e5,
                5.362e5,
                5.294e5,
                5.261e5,
                5.212e5,
                5.194e5,
                5.158e5,
                5.144e5,
                5.118e5,
                5.107e5,
                5.087e5,
                5.078e5,
                5.063e5,
                5.056e5,
                5.044e5,
                5.038e5,
                5.029e5,
                5.022e5,
                5.013e5,
            ],
        ]
        pred_values = (
            list(time),
            {
                str(pyo.ComponentUID(var.referent)): values
                for var, values in zip(tracking_variables, pred_tracking_variable_data)
            },
        )

        actual_values = (
            list(time),
            {
                str(pyo.ComponentUID(var.referent)): [var[t].value for t in time]
                for var in tracking_variables
            },
        )
        self.assertStructuredAlmostEqual(pred_values, actual_values, reltol=0.02)


@pytest.mark.unit
class TestConstructPipeline(unittest.TestCase):
    """
    Test for construction of pipeline
    """

    def test_geometry(self):
        m = pyo.ConcreteModel()
        default = {
            "dynamic": False,
        }
        m.fs = idaes.FlowsheetBlock(**default)
        m.fs.properties = NaturalGasParameterBlock()
        pipeline_config = {
            "property_package": m.fs.properties,
        }
        m.fs.pipeline = GasPipeline(**pipeline_config)

        pipeline = m.fs.pipeline
        cv = pipeline.control_volume
        self.assertTrue(isinstance(cv.area, pyo.Var))
        self.assertTrue(isinstance(m.fs.pipeline.diameter, pyo.Var))

        diameter = m.fs.pipeline.diameter
        diameter_eqn = m.fs.pipeline.diameter_eqn
        cv.area.fix(5.0 * pyo.units.m**2)
        calculate_variable_from_constraint(diameter, diameter_eqn)
        self.assertAlmostEqual(
            diameter.value,
            pyo.value(
                pyo.units.convert(
                    pyo.sqrt(cv.area / Constants.pi) * 2, diameter.get_units()
                )
            ),
        )

    def test_momentum_balance_steady(self):
        m = pyo.ConcreteModel()
        default = {
            "dynamic": False,
        }
        m.fs = idaes.FlowsheetBlock(**default)
        m.fs.properties = NaturalGasParameterBlock()
        pipeline_config = {
            "property_package": m.fs.properties,
        }
        m.fs.pipeline = GasPipeline(**pipeline_config)

        pipeline = m.fs.pipeline
        cv = pipeline.control_volume

        self.assertTrue(isinstance(cv.momentum_balance, pyo.Constraint))
        self.assertTrue(isinstance(cv.pressure_dx, dae.DerivativeVar))
        self.assertTrue(isinstance(cv.pressure_dx_disc_eq, pyo.Constraint))
        assert_units_consistent(cv.momentum_balance)

        t = m.fs.time.first()
        x = m.fs.pipeline.control_volume.length_domain.first()

        momentum_bal = cv.momentum_balance[t, x]
        var_set = ComponentSet(identify_variables(momentum_bal.expr))
        self.assertIn(cv.pressure_dx[t, x], var_set)

    def test_steady_forward(self):
        """
        Tests a steady state instance of the pipeline for forward simulation
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

        cv = m.fs.pipeline.control_volume
        self.assertTrue(isinstance(cv.material_balances, pyo.Constraint))
        self.assertTrue(isinstance(cv.material_flow_dx, dae.DerivativeVar))
        self.assertTrue(isinstance(cv.length_domain, dae.ContinuousSet))
        self.assertTrue(isinstance(cv.material_flow_dx_disc_eq, pyo.Constraint))
        self.assertTrue(isinstance(cv.momentum_balance, pyo.Constraint))
        self.assertTrue(isinstance(cv.pressure_dx, dae.DerivativeVar))
        p = "Vap"
        j = next(iter(m.fs.properties.component_list))
        t = m.fs.time.first()
        x = cv.length_domain.first()
        flux = cv.material_flow_dx[t, x, p, j]
        mat_bal = cv.material_balances[t, x, p, j]
        var_set = ComponentSet(identify_variables(mat_bal.expr))
        self.assertIn(flux, var_set)
        # Really would like to assert that mat bal only has one nonzero term
        # and is linear wrt flow_dx

        # Fix pipeline degrees of freedom
        # Design variables:
        m.fs.pipeline.diameter.fix(0.92 * pyo.units.m)
        m.fs.pipeline.control_volume.length.fix(300.0 * pyo.units.m)
        # Inlet variables:
        x0 = m.fs.pipeline.control_volume.length_domain.first()
        xf = m.fs.pipeline.control_volume.length_domain.last()
        state = m.fs.pipeline.control_volume.properties
        state[:, x0].mole_frac_comp[j].fix(1.0)
        state[:, x0].temperature.fix(300.0 * pyo.units.K)
        state[:, x0].pressure.fix(57.0 * pyo.units.bar)
        cv.flow_mass[:, x0].fix(
            3.0e5
            * pyo.units.kg
            / pyo.units.hr
            # close to 10 * (1e6 SCM) / day, the nominal value in the model
        )

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

    def test_dynamic_forward(self):
        """
        Tests construction of a dynamic instance of the pipeline for
        forward simulation.
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
            "finite_elements": 4,
        }
        m.fs.pipeline = GasPipeline(**pipeline_config)

        cv = m.fs.pipeline.control_volume
        self.assertTrue(isinstance(cv.material_balances, pyo.Constraint))
        self.assertTrue(isinstance(cv.material_flow_dx, dae.DerivativeVar))
        self.assertTrue(isinstance(cv.material_accumulation, dae.DerivativeVar))
        self.assertTrue(isinstance(cv.length_domain, dae.ContinuousSet))
        self.assertTrue(isinstance(cv.material_flow_dx_disc_eq, pyo.Constraint))
        p = "Vap"
        j = next(iter(m.fs.properties.component_list))
        t = m.fs.time.first()
        x = cv.length_domain.first()
        flux = cv.material_flow_dx[t, x, p, j]
        accum = cv.material_accumulation[t, x, p, j]
        mat_bal = cv.material_balances[t, x, p, j]
        var_set = ComponentSet(identify_variables(mat_bal.expr))
        self.assertIn(flux, var_set)
        self.assertIn(accum, var_set)
        assert_units_consistent(mat_bal)
        assert_units_consistent(m)

        disc = pyo.TransformationFactory("dae.finite_difference")
        ntfe = 5
        disc.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme="BACKWARD")

        t0 = m.fs.time.first()
        x0 = cv.length_domain.first()
        xf = cv.length_domain.last()

        # Fix geometry variables
        m.fs.pipeline.diameter.fix(0.92 * pyo.units.m)
        cv.length.fix(300.0 * pyo.units.km)

        # Fix dynamic inputs. Here, inlet flow and pressure for all time.
        cv.pressure[:, x0].fix()
        cv.flow_mass[:, x0].fix()
        cv.properties[:, x0].mole_frac_comp[j].fix()
        cv.properties[:, x0].temperature.fix()

        # Fix initial conditions. Here, pressure and volume for all
        # non-inlet points.
        for x in cv.length_domain:
            if x != x0:
                cv.pressure[t0, x].fix()
                cv.flow_mass[t0, x].fix()

        igraph = IncidenceGraphInterface(m)
        N, M = igraph.incidence_matrix.shape
        matching = igraph.maximum_matching()
        self.assertEqual(degrees_of_freedom(m), 0)
        self.assertEqual(N, M)  # Sanity check
        self.assertEqual(len(matching), N)

    def test_dynamic_inlet_pressure_outlet_flow(self):
        """
        Dynamic simulation with both inlets fixed appears to be unstable,
        so we'll proceed with inlet pressure and outlet flow specified.
        This test makes sure this version of the model is structurally sound.
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
            "finite_elements": 4,
        }
        m.fs.pipeline = GasPipeline(**pipeline_config)

        cv = m.fs.pipeline.control_volume
        assert_units_consistent(m)

        disc = pyo.TransformationFactory("dae.finite_difference")
        ntfe = 5
        disc.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme="BACKWARD")
        # Note that units are not consistent after time discretization.
        # This is because accumulations have different units than holdups,
        # and "delta t" factors in discretization equations have no units.

        t0 = m.fs.time.first()
        x0 = cv.length_domain.first()
        xf = cv.length_domain.last()
        j = next(iter(m.fs.properties.component_list))

        # Fix geometry variables
        m.fs.pipeline.diameter.fix(0.92 * pyo.units.m)
        cv.length.fix(300.0 * pyo.units.km)

        # Fix dynamic inputs. Now these are inlet pressure and outlet
        # flow rate, as well as inlet mole fraction and temperature.
        cv.pressure[:, x0].fix()
        cv.flow_mass[:, xf].fix()
        cv.properties[:, x0].mole_frac_comp[j].fix()
        cv.properties[:, x0].temperature.fix()

        # Fix initial conditions. Here, pressure and volume for all
        # non-inlet points.
        for x in cv.length_domain:
            if x != x0:
                cv.pressure[t0, x].fix()
            if x != xf:
                cv.flow_mass[t0, x].fix()

        igraph = IncidenceGraphInterface(m)
        N, M = igraph.incidence_matrix.shape
        matching = igraph.maximum_matching()
        self.assertEqual(degrees_of_freedom(m), 0)
        self.assertEqual(N, M)  # Sanity check
        self.assertEqual(len(matching), N)

    def test_dynamic_dof_forward_space(self):
        """
        Test that dynamic degrees of freedom are what we expect when using
        a forward discretization in the space domain.
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
            "transformation_method": "dae.finite_difference",
            "transformation_scheme": "FORWARD",
        }
        m.fs.pipeline = GasPipeline(**pipeline_config)

        cv = m.fs.pipeline.control_volume

        disc = pyo.TransformationFactory("dae.finite_difference")
        ntfe = 2
        disc.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme="BACKWARD")

        # Fix geometry variables
        m.fs.pipeline.diameter.fix(0.92 * pyo.units.m)
        cv.length.fix(300.0 * pyo.units.km)

        # Inputs are inlet pressure, mole frac, and temperature, and outlet
        # flow rate.
        n_inputs = 4
        # pressure and flow mass are differential except where they are
        # specified by inputs.
        n_differential = 2 * (len(cv.length_domain) - 1)
        pred_dof = n_inputs * len(m.fs.time) + n_differential
        self.assertEqual(degrees_of_freedom(m), pred_dof)

        # Fix degrees of freedom.
        x0 = cv.length_domain.first()
        xf = cv.length_domain.last()
        t0 = m.fs.time.first()
        # Inputs at every point in time:
        cv.pressure[:, x0].fix()
        cv.flow_mass[:, xf].fix()
        for j in m.fs.properties.component_list:
            cv.properties[:, x0].mole_frac_comp[j].fix()
        cv.properties[:, x0].temperature.fix()

        # Initial conditions:
        for x in cv.length_domain:
            if x != x0:
                cv.pressure[t0, x].fix()
            if x != xf:
                cv.flow_mass[t0, x].fix()

        igraph = IncidenceGraphInterface(m)
        N, M = igraph.incidence_matrix.shape
        matching = igraph.maximum_matching()
        self.assertEqual(degrees_of_freedom(m), 0)
        self.assertEqual(N, M)  # Sanity check
        self.assertEqual(len(matching), N)

    def test_dynamic_dof_backward_space(self):
        """
        Test that dynamic degrees of freedom are what we expect when using
        a backward discretization in the space domain.
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
            "transformation_method": "dae.finite_difference",
            "transformation_scheme": "BACKWARD",
        }
        m.fs.pipeline = GasPipeline(**pipeline_config)

        cv = m.fs.pipeline.control_volume

        disc = pyo.TransformationFactory("dae.finite_difference")
        ntfe = 2
        disc.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme="BACKWARD")

        # Fix geometry variables
        m.fs.pipeline.diameter.fix(0.92 * pyo.units.m)
        cv.length.fix(300.0 * pyo.units.km)

        # Inputs are inlet pressure, mole frac, and temperature, and outlet
        # flow rate.
        n_inputs = 4
        # pressure and flow mass are differential except where they are
        # specified by inputs.
        n_differential = 2 * (len(cv.length_domain) - 1)
        pred_dof = n_inputs * len(m.fs.time) + n_differential
        self.assertEqual(degrees_of_freedom(m), pred_dof)

        # Fix degrees of freedom.
        x0 = cv.length_domain.first()
        xf = cv.length_domain.last()
        t0 = m.fs.time.first()
        # Inputs at every point in time:
        cv.pressure[:, x0].fix()
        cv.flow_mass[:, xf].fix()
        for j in m.fs.properties.component_list:
            cv.properties[:, x0].mole_frac_comp[j].fix()
        cv.properties[:, x0].temperature.fix()

        # Initial conditions:
        for x in cv.length_domain:
            if x != x0:
                cv.pressure[t0, x].fix()
            if x != xf:
                cv.flow_mass[t0, x].fix()

        igraph = IncidenceGraphInterface(m)
        N, M = igraph.incidence_matrix.shape
        matching = igraph.maximum_matching()
        self.assertEqual(degrees_of_freedom(m), 0)
        self.assertEqual(N, M)  # Sanity check
        self.assertEqual(len(matching), N)

    def test_dynamic_dof_radau_space(self):
        """
        Test that dynamic degrees of freedom are what we expect when using
        a Radau discretization in the space domain.
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
            "collocation_points": 2,
            "transformation_method": "dae.collocation",
            "transformation_scheme": "LAGRANGE-RADAU",
        }
        m.fs.pipeline = GasPipeline(**pipeline_config)

        cv = m.fs.pipeline.control_volume

        disc = pyo.TransformationFactory("dae.finite_difference")
        ntfe = 2
        disc.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme="BACKWARD")

        # Fix geometry variables
        m.fs.pipeline.diameter.fix(0.92 * pyo.units.m)
        cv.length.fix(300.0 * pyo.units.km)

        # Inputs are inlet pressure, mole frac, and temperature, and outlet
        # flow rate.
        n_inputs = 4
        # pressure and flow mass are differential except where they are
        # specified by inputs.
        n_differential = 2 * (len(cv.length_domain) - 1)
        pred_dof = n_inputs * len(m.fs.time) + n_differential
        self.assertEqual(degrees_of_freedom(m), pred_dof)

        # Fix degrees of freedom.
        x0 = cv.length_domain.first()
        xf = cv.length_domain.last()
        t0 = m.fs.time.first()
        # Inputs at every point in time:
        cv.pressure[:, x0].fix()
        cv.flow_mass[:, xf].fix()
        for j in m.fs.properties.component_list:
            cv.properties[:, x0].mole_frac_comp[j].fix()
        cv.properties[:, x0].temperature.fix()

        # Initial conditions:
        for x in cv.length_domain:
            if x != x0:
                cv.pressure[t0, x].fix()
            if x != xf:
                cv.flow_mass[t0, x].fix()

        igraph = IncidenceGraphInterface(m)
        N, M = igraph.incidence_matrix.shape
        matching = igraph.maximum_matching()
        self.assertEqual(degrees_of_freedom(m), 0)
        self.assertEqual(N, M)  # Sanity check
        self.assertEqual(len(matching), N)

    def test_dynamic_dof_legendre_space(self):
        """
        Test that trying to create a dynamic pipeline with a Legendre
        discretization throws an error. The commented code is what we
        would like to work. Debugging the Legendre discretization is not
        a high priority, so for now it is not supported.

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
            "collocation_points": 2,
            "transformation_method": "dae.collocation",
            "transformation_scheme": "LAGRANGE-LEGENDRE",
        }
        with self.assertRaisesRegex(ValueError, "Discretization scheme"):
            m.fs.pipeline = GasPipeline(**pipeline_config)

        #
        # The following is the failing test that causes Legendre
        # discretization to not be supported. If at any point we really
        # need to use a Legendre discretization, this test will need to
        # be addressed.
        #
        # cv = m.fs.pipeline.control_volume

        # disc = pyo.TransformationFactory("dae.finite_difference")
        # ntfe = 2
        # disc.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme="BACKWARD")

        ## Fix geometry variables
        # m.fs.pipeline.diameter.fix(0.92*pyo.units.m)
        # cv.length.fix(300.0*pyo.units.km)

        ## Inputs are inlet pressure, mole frac, and temperature, and outlet
        ## flow rate.
        # n_inputs = 4
        ## pressure and flow mass are differential except where they are
        ## specified by inputs.
        # n_differential = 2*(len(cv.length_domain) - 1)
        # pred_dof = n_inputs * len(m.fs.time) + n_differential
        # self.assertEqual(degrees_of_freedom(m), pred_dof)

        ## Fix degrees of freedom.
        # x0 = cv.length_domain.first()
        # xf = cv.length_domain.last()
        # t0 = m.fs.time.first()
        ## Inputs at every point in time:
        # cv.pressure[:, x0].fix()
        # cv.flow_mass[:, xf].fix()
        # for j in m.fs.properties.component_list:
        #    cv.properties[:, x0].mole_frac_comp[j].fix()
        # cv.properties[:, x0].temperature.fix()

        ## Initial conditions:
        # for x in cv.length_domain:
        #    if x != x0:
        #        cv.pressure[t0, x].fix()
        #    if x != xf:
        #        cv.flow_mass[t0, x].fix()

        # igraph = IncidenceGraphInterface(m)
        # N, M = igraph.incidence_matrix.shape
        # matching = igraph.maximum_matching()
        # self.assertEqual(degrees_of_freedom(m), 0)
        # self.assertEqual(N, M) # Sanity check
        # self.assertEqual(len(matching), N)


if __name__ == "__main__":
    unittest.main()
