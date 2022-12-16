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
import pytest

import pyomo.common.unittest as unittest
from pyomo.common.collections import ComponentMap
import pyomo.environ as pyo

from pyomo.contrib.incidence_analysis import (
    IncidenceGraphInterface,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.util.subsystems import ParamSweeper

import idaes.core as idaes
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models_extra.gas_distribution.properties.natural_gas import (
    NaturalGasParameterBlock,
)
from idaes.models_extra.gas_distribution.unit_models.compressor import (
    IsothermalCompressor,
)

"""
Test for the simple compressor.
"""


@pytest.mark.component
class TestCompressorValues(unittest.TestCase):
    """
    This class tests whether the compressor correctly computes power
    requirement.
    """

    # Here we record predicted values for power requirement calculated
    # from each set of inputs.
    # These values were calculated with a pure-Pyomo simulation of the
    # compressor.
    # This dict maps pressure (bar), boost pressure (bar), and flow
    # rate (SCM/hr) to power requirement (kW).
    # Flow rate will later be converted to kmol/hr.
    expected_power = {
        (30.0, 2.0, 30.0e4): 560.0079273449446,
        (30.0, 2.0, 40.0e4): 746.6772364599237,
        (30.0, 2.0, 50.0e4): 933.3465455749101,
        (30.0, 2.0, 60.0e4): 1120.0158546898892,
        (30.0, 5.0, 30.0e4): 1350.2398277565735,
        (30.0, 5.0, 40.0e4): 1800.3197703420956,
        (30.0, 5.0, 50.0e4): 2250.399712927625,
        (30.0, 5.0, 60.0e4): 2700.479655513147,
        (30.0, 10.0, 30.0e4): 2555.623582229353,
        (30.0, 10.0, 40.0e4): 3407.498109639142,
        (30.0, 10.0, 50.0e4): 4259.372637048931,
        (30.0, 10.0, 60.0e4): 5111.247164458706,
        (40.0, 2.0, 30.0e4): 422.65881875472405,
        (40.0, 2.0, 40.0e4): 563.5450916729678,
        (40.0, 2.0, 50.0e4): 704.4313645912043,
        (40.0, 2.0, 60.0e4): 845.3176375094481,
        (40.0, 5.0, 30.0e4): 1027.7481241674832,
        (40.0, 5.0, 40.0e4): 1370.3308322233133,
        (40.0, 5.0, 50.0e4): 1712.9135402791435,
        (40.0, 5.0, 60.0e4): 2055.4962483349664,
        (40.0, 10.0, 30.0e4): 1968.8273762191675,
        (40.0, 10.0, 40.0e4): 2625.1031682922257,
        (40.0, 10.0, 50.0e4): 3281.378960365284,
        (40.0, 10.0, 60.0e4): 3937.654752438335,
        (50.0, 2.0, 30.0e4): 339.4199310524127,
        (50.0, 2.0, 40.0e4): 452.5599080698812,
        (50.0, 2.0, 50.0e4): 565.6998850873497,
        (50.0, 2.0, 60.0e4): 678.8398621048254,
        (50.0, 5.0, 30.0e4): 829.6932850445155,
        (50.0, 5.0, 40.0e4): 1106.2577133926825,
        (50.0, 5.0, 50.0e4): 1382.8221417408495,
        (50.0, 5.0, 60.0e4): 1659.386570089031,
        (50.0, 10.0, 30.0e4): 1601.739808280763,
        (50.0, 10.0, 40.0e4): 2135.653077707684,
        (50.0, 10.0, 50.0e4): 2669.566347134605,
        (50.0, 10.0, 60.0e4): 3203.479616561526,
        (60.0, 2.0, 30.0e4): 283.57512892647355,
        (60.0, 2.0, 40.0e4): 378.1001719019623,
        (60.0, 2.0, 50.0e4): 472.62521487745107,
        (60.0, 2.0, 60.0e4): 567.1502578529471,
        (60.0, 5.0, 30.0e4): 695.670690390034,
        (60.0, 5.0, 40.0e4): 927.5609205200453,
        (60.0, 5.0, 50.0e4): 1159.4511506500567,
        (60.0, 5.0, 60.0e4): 1391.341380780068,
        (60.0, 10.0, 30.0e4): 1350.2398277565735,
        (60.0, 10.0, 40.0e4): 1800.3197703420956,
        (60.0, 10.0, 50.0e4): 2250.399712927625,
        (60.0, 10.0, 60.0e4): 2700.479655513147,
    }

    def test_calculate_power(self):
        m = pyo.ConcreteModel()
        default = {
            "dynamic": False,
        }
        m.fs = idaes.FlowsheetBlock(**default)
        m.fs.properties = NaturalGasParameterBlock()
        compressor_config = {
            "property_package": m.fs.properties,
        }
        m.fs.compressor = IsothermalCompressor(**compressor_config)

        time = m.fs.time
        t0 = m.fs.time.first()

        # Preliminary constants and units
        # TODO: Should this nominal density constant live on the model?
        # Probably, even though it's not really necessary.
        density = 0.72 * pyo.units.kg / pyo.units.m**3
        kghr = pyo.units.kg / pyo.units.hr
        K = pyo.units.K
        bar = pyo.units.bar
        mw = m.fs.compressor.inlet_state[t0].mw
        scmhr = pyo.units.m**3 / pyo.units.hr
        kW = pyo.units.kW

        pressure_list = [30.0 * bar, 40.0 * bar, 50.0 * bar, 60.0 * bar]
        boost_pressure_list = [2.0 * bar, 5.0 * bar, 10.0 * bar]
        flow_list = [30.0e4 * scmhr, 40.0e4 * scmhr, 50.0e4 * scmhr, 60.0e4 * scmhr]
        n_scen = len(pressure_list) * len(boost_pressure_list) * len(flow_list)
        j = next(iter(m.fs.properties.component_list))
        inlet_values = {
            "temperature": [],
            "flow_mol": [],
            "pressure": [],
            "mole_frac_comp[%s]" % j: [],
        }
        input_values = {
            "boost_pressure[%s]" % t0: [],
        }
        target_values = {
            "power[%s]" % t0: [],
        }
        for p, dp, f in itertools.product(
            pressure_list,
            boost_pressure_list,
            flow_list,
        ):

            # Redundant unit conversion, for sanity.
            p_val = pyo.value(pyo.units.convert(p, bar))
            dp_val = pyo.value(pyo.units.convert(dp, bar))
            f_val = pyo.value(pyo.units.convert(f, scmhr))
            target_values["power[%s]" % t0].append(
                self.expected_power[p_val, dp_val, f_val] * kW
            )

            # IDAES generic property package requires kmol/hr.
            # I would prefer kg/hr...
            # SCM/hr -> kg/hr -> kmol/hr
            f = f * density / mw
            inlet_values["temperature"].append(293.15 * K)
            inlet_values["flow_mol"].append(f)
            inlet_values["pressure"].append(p)
            inlet_values["mole_frac_comp[%s]" % j].append(1.0)
            input_values["boost_pressure[%s]" % t0].append(dp)
        state = m.fs.compressor.inlet_state[t0]
        inlet_values = ComponentMap(
            (state.find_component(key), val) for key, val in inlet_values.items()
        )
        input_values = ComponentMap(
            (m.fs.compressor.find_component(key), val)
            for key, val in input_values.items()
        )
        input_values.update(inlet_values)
        target_values = ComponentMap(
            (m.fs.compressor.find_component(key), val)
            for key, val in target_values.items()
        )

        ipopt = pyo.SolverFactory("ipopt")
        param_sweeper = ParamSweeper(
            n_scen,
            input_values,
            output_values=target_values,
            to_fix=input_values,
        )
        with param_sweeper:
            self.assertEqual(degrees_of_freedom(m), 0)
            for inputs, target in param_sweeper:
                ipopt.solve(m, tee=False)

                # Sanity check that inputs have been properly set.
                for var, val in inputs.items():
                    val = pyo.value(pyo.units.convert(val, var.get_units()))
                    self.assertEqual(val, var.value)

                # Check that calculated power requirements are close to what
                # we expect. Values are different due to discrepancy in cv
                # between paper and IDAES.
                for var, val in target.items():
                    val = pyo.value(pyo.units.convert(val, var.get_units()))
                    # self.assertAlmostEqual(val, var.value, delta=1e2)
                    # self.assertAlmostEqual(val, var.value, reltol=0.1)
                    # Forget the relative tolerance argument for
                    # assertAlmostEqual, so using pytest.approx for now...
                    self.assertEqual(val, pytest.approx(var.value, rel=0.1))
                    # We use a pretty generous margin here because the
                    # property package's heat capacity cv is different
                    # from that used in the paper. (Paper does not use
                    # exactly (cp - R))


@pytest.mark.unit
class TestSimpleCompressor(unittest.TestCase):
    def test_compressor_methane(self):
        """
        Test constructing the compressor with a simple methane
        property package.
        """
        m = pyo.ConcreteModel()
        default = {
            "dynamic": False,
        }
        m.fs = idaes.FlowsheetBlock(**default)
        m.fs.properties = NaturalGasParameterBlock()
        compressor_config = {
            "property_package": m.fs.properties,
        }
        m.fs.compressor = IsothermalCompressor(**compressor_config)

        m.fs.compressor.inlet_state[0].temperature.fix(300.0 * pyo.units.K)
        m.fs.compressor.inlet_state[0].pressure.fix(20.0 * pyo.units.bar)
        m.fs.compressor.boost_pressure[0].fix(10.0 * pyo.units.bar)
        m.fs.compressor.inlet_state[0].flow_mol.fix(1000.0)
        # Fix mole_frac as defined_state is True at inlet
        j = next(iter(m.fs.properties.component_list))
        m.fs.compressor.inlet_state[0].mole_frac_comp[j].fix(1.0)

        # This test asserts:
        # (a) consistent units
        # (b) zero degrees of freedom
        # (c) a perfect matching between constraints and variables
        assert_units_consistent(m)

        self.assertEqual(degrees_of_freedom(m), 0)
        igraph = IncidenceGraphInterface(m)
        N, M = igraph.incidence_matrix.shape
        self.assertEqual(N, M)
        matching = igraph.maximum_matching()
        self.assertEqual(len(matching), N)


if __name__ == "__main__":
    unittest.main()
