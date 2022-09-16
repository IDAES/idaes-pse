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
This module contains heat utility minimization functions for use in
IDAES models.
"""

__author__ = "Alejandro Garciadiego"

import pytest
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Objective,
    units as pyunits,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver

from idaes.core import VaporPhase, Component
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
import idaes.models.properties.modular_properties.pure.RPP4 as RPP4

from idaes.models.unit_models import Heater
from idaes.core import FlowsheetBlock

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from idaes.core.util.utility_minimization import (
    min_utility,
    PinchDataClass,
    heat_ex_data,
    gen_curves,
)


# -----------------------------------------------------------------------------
# [1] systematic Methods of Chemical Process Design (1997)
#     Chemical Engineering Series - L. T. Biegler, I. E. Grossmann,
#     A. W. Westerberg, page 529, Example 16.1
# Get default solver for testing
solver = get_solver()


@pytest.mark.unit
def test_PinchDataClass():

    PD = PinchDataClass(100, 80)
    PD.initQAh = {"fs.H101": 540756, "fs.H102": 719998, "fs.H103": 0, "fs.H104": 59787}
    PD.initQAc = {
        "fs.H101": 477515,
        "fs.H102": 554999,
        "fs.H103": 29875,
        "fs.H104": 119502,
    }

    assert PD.initQs == 100
    assert PD.initQw == 80
    assert PD.initQAh == {
        "fs.H101": 540756,
        "fs.H102": 719998,
        "fs.H103": 0,
        "fs.H104": 59787,
    }
    assert PD.initQAc == {
        "fs.H101": 477515,
        "fs.H102": 554999,
        "fs.H103": 29875,
        "fs.H104": 119502,
    }


# Author: Alejandro Garciadiego
class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        configuration = {
            # Specifying components
            "components": {
                "A": {
                    "type": Component,
                    "enth_mol_ig_comp": RPP4,
                    "entr_mol_ig_comp": RPP4,
                    "parameter_data": {
                        "mw": (2.016e-3, pyunits.kg / pyunits.mol),  # [1]
                        "pressure_crit": (12.9e5, pyunits.Pa),  # [1]
                        "temperature_crit": (33.2, pyunits.K),  # [1]
                        "omega": -0.218,
                        "cp_mol_ig_comp_coeff": {
                            "A": (1000000, pyunits.J / pyunits.mol / pyunits.K),
                            "B": (0, pyunits.J / pyunits.mol / pyunits.K**2),
                            "C": (0, pyunits.J / pyunits.mol / pyunits.K**3),
                            "D": (0, pyunits.J / pyunits.mol / pyunits.K**4),
                        },
                        "entr_mol_form_vap_comp_ref": (
                            0,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),
                        "enth_mol_form_vap_comp_ref": (0.0, pyunits.J / pyunits.mol),
                    },
                }
            },
            # Specifying phases
            "phases": {"Vap": {"type": VaporPhase, "equation_of_state": Ideal}},
            # Set base units of measurement
            "base_units": {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
            # Specifying state definition
            "state_definition": FTPx,
            "state_bounds": {
                "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (273.15, 300, 450, pyunits.K),
                "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
            },
            "pressure_ref": (1e5, pyunits.Pa),
            "temperature_ref": (300, pyunits.K),
        }

        # Create the ConcreteModel and the FlowsheetBlock, and attach the
        # flowsheet block to it.
        model = ConcreteModel()

        model.fs = FlowsheetBlock(dynamic=False)

        # Add properties parameter blocks to the flowsheet with specifications
        model.fs.props = GenericParameterBlock(**configuration)

        # Create an instance of the units, attaching them to the flowsheet
        # Specify that the property package to be used with with eash unit.
        model.fs.H101 = Heater(
            property_package=model.fs.props,
            has_pressure_change=False,
            has_phase_equilibrium=False,
        )

        model.fs.H102 = Heater(
            property_package=model.fs.props,
            has_pressure_change=False,
            has_phase_equilibrium=False,
        )

        model.fs.H103 = Heater(
            property_package=model.fs.props,
            has_pressure_change=False,
            has_phase_equilibrium=False,
        )

        model.fs.H104 = Heater(
            property_package=model.fs.props,
            has_pressure_change=False,
            has_phase_equilibrium=False,
        )

        model.fs.H101.inlet.mole_frac_comp[0, "A"].fix(1)
        model.fs.H101.inlet.flow_mol[0].fix(1.5)
        model.fs.H101.inlet.temperature.fix(160)
        model.fs.H101.inlet.pressure.fix(100000)

        model.fs.H101.outlet.temperature.fix(400)

        model.fs.H102.inlet.mole_frac_comp[0, "A"].fix(1)
        model.fs.H102.inlet.flow_mol[0].fix(1.3)
        model.fs.H102.inlet.temperature.fix(100)
        model.fs.H102.inlet.pressure.fix(100000)

        model.fs.H102.outlet.temperature.fix(250)

        model.fs.H103.inlet.mole_frac_comp[0, "A"].fix(1)
        model.fs.H103.inlet.flow_mol[0].fix(1)
        model.fs.H103.inlet.temperature.fix(400)
        model.fs.H103.inlet.pressure.fix(100000)

        model.fs.H103.outlet.temperature.fix(120)

        model.fs.H104.inlet.mole_frac_comp[0, "A"].fix(1)
        model.fs.H104.inlet.flow_mol[0].fix(2)
        model.fs.H104.inlet.temperature.fix(340)
        model.fs.H104.inlet.pressure.fix(100000)

        model.fs.H104.outlet.temperature.fix(120)

        assert degrees_of_freedom(model) == 0

        return model

    @pytest.mark.unit
    def test_dof(self, model):

        # add heat integration constraints
        min_utility(
            model.fs, [model.fs.H101, model.fs.H102], [model.fs.H103, model.fs.H104], 20
        )

        assert degrees_of_freedom(model) == 1

        model.fs.objective = Objective(expr=model.fs.Qs + model.fs.Qw)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.unit
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.unit
    def test_solution(self, model):
        assert model.fs.QAh["fs.H101"].value == pytest.approx(540, abs=1e-0)
        assert model.fs.QAh["fs.H102"].value == pytest.approx(719, abs=1e-0)
        assert model.fs.QAh["fs.H103"].value == pytest.approx(0, abs=1e-0)
        assert model.fs.QAh["fs.H104"].value == pytest.approx(60, abs=1e-0)
        assert model.fs.QAc["fs.H101"].value == pytest.approx(476, abs=1e-0)
        assert model.fs.QAc["fs.H102"].value == pytest.approx(554, abs=1e-0)
        assert model.fs.QAc["fs.H103"].value == pytest.approx(30, abs=1e-0)
        assert model.fs.QAc["fs.H104"].value == pytest.approx(120, abs=1e-0)
        assert model.fs.Qs.value == pytest.approx(60, abs=1e-0)
        assert model.fs.Qw.value == pytest.approx(225, abs=1e-0)

    @pytest.mark.unit
    def test_heat_ex_data(self, model):

        CD = heat_ex_data(
            model.fs, [model.fs.H101, model.fs.H102], [model.fs.H103, model.fs.H104]
        )

        assert CD.Cooling_Tin[0] == 401
        assert CD.Cooling_Tin[1] == 341
        assert CD.Cooling_Tout[0] == 120
        assert CD.Cooling_Tout[1] == 120
        assert CD.Cooling_Q[0] == -280
        assert CD.Cooling_Q[1] == -440
        assert CD.Heating_Tin[0] == 160
        assert CD.Heating_Tin[1] == 100
        assert CD.Heating_Tout[0] == 401
        assert CD.Heating_Tout[1] == 251
        assert CD.Heating_Q[0] == 360
        assert CD.Heating_Q[1] == 195

    @pytest.mark.unit
    def test_curve_data(self):

        T_test, Q_test = gen_curves([401, 341], [120, 120], [-280, -440])
        Ttest = [120, 120, 341, 401]
        for i in range(len(Ttest)):
            assert T_test[i] == Ttest[i]
        Qtest = [0, 0, -660.213, -720]
        for i in range(len(Qtest)):
            assert Q_test[i] == pytest.approx(Qtest[i], abs=1e-3)
