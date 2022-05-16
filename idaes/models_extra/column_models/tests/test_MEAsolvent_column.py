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
Tests for MEA solvent column model
Author: Anuja Deshpande, Andrew Lee
"""
# Import Python libraries
import copy
import pytest

# Import Pyomo libraries
from pyomo.environ import (ConcreteModel,
                           value,
                           Constraint,
                           check_optimal_termination,
                           units as pyunits)

# Import IDAES Libraries
from idaes.core import FlowsheetBlock
from idaes.models_extra.column_models.MEAsolvent_column import MEAColumn
from idaes.models.properties.modular_properties.base.generic_property import (
        GenericParameterBlock)
from idaes.models_extra.column_models.properties.MEA_vapor \
    import flue_gas as vaporconfig
from idaes.models_extra.column_models.properties.MEA_solvent \
    import configuration as liquidconfig

from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver

import idaes.logger as idaeslog


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()
solver.options["bound_push"] = 1e-10

class TestAbsorber:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
    
        # Set up molar flow bounds
        vaporconfig_copy = copy.deepcopy(vaporconfig)
        liquidconfig_copy = copy.deepcopy(liquidconfig)
        vaporconfig_copy["state_bounds"]["flow_mol"] = (0, 1, 1e6, pyunits.mol / pyunits.s)
        liquidconfig_copy["state_bounds"]["flow_mol"] = (0, 1, 1e6, pyunits.mol / pyunits.s)
    
        # Set up property package
        m.fs.vapor_properties = GenericParameterBlock(default=vaporconfig_copy)
        m.fs.liquid_properties = GenericParameterBlock(default=liquidconfig_copy)
    
        # Number of finite elements and finite element list in the spatial domain
        x_nfe = 40
    
        # Create an instance of the column in the flowsheet
        m.fs.unit = MEAColumn(default={
            "finite_elements": x_nfe,
            "vapor_phase": {
                "property_package": m.fs.vapor_properties},
            "liquid_phase":{
                "property_package": m.fs.liquid_properties}})
        
        # Fix column design variables
        # Absorber diameter 
        m.fs.unit.diameter_column.fix(18)
        
        # Absorber length
        m.fs.unit.length_column.fix(21.6) # meter
        
        # Fix operating conditions
        
        # Flue gas inlet
        m.fs.unit.vapor_inlet.flow_mol.fix(19000) # mol/sec
        m.fs.unit.vapor_inlet.temperature.fix(313.15) # K
        m.fs.unit.vapor_inlet.pressure.fix(100000) # Pa
        m.fs.unit.vapor_inlet.mole_frac_comp[0,"CO2"].fix(0.041)
        m.fs.unit.vapor_inlet.mole_frac_comp[0,"H2O"].fix(0.074)
        m.fs.unit.vapor_inlet.mole_frac_comp[0,"N2"].fix(0.763)
        m.fs.unit.vapor_inlet.mole_frac_comp[0,"O2"].fix(0.122)
        
        # Lean solvent inlet
        m.fs.unit.liquid_inlet.flow_mol.fix(24000) # mol/sec
        m.fs.unit.liquid_inlet.temperature.fix(313.15) # K
        m.fs.unit.liquid_inlet.pressure.fix(100000) # Pa
        m.fs.unit.liquid_inlet.mole_frac_comp[0,"CO2"].fix(0.022)
        m.fs.unit.liquid_inlet.mole_frac_comp[0,"H2O"].fix(0.868)
        m.fs.unit.liquid_inlet.mole_frac_comp[0,"MEA"].fix(0.11)

        return m

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model) == 0
        
    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model, outlvl=idaeslog.DEBUG)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        res = solver.solve(model)

        # Solver status/condition
        assert check_optimal_termination(res)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(19281.692273, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.flow_mol[0])
        assert pytest.approx(0.00135604, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"])
        assert pytest.approx(0.1265732, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "H2O"])
        assert pytest.approx(0.75185309, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "N2"])
        assert pytest.approx(0.12021766, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "O2"])
        assert pytest.approx(100000.0, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.pressure[0])
        assert pytest.approx(336.10714, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0])

        assert pytest.approx(23718.307726, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.flow_mol[0])
        assert pytest.approx(0.05400272, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "CO2"])
        assert pytest.approx(0.83469085, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "H2O"])
        assert pytest.approx(0.11130642, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "MEA"])
        assert pytest.approx(100000.0, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.pressure[0])
        assert pytest.approx(315.34352506, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.temperature[0])

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        vap_in = model.fs.unit.vapor_phase.properties[0, 0]
        vap_out = model.fs.unit.vapor_phase.properties[0, 1]
        liq_in = model.fs.unit.liquid_phase.properties[0, 1]
        liq_out = model.fs.unit.liquid_phase.properties[0, 0]

        # Material conservation
        for j in ["CO2", "H2O"]:
            assert 1e-6 >= abs(value(
                vap_in.get_material_flow_terms("Vap", j) +
                liq_in.get_material_flow_terms("Liq", j) -
                vap_out.get_material_flow_terms("Vap", j) -
                liq_out.get_material_flow_terms("Liq", j)))
        for j in ["N2", "O2"]:
            assert 1e-6 >= abs(value(
                vap_in.get_material_flow_terms("Vap", j) -
                vap_out.get_material_flow_terms("Vap", j)))
        for j in ["MEA"]:
            assert 1e-6 >= abs(value(
                liq_in.get_material_flow_terms("Liq", j) -
                liq_out.get_material_flow_terms("Liq", j)))

        # Energy conservation
        assert 1e-6 >= abs(value(
            vap_in.get_enthalpy_flow_terms("Vap") +
            liq_in.get_enthalpy_flow_terms("Liq") -
            vap_out.get_enthalpy_flow_terms("Vap") -
            liq_out.get_enthalpy_flow_terms("Liq")))
