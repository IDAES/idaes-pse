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
Author: Paul Akula, Anuja Deshpande, Andrew Lee
"""
# Import Python libraries
import pytest

# Import Pyomo libraries
from pyomo.environ import ConcreteModel, value, Param, TransformationFactory, \
    check_optimal_termination, units as pyunits

# Import IDAES Libraries
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models.column_models.solvent_column \
    import PackedColumn
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.MEA_vapor \
    import flue_gas as vaporconfig_absorber
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.MEA_vapor \
    import wet_co2 as vaporconfig_stripper
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.MEA_solvent \
    import configuration as liquidconfig

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver


# -----------------------------------------------------------------------------
solver = get_solver()


class TestColumn(object):
    """
    Tests for the column model.

    All inputs for state variables are in SI units:
    -Flowrate: mol/s
    -Temperature: K
    -Pressure: Pa
    """
    @pytest.fixture(scope="class")
    def model_absorber_steady_state(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up property package
        m.fs.vapor_properties = GenericParameterBlock(default=vaporconfig_absorber)
        m.fs.liquid_properties = GenericParameterBlock(default=liquidconfig)
        
        # Number of finite elements and finite element list in the spatial domain
        x_nfe = 10
        x_nfe_list = [i / x_nfe for i in range(x_nfe + 1)]
        
        # Create an instance of the column in the flowsheet
        m.fs.unit = PackedColumn(default={
            "finite_elements": x_nfe,
            "length_domain_set": x_nfe_list,
            "transformation_method": "dae.finite_difference",
            "vapor_side": {
                "transformation_scheme": "BACKWARD",
                "property_package": m.fs.vapor_properties,
                "has_pressure_change": False},
            "liquid_side":
            {
                "transformation_scheme": "FORWARD",
                "property_package": m.fs.liquid_properties
            }})
            
        # Fix column design variables
        m.fs.unit.diameter_column.fix(0.64135)
        m.fs.unit.length_column.fix(18.15)
        
        # Fix operating conditions 
        for t in m.fs.time:
            # Flue gas
            m.fs.unit.vapor_inlet.flow_mol[t].fix(21.48)
            m.fs.unit.vapor_inlet.temperature[t].fix(317.88)
            m.fs.unit.vapor_inlet.pressure[t].fix(107650)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "CO2"].fix(0.11453)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "H2O"].fix(0.08526)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "N2"].fix(0.73821)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "O2"].fix(0.06200)
            # Solvent liquid
            m.fs.unit.liquid_inlet.flow_mol[t].fix(37.55)
            m.fs.unit.liquid_inlet.temperature[t].fix(319.87)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "CO2"].fix(0.00963)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "H2O"].fix(0.87435)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "MEA"].fix(0.11602)
            
        # Fix vapor phase mass transfer coefficient values
        k_v_values = [[0, 0],[2.837e-05, 3.728e-05],[ 2.862e-05, 3.757e-05],
                      [ 2.891e-05, 3.788e-05],[ 2.924e-05, 3.825e-05],
                      [ 2.965e-05, 3.87e-05],[ 3.018e-05, 3.929e-05],
                      [ 3.092e-05, 4.011e-05],[ 3.195e-05, 4.126e-05],
                      [ 3.305e-05, 4.251e-05],[ 3.18e-05, 4.121e-05]]
        
        for t in m.fs.time:
            for i,x in enumerate(m.fs.unit.vapor_phase.length_domain):
                for j,comp in enumerate(['CO2','H2O']):
                    if x == m.fs.unit.vapor_phase.length_domain.first():
                        m.fs.unit.k_v[t, x, comp].fix(0.001)
                    else:
                        m.fs.unit.k_v[t, x, comp].fix(k_v_values[i][j])
        
        # Fix liquid phase mass transfer coefficient values        
        k_l_values = [9.613e-05, 9.861e-05, 0.0001012, 0.000104, 
                      0.0001072, 0.0001111, 0.0001159, 0.0001222, 
                      0.0001294, 0.0001311, 0.001]
        
        for t in m.fs.time:
            for i,x in enumerate(m.fs.unit.liquid_phase.length_domain):
                for j,comp in enumerate(['CO2']):
                    if x == m.fs.unit.liquid_phase.length_domain.last():
                        m.fs.unit.k_l[t, x, comp].fix(0.001)
                    else:
                        m.fs.unit.k_l[t, x, comp].fix(k_l_values[i])
        
        # Fix vapor phase heat transfer coefficient values
        h_v_values = [100, 102.3, 103.1, 103.9, 104.9, 106.1, 107.6, 109.7, 
                      112.6, 115.5, 111.8]
        for t in m.fs.time:
            for i,x in enumerate(m.fs.unit.vapor_phase.length_domain):
                if x == m.fs.unit.vapor_phase.length_domain.first():
                    m.fs.unit.h_v[t, x].fix(100)
                else:
                    m.fs.unit.h_v[t, x].fix(h_v_values[i])
        
        # Fix interfacial area values
        interfacial_area_values = [0, 198.2, 198.5, 198.8, 199.2, 199.6, 
                                   200.2, 201, 202.2, 203.3, 201.3]
        
        for t in m.fs.time:
            for i,x in enumerate(m.fs.unit.vapor_phase.length_domain):
                if x == m.fs.unit.vapor_phase.length_domain.first():
                    m.fs.unit.area_interfacial[t, x].fix(0)
                else:
                    m.fs.unit.area_interfacial[t, x].fix(interfacial_area_values[i])
        
        # Fix enhancement factor values       
        enhancement_factor_values = [11.81960366, 13.21436568, 14.8235168,
                                     16.80737692, 19.43845149, 23.23126553,
                                     29.47937877, 41.78076923, 74.63068006,
                                     188.3501144, 10]
        
        for t in m.fs.time:
            for i,x in enumerate(m.fs.unit.liquid_phase.length_domain):
                if x == m.fs.unit.liquid_phase.length_domain.last():
                   m.fs.unit.enhancement_factor[t, x].fix(10)
                else:
                    m.fs.unit.enhancement_factor[t, x].fix(enhancement_factor_values[i])
                    
        return(m)
    
    @pytest.mark.build
    @pytest.mark.unit
    def test_steady_state_absorber_build(self, model_absorber_steady_state):

        assert model_absorber_steady_state.fs.unit.config.dynamic is False
        assert model_absorber_steady_state.fs.unit.config.liquid_side.transformation_scheme ==\
            'FORWARD'
        assert model_absorber_steady_state.fs.unit.config.vapor_side.transformation_scheme == \
            'BACKWARD'

        assert hasattr(model_absorber_steady_state.fs.unit, "vapor_inlet")
        assert hasattr(model_absorber_steady_state.fs.unit, "vapor_outlet")
        assert hasattr(model_absorber_steady_state.fs.unit, "liquid_inlet")
        assert hasattr(model_absorber_steady_state.fs.unit, "liquid_outlet")
    
    @pytest.mark.unit
    def test_dof_absorber(self, model_absorber_steady_state):
        assert degrees_of_freedom(model_absorber_steady_state) == 0
        
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_solve_absorber(self, model_absorber_steady_state):
        initialization_tester(model_absorber_steady_state)
    
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve_absorber(self, model_absorber_steady_state):
        results=solver.solve(model_absorber_steady_state)
        
        # Solver status and condition
        assert check_optimal_termination(results)
        
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_absorber_conservation(self, model_absorber_steady_state):
        vap_comp = model_absorber_steady_state.fs.unit.config.vapor_side.property_package.component_list
        liq_comp = model_absorber_steady_state.fs.unit.config.liquid_side.property_package.component_list
        equilibrium_comp = vap_comp & liq_comp
        solvent_comp_list = \
            model_absorber_steady_state.fs.unit.config.liquid_side.property_package.solvent_set
        solute_comp_list = model_absorber_steady_state.fs.unit.config.liquid_side.property_package.solute_set
        
        # Mass conservation test
        
        vap_in = model_absorber_steady_state.fs.unit.vapor_phase.properties[0, 0]
        vap_out = model_absorber_steady_state.fs.unit.vapor_phase.properties[0, 1]
        liq_in = model_absorber_steady_state.fs.unit.liquid_phase.properties[0, 1]
        liq_out = model_absorber_steady_state.fs.unit.liquid_phase.properties[0, 0]
        
        # Material conservation
        for j in liq_comp:
            if j in equilibrium_comp:
                assert 1e-6 >= abs(value(
                    vap_in.get_material_flow_terms("Vap", j) +
                    liq_in.get_material_flow_terms("Liq", j) -
                    vap_out.get_material_flow_terms("Vap", j) -
                    liq_out.get_material_flow_terms("Liq", j)))
            elif j in solvent_comp_list:
                assert 1e-6 >= abs(value(
                liq_in.get_material_flow_terms("Liq", j) -
                liq_out.get_material_flow_terms("Liq", j)))
    
        for j in vap_comp:
            if j not in equilibrium_comp:
                assert 1e-6 >= abs(value(
                    vap_in.get_material_flow_terms("Vap", j) -
                    vap_out.get_material_flow_terms("Vap", j)))
                
        # Energy conservation
        assert 1e-6 >= abs(value(
            vap_in.get_enthalpy_flow_terms("Vap") +
            liq_in.get_enthalpy_flow_terms("Liq") -
            vap_out.get_enthalpy_flow_terms("Vap") -
            liq_out.get_enthalpy_flow_terms("Liq")))
        
        
    @pytest.fixture(scope="class")
    def model_stripper_steady_state(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up property package
        m.fs.vapor_properties = GenericParameterBlock(default=vaporconfig_stripper)
        m.fs.liquid_properties = GenericParameterBlock(default=liquidconfig)
        
        # Number of finite elements and finite element list in the spatial domain
        x_nfe = 10
        x_nfe_list = [i / x_nfe for i in range(x_nfe + 1)]
        
        # Create an instance of the column in the flowsheet
        m.fs.unit = PackedColumn(default={
            "finite_elements": x_nfe,
            "length_domain_set": x_nfe_list,
            "transformation_method": "dae.finite_difference",
            "vapor_side": {
                "transformation_scheme": "BACKWARD",
                "property_package": m.fs.vapor_properties,
                "has_pressure_change": False},
            "liquid_side":
            {
                "transformation_scheme": "FORWARD",
                "property_package": m.fs.liquid_properties
            }})
            
        # Fix column design variables
        m.fs.unit.diameter_column.fix(0.64135)
        m.fs.unit.length_column.fix(12.1)
        
        # Fix operating conditions         
        for t in m.fs.time:
            # Flue gas
            m.fs.unit.vapor_inlet.flow_mol[t].fix(17.496)
            m.fs.unit.vapor_inlet.temperature[t].fix(396.6)
            m.fs.unit.vapor_inlet.pressure[t].fix(183430)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "CO2"].fix(0.0145)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "H2O"].fix(0.9855)
            # Solvent liquid
            m.fs.unit.liquid_inlet.flow_mol[t].fix(84.48)
            m.fs.unit.liquid_inlet.temperature[t].fix(382.15)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "CO2"].fix(0.0331)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "H2O"].fix(0.8547)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "MEA"].fix(0.1122)
            
        # Fix vapor phase mass transfer coefficient values
        k_v_values = [[0, 0],[2.837e-05, 3.728e-05],[ 2.862e-05, 3.757e-05],
                      [ 2.891e-05, 3.788e-05],[ 2.924e-05, 3.825e-05],
                      [ 2.965e-05, 3.87e-05],[ 3.018e-05, 3.929e-05],
                      [ 3.092e-05, 4.011e-05],[ 3.195e-05, 4.126e-05],
                      [ 3.305e-05, 4.251e-05],[ 3.18e-05, 4.121e-05]]
        
        for t in m.fs.time:
            for i,x in enumerate(m.fs.unit.vapor_phase.length_domain):
                for j,comp in enumerate(['CO2','H2O']):
                    if x == m.fs.unit.vapor_phase.length_domain.first():
                        m.fs.unit.k_v[t, x, comp].fix(0.001)
                    else:
                        m.fs.unit.k_v[t, x, comp].fix(k_v_values[i][j])
        
        # Fix liquid phase mass transfer coefficient values        
        k_l_values = [9.613e-05, 9.861e-05, 0.0001012, 0.000104, 
                      0.0001072, 0.0001111, 0.0001159, 0.0001222, 
                      0.0001294, 0.0001311, 0.001]
        
        for t in m.fs.time:
            for i,x in enumerate(m.fs.unit.liquid_phase.length_domain):
                for j,comp in enumerate(['CO2']):
                    if x == m.fs.unit.liquid_phase.length_domain.last():
                        m.fs.unit.k_l[t, x, comp].fix(0.001)
                    else:
                        m.fs.unit.k_l[t, x, comp].fix(k_l_values[i])
        
        # Fix vapor phase heat transfer coefficient values
        h_v_values = [100, 102.3, 103.1, 103.9, 104.9, 106.1, 107.6, 109.7, 
                      112.6, 115.5, 111.8]
        for t in m.fs.time:
            for i,x in enumerate(m.fs.unit.vapor_phase.length_domain):
                if x == m.fs.unit.vapor_phase.length_domain.first():
                    m.fs.unit.h_v[t, x].fix(100)
                else:
                    m.fs.unit.h_v[t, x].fix(h_v_values[i])
        
        # Fix interfacial area values
        interfacial_area_values = [0, 198.2, 198.5, 198.8, 199.2, 199.6, 
                                   200.2, 201, 202.2, 203.3, 201.3]
        
        for t in m.fs.time:
            for i,x in enumerate(m.fs.unit.vapor_phase.length_domain):
                if x == m.fs.unit.vapor_phase.length_domain.first():
                    m.fs.unit.area_interfacial[t, x].fix(0)
                else:
                    m.fs.unit.area_interfacial[t, x].fix(interfacial_area_values[i])
        
        # Fix enhancement factor values       
        enhancement_factor_values = [11.81960366, 13.21436568, 14.8235168,
                                     16.80737692, 19.43845149, 23.23126553,
                                     29.47937877, 41.78076923, 74.63068006,
                                     188.3501144, 10]
        
        for t in m.fs.time:
            for i,x in enumerate(m.fs.unit.liquid_phase.length_domain):
                if x == m.fs.unit.liquid_phase.length_domain.last():
                   m.fs.unit.enhancement_factor[t, x].fix(10)
                else:
                    m.fs.unit.enhancement_factor[t, x].fix(enhancement_factor_values[i])
                    
        return(m)
    
    @pytest.mark.build
    @pytest.mark.unit
    def test_steady_state_stripper_build(self, model_stripper_steady_state):

        assert model_stripper_steady_state.fs.unit.config.dynamic is False
        assert model_stripper_steady_state.fs.unit.config.liquid_side.transformation_scheme ==\
            'FORWARD'
        assert model_stripper_steady_state.fs.unit.config.vapor_side.transformation_scheme == \
            'BACKWARD'

        assert hasattr(model_stripper_steady_state.fs.unit, "vapor_inlet")
        assert hasattr(model_stripper_steady_state.fs.unit, "vapor_outlet")
        assert hasattr(model_stripper_steady_state.fs.unit, "liquid_inlet")
        assert hasattr(model_stripper_steady_state.fs.unit, "liquid_outlet")
    
    @pytest.mark.unit
    def test_dof_stripper(self, model_stripper_steady_state):
        assert degrees_of_freedom(model_stripper_steady_state) == 0
        
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_solve_stripper(self, model_stripper_steady_state):
        initialization_tester(model_stripper_steady_state)
     
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve_stripper(self, model_stripper_steady_state):
        results=solver.solve(model_stripper_steady_state)
        
        # Solver status and condition
        assert check_optimal_termination(results)
        
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_stripper_conservation(self, model_stripper_steady_state):
        
        vap_comp = model_stripper_steady_state.fs.unit.config.vapor_side.property_package.component_list
        liq_comp = model_stripper_steady_state.fs.unit.config.liquid_side.property_package.component_list
        equilibrium_comp = vap_comp & liq_comp
        solvent_comp_list = \
            model_stripper_steady_state.fs.unit.config.liquid_side.property_package.solvent_set
        solute_comp_list = model_stripper_steady_state.fs.unit.config.liquid_side.property_package.solute_set
        
        # Mass conservation test
        
        vap_in = model_stripper_steady_state.fs.unit.vapor_phase.properties[0, 0]
        vap_out = model_stripper_steady_state.fs.unit.vapor_phase.properties[0, 1]
        liq_in = model_stripper_steady_state.fs.unit.liquid_phase.properties[0, 1]
        liq_out = model_stripper_steady_state.fs.unit.liquid_phase.properties[0, 0]
        
        # Material conservation
        for j in liq_comp:
            if j in equilibrium_comp:
                assert 1e-6 >= abs(value(
                    vap_in.get_material_flow_terms("Vap", j) +
                    liq_in.get_material_flow_terms("Liq", j) -
                    vap_out.get_material_flow_terms("Vap", j) -
                    liq_out.get_material_flow_terms("Liq", j)))
            elif j in solvent_comp_list:
                assert 1e-6 >= abs(value(
                liq_in.get_material_flow_terms("Liq", j) -
                liq_out.get_material_flow_terms("Liq", j)))
    
        for j in vap_comp:
            if j not in equilibrium_comp:
                assert 1e-6 >= abs(value(
                    vap_in.get_material_flow_terms("Vap", j) -
                    vap_out.get_material_flow_terms("Vap", j)))
                
        # Energy conservation
        assert 1e-6 >= abs(value(
            vap_in.get_enthalpy_flow_terms("Vap") +
            liq_in.get_enthalpy_flow_terms("Liq") -
            vap_out.get_enthalpy_flow_terms("Vap") -
            liq_out.get_enthalpy_flow_terms("Liq")))