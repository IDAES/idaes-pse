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
Author: Paul Akula, Anuja Deshpande
"""
# Import Python libraries
import pytest

# Import Pyomo libraries
from pyomo.environ import ConcreteModel, value, Param, TransformationFactory,\
    TerminationCondition, units as pyunits

# Import IDAES Libraries
from idaes.core import FlowsheetBlock
from idaes.power_generation.carbon_capture.mea_solvent_system.unit_models.column \
    import PackedColumn, ProcessType
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.vapor_prop \
    import VaporParameterBlock
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.liquid_prop \
    import LiquidParameterBlock

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver


# -----------------------------------------------------------------------------
solver = get_solver()


class TestColumn:
    """
    Tests for the column model.

    All inputs for state variables are in SI units:
    -Flowrate: mol/s
    -Temperature: K
    -Pressure: Pa
    """

    @pytest.mark.integration
    def test_steady_state_initialization(self):
        # Spacial domain finite elemets and finite element list
        x_nfe = 10
        x_nfe_list = [i / x_nfe for i in range(x_nfe + 1)]

        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up property package
        m.fs.vapor_properties = VaporParameterBlock()
        m.fs.liquid_properties = LiquidParameterBlock()

        # Create instance of column on flowsheet
        m.fs.unit = PackedColumn(default={
            "finite_elements": x_nfe,
            "length_domain_set": x_nfe_list,
            "transformation_method": "dae.finite_difference",
            "vapor_side": {
                "transformation_scheme": "BACKWARD",
                "property_package": m.fs.vapor_properties,
                "has_pressure_change": False,
                "pressure_drop_type": None},
            "liquid_side":
            {
                "transformation_scheme": "FORWARD",
                "property_package": m.fs.liquid_properties
            }})

        # Fix  input variables
        m.fs.unit.diameter_column.fix(0.64135)
        m.fs.unit.length_column.fix(18.15)
        for t in m.fs.time:
            # Vapor
            m.fs.unit.vapor_inlet.flow_mol[t].fix(21.48)
            m.fs.unit.vapor_inlet.temperature[t].fix(317.88)
            m.fs.unit.vapor_inlet.pressure[t].fix(107650)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "CO2"].fix(0.11453)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "H2O"].fix(0.08526)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "N2"].fix(0.73821)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "O2"].fix(0.06200)
            # Liquid
            m.fs.unit.liquid_inlet.flow_mol[t].fix(37.55)
            m.fs.unit.liquid_inlet.temperature[t].fix(319.87)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "CO2"].fix(0.00963)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "H2O"].fix(0.87435)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "MEA"].fix(0.11602)
        initialization_tester(m)

    @pytest.mark.integration
    def test_dynamic_initialization(self):
        # Spacial domain finite elemets and finite element list
        x_nfe = 10
        x_nfe_list = [i / x_nfe for i in range(x_nfe + 1)]

        # Time horizon
        t_nfe = 2
        time_set = [0, 4]

        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": True,
                                       'time_units': pyunits.s,
                                       "time_set": time_set})

        # Set up property package
        m.fs.vapor_properties = VaporParameterBlock()
        m.fs.liquid_properties = LiquidParameterBlock()

        # Create instance of column on flowsheet
        m.fs.unit = PackedColumn(default={
            "finite_elements": x_nfe,
            "length_domain_set": x_nfe_list,
            "transformation_method": "dae.finite_difference",
            "vapor_side": {
                "transformation_scheme": "BACKWARD",
                "property_package": m.fs.vapor_properties,
                "has_pressure_change": False,
                "pressure_drop_type": None},
            "liquid_side":
            {
                "transformation_scheme": "FORWARD",
                "property_package": m.fs.liquid_properties
            }})

        # Time discretization
        discretizer = TransformationFactory('dae.finite_difference')
        discretizer.apply_to(m.fs, wrt=m.fs.time, nfe=t_nfe, scheme='BACKWARD')

        # Fix  input variables
        m.fs.unit.diameter_column.fix(0.64135)
        m.fs.unit.length_column.fix(18.15)
        for t in m.fs.time:
            # Vapor
            m.fs.unit.vapor_inlet.flow_mol[t].fix(21.48)
            m.fs.unit.vapor_inlet.temperature[t].fix(317.88)
            m.fs.unit.vapor_inlet.pressure[t].fix(107650)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "CO2"].fix(0.11453)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "H2O"].fix(0.08526)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "N2"].fix(0.73821)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "O2"].fix(0.06200)
            # Liquid
            m.fs.unit.liquid_inlet.flow_mol[t].fix(37.55)
            m.fs.unit.liquid_inlet.temperature[t].fix(319.87)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "CO2"].fix(0.00963)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "H2O"].fix(0.87435)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "MEA"].fix(0.11602)
        initialization_tester(m)

    @pytest.fixture(scope="module",
                    params=[ProcessType.absorber, ProcessType.stripper])
    def column_model_ss(self, request):
        """ Setup for steady-state column"""

        Fpar = request.param

        # Spacial domain finite elemets and finite element list
        x_nfe = 10
        x_nfe_list = [i / x_nfe for i in range(x_nfe + 1)]

        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up property package
        m.fs.vapor_properties = VaporParameterBlock(
            default={'process_type': Fpar})
        m.fs.liquid_properties = LiquidParameterBlock(
            default={'process_type': Fpar})

        if Fpar == ProcessType.absorber:
            col_pressure = 107650
        elif Fpar == ProcessType.stripper:
            col_pressure = 183430

        # Create instance of column on flowsheet
        m.fs.unit = PackedColumn(default={
            "process_type": Fpar,
            "finite_elements": x_nfe,
            "length_domain_set": x_nfe_list,
            "column_pressure": col_pressure,
            "transformation_method": "dae.finite_difference",
            "vapor_side": {
                "transformation_scheme": "BACKWARD",
                "property_package": m.fs.vapor_properties,
                "has_pressure_change": False,
                "pressure_drop_type": None},
            "liquid_side":
            {
                "transformation_scheme": "FORWARD",
                "property_package": m.fs.liquid_properties
            }})

        # Fix  input variables
        m.fs.unit.diameter_column.fix(0.64135)
        if Fpar == ProcessType.absorber:
            m.fs.unit.length_column.fix(18.15)
            # Stream expected exit temperature
            m.fs.unit.vapor_exit_temp = Param(initialize=346.17684550190495)
            m.fs.unit.liquid_exit_temp = Param(initialize=322.95279074281416)
            for t in m.fs.time:
                # Vapor
                m.fs.unit.vapor_inlet.flow_mol[t].fix(21.48)
                m.fs.unit.vapor_inlet.temperature[t].fix(317.88)
                m.fs.unit.vapor_inlet.pressure[t].fix(col_pressure)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "CO2"].fix(0.11453)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "H2O"].fix(0.08526)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "N2"].fix(0.73821)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "O2"].fix(0.06200)
                # Liquid
                m.fs.unit.liquid_inlet.flow_mol[t].fix(37.55)
                m.fs.unit.liquid_inlet.temperature[t].fix(319.87)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "CO2"].fix(0.00963)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "H2O"].fix(0.87435)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "MEA"].fix(0.11602)

            # Initialize column
            m.fs.unit.initialize()
        elif Fpar == ProcessType.stripper:
            m.fs.unit.length_column.fix(12.1)
            # Stream expected  exit temperature
            m.fs.unit.vapor_exit_temp = Param(initialize=396.51574931626186)
            m.fs.unit.liquid_exit_temp = Param(initialize=393.86554487492464)
            for t in m.fs.time:
                # Vapor
                m.fs.unit.vapor_inlet.flow_mol[t].fix(17.496)
                m.fs.unit.vapor_inlet.temperature[t].fix(396.6)
                m.fs.unit.vapor_inlet.pressure[t].fix(col_pressure)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "CO2"].fix(0.0145)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "H2O"].fix(0.9855)
                # Liquid
                m.fs.unit.liquid_inlet.flow_mol[t].fix(84.48)
                m.fs.unit.liquid_inlet.temperature[t].fix(382.15)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "CO2"].fix(0.0331)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "H2O"].fix(0.8547)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "MEA"].fix(0.1122)

            # initialize column
            m.fs.unit.initialize(
                homotopy_steps_h=[0.1, 0.2, 0.4, 0.6, 0.8, 1])
        return m

    @pytest.fixture(scope="module",
                    params=[ProcessType.absorber, ProcessType.stripper])
    def column_model_dyn(self, request):
        """ Setup for dynamic column"""

        Fpar = request.param
        # Spacial domain finite elemets and finite element list
        x_nfe = 10
        x_nfe_list = [i / x_nfe for i in range(x_nfe + 1)]

        # Time horizon
        t_nfe = 2
        time_set = [0, 4]

        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": True,
                                       'time_units': pyunits.s,
                                       "time_set": time_set})
        # Set up property package

        m.fs.vapor_properties = VaporParameterBlock(
            default={'process_type': Fpar})
        m.fs.liquid_properties = LiquidParameterBlock(
            default={'process_type': Fpar})

        if Fpar == ProcessType.absorber:
            col_pressure = 107650
        elif Fpar == ProcessType.stripper:
            col_pressure = 183430

        # Create instance of column  on flowsheet
        m.fs.unit = PackedColumn(default={
            "process_type": Fpar,
            "finite_elements": x_nfe,
            "length_domain_set": x_nfe_list,
            "column_pressure": col_pressure,
            "transformation_method": "dae.finite_difference",
            "vapor_side": {
                "transformation_scheme": "BACKWARD",
                "property_package": m.fs.vapor_properties,
                "has_pressure_change": False,
                "pressure_drop_type": None},
            "liquid_side":
            {
                "transformation_scheme": "FORWARD",
                "property_package": m.fs.liquid_properties
            }})

        # Time discretization
        discretizer = TransformationFactory('dae.finite_difference')
        discretizer.apply_to(m.fs, wrt=m.fs.time, nfe=t_nfe, scheme='BACKWARD')

        # Fix inputs variables
        m.fs.unit.diameter_column.fix(0.64135)
        if Fpar == ProcessType.absorber:
            m.fs.unit.length_column.fix(18.15)
            # Rich loading at the bottom of column @ final time
            m.fs.unit.loading = Param(initialize=0.4927155969073804)
            for t in m.fs.time:
                # Vapor
                m.fs.unit.vapor_inlet.flow_mol[t].fix(21.48)
                m.fs.unit.vapor_inlet.temperature[t].fix(317.88)
                m.fs.unit.vapor_inlet.pressure[t].fix(col_pressure)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "CO2"].fix(0.11453)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "H2O"].fix(0.08526)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "N2"].fix(0.73821)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "O2"].fix(0.06200)
                # Liquid
                m.fs.unit.liquid_inlet.flow_mol[t].fix(37.55)
                m.fs.unit.liquid_inlet.temperature[t].fix(319.87)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "CO2"].fix(0.00963)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "H2O"].fix(0.87435)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "MEA"].fix(0.11602)

            # Initialize column
            m.fs.unit.initialize()
        elif Fpar == ProcessType.stripper:
            m.fs.unit.length_column.fix(12.1)
            # Lean loading at the bottom of column @ final time
            m.fs.unit.loading = Param(initialize=0.17982818165156983)
            for t in m.fs.time:
                # Vapor
                m.fs.unit.vapor_inlet.flow_mol[t].fix(17.496)
                m.fs.unit.vapor_inlet.temperature[t].fix(396.6)
                m.fs.unit.vapor_inlet.pressure[t].fix(col_pressure)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "CO2"].fix(0.0145)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "H2O"].fix(0.9855)
                # Liquid
                m.fs.unit.liquid_inlet.flow_mol[t].fix(84.48)
                m.fs.unit.liquid_inlet.temperature[t].fix(382.15)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "CO2"].fix(0.0331)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "H2O"].fix(0.8547)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "MEA"].fix(0.1122)

            # Initialize column
            m.fs.unit.initialize(
                homotopy_steps_h=[
                    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
        return m

    # -------------------------------------------------------------------------
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_steady_state_column_build(self, column_model_ss):

        assert column_model_ss.fs.unit.config.dynamic is False
        assert column_model_ss.fs.unit.config.liquid_side.transformation_scheme ==\
            'FORWARD'
        assert column_model_ss.fs.unit.config.vapor_side.transformation_scheme == \
            'BACKWARD'
        assert degrees_of_freedom(column_model_ss) == 0

        assert hasattr(column_model_ss.fs.unit, "vapor_inlet")
        assert hasattr(column_model_ss.fs.unit, "vapor_outlet")
        assert hasattr(column_model_ss.fs.unit, "liquid_inlet")
        assert hasattr(column_model_ss.fs.unit, "liquid_outlet")

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_dynamic_column_build(self, column_model_dyn):

        assert column_model_dyn.fs.unit.config.dynamic is True
        assert column_model_dyn.fs.unit.config.liquid_side.transformation_scheme ==\
            'FORWARD'
        assert column_model_dyn.fs.unit.config.vapor_side.transformation_scheme == \
            'BACKWARD'
        assert degrees_of_freedom(column_model_dyn) == 0
        assert hasattr(column_model_dyn.fs.unit, "vapor_inlet")
        assert hasattr(column_model_dyn.fs.unit, "vapor_outlet")
        assert hasattr(column_model_dyn.fs.unit, "liquid_inlet")
        assert hasattr(column_model_dyn.fs.unit, "liquid_outlet")

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_steady_state_column(self, column_model_ss):

        m = column_model_ss
        res = solver.solve(m)

        # Solver status/condition
        assert res.solver.termination_condition == TerminationCondition.optimal

        # Outlet Stream Condition Testing
        assert m.fs.unit.vapor_phase.properties[0, 1].temperature.value ==\
            pytest.approx(value(m.fs.unit.vapor_exit_temp), abs=1e-1)

        assert m.fs.unit.liquid_phase.properties[0, 0].temperature.value ==\
            pytest.approx(value(m.fs.unit.liquid_exit_temp), abs=1e-1)

        # material conservation check
        assert value(m.fs.unit.vapor_inlet.flow_mol[0] *
                     m.fs.unit.vapor_phase.properties[0, 0].mw -
                     m.fs.unit.vapor_outlet.flow_mol[0] *
                     m.fs.unit.vapor_phase.properties[0, 1].mw +
                     m.fs.unit.liquid_inlet.flow_mol[0] *
                     m.fs.unit.liquid_phase.properties[0, 1].mw -
                     m.fs.unit.liquid_outlet.flow_mol[0] *
                     m.fs.unit.liquid_phase.properties[0, 0].mw) ==\
            pytest.approx(0, abs=1e-2)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_dynamic_column(self, column_model_dyn):

        m = column_model_dyn
        res = solver.solve(m)

        # Solver status/condition
        assert res.solver.termination_condition == TerminationCondition.optimal

        # Performance Condition Testing at final time
        assert value(m.fs.unit.liquid_phase.properties[
                     m.fs.time.last(), 0].mole_frac_comp['CO2'] /
                     m.fs.unit.liquid_phase.properties[
                     m.fs.time.last(), 0].mole_frac_comp['MEA']) == \
            pytest.approx(value(m.fs.unit.loading), abs=1e-1)
