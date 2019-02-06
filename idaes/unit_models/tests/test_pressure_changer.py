##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Tests for Pressure Changer unit model.

Author: Andrew Lee, Emmanuel Ogbe
"""
import pytest
from pyomo.environ import ConcreteModel, SolverFactory
from idaes.core import FlowsheetBlock, declare_process_block_class
from idaes.unit_models.pressure_changer import PressureChanger
from idaes.ui.report import degrees_of_freedom

# Import property package for testing
from idaes.property_models import iapws95_ph as pp


# -----------------------------------------------------------------------------
# General test classes
@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlock):
    def build(self):
        super(_Flowsheet, self).build()


if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6, 'bound_push': 1e-8}
else:
    solver = None


def test_build_pc():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.props = pp.Iapws95ParameterBlock()
    m.fs.pc = PressureChanger(default={"property_package": m.fs.props})

    assert hasattr(m.fs.pc, "inlet")
    assert hasattr(m.fs.pc, "outlet")
    assert len(m.fs.pc.inlet.vars) == 3
    assert len(m.fs.pc.outlet.vars) == 3


def test_set_geometry_include_holdup_true():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.props = pp.Iapws95ParameterBlock()
    m.fs.pc = PressureChanger(default={"property_package": m.fs.props,
                                       "has_holdup": True})

    assert hasattr(m.fs.pc, "volume")
    assert hasattr(m.fs.pc.control_volume, "material_holdup")


def test_make_performance():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.props = pp.Iapws95ParameterBlock()
    m.fs.pc = PressureChanger(default={"property_package": m.fs.props})

    assert hasattr(m.fs.pc, "work_mechanical")
    assert m.fs.pc.work_mechanical == m.fs.pc.control_volume.work
    assert hasattr(m.fs.pc, "deltaP")
    assert hasattr(m.fs.pc, "ratioP")
    assert hasattr(m.fs.pc, "ratioP_calculation")


def test_make_isothermal():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.props = pp.Iapws95ParameterBlock()

    m.fs.pc = PressureChanger(default={
            "property_package": m.fs.props,
            "thermodynamic_assumption": 'isothermal'})

    assert hasattr(m.fs.pc, "isothermal")


def test_make_pump():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.props = pp.Iapws95ParameterBlock()

    m.fs.pc = PressureChanger(default={
            "property_package": m.fs.props,
            "thermodynamic_assumption": 'pump'})

    assert hasattr(m.fs.pc, "work_fluid")
    assert hasattr(m.fs.pc, "efficiency_pump")
    assert hasattr(m.fs.pc, "fluid_work_calculation")
    assert hasattr(m.fs.pc, "actual_work")


def test_make_isentropic():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.props = pp.Iapws95ParameterBlock()

    m.fs.pc = PressureChanger(default={
            "property_package": m.fs.props,
            "thermodynamic_assumption": 'isentropic'})

    assert hasattr(m.fs.pc, "efficiency_isentropic")
    assert hasattr(m.fs.pc, "work_isentropic")
    assert hasattr(m.fs.pc, "isentropic_pressure")
    assert hasattr(m.fs.pc, "isentropic_material")
    assert hasattr(m.fs.pc, "isentropic")
    assert hasattr(m.fs.pc, "isentropic_energy_balance")
    assert hasattr(m.fs.pc, "actual_work")


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialization_isothermal():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.props = pp.Iapws95ParameterBlock()

    m.fs.pc = PressureChanger(default={
            "property_package": m.fs.props,
            "thermodynamic_assumption": 'isothermal'})

    m.fs.pc.deltaP.fix(-1e3)
    m.fs.pc.inlet.flow_mol.fix(27.5e3)
    m.fs.pc.inlet.enth_mol.fix(4000)
    m.fs.pc.inlet.pressure.fix(2e6)

    assert degrees_of_freedom(m) == 0

    init_state = {
        "flow_mol": 27.5e3,
        "pressure": 2e6,
        "enth_mol": 4000
    }

    m.fs.pc.initialize(state_args=init_state, outlvl=5)

    assert (pytest.approx(27500.0, abs=1e-2) ==
            m.fs.pc.outlet.flow_mol[0].value)
    assert (pytest.approx(3999.984582673592, abs=1e-2) ==
            m.fs.pc.outlet.enth_mol[0].value)
    assert (pytest.approx(1999000.0, abs=1e-2) ==
            m.fs.pc.outlet.pressure[0].value)

    solver.solve(m)


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialization_pump():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.props = pp.Iapws95ParameterBlock()

    m.fs.pc = PressureChanger(default={
            "property_package": m.fs.props,
            "thermodynamic_assumption": 'pump'})

    m.fs.pc.inlet.flow_mol.fix(27.5e3)
    m.fs.pc.inlet.enth_mol.fix(4000)
    m.fs.pc.inlet.pressure.fix(2e6)
    m.fs.pc.deltaP.fix(-1e3)
    m.fs.pc.efficiency_pump.fix(0.9)

    assert degrees_of_freedom(m) == 0

    init_state = {
        "flow_mol": 27.5e3,
        "pressure": 2e6,
        "enth_mol": 4000
    }

    m.fs.pc.initialize(state_args=init_state, outlvl=5,
                       optarg={'tol': 1e-6})

    assert (pytest.approx(27.5e3, abs=1e-2) ==
            m.fs.pc.outlet.flow_mol[0].value)
    assert (pytest.approx(3999.979732728688, abs=1e-2) ==
            m.fs.pc.outlet.enth_mol[0].value)
    assert (pytest.approx(1999000.0, abs=1e-2) ==
            m.fs.pc.outlet.pressure[0].value)

    solver.solve(m)


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialization_adiabatic():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.props = pp.Iapws95ParameterBlock()

    m.fs.pc = PressureChanger(default={
            "property_package": m.fs.props,
            "thermodynamic_assumption": 'adiabatic'})

    init_state = {
        "flow_mol": 27.5e3,
        "pressure": 2e6,
        "enth_mol": 4000
    }

    m.fs.pc.inlet.flow_mol.fix(27.5e3)
    m.fs.pc.inlet.enth_mol.fix(4000)
    m.fs.pc.inlet.pressure.fix(2e6)
    m.fs.pc.deltaP.fix(-1e3)

    assert degrees_of_freedom(m) == 0

    m.fs.pc.initialize(state_args=init_state, outlvl=5,
                       optarg={'tol': 1e-6})

    assert (pytest.approx(27.5e3, abs=1e-2) ==
            m.fs.pc.outlet.flow_mol[0].value)
    assert (pytest.approx(4000, abs=1e-2) ==
            m.fs.pc.outlet.enth_mol[0].value)
    assert (pytest.approx(1999000.0, abs=1e-2) ==
            m.fs.pc.outlet.pressure[0].value)

    solver.solve(m)


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialization_isentropic():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.props = pp.Iapws95ParameterBlock()

    m.fs.pc = PressureChanger(default={
            "property_package": m.fs.props,
            "thermodynamic_assumption": 'isentropic'})

    init_state = {
        "flow_mol": 27.5e3,
        "pressure": 2e6,
        "enth_mol": 4000
    }

    m.fs.pc.inlet.flow_mol.fix(27.5e3)
    m.fs.pc.inlet.enth_mol.fix(4000)
    m.fs.pc.inlet.pressure.fix(2e6)
    m.fs.pc.deltaP.fix(-1e3)
    m.fs.pc.efficiency_isentropic.fix(0.83)

    assert degrees_of_freedom(m) == 0

    m.fs.pc.initialize(state_args=init_state, outlvl=5,
                       optarg={'tol': 1e-6})

    assert (pytest.approx(27.5e3, abs=1e-2) ==
            m.fs.pc.outlet.flow_mol[0].value)
    assert (pytest.approx(3999.979732728688, abs=1e-2) ==
            m.fs.pc.outlet.enth_mol[0].value)
    assert (pytest.approx(1999000.0, abs=1e-2) ==
            m.fs.pc.outlet.pressure[0].value)

    solver.solve(m)
