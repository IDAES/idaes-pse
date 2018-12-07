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
from idaes.core import FlowsheetBlockData, declare_process_block_class, \
                        ControlVolume0D
from idaes.unit_models.pressure_changer import PressureChanger, PressureChangerData
from idaes.ui.report import degrees_of_freedom

# Import property package for testing
from idaes.property_models import iapws95_ph as fe

# -----------------------------------------------------------------------------
# General test classes
@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(_Flowsheet, self).build()


@declare_process_block_class("PressureChangerFrame")
class _PressureChanger(PressureChangerData):
    def build(self):
        pass


@pytest.fixture()
def m():
    """
    Build a simple flowsheet model for build tests from.
    """
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.props = fe.Iapws95ParameterBlock()
    m.fs.pc = PressureChangerFrame(property_package=m.fs.props)
    return m


# -----------------------------------------------------------------------------
# Test class methods
def test_default_config_block(m):
    assert len(m.fs.pc.config) == 20
    assert m.fs.pc.config.dynamic == 'use_parent_value'
    assert not m.fs.pc.config.include_holdup
    assert m.fs.pc.config.material_balance_type == 'component_phase'
    assert m.fs.pc.config.energy_balance_type == 'total'
    assert m.fs.pc.config.momentum_balance_type == 'total'
    assert not m.fs.pc.config.has_rate_reactions
    assert not m.fs.pc.config.has_equilibrium_reactions
    assert m.fs.pc.config.has_phase_equilibrium
    assert not m.fs.pc.config.has_mass_transfer
    assert not m.fs.pc.config.has_heat_transfer
    assert m.fs.pc.config.has_work_transfer
    assert m.fs.pc.config.has_pressure_change
    assert m.fs.pc.config.property_package == m.fs.props
    assert m.fs.pc.config.property_package_args == {}
    assert m.fs.pc.config.inlet_list is None
    assert m.fs.pc.config.num_inlets is None
    assert m.fs.pc.config.outlet_list is None
    assert m.fs.pc.config.num_outlets is None

    assert m.fs.pc.config.compressor
    m.fs.pc.config.compressor = True
    m.fs.pc.config.compressor = False
    with pytest.raises(ValueError):
        m.fs.pc.config.compressor = 'foo'
    with pytest.raises(ValueError):
        m.fs.pc.config.compressor = 10

    assert m.fs.pc.config.thermodynamic_assumption == 'isothermal'
    m.fs.pc.config.thermodynamic_assumption = 'isothermal'
    m.fs.pc.config.thermodynamic_assumption = 'isentropic'
    m.fs.pc.config.thermodynamic_assumption = 'pump'
    with pytest.raises(ValueError):
        m.fs.pc.config.thermodynamic_assumption = 'foo'
    with pytest.raises(ValueError):
        m.fs.pc.config.thermodynamic_assumption = True
    with pytest.raises(ValueError):
        m.fs.pc.config.thermodynamic_assumption = 10


def test_super_build(m):
    super(PressureChangerData, m.fs.pc).build()
    assert not m.fs.pc.config.dynamic


def test_adding_holdup(m):
    super(PressureChangerData, m.fs.pc).build()
    m.fs.pc.control_volume = ControlVolume0D()

    assert not m.fs.pc.control_volume.config.include_holdup
    assert m.fs.pc.control_volume.config.material_balance_type == 'component_phase'
    assert m.fs.pc.control_volume.config.energy_balance_type == 'total'
    assert m.fs.pc.control_volume.config.momentum_balance_type == 'total'
    assert not m.fs.pc.control_volume.config.has_rate_reactions
    assert not m.fs.pc.control_volume.config.has_equilibrium_reactions
    # Properties do not have phase equilibrium, so this should be False
    assert not m.fs.pc.control_volume.config.has_phase_equilibrium
    assert not m.fs.pc.control_volume.config.has_mass_transfer
    assert not m.fs.pc.control_volume.config.has_heat_transfer
    assert m.fs.pc.control_volume.config.has_work_transfer
    assert m.fs.pc.control_volume.config.has_pressure_change
    assert m.fs.pc.control_volume.config.property_package == m.fs.props
    assert m.fs.pc.control_volume.config.property_package_args == {}


def test_set_geometry_include_holdup_false(m):
    super(PressureChangerData, m.fs.pc).build()
    m.fs.pc.control_volume = ControlVolume0D()
    m.fs.pc.set_geometry()

    assert not hasattr(m.fs.pc, "volume")


def test_set_geometry_include_holdup_true(m):
    m.fs.pc.config.include_holdup = True
    super(PressureChangerData, m.fs.pc).build()
    m.fs.pc.control_volume = ControlVolume0D()
    m.fs.pc.set_geometry()

    assert hasattr(m.fs.pc, "volume")
    assert hasattr(m.fs.pc.control_volume, "volume")
    assert m.fs.pc.volume == m.fs.pc.control_volume.volume


def test_make_performance(m):
    super(PressureChangerData, m.fs.pc).build()
    m.fs.pc.control_volume = ControlVolume0D()
    m.fs.pc.set_geometry()
    m.fs.pc.add_performance()

    assert hasattr(m.fs.pc, "work_mechanical")
    assert m.fs.pc.work_mechanical == m.fs.pc.control_volume.work
    assert hasattr(m.fs.pc, "deltaP")
    assert hasattr(m.fs.pc, "ratioP")
    assert hasattr(m.fs.pc, "ratioP_calculation")


def test_make_pump(m):
    super(PressureChangerData, m.fs.pc).build()
    m.fs.pc.control_volume = ControlVolume0D()
    m.fs.pc.set_geometry()
    m.fs.pc.make_performance()
    m.fs.pc.make_pump()

    assert hasattr(m.fs.pc, "work_fluid")
    assert hasattr(m.fs.pc, "efficiency_pump")
    assert hasattr(m.fs.pc, "fluid_work_calculation")
    assert hasattr(m.fs.pc, "actual_work")


def test_make_isothermal(m):
    super(PressureChangerData, m.fs.pc).build()
    m.fs.pc.control_volume = COntrolVolume0D()
    m.fs.pc.set_geometry()
    m.fs.pc.add_performance()
    m.fs.pc.add_isothermal()

    assert hasattr(m.fs.pc, "isothermal")


def test_make_isentropic(m):
    super(PressureChangerData, m.fs.pc).build()
    m.fs.pc.control_volume = ControlVolume0D()
    m.fs.pc.set_geometry()
    m.fs.pc.add_performance()
    m.fs.pc.add_isentropic()

    assert hasattr(m.fs.pc, "efficiency_isentropic")
    assert hasattr(m.fs.pc, "work_isentropic")
    assert hasattr(m.fs.pc, "properties_isentropic")
    assert hasattr(m.fs.pc, "isentropic_pressure")
    assert hasattr(m.fs.pc, "isentropic_material")
    assert hasattr(m.fs.pc, "isentropic")
    assert hasattr(m.fs.pc, "isentropic_energy_balance")
    assert hasattr(m.fs.pc, "actual_work")


# -----------------------------------------------------------------------------
# Test build method for each case - use compressor = False for all case
def test_build_isothermal():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.props = fe.PropertyParameterBlock()
    m.fs.pc = PressureChanger(property_package=m.fs.props,
                              compressor=False,
                              thermodynamic_assumption='isothermal')

    assert hasattr(m.fs.pc, "isothermal")


def test_build_pump():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.props = fe.PropertyParameterBlock()
    m.fs.pc = PressureChanger(property_package=m.fs.props,
                              compressor=False,
                              thermodynamic_assumption='pump')

    assert hasattr(m.fs.pc, "efficiency_pump")


def test_build_isentropic():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.props = fe.PropertyParameterBlock()
    m.fs.pc = PressureChanger(property_package=m.fs.props,
                              compressor=False,
                              thermodynamic_assumption='isentropic')

    assert hasattr(m.fs.pc, "efficiency_isentropic")


# -----------------------------------------------------------------------------
# Test post_transform_build and connections
@pytest.fixture()
def m2():
    """
    Build a simple flowsheet model for build tests from.
    """
    m2 = ConcreteModel()
    m2.fs = Flowsheet(dynamic=False)
    m2.fs.props = fe.PropertyParameterBlock()
    m2.fs.pc = PressureChanger(property_package=m2.fs.props)
    return m2


def test_post_transform_build_inlet_list(m2):
    m2.fs.pc.config.inlet_list = ["a", "b"]

    assert hasattr(m2.fs.pc, "inlet")
    assert hasattr(m2.fs.pc, "outlet")

    # Test that inlet and outlet have correct indexing sets
    for k in m2.fs.pc.inlet:
        assert k in [(0, "a"), (0, "b")]
    for k in m2.fs.pc.outlet:
        assert k in [0]


def test_post_transform_build_num_inlets(m2):
    m2.fs.pc.config.num_inlets = 2

    assert hasattr(m2.fs.pc, "inlet")
    assert hasattr(m2.fs.pc, "outlet")

    # Test that inlet and outlet have correct indexing sets
    for k in m2.fs.pc.inlet:
        assert k in [(0, "1"), (0, "2")]
    for k in m2.fs.pc.outlet:
        assert k in [0]


def test_post_transform_build_outlet_list(m2):
    m2.fs.pc.config.outlet_list = ["a", "b"]

    assert hasattr(m2.fs.pc, "inlet")
    assert hasattr(m2.fs.pc, "outlet")

    # Test that inlet and outlet have correct indexing sets
    for k in m2.fs.pc.inlet:
        assert k in [0]
    for k in m2.fs.pc.outlet:
        assert k in [(0, "a"), (0, "b")]


def test_post_transform_build_num_outlet2(m2):
    m2.fs.pc.config.num_outlets = 2

    assert hasattr(m2.fs.pc, "inlet")
    assert hasattr(m2.fs.pc, "outlet")

    # Test that inlet and outlet have correct indexing sets
    for k in m2.fs.pc.inlet:
        assert k in [0]
    for k in m2.fs.pc.outlet:
        assert k in [(0, "1"), (0, "2")]


# -----------------------------------------------------------------------------
# Test model_checks
def test_model_check_compressor_true():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.props = fe.PropertyParameterBlock()

    m.fs.pc = PressureChanger(property_package=m.fs.props,
                              compressor=True)
    m.fs.pc.model_check()


def test_model_check_compressor_false():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.props = fe.PropertyParameterBlock()

    m.fs.pc = PressureChanger(property_package=m.fs.props,
                              compressor=False)
    m.fs.pc.model_check()


def test_model_check_isentropic():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.props = fe.PropertyParameterBlock()

    m.fs.pc = PressureChanger(property_package=m.fs.props,
                              thermodynamic_assumption='isentropic')
    m.fs.pc.model_check()



# -----------------------------------------------------------------------------
# Test initialization
# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialization_isothermal():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.props = fe.PropertyParameterBlock()

    m.fs.pc = PressureChanger(property_package=m.fs.props,
                              thermodynamic_assumption='isothermal')

    m.fs.pc.inlet[:].vars["flow_mol_comp"]["H2O"].fix(27.5e3)
    m.fs.pc.inlet[:].vars["temperature"].fix(866.5)
    m.fs.pc.inlet[:].vars["pressure"].fix(2.97e7)

    m.fs.pc.deltaP.fix(-1e7)

    assert degrees_of_freedom(m) == 0

    m.fs.pc.initialize()

    assert (pytest.approx(27500.0, abs=1e-2) ==
            m.fs.pc.outlet[0].vars["flow_mol_comp"]["H2O"].value)
    assert (pytest.approx(866.5, abs=1e-2) ==
            m.fs.pc.outlet[0].vars["temperature"].value)
    assert (pytest.approx(19700000.0, abs=1e-2) ==
            m.fs.pc.outlet[0].vars["pressure"].value)
    assert (pytest.approx(57753341.62230536, abs=1e-2) ==
            m.fs.pc.work_mechanical[0].value)


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialization_pump():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.props = fe.PropertyParameterBlock()

    m.fs.pc = PressureChanger(property_package=m.fs.props,
                              thermodynamic_assumption='pump')

    m.fs.pc.inlet[:].vars["flow_mol_comp"]["H2O"].fix(27.5e3)
    m.fs.pc.inlet[:].vars["temperature"].fix(866.5)
    m.fs.pc.inlet[:].vars["pressure"].fix(2.97e7)

    m.fs.pc.deltaP.fix(-1e7)
    m.fs.pc.efficiency_pump.fix(0.9)

    assert degrees_of_freedom(m) == 0

    m.fs.pc.initialize()

    assert (pytest.approx(27500.0, abs=1e-2) ==
            m.fs.pc.outlet[0].vars["flow_mol_comp"]["H2O"].value)
    assert (pytest.approx(775.658, abs=1e-2) ==
            m.fs.pc.outlet[0].vars["temperature"].value)
    assert (pytest.approx(19700000.0, abs=1e-2) ==
            m.fs.pc.outlet[0].vars["pressure"].value)
    assert (pytest.approx(-79770875.82, abs=1e-2) ==
            m.fs.pc.work_mechanical[0].value)


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialization_isentropic():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.props = fe.PropertyParameterBlock()

    m.fs.pc = PressureChanger(property_package=m.fs.props,
                              thermodynamic_assumption='isentropic')
    m.fs.pc.inlet[:].vars["flow_mol_comp"]["H2O"].fix(27.5e3)
    m.fs.pc.inlet[:].vars["temperature"].fix(866.5)
    m.fs.pc.inlet[:].vars["pressure"].fix(2.97e7)

    m.fs.pc.deltaP.fix(-1e7)
    m.fs.pc.efficiency_isentropic.fix(0.83)

    assert degrees_of_freedom(m) == 0

    m.fs.pc.initialize()

    assert (pytest.approx(27500.0, abs=1e-2) ==
            m.fs.pc.outlet[0].vars["flow_mol_comp"]["H2O"].value)
    assert (pytest.approx(745.582, abs=1e-2) ==
            m.fs.pc.outlet[0].vars["temperature"].value)
    assert (pytest.approx(19700000.0, abs=1e-2) ==
            m.fs.pc.outlet[0].vars["pressure"].value)
    assert (pytest.approx(-136838091.58, abs=1e-2) ==
            m.fs.pc.work_mechanical[0].value)
