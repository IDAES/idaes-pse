##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Pytest for Compression Model

"""
import pytest
# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.util.check_units import assert_units_consistent


# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock)
from idaes.core.util.model_statistics import degrees_of_freedom

# Import Unit Model Modules
from idaes.power_generation.carbon_capture.compression_system.compressor \
      import (CompressionStage, VaneDiffuserType, ImpellerType)

# Import Property Package
from idaes.power_generation.properties.natural_gas_PR import get_prop

from idaes.core.util.testing import get_default_solver, initialization_tester

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


@pytest.fixture(scope="module")
def build_unit():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()
    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})
    NG_config = get_prop(
        components=["H2O", 'CO2'])
    m.fs.NG_props = GenericParameterBlock(default=NG_config)

    m.fs.unit = CompressionStage(
        default={"property_package": m.fs.NG_props,
                 "impeller_type": ImpellerType.open_impeller,
                 "vane_diffuser_type": VaneDiffuserType.vane_diffuser})

    # Set the compressor inlet conditions and an initial flow guess
    p = 1.12308 * 1e5  # Pa
    t = 38.4735 + 273.15  # K
    fin = 1698.37  # mol/s

    m.fs.unit.inlet.flow_mol[:].fix(fin)
    m.fs.unit.inlet.temperature.fix(t)
    m.fs.unit.inlet.pressure[:].fix(p)
    m.fs.unit.inlet.mole_frac_comp[:, "H2O"].fix(0.061396)
    m.fs.unit.inlet.mole_frac_comp[:, "CO2"].fix(0.938604)

    # inlet specifications
    m.fs.unit.U2.fix(315.3)
    m.fs.unit.outlet.pressure[:].fix(2.73423*1e5)

    # fix compressor specification
    m.fs.unit.speed_of_sound.fix(279.214)
    m.fs.unit.r2.fix(0.67654)
    m.fs.unit.z_s.fix(0.9952)
    m.fs.unit.z_d1.fix(0.97373)
    m.fs.unit.efficiency_mech.fix(0.97)
    m.fs.unit.eff_drive.fix(1.0)

    return m


@pytest.mark.unit
def test_basic_build(build_unit):
    """Make a model and make sure it doesn't throw exception"""
    m = build_unit
    assert degrees_of_freedom(m) == 0
    # Check unit config arguments
    assert len(m.fs.unit.config) == 16
    assert m.fs.unit.config.thermodynamic_assumption
    assert m.fs.unit.config.property_package is m.fs.NG_props


@pytest.mark.integration
def test_units(build_unit):
    assert_units_consistent(build_unit)


@pytest.mark.component
def test_initialize(build_unit):
    initialization_tester(build_unit, dof=0)


@pytest.mark.component
def test_run(build_unit):
    m = build_unit

    optarg = {"tol": 1e-6,
              "linear_solver": "ma27",
              "max_iter": 50}
    solver.options = optarg
    # solve model
    results = solver.solve(m, tee=True)
    # Check for optimal solution
    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok
    assert degrees_of_freedom(m) == 0

    # work
    assert (pytest.approx(5092763.0, abs=0.1) ==
            pyo.value(m.fs.unit.work_mechanical[0]))
    # pressure change
    assert (pytest.approx(161115.0, abs=0.1) ==
            pyo.value(m.fs.unit.deltaP[0]))
    # mass balance
    assert (pytest.approx(0, abs=1e-2) ==
            pyo.value(m.fs.unit.inlet.flow_mol[0]
                      - m.fs.unit.outlet.flow_mol[0]
                      ))
