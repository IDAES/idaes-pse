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
Pytest for Compression Model

"""
import pytest

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.util.check_units import assert_units_consistent


# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

# Import Unit Model Modules
import idaes.models.properties.swco2 as swco2
from idaes.models_extra.carbon_capture.co2_compressor import (
    CompressionStage,
    VaneDiffuserType,
    ImpellerType,
)

from idaes.core.solvers import get_solver
from idaes.core.util.testing import initialization_tester

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


@pytest.fixture(scope="module")
def build_unit():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()
    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_co2 = swco2.SWCO2ParameterBlock(
        phase_presentation=swco2.PhaseType.G
    )
    m.fs.unit = CompressionStage(
        property_package=m.fs.properties_co2,
        impeller_type=ImpellerType.open_impeller,
        vane_diffuser_type=VaneDiffuserType.vane_diffuser,
    )

    # Set the compressor inlet conditions and an initial flow guess
    p = 1.0951 * 1e5  # Pa
    t = 40.0113 + 273.15  # K #
    fin = 1159.44  # mol/s #

    hin_co2 = swco2.htpx(T=t * pyo.units.K, P=p * pyo.units.Pa)

    m.fs.unit.inlet.flow_mol[:].fix(fin)
    m.fs.unit.inlet.enth_mol[:].fix(hin_co2)
    m.fs.unit.inlet.pressure[:].fix(p)

    # inlet specifications
    m.fs.unit.U2.fix(315.3)
    m.fs.unit.outlet.pressure[:].fix(2.53161 * 1e5)

    # fix compressor specification
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
    assert m.fs.unit.config.property_package is m.fs.properties_co2


@pytest.mark.integration
def test_units(build_unit):
    assert_units_consistent(build_unit)


@pytest.mark.component
def test_initialize(build_unit):
    initialization_tester(build_unit, dof=0)


@pytest.mark.component
def test_run(build_unit):
    m = build_unit

    optarg = {"tol": 1e-7, "linear_solver": "ma27", "max_iter": 50}
    solver.options = optarg
    # solve model
    results = solver.solve(m, tee=True)
    # Check for optimal solution
    assert pyo.check_optimal_termination(results)
    assert degrees_of_freedom(m) == 0

    # energy balance
    assert pytest.approx(0, abs=1e-3) == pyo.value(
        m.fs.unit.inlet.flow_mol[0] * m.fs.unit.inlet.enth_mol[0]
        - m.fs.unit.outlet.flow_mol[0] * m.fs.unit.outlet.enth_mol[0]
        + m.fs.unit.work_mechanical[0]
    )
    # pressure change
    assert pytest.approx(143651.0, abs=0.1) == pyo.value(m.fs.unit.deltaP[0])
    # mass balance
    assert pytest.approx(0, abs=1e-2) == pyo.value(
        m.fs.unit.inlet.flow_mol[0] - m.fs.unit.outlet.flow_mol[0]
    )
