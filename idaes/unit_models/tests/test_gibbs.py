##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Tests for ControlVolumeBlockData.

Author: Andrew Lee
"""
import pytest

from pyomo.environ import ConcreteModel, SolverFactory

from idaes.core import FlowsheetBlock
from idaes.unit_models.gibbs_reactor import GibbsReactor
from idaes.property_models.examples.methane_combustion_ideal import (
    MethaneCombustionParameterBlock)
from idaes.core.util.model_statistics import degrees_of_freedom


# -----------------------------------------------------------------------------
# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None


# -----------------------------------------------------------------------------
def test_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = MethaneCombustionParameterBlock()

    m.fs.gibbs = GibbsReactor(default={"property_package": m.fs.properties,
                                       "has_heat_transfer": True})

    assert hasattr(m.fs.gibbs, "inlet")
    assert len(m.fs.gibbs.inlet.vars) == 3
    assert hasattr(m.fs.gibbs.inlet, "flow_mol_comp")
    assert hasattr(m.fs.gibbs.inlet, "temperature")
    assert hasattr(m.fs.gibbs.inlet, "pressure")

    assert hasattr(m.fs.gibbs, "outlet")
    assert len(m.fs.gibbs.outlet.vars) == 3
    assert hasattr(m.fs.gibbs.outlet, "flow_mol_comp")
    assert hasattr(m.fs.gibbs.outlet, "temperature")
    assert hasattr(m.fs.gibbs.outlet, "pressure")

    assert hasattr(m.fs.gibbs, "gibbs_minimization")
    assert hasattr(m.fs.gibbs.control_volume, "heat")
    assert hasattr(m.fs.gibbs, "heat_duty")


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize_temperature():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = MethaneCombustionParameterBlock()

    m.fs.gibbs = GibbsReactor(default={"property_package": m.fs.properties,
                                       "has_heat_transfer": True})

    m.fs.gibbs.inlet.flow_mol_comp[0, "H2"].fix(10.0)
    m.fs.gibbs.inlet.flow_mol_comp[0, "N2"].fix(150.0)
    m.fs.gibbs.inlet.flow_mol_comp[0, "O2"].fix(40.0)
    m.fs.gibbs.inlet.flow_mol_comp[0, "CO2"].fix(1e-5)
    m.fs.gibbs.inlet.flow_mol_comp[0, "CH4"].fix(30.0)
    m.fs.gibbs.inlet.flow_mol_comp[0, "CO"].fix(1e-5)
    m.fs.gibbs.inlet.flow_mol_comp[0, "H2O"].fix(1e-5)
    m.fs.gibbs.inlet.flow_mol_comp[0, "NH3"].fix(1e-5)
    m.fs.gibbs.inlet.temperature[0].fix(1500.0)
    m.fs.gibbs.inlet.pressure[0].fix(101325.0)

    m.fs.gibbs.outlet.temperature[0].fix(2844.38)

    assert degrees_of_freedom(m) == 0

    m.fs.gibbs.initialize(outlvl=5,
                          optarg={'tol': 1e-6},
                          state_args={'temperature': 2844.38,
                                      'pressure': 101325.0,
                                      'flow_mol_comp': {'CH4': 1e-5,
                                                        'CO': 23.0,
                                                        'CO2': 7.05,
                                                        'H2': 29.0,
                                                        'H2O': 41.0,
                                                        'N2': 150.0,
                                                        'NH3': 1e-5,
                                                        'O2': 1.0}})

    assert (pytest.approx(0.0, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "CH4"].value)
    assert (pytest.approx(22.614, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "CO"].value)
    assert (pytest.approx(7.386, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "CO2"].value)
    assert (pytest.approx(28.806, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "H2"].value)
    assert (pytest.approx(41.194, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "H2O"].value)
    assert (pytest.approx(150.0, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "N2"].value)
    assert (pytest.approx(0.0, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "NH3"].value)
    assert (pytest.approx(0.710, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "O2"].value)
    assert (pytest.approx(161882.3, abs=1e-2) ==
            m.fs.gibbs.heat_duty[0].value)
    assert (pytest.approx(101325.0, abs=1e-2) ==
            m.fs.gibbs.outlet.pressure[0].value)


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize_heat_duty():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = MethaneCombustionParameterBlock()

    m.fs.gibbs = GibbsReactor(default={"property_package": m.fs.properties,
                                       "has_heat_transfer": True})

    m.fs.gibbs.inlet.flow_mol_comp[0, "H2"].fix(10.0)
    m.fs.gibbs.inlet.flow_mol_comp[0, "N2"].fix(150.0)
    m.fs.gibbs.inlet.flow_mol_comp[0, "O2"].fix(40.0)
    m.fs.gibbs.inlet.flow_mol_comp[0, "CO2"].fix(1e-5)
    m.fs.gibbs.inlet.flow_mol_comp[0, "CH4"].fix(30.0)
    m.fs.gibbs.inlet.flow_mol_comp[0, "CO"].fix(1e-5)
    m.fs.gibbs.inlet.flow_mol_comp[0, "H2O"].fix(1e-5)
    m.fs.gibbs.inlet.flow_mol_comp[0, "NH3"].fix(1e-5)
    m.fs.gibbs.inlet.temperature[0].fix(1500.0)
    m.fs.gibbs.inlet.pressure[0].fix(101325.0)

    m.fs.gibbs.heat_duty.fix(161882.303661)

    assert degrees_of_freedom(m) == 0

    m.fs.gibbs.initialize(outlvl=5,
                          optarg={'tol': 1e-6},
                          state_args={'temperature': 2844.38,
                                      'pressure': 101325.0,
                                      'flow_mol_comp': {'CH4': 1e-5,
                                                        'CO': 23.0,
                                                        'CO2': 7.05,
                                                        'H2': 29.0,
                                                        'H2O': 41.0,
                                                        'N2': 150.0,
                                                        'NH3': 1e-5,
                                                        'O2': 1.0}})

    assert (pytest.approx(0.0, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "CH4"].value)
    assert (pytest.approx(22.614, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "CO"].value)
    assert (pytest.approx(7.386, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "CO2"].value)
    assert (pytest.approx(28.806, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "H2"].value)
    assert (pytest.approx(41.194, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "H2O"].value)
    assert (pytest.approx(150.0, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "N2"].value)
    assert (pytest.approx(0.0, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "NH3"].value)
    assert (pytest.approx(0.710, abs=1e-2) ==
            m.fs.gibbs.outlet.flow_mol_comp[0, "O2"].value)
    assert (pytest.approx(2844.38, abs=1e-2) ==
            m.fs.gibbs.outlet.temperature[0].value)
    assert (pytest.approx(101325.0, abs=1e-2) ==
            m.fs.gibbs.outlet.pressure[0].value)
