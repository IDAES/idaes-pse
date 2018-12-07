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
Tests for ControlVolumeBase.

Author: Emmanuel Ogbe
"""
import pytest

from pyomo.environ import ConcreteModel, SolverFactory

from idaes.core import FlowsheetBlock
from idaes.unit_models.pressure_changer import PressureChanger
from idaes.property_models.iapws95.iapws95_wrap_ph import (
                        Iapws95ParameterBlockData)
from idaes.ui.report import degrees_of_freedom


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

    m.fs.properties = Iapws95ParameterBlockData()

    m.fs.presschanger = PressureChanger(default={"property_package": m.fs.properties,
                                       "has_heat_transfer": True})

    assert hasattr(m.fs.presschanger, "inlet")
    assert len(m.fs.presschanger.inlet[0].vars) == 3
    assert hasattr(m.fs.presschanger.inlet[0], "flow_mol_comp")
    assert hasattr(m.fs.presschanger.inlet[0], "temperature")
    assert hasattr(m.fs.presschanger.inlet[0], "pressure")

    assert hasattr(m.fs.presschanger, "outlet")
    assert len(m.fs.presschanger.outlet[0].vars) == 3
    assert hasattr(m.fs.presschanger.outlet[0], "flow_mol_comp")
    assert hasattr(m.fs.presschanger.outlet[0], "temperature")
    assert hasattr(m.fs.presschanger.outlet[0], "pressure")

    assert hasattr(m.fs.presschanger, "gibbs_minimization")
    assert hasattr(m.fs.presschanger.control_volume, "heat")
    assert hasattr(m.fs.presschanger, "heat_duty")
