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

