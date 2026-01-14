#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Simple Flash flowsheet for use in testing.
"""

from pyomo.environ import ConcreteModel, SolverFactory
from idaes.core import FlowsheetBlock

# Import idaes logger to set output levels
import idaes.logger as idaeslog
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.unit_models import Flash
from ..fsrunner import FlowsheetRunner


FS = FlowsheetRunner()

# # Flash Unit Model
#
# Author: Jaffer Ghouse
# Maintainer: Andrew Lee
# Updated: 2023-06-01


@FS.step("build")
def build_model(ctx):
    """Build the model."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BTXParameterBlock(
        valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal", state_vars="FTPz"
    )
    m.fs.flash = Flash(property_package=m.fs.properties)
    # assert degrees_of_freedom(m) == 7
    ctx.model = m


@FS.step("set_operating_conditions")
def set_operating_conditions(ctx):
    """Set operating conditions."""
    m = ctx.model
    m.fs.flash.inlet.flow_mol.fix(1)
    m.fs.flash.inlet.temperature.fix(368)
    m.fs.flash.inlet.pressure.fix(101325)
    m.fs.flash.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.flash.inlet.mole_frac_comp[0, "toluene"].fix(0.5)
    m.fs.flash.heat_duty.fix(0)
    m.fs.flash.deltaP.fix(0)


@FS.step("initialize")
def init_model(ctx):
    """ "Initialize the model."""
    m = ctx.model
    m.fs.flash.initialize(outlvl=idaeslog.INFO)


@FS.step("set_solver")
def set_solver(ctx):
    """Set the solver."""
    ctx.solver = SolverFactory("ipopt")


@FS.step("solve_initial")
def solve(ctx):
    """Perform the initial model solve."""
    ctx["status"] = ctx.solver.solve(ctx.model, tee=ctx["tee"])
