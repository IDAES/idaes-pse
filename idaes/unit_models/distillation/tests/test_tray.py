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
Tests for tray mixer unit model.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import (ConcreteModel, TerminationCondition,
                           SolverStatus, value)

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
                        MomentumBalanceType)
from idaes.unit_models.distillation import Tray
from idaes.property_models.activity_coeff_models.BTX_activity_coeff_VLE \
    import BTXParameterBlock
from idaes.core.util.model_statistics import degrees_of_freedom, \
    number_variables, number_total_constraints
from idaes.core.util.testing import get_default_solver


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = BTXParameterBlock(default={"valid_phase":
                                             ('Liq', 'Vap'),
                                             "activity_coeff_model":
                                             "Ideal"})

###############################################################################

m.fs.tray = Tray(default={"property_package": m.fs.properties,
                          "is_feed_tray": True,
                          "has_liquid_side_draw": False,
                          "has_vapor_side_draw": False,
                          "has_heat_transfer": False,
                          "has_pressure_change": True})

# Set inputs
m.fs.tray.feed.flow_mol.fix(1)
m.fs.tray.feed.temperature.fix(369)
m.fs.tray.feed.pressure.fix(101325)
m.fs.tray.feed.mole_frac_comp[0, "benzene"].fix(0.5)
m.fs.tray.feed.mole_frac_comp[0, "toluene"].fix(0.5)

m.fs.tray.liq_in.flow_mol.fix(1)
m.fs.tray.liq_in.temperature.fix(369)
m.fs.tray.liq_in.pressure.fix(101325)
m.fs.tray.liq_in.mole_frac_comp[0, "benzene"].fix(0.5)
m.fs.tray.liq_in.mole_frac_comp[0, "toluene"].fix(0.5)

m.fs.tray.vap_in.flow_mol.fix(1)
m.fs.tray.vap_in.temperature.fix(372)
m.fs.tray.vap_in.pressure.fix(101325)
m.fs.tray.vap_in.mole_frac_comp[0, "benzene"].fix(0.5)
m.fs.tray.vap_in.mole_frac_comp[0, "toluene"].fix(0.5)

m.fs.tray.deltaP.fix(0)
# m.fs.tray.heat_duty.fix(0)
# m.fs.tray.liq_side_sf.fix(0.5)
# m.fs.tray.vap_side_sf.fix(0.5)
m.fs.tray.initialize(outlvl=2, solver=solver)

solve_status = solver.solve(m, tee=True)

m.fs.tray.properties_out[0].flow_mol.display()
m.fs.tray.properties_out[0].flow_mol_phase.display()
m.fs.tray.properties_out[0].mole_frac_comp.display()
m.fs.tray.properties_out[0].mole_frac_phase_comp.display()
m.fs.tray.properties_out[0].temperature.display()

m.fs.tray.liq_in.display()
m.fs.tray.vap_in.display()

m.fs.tray.vap_out.display()
m.fs.tray.liq_out.display()
