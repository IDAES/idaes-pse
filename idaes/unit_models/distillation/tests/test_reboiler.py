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
Tests for Reboiler unit model.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import (ConcreteModel, TerminationCondition,
                           SolverStatus, value)

from idaes.core import FlowsheetBlock
from idaes.unit_models.distillation import Reboiler
from idaes.property_models.activity_coeff_models.BTX_activity_coeff_VLE \
    import BTXParameterBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import get_default_solver

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = BTXParameterBlock(default={"valid_phase":
                                             ('Liq', 'Vap'),
                                             "activity_coeff_model":
                                             "Ideal"})
m.fs.properties_2 = BTXParameterBlock(default={"valid_phase":
                                               ('Liq', 'Vap'),
                                               "activity_coeff_model":
                                               "Ideal",
                                               "state_vars":
                                               "FcTP"})

###############################################################################
# total reboiler with FTPz
m.fs.R101 = Reboiler(default={"property_package": m.fs.properties})

# Fix the partial reboiler variables
m.fs.R101.boilup_ratio.fix(1)
m.fs.R101.deltaP.fix(0)

# Fix the inputs (typically this will be the outlet vapor from the top tray)
m.fs.R101.inlet.flow_mol.fix(1)
m.fs.R101.inlet.temperature.fix(362)
m.fs.R101.inlet.pressure.fix(101325)
m.fs.R101.inlet.mole_frac[0, "benzene"].fix(0.5)
m.fs.R101.inlet.mole_frac[0, "toluene"].fix(0.5)

print("-----------------------------------------------------------------------")
print("Total Reboiler - FTPz")
print("The degrees of freedom is ", degrees_of_freedom(m))


solver = get_default_solver()
m.fs.R101.initialize(solver=solver, outlvl=2)

m.fs.R101.bottoms.display()
m.fs.R101.vapor_reboil.display()
m.fs.R101.heat_duty.display()

print("-----------------------------------------------------------------------")

###############################################################################
# total condenser with FcTP
m.fs.R101_FcTP = Reboiler(default={"property_package": m.fs.properties_2})

# Fix the partial condenser variables
m.fs.R101_FcTP.boilup_ratio.fix(1)
m.fs.R101_FcTP.deltaP.fix(0)

# Fix the inputs (typically this will be the outlet vapor from the top tray)
m.fs.R101_FcTP.inlet.flow_mol_comp[0, "benzene"].fix(0.5)
m.fs.R101_FcTP.inlet.flow_mol_comp[0, "toluene"].fix(0.5)
m.fs.R101_FcTP.inlet.temperature.fix(362)
m.fs.R101_FcTP.inlet.pressure.fix(101325)


print("-----------------------------------------------------------------------")
print("Reboiler - FcTP")
print("The degrees of freedom is ", degrees_of_freedom(m.fs.R101_FcTP))


solver = get_default_solver()
m.fs.R101_FcTP.initialize(solver=solver, outlvl=2)

m.fs.R101_FcTP.bottoms.display()
m.fs.R101_FcTP.vapor_reboil.display()
m.fs.R101_FcTP.heat_duty.display()

print("-----------------------------------------------------------------------")
