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
Tests for Condenser unit model.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import (ConcreteModel, TerminationCondition,
                           SolverStatus, value)
from pyomo.common.config import ConfigBlock

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
                        MomentumBalanceType, useDefault)
from idaes.unit_models.distillation import Condenser
from idaes.unit_models.distillation.condenser import CondenserType
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
# total condenser with FTPz
m.fs.C101_total = Condenser(default={"property_package": m.fs.properties,
                                     "condenser_type": CondenserType.totalCondenser})

# Fix the partial condenser variables
m.fs.C101_total.reflux_ratio.fix(1)
m.fs.C101_total.deltaP.fix(0)

# Fix the inputs (typically this will be the outlet vapor from the top tray)
m.fs.C101_total.inlet.flow_mol.fix(1)
m.fs.C101_total.inlet.temperature.fix(368)
m.fs.C101_total.inlet.pressure.fix(101325)
m.fs.C101_total.inlet.mole_frac[0, "benzene"].fix(0.5)
m.fs.C101_total.inlet.mole_frac[0, "toluene"].fix(0.5)

print("-----------------------------------------------------------------------")
print("Total Condenser - FTPz")
print("The degrees of freedom is ", degrees_of_freedom(m))


solver = get_default_solver()
m.fs.C101_total.initialize(solver=solver, outlvl=2)

m.fs.C101_total.reflux.display()
m.fs.C101_total.distillate.display()
m.fs.C101_total.heat_duty.display()

print("-----------------------------------------------------------------------")

###############################################################################
# Partial Condenser with FTPz
m.fs.C101_partial = Condenser(default={"property_package": m.fs.properties,
                              "condenser_type": CondenserType.partialCondenser})

# Fix the partial condenser variables
m.fs.C101_partial.reflux_ratio.fix(0.5)
m.fs.C101_partial.deltaP.fix(0)
m.fs.C101_partial.vapor_outlet.temperature.fix(367.5)

# Fix the inputs (typically this will be the outlet vapor from the top tray)
m.fs.C101_partial.inlet.flow_mol.fix(1)
m.fs.C101_partial.inlet.temperature.fix(368)
m.fs.C101_partial.inlet.pressure.fix(101325)
m.fs.C101_partial.inlet.mole_frac[0, "benzene"].fix(0.5)
m.fs.C101_partial.inlet.mole_frac[0, "toluene"].fix(0.5)

print("Partial Condenser - FTPz")
print("The degrees of freedom is ", degrees_of_freedom(m))

m.fs.C101_partial.initialize(solver=solver, outlvl=2)

m.fs.C101_partial.reflux.display()
m.fs.C101_partial.distillate.display()
m.fs.C101_partial.vapor_outlet.display()
m.fs.C101_partial.heat_duty.display()

###############################################################################
# total condenser with FcTP
m.fs.C101_total_FcTP = Condenser(default={"property_package": m.fs.properties_2,
                                          "condenser_type":
                                          CondenserType.totalCondenser})

# Fix the partial condenser variables
m.fs.C101_total_FcTP.reflux_ratio.fix(1)
m.fs.C101_total_FcTP.deltaP.fix(0)

# Fix the inputs (typically this will be the outlet vapor from the top tray)
m.fs.C101_total_FcTP.inlet.flow_mol_comp[0, "benzene"].fix(0.5)
m.fs.C101_total_FcTP.inlet.flow_mol_comp[0, "toluene"].fix(0.5)
m.fs.C101_total_FcTP.inlet.temperature.fix(368)
m.fs.C101_total_FcTP.inlet.pressure.fix(101325)


print("-----------------------------------------------------------------------")
print("Total Condenser - FcTP")
print("The degrees of freedom is ", degrees_of_freedom(m.fs.C101_total_FcTP))


solver = get_default_solver()
m.fs.C101_total_FcTP.initialize(solver=solver, outlvl=2)

m.fs.C101_total_FcTP.reflux.display()
m.fs.C101_total_FcTP.distillate.display()
m.fs.C101_total_FcTP.heat_duty.display()

print("-----------------------------------------------------------------------")

# Partial Condenser with FcTP
m.fs.C101_partial_FcTP = Condenser(default={"property_package": m.fs.properties_2,
                                            "condenser_type":
                                            CondenserType.partialCondenser})

# Fix the partial condenser variables
m.fs.C101_partial_FcTP.reflux_ratio.fix(0.5)
m.fs.C101_partial_FcTP.deltaP.fix(0)
m.fs.C101_partial_FcTP.vapor_outlet.temperature.fix(367.5)

# Fix the inputs (typically this will be the outlet vapor from the top tray)
m.fs.C101_partial_FcTP.inlet.flow_mol_comp[0, "benzene"].fix(0.5)
m.fs.C101_partial_FcTP.inlet.flow_mol_comp[0, "toluene"].fix(0.5)
m.fs.C101_partial_FcTP.inlet.temperature.fix(368)
m.fs.C101_partial_FcTP.inlet.pressure.fix(101325)


print("Partial Condenser - FcTP")
print("The degrees of freedom is ", degrees_of_freedom(m.fs.C101_partial_FcTP))

m.fs.C101_partial_FcTP.initialize(solver=solver, outlvl=2)

m.fs.C101_partial_FcTP.reflux.display()
m.fs.C101_partial_FcTP.distillate.display()
m.fs.C101_partial_FcTP.vapor_outlet.display()
m.fs.C101_partial_FcTP.heat_duty.display()
