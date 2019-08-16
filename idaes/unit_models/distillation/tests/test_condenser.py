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
Tests for Heat Exchanger 1D unit model.

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
m.fs.C101 = Condenser(default={"property_package": m.fs.properties,
                      "condenser_type": CondenserType.totalCondenser})

# Fix the condenser variables
m.fs.C101.reflux_ratio.fix(1)
m.fs.C101.deltaP.fix(0)

# Fix the inputs (typically this will be the outlet vapor from the top tray)
m.fs.C101.inlet.flow_mol.fix(1)
m.fs.C101.inlet.temperature.fix(368)
m.fs.C101.inlet.pressure.fix(101325)
m.fs.C101.inlet.mole_frac[0, "benzene"].fix(0.5)
m.fs.C101.inlet.mole_frac[0, "toluene"].fix(0.5)

print("The degrees of freedom is ", degrees_of_freedom(m))


solver = get_default_solver()
m.fs.C101.initialize(solver=solver, outlvl=2)

m.fs.C101.reflux.display()
m.fs.C101.distillate.display()
m.fs.C101.heat_duty.display()
