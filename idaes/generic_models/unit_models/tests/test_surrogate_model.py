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
Generic template for a surrogate unit model.
"""

import pytest
from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           units,
                           value,
                           Var,
                           Constraint)
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)
from idaes.generic_models.unit_models.surrogate_model import SurrogateModel
from idaes.generic_models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock)
from idaes.generic_models.properties.examples.saponification_reactions import (
    SaponificationReactionParameterBlock)
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util.testing import (get_default_solver,
                                     PhysicalParameterTestBlock,
                                     ReactionParameterTestBlock,
                                     initialization_tester)
from pyomo.util.check_units import (assert_units_consistent,
                                    assert_units_equivalent)

__author__ = "Jaffer Ghouse"


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


# -----------------------------------------------------------------------------
m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.bfb_surrogate = SurrogateModel()
m.fs.bfb_surrogate.x_1 = Var(initialize=0.5)
m.fs.bfb_surrogate.x_2 = Var(initialize=0.5)
m.fs.bfb_surrogate.y_1 = Var(initialize=0.5)
m.fs.bfb_surrogate.y_2 = Var(initialize=0.5)
m.fs.bfb_surrogate.c1 = Constraint(
    expr=m.fs.bfb_surrogate.y_1 == m.fs.bfb_surrogate.x_1)
m.fs.bfb_surrogate.c2 = Constraint(
    expr=m.fs.bfb_surrogate.y_2 == m.fs.bfb_surrogate.x_2)

inlet_dict = {"input_1": m.fs.bfb_surrogate.x_1,
              "input_2": m.fs.bfb_surrogate.x_2}
outlet_dict = {"output_1": m.fs.bfb_surrogate.y_1,
               "output_2": m.fs.bfb_surrogate.y_2}

m.fs.bfb_surrogate._add_ports(name="inlet",
                              member_list=inlet_dict)
m.fs.bfb_surrogate._add_ports(name="outlet",
                              member_list=outlet_dict)
