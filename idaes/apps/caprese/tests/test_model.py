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
"""

import pyomo.environ as aml
import pyomo.dae as dae
import pyomo.network as pyn
from pyomo.common.collections import ComponentSet
from pyomo.core.expr.visitor import identify_variables
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
        MomentumBalanceType)
from idaes.core.util.model_statistics import (degrees_of_freedom,
        activated_equalities_generator, unfixed_variables_generator)
from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util.exceptions import ConfigurationError
from idaes.apps.caprese.tests.test_simple_model import (
        make_model, 
        make_small_model,
        initialize_t0,
        copy_values_forward,
        )
from idaes.apps.caprese.model import DynamicBlock
import idaes.logger as idaeslog
import random
import pytest

__author__ = "Robert Parker"

class TestDynamicBlock(object):

    def make_block(self, sample_time=0.5, horizon=1, nfe=2):
        model = make_model(horizon=horizon, nfe=nfe)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in[t0]]
        block = DynamicBlock(model=model, time=time, inputs=inputs)
        return block

    def test_constructor(self):
        block = self.make_block()
