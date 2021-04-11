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
Test categorize_dae_variables_and_constraints function.
"""

import pytest
import pyomo.environ as pyo
from pyomo.network import Arc
from pyomo.common.collections import ComponentSet
from pyomo.dae.flatten import flatten_dae_components

from idaes.apps.caprese.categorize import (
        categorize_dae_variables_and_constraints,
        )
from idaes.apps.caprese.common.config import VariableCategory as VC
from idaes.apps.caprese.common.config import ConstraintCategory as VC

__author__ = "Robert Parker"

@pytest.mark.unit
def test_categorize_deriv():
    """ The simplest test. Identify a differential and a derivative var.
    """
    pass
