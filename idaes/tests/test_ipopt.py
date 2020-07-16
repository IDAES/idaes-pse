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
from pyomo.environ import SolverFactory
import pytest


@pytest.mark.unit
def test_ipopt_available():
    """
    Tries to set-up the IPOPT and returns exception if not available
    """
    if not SolverFactory('ipopt').available():
        raise Exception(
            "Could not find IPOPT. Users are strongly encouraged to have a "
            "version of IPOPT available, as it is the default solver assumed "
            "by many IDAES examples and tests. See the IDAES install "
            "documentation for instructions on how to get IPOPT.")
