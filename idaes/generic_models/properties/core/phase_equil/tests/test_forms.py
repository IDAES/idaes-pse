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
Library of common forms for phase equilibrium constraints
"""
from pyomo.environ import ConcreteModel, Var

from idaes.generic_models.properties.core.generic.generic_property import \
    GenericParameterBlock
from idaes.generic_models.properties.core.state_definitions import FTPx

from idaes.generic_models.properties.core.phase_equil.forms import *


# Dummy EoS to use for fugacity calls
class DummyEoS(object):
    def common(self):
        pass

    def fug_phase_comp(b, p, j):
        return b.x[p, j]


def test_fugacity():
    m = ConcreteModel()

    # Add a dummy var for use in constructing expressions
    m.x = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    # Create a dummy parameter block
    m.params = GenericParameterBlock(default={
        "components": {"H2O": {"temperature_crit": 647.3,
                               "phase_equilibrium_form": fugacity}},
        "phases": {"Liq": {"equation_of_state": DummyEoS},
                   "Vap": {"equation_of_state": DummyEoS}},
        "state_definition": FTPx,
        "pressure_ref": 1e5,
        "temperature_ref": 300})

    assert str(fugacity(m, "Vap", "Liq", "H2O")) == str(
        m.x["Vap", "H2O"] == m.x["Liq", "H2O"])
