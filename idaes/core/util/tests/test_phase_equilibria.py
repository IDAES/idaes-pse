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
This module contains phase equilibria utility functions for use in IDAES models.
"""

import pytest

from pyomo.environ import ConcreteModel, Expression, Set, Block, Var
# Import Pyomo units
from pyomo.environ import units as pyunits
from idaes.core import LiquidPhase, VaporPhase, Component

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.generic_models.properties.core.phase_equil import smooth_VLE
from idaes.generic_models.properties.core.phase_equil.bubble_dew import \
        IdealBubbleDew
from idaes.generic_models.properties.core.phase_equil.forms import fugacity
import idaes.generic_models.properties.core.pure.NIST as NIST
from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)

from idaes.core.util.phase_equilibria import (TXYDataClass,Txy_data)

# Author: Alejandro Garciadiego
@pytest.mark.unit
def test_Txy_data():

    TD = TXYDataClass('benzene','toluene',kg / m / s ** 2,K,101325)
    TD.TBubb = [353.3205, 365.3478, 383.8817]
    TD.TDew = [353.3237, 372.0203, 383.8845]
    TD.x = [0.9999, 0.5, 9.9999]

    assert TD.Component_1 == 'benzene'
    assert TD.Component_2 == 'toluene'
    assert TD.Punits == pyunits.kg / pyunits.m / pyunits.s ** 2
    assert TD.Tunits == pyunits.K
    assert TD.P == 101325
    assert TD.TBubb == [pytest.approx(353.3205, abs=1e-4), pytest.approx(365.3478, abs=1e-4), pytest.approx(383.8817, abs=1e-4)]
    assert TD.TDew == [pytest.approx(353.3237, abs=1e-4), pytest.approx(372.0203, abs=1e-4), pytest.approx(383.8845, abs=1e-4)]
    assert TD.x == [pytest.approx(0.9999, abs=1e-4), pytest.approx(0.5, abs=1e-4), pytest.approx(9.9999, abs=1e-4)]

# Author: Alejandro Garciadiego
@pytest.mark.unit
def test_Txy_data():
    configuration = {
        # Specifying components
        "components": {
            'benzene': {"type": Component,
                        "pressure_sat_comp": NIST,
                        "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                        "parameter_data": {
                            "mw": (78.1136E-3, pyunits.kg/pyunits.mol),  # [1]
                            "pressure_crit": (48.9e5, pyunits.Pa),  # [1]
                            "temperature_crit": (562.2, pyunits.K),  # [1]
                            "pressure_sat_comp_coeff":{"A": (4.72583, None),  # [2]
                                                        "B": (1660.652, pyunits.K),
                                                        "C": (-1.461, pyunits.K)}}},
            'toluene': {"type": Component,
                        "pressure_sat_comp": NIST,
                        "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                        "parameter_data": {
                            "mw": (92.1405E-3, pyunits.kg/pyunits.mol),  # [1]
                            "pressure_crit": (41e5, pyunits.Pa),  # [1]
                            "temperature_crit": (591.8, pyunits.K),  # [1]
                            "pressure_sat_comp_coeff":{"A": (4.07827, None),  # [2]
                                                        "B": (1343.943, pyunits.K),
                                                        "C": (-53.773, pyunits.K)}}}},

        # Specifying phases
        "phases":  {'Liq': {"type": LiquidPhase,
                            "equation_of_state": Ideal},
                    'Vap': {"type": VaporPhase,
                            "equation_of_state": Ideal}},

        # Set base units of measurement
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K},

        # Specifying state definition
        "state_definition": FTPx,
        "state_bounds": {"flow_mol": (0, 100, 1000, pyunits.mol/pyunits.s),
                         "temperature": (273.15, 300, 450, pyunits.K),
                         "pressure": (5e4, 1e5, 1e6, pyunits.Pa)},
        "pressure_ref": (1e5, pyunits.Pa),
        "temperature_ref": (300, pyunits.K),

        # Defining phase equilibria
        "phases_in_equilibrium": [("Vap", "Liq")],
        "phase_equilibrium_state": {("Vap", "Liq"): smooth_VLE},
        "bubble_dew_method": IdealBubbleDew}

    model = ConcreteModel()

    model.params = GenericParameterBlock(default=configuration)

    TD = Txy_data('benzene','toluene', 101325, 3, model)

    assert TD.Component_1 == 'benzene'
    assert TD.Component_2 == 'toluene'
    assert TD.Punits == pyunits.kg / pyunits.m / pyunits.s ** 2
    assert TD.Tunits == pyunits.K
    assert TD.P == 101325
    assert TD.TBubb == [pytest.approx(353.3853, abs=1e-4), pytest.approx(365.2127, abs=1e-4), pytest.approx(383.5312, abs=1e-4)]
    assert TD.TDew == [pytest.approx(353.5429, abs=1e-2), pytest.approx(371.8702, abs=1e-4), pytest.approx(383.6709, abs=1e-4)]
    assert TD.x == [pytest.approx(0.995, abs=1e-4), pytest.approx(0.5, abs=1e-4), pytest.approx(0.005, abs=1e-4)]
