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
Tests for eNRTL methods

Author: Andrew Lee
"""
import pytest

from pyomo.environ import (ConcreteModel,
                           value,
                           Var,
                           units as pyunits)

from idaes.core import (AqueousPhase,
                        Solvent,
                        Solute,
                        Apparent,
                        Anion,
                        Cation)
from idaes.generic_models.properties.core.eos.enrtl import ENRTL
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog


configuration = {
    "components": {
        "H2O": {"type": Solvent},
        "C6H12": {"type": Solute},
        "NaCl": {"type": Apparent},
        "HCl": {"type": Apparent},
        "NaOH": {"type": Apparent},
        "Na+": {"type": Cation,
                "charge": +1},
        "H+": {"type": Cation,
               "charge": +1},
        "Cl-": {"type": Anion,
                "charge": -1},
        "OH-": {"type": Anion,
                "charge": -1}},
    "phases": {
        "Liq": {"type": AqueousPhase,
                "equation_of_state": ENRTL}},
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "state_definition": FTPx,
    "pressure_ref": 1e5,
    "temperature_ref": 300}


class TestParameters(object):
    @pytest.mark.unit
    def test_parameters_no_assignment(self):
        m = ConcreteModel()

        m.params = GenericParameterBlock(default=configuration)

        assert isinstance(m.params.Liq.alpha, Var)
        assert len(m.params.Liq.alpha) == 15
        for (i, j) in m.params.Liq.alpha:
            if i != j:
                assert (j, i) not in m.params.Liq.alpha
            if (i, j) in [
                    ("C6H12", "C6H12"), ("H2O", "H2O"), ("H2O", "C6H12")]:
                assert m.params.Liq.alpha[(i, j)].value == 0.3
                assert m.params.Liq.alpha[(i, j)].fixed
            else:
                assert m.params.Liq.alpha[(i, j)].value == 0.2
                assert m.params.Liq.alpha[(i, j)].fixed

    @pytest.mark.unit
    def test_parameters_assignment(self):
        test_config = dict(configuration)
        test_config["parameter_data"] = {}
        test_config["parameter_data"]["Liq_alpha"] = {}
        test_config["parameter_data"]["Liq_alpha"][("H2O", "NaCl")] = 0.6

        m = ConcreteModel()

        m.params = GenericParameterBlock(default=test_config)

        assert isinstance(m.params.Liq.alpha, Var)
        assert len(m.params.Liq.alpha) == 15
        for (i, j) in m.params.Liq.alpha:
            if i != j:
                assert (j, i) not in m.params.Liq.alpha
            if (i, j) == ("H2O", "NaCl"):
                assert m.params.Liq.alpha[(i, j)].value == 0.6
                assert m.params.Liq.alpha[(i, j)].fixed
            elif (i, j) in [
                    ("C6H12", "C6H12"), ("H2O", "H2O"), ("H2O", "C6H12")]:
                assert m.params.Liq.alpha[(i, j)].value == 0.3
                assert m.params.Liq.alpha[(i, j)].fixed
            else:
                assert m.params.Liq.alpha[(i, j)].value == 0.2
                assert m.params.Liq.alpha[(i, j)].fixed

    @pytest.mark.unit
    def test_parameters_unsymmetric_alpha(self):
        test_config = dict(configuration)
        test_config["parameter_data"] = {}
        test_config["parameter_data"]["Liq_alpha"] = {}
        test_config["parameter_data"]["Liq_alpha"][("H2O", "NaCl")] = 0.6
        test_config["parameter_data"]["Liq_alpha"][("NaCl", "H2O")] = 0.8

        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params.Liq eNRTL alpha parameter assigned "
                           "non-symmetric value for pair \('H2O', 'NaCl'\). "
                           "Please assign only one value for component pair."):
            m.params = GenericParameterBlock(default=test_config)

    @pytest.mark.unit
    def test_parameters_alpha_symmetry_duplicate(self, caplog):
        caplog.set_level(
            idaeslog.INFO,
            logger=("idaes.generic_models.properties.core."
                    "generic.generic_property"))

        test_config = dict(configuration)
        test_config["parameter_data"] = {}
        test_config["parameter_data"]["Liq_alpha"] = {}
        test_config["parameter_data"]["Liq_alpha"][("H2O", "NaCl")] = 0.6
        test_config["parameter_data"]["Liq_alpha"][("NaCl", "H2O")] = 0.6

        m = ConcreteModel()

        m.params = GenericParameterBlock(default=test_config)

        assert ("eNRTL alpha value provided for both ('H2O', 'NaCl') and "
                "('NaCl', 'H2O'). It is only necessary to provide a "
                "value for one of these due to symmetry." in caplog.text)

    @pytest.mark.unit
    def test_parameters_alpha_unused_parameter(self):
        test_config = dict(configuration)
        test_config["parameter_data"] = {}
        test_config["parameter_data"]["Liq_alpha"] = {}
        test_config["parameter_data"]["Liq_alpha"][("H2O", "Na+")] = 0.6

        m = ConcreteModel()

        # TODO: Can't get regex to match tuple for some reason
        # For now, just chck the start of the expected string.
        with pytest.raises(ConfigurationError,
                           match="params.Liq eNRTL alpha parameter provided "
                           "for invalid component pair "):
            m.params = GenericParameterBlock(default=test_config)
