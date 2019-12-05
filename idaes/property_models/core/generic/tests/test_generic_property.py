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
Tests for generic property package core code

Author: Andrew Lee
"""
import pytest
from sys import modules

from pyomo.environ import Block, ConcreteModel, Set
from pyomo.common.config import ConfigBlock, ConfigValue

from idaes.property_models.core.generic.generic_property import (
        GenericPropertyPackageError,
        get_method,
        GenericParameterData,
        GenericStateBlock)

from idaes.core import declare_process_block_class
from idaes.core.util.exceptions import ConfigurationError, PropertyPackageError
from idaes.core.util.misc import add_object_reference


# -----------------------------------------------------------------------------
supported_properties = ["phase_component_list",
                        "state_definition",
                        "state_bounds",
                        "phase_equilibrium_formulation",
                        "phase_equilibrium_dict",
                        "equation_of_state",
                        "bubble_temperature",
                        "dew_temperature",
                        "bubble_pressure",
                        "dew_pressure",
                        "dens_mol_liq_comp",
                        "enth_mol_liq_comp",
                        "enth_mol_ig_comp",
                        "entr_mol_liq_comp",
                        "entr_mol_ig_comp",
                        "pressure_sat_comp"]


# -----------------------------------------------------------------------------
def test_GenericPropertyPackageError():
    with pytest.raises(PropertyPackageError):
        raise GenericPropertyPackageError("block", "prop")


# Dummy method for testing get_method
def dummy_option(b):
    return b.dummy_response


class TestGetMethod(object):
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = Block()

        # Add necessary parameters to parameter block
        m.params.config = ConfigBlock()
        m.params.config.declare("dummy_option", ConfigValue(default=None))

        # Create a dummy state block
        m.props = Block([1])
        m.props[1].config = ConfigBlock()
        add_object_reference(m.props[1], "_params", m.params)

        m.props[1].dummy_response = "foo"

        return m

    def test_None(self, frame):
        with pytest.raises(GenericPropertyPackageError):
            get_method(frame.props[1], "dummy_option")

    def test_invlaid_option(self, frame):
        with pytest.raises(AttributeError):
            get_method(frame.props[1], "foo")

    def test_method(self, frame):
        # Test that get_method works when provided with the method directly
        frame.params.config.dummy_option = dummy_option

        assert get_method(frame.props[1], "dummy_option") is dummy_option
        assert get_method(frame.props[1], "dummy_option")(frame.props[1]) is \
            frame.props[1].dummy_response

    def test_module(self, frame):
        # Test that get_method works when pointed to a module with the method
        frame.params.config.dummy_option = modules[__name__]

        assert get_method(frame.props[1], "dummy_option") is dummy_option
        assert get_method(frame.props[1], "dummy_option")(frame.props[1]) is \
            frame.props[1].dummy_response

    def test_not_method_or_module(self, frame):
        # Test that get_method works when provided with the method directly
        frame.params.config.dummy_option = "foo"

        with pytest.raises(ConfigurationError):
            get_method(frame.props[1], "dummy_option")


# -----------------------------------------------------------------------------
@declare_process_block_class("DummyParameterBlock")
class DummyParameterData(GenericParameterData):
    def configure(self):
        self.config.component_list = ["a", "b", "c"]
        self.config.phase_list = [1, 2]

    def parameters(self):
        self.parameters_set = True


class TestGenericParameterBlock(object):
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = DummyParameterBlock()

        return m

    def test_config(self, frame):
        assert frame.params.config.component_list == ["a", "b", "c"]
        assert frame.params.config.phase_list == [1, 2]

        # Tesrt number of config arguments. Note 1 inherited argument
        assert len(frame.params.config) == len(supported_properties) + 2 + 1

        for i in supported_properties:
            assert frame.params.config[i] is None

    def test_build(self, frame):
        # Build core components
        assert frame.params.state_block_class is GenericStateBlock

        assert isinstance(frame.params.component_list, Set)
        assert len(frame.params.component_list) == len(
                frame.params.config.component_list)
        for j in frame.params.component_list:
            assert j in frame.params.config.component_list

        assert isinstance(frame.params.phase_list, Set)
        assert len(frame.params.phase_list) == len(
                frame.params.config.phase_list)
        for p in frame.params.phase_list:
            assert p in frame.params.config.phase_list

#
#        # Validate that user provided either both a phase equilibrium
#        # formulation and a dict of phase equilibria or neither
#        if ((self.config.phase_equilibrium_formulation is not None) ^
#                (self.config.phase_equilibrium_dict is not None)):
#            raise ConfigurationError(
#                    "{} Generic Property Package provided with only one of "
#                    "phase_equilibrium_formulation and phase_equilibrium_dict."
#                    " Either both of these arguments need to be provided or "
#                    "neither.".format(self.name))
#
#        # Build phase equilibrium list
#        if self.config.phase_equilibrium_dict is not None:
#            self.phase_equilibrium_list = self.config.phase_equilibrium_dict
#
#            pe_set = []
#            for k in self.config.phase_equilibrium_dict.keys():
#                pe_set.append(k)
#            self.phase_equilibrium_idx = Set(initialize=pe_set,
#                                             ordered=True)
