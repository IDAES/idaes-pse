#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for column model utitliy functions

Author: Alejandro Garciadiego, Andrew Lee
"""
import pytest
from types import MethodType

from pyomo.environ import ConcreteModel, Block, Var, Expression, Reference, Set, value
from pyomo.network import Port

from idaes.models_extra.column_models.util import make_phase_split

from idaes.core.util.exceptions import (
    PropertyPackageError,
    PropertyNotSupportedError,
    PropertyPackageError,
)

# Mark module as all unit tests
pytestmark = pytest.mark.unit


@pytest.fixture
def m():
    # Mock up something that looks like a tray model
    m = ConcreteModel()
    m.time = Set(initialize=[1, 2, 3])

    def flowsheet(m):
        return m

    m.flowsheet = MethodType(flowsheet, m)

    m.properties_out = Block(m.time)

    m.port = Port()

    return m


def test_make_phase_split_general(m):
    # Test for variables without special indices
    # Mock define_port_members method
    def define_port_members(m):
        return {"temperature": m.temperature, "indexed": m.indexed}

    # Add variables and define_port_members method to dummy state block
    for i in m.time:
        m.properties_out[i].temperature = Var()
        m.properties_out[i].indexed = Var(["a", "b"])

        m.properties_out[i].define_port_members = MethodType(
            define_port_members, m.properties_out[i]
        )

    # Call make_phase_split function
    make_phase_split(m, m.port)

    # Check for expected results
    for i in m.time:
        assert m.port.temperature[i] is m.properties_out[i].temperature
        for j in ["a", "b"]:
            assert m.port.indexed[i, j] is m.properties_out[i].indexed[j]


def test_make_phase_split_frac(m):
    # Mock define_port_members method
    def define_port_members(m):
        return {
            "flow_mass_phase": m.flow_mass_phase,
            "mass_frac_phase_comp": m.mass_frac_phase_comp,
        }

    def _liquid_set(m):
        return ["L", "L1"]

    def _vapor_set(m):
        return ["V", "V1"]

    m.components = Set(initialize=["a1", "b1"])
    m.phases = Set(initialize=["L", "V"])

    m.config.property_package = DummyData()
    m.config.property_package.component_list = m.components
    m.config.property_package.phase_list = m.phases

    m._liquid_set = MethodType(_liquid_set, m)
    m._vapor_set = MethodType(_vapor_set, m)

    # Add variables and define_port_members method to dummy state block
    for i in m.time:
        m.properties_out[i].flow_mass_phase = Var(m.phases)
        m.properties_out[i].mass_frac_phase_comp = Var(m.phases, m.components)

        m.properties_out[i].define_port_members = MethodType(
            define_port_members, m.properties_out[i]
        )

    # Call make_phase_split function
    make_phase_split(m, m.port)

    # Check for expected results
    for i in m.time:
        for p in m.config.property_package.phase_list:
            for j in m.config.property_package.component_list:
                assert (
                    m.port.flow_mass_phase[i, p, j] is m.e_flow_mass_phaseport[i, p, j]
                )
                assert value(m.port.flow_mass_phase[i, p, j]) == 1e-8
                assert (
                    m.port.mass_frac_phase_comp[i, p, j]
                    is m.properties_out[i].mass_frac_phase_comp[p, j]
                )


def test_make_phase_mole_frac(m):
    # Mock define_port_members method
    def define_port_members(m):
        return {
            "flow_mol_phase_comp": m.flow_mol_phase_comp,
            "mole_frac_comp": m.mole_frac_comp,
        }

    def _liquid_set(m):
        return ["L", "L1"]

    def _vapor_set(m):
        return ["V", "V1"]

    m.components = Set(initialize=["a1", "b1"])
    m.phases = Set(initialize=["L", "V"])

    m.config.property_package = DummyData()
    m.config.property_package.component_list = m.components
    m.config.property_package.phase_list = m.phases

    m._liquid_set = MethodType(_liquid_set, m)
    m._vapor_set = MethodType(_vapor_set, m)

    # Add variables and define_port_members method to dummy state block
    for i in m.time:
        m.properties_out[i].flow_mol_phase_comp = Var(m.phases, m.components)
        m.properties_out[i].mole_frac_comp = Var(m.components)

        m.properties_out[i].define_port_members = MethodType(
            define_port_members, m.properties_out[i]
        )

    # Call make_phase_split function
    make_phase_split(m, m.port, m.phases)

    # Check for expected results
    for i in m.time:
        for p in m.config.property_package.phase_list:
            for j in m.config.property_package.component_list:
                m.port.flow_mol_phase_comp[i, p, j] is m.e_flow_mol_phase_compport[
                    i, p, j
                ]
                assert value(m.port.flow_mol_phase_comp[i, p, j]) == 1e-8
                assert m.port.mole_frac_comp[i, j] is m.e_mole_frac_port[i, j]


def test_make_phase_split_mass_flow(m):
    # Mock define_port_members method
    def define_port_members(m):
        return {
            # "flow_mass": m.flow_mass,
            "mass_frac_comp": m.mass_frac_comp,
            "flow_mass_phase": m.flow_mass_phase,
            "mass_frac_phase_comp": m.mass_frac_phase_comp,
        }

    def _liquid_set(m):
        return ["a1", "b1"]

    def _vapor_set(m):
        return ["a1", "b1"]

    m.components = Set(initialize=["a1", "b1"])
    m.phases = Set(initialize=["L", "V"])

    m.config.property_package = DummyData()
    m.config.property_package.component_list = m.components
    m.config.property_package.phase_list = m.phases

    m._liquid_set = MethodType(_liquid_set, m)

    m._vapor_set = MethodType(_vapor_set, m)

    # Add variables and define_port_members method to dummy state block
    for i in m.time:
        # m.properties_out[i].flow_mass = Var()
        m.properties_out[i].mass_frac_comp = Var(m.components)
        m.properties_out[i].mass_frac_phase_comp = Var(m.phases, m.components)
        m.properties_out[i].flow_mass_phase = Var(m.phases)

        m.properties_out[i].define_port_members = MethodType(
            define_port_members, m.properties_out[i]
        )

    # Call make_phase_split function
    make_phase_split(m, m.port, m.phases, side_sf=1)

    # Check for expected results
    for i in m.time:
        for p in m.config.property_package.phase_list:
            for j in m.config.property_package.component_list:
                assert (
                    m.port.mass_frac_phase_comp[i, p, j]
                    is m._mass_frac_phase_comp_port_ref[i, p, j]
                )
                assert m.port.mass_frac_comp[i, j] is m.e_mole_frac_port[i, j]


def test_make_phase_split_mass(m):
    # Mock define_port_members method
    def define_port_members(m):
        return {
            # "flow_mass": m.flow_mass,
            "mass_frac_comp": m.mass_frac_comp,
            "flow_mass_phase_comp": m.flow_mass_phase_comp,
        }

    def _liquid_set(m):
        return ["a1", "b1"]

    def _vapor_set(m):
        return ["a1", "b1"]

    m.components = Set(initialize=["a1", "b1"])
    m.phases = Set(initialize=["L", "V"])

    m.config.property_package = DummyData()
    m.config.property_package.component_list = m.components
    m.config.property_package.phase_list = m.phases

    m._liquid_set = MethodType(_liquid_set, m)

    m._vapor_set = MethodType(_vapor_set, m)

    # Add variables and define_port_members method to dummy state block
    for i in m.time:
        # m.properties_out[i].flow_mass = Var()
        m.properties_out[i].mass_frac_comp = Var(m.components)
        m.properties_out[i].flow_mass_phase_comp = Var(m.phases, m.components)

        m.properties_out[i].define_port_members = MethodType(
            define_port_members, m.properties_out[i]
        )

    # Call make_phase_split function
    make_phase_split(m, m.port, m.phases, side_sf=1)

    # Check for expected results
    for i in m.time:
        for p in m.config.property_package.phase_list:
            for j in m.config.property_package.component_list:
                assert m.port.mass_frac_comp[i, j] is m.e_mole_frac_port[i, j]


def test_make_phase_split_enth(m):
    # Mock define_port_members method
    def define_port_members(m):
        return {
            "enth_mol_phase": m.enth_mol_phase,
            "enth_mol_phase_comp": m.enth_mol_phase_comp,
        }

    m.components = Set(initialize=["a1", "b1"])
    m.phases = Set(initialize=["L", "V"])

    m.config.property_package = DummyData()
    m.config.property_package.component_list = m.components
    m.config.property_package.phase_list = m.phases

    # Add variables and define_port_members method to dummy state block
    for i in m.time:
        m.properties_out[i].enth_mol_phase = Var(m.phases)
        m.properties_out[i].enth_mol_phase_comp = Var(m.phases, m.components)

        m.properties_out[i].define_port_members = MethodType(
            define_port_members, m.properties_out[i]
        )

    # Call make_phase_split function
    make_phase_split(m, m.port)

    # Check for expected results
    for i in m.time:
        for p in m.config.property_package.phase_list:
            assert m.port.enth_mol_phase[i, p] is m.properties_out[i].enth_mol_phase[p]
            for j in m.config.property_package.component_list:
                assert (
                    m.port.enth_mol_phase_comp[i, p, j]
                    is m.properties_out[i].enth_mol_phase_comp[p, j]
                )


def test_make_phase_split_enth_no_phase(m):
    # Mock define_port_members method
    def define_port_members(m):
        return {
            "enth_mass": m.enth_mass,
            "enth_mass_phase": m.enth_mass_phase,
        }

    # m.components = Set(initialize=["a1", "b1"])
    m.phases = Set(initialize=["L", "V"])

    m.config.property_package = DummyData()
    # m.config.property_package.component_list = m.components
    m.config.property_package.phase_list = m.phases

    # Add variables and define_port_members method to dummy state block
    for i in m.time:
        m.properties_out[i].enth_mass = Var()
        m.properties_out[i].enth_mass_phase = Var(m.phases)

        m.properties_out[i].define_port_members = MethodType(
            define_port_members, m.properties_out[i]
        )

    # Call make_phase_split function
    make_phase_split(m, m.port, m.phases)
    # Check for expected results
    for i in m.time:
        assert m.port.enth_mass[i] is m.e_enth_port[i]


def test_make_phase_non_mol_error(m):
    # Mock define_port_members method
    def define_port_members(m):
        return {
            "frac_comp": m.frac_comp,
        }

    m.components = Set(initialize=["a1", "b1"])
    m.phases = Set(initialize=["L", "V"])

    m.config.property_package = DummyData()
    m.config.property_package.component_list = m.components
    m.config.property_package.phase_list = m.phases

    # Add variables and define_port_members method to dummy state block
    for i in m.time:
        m.properties_out[i].frac_comp = Var(m.phases, m.components)

        m.properties_out[i].define_port_members = MethodType(
            define_port_members, m.properties_out[i]
        )

    # Check for expected results
    with pytest.raises(
        PropertyNotSupportedError,
        match="No mass frac or mole frac variables encountered "
        " while building ports. "
        "phase_frac as a state variable is not "
        "supported with distillation unit models.",
    ):
        make_phase_split(m, m.port)


def test_make_phase_mol_error(m):
    # Mock define_port_members method
    def define_port_members(m):
        return {
            "mole_frac_comp": m.mole_frac_comp,
        }

    m.components = Set(initialize=["a1", "b1"])
    m.phases = Set(initialize=["L", "V"])

    m.config.property_package = DummyData()
    m.config.property_package.component_list = m.components
    m.config.property_package.phase_list = m.phases

    # Add variables and define_port_members method to dummy state block
    for i in m.time:
        m.properties_out[i].mole_frac_comp = Var(m.phases, m.components)

        m.properties_out[i].define_port_members = MethodType(
            define_port_members, m.properties_out[i]
        )

    # Check for expected results
    with pytest.raises(
        PropertyNotSupportedError,
        match="No mole_frac_phase_comp or flow_mol_phase or"
        " flow_mol_phase_comp variables encountered "
        "while building ports. ",
    ):
        make_phase_split(m, m.port)


def test_make_phase_enth_error(m):
    # Mock define_port_members method
    def define_port_members(m):
        return {
            "enth_mol": m.enth_mol,
        }

    m.phases = Set(initialize=["L", "V"])

    m.config.property_package = DummyData()

    # Add variables and define_port_members method to dummy state block
    for i in m.time:
        m.properties_out[i].enth_mol = Var(m.phases)

        m.properties_out[i].define_port_members = MethodType(
            define_port_members, m.properties_out[i]
        )

    # Check for expected results
    with pytest.raises(
        PropertyPackageError,
        match="Enthalpy is indexed but the variable "
        "name does not reflect the presence of an index. "
        "Please follow the naming convention outlined "
        "in the documentation for state variables.",
    ):
        make_phase_split(m, m.port)


def test_make_phase_mass_error(m):
    # Mock define_port_members method
    def define_port_members(m):
        return {
            "mass_frac_comp": m.mass_frac_comp,
        }

    m.components = Set(initialize=["a1", "b1"])
    m.phases = Set(initialize=["L", "V"])

    m.config.property_package = DummyData()
    m.config.property_package.component_list = m.components
    m.config.property_package.phase_list = m.phases

    # Add variables and define_port_members method to dummy state block
    for i in m.time:
        m.properties_out[i].mass_frac_comp = Var(m.phases, m.components)

        m.properties_out[i].define_port_members = MethodType(
            define_port_members, m.properties_out[i]
        )

    # Check for expected results
    with pytest.raises(
        PropertyNotSupportedError,
        match="No mass_frac_phase_comp or flow_mass_phase or"
        " flow_mass_phase_comp variables encountered "
        "while building ports.",
    ):
        make_phase_split(m, m.port)


class DummyData:
    """
    Class containing necessary data to generate component list
    """

    def __init__(self):
        self.config = self.config()

    class config:
        def __init__(self):
            self.property_package = self.property_package()

        class property_package:
            def component_list(self, component_list):
                self.component_list = component_list

            def phase_list(self, phase_list):
                self.phase_list = phase_list
