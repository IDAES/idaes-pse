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
General purpose separator block for IDAES models
"""
from __future__ import absolute_import  # disable implicit relative imports
from __future__ import division, print_function

import logging

from pyomo.environ import (Constraint,
                           Expression,
                           Reference,
                           Set,
                           SolverFactory,
                           TerminationCondition,
                           Var,
                           value)
from pyomo.network import Port
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyutilib.enum import Enum

from idaes.core import (declare_process_block_class,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import (is_physical_parameter_block,
                                    is_state_block,
                                    list_of_strings)
from idaes.core.util.exceptions import (BurntToast,
                                        ConfigurationError)

__author__ = "Andrew Lee"


# Set up logger
_log = logging.getLogger(__name__)


# Enumerate options for balances
SplittingType = Enum(
    'totalFlow',
    'phaseFlow',
    'componentFlow',
    'phaseComponentFlow')

EnergySplittingType = Enum(
    'equal_temperature',
    'equal_molar_enthalpy')


@declare_process_block_class("Separator")
class SeparatorData(UnitModelBlockData):
    """
    This is a general purpose model for a Separator block with the IDAES
    modeling framework. This block can be used either as a stand-alone
    Separator unit operation, or as a sub-model within another unit operation.

    This model creates a number of StateBlocks to represent the outgoing
    streams, then writes a set of phase-component material balances, an
    overall enthalpy balance (2 options), and a momentum balance (2 options)
    linked to a mixed-state StateBlock. The mixed-state StateBlock can either
    be specified by the user (allowing use as a sub-model), or created by the
    Separator.

    When being used as a sub-model, Separator should only be used when a
    set of new StateBlocks are required for the streams to be separated. It
    should not be used to separate streams to go to mutiple ControlVolumes in a
    single unit model - in these cases the unit model developer should write
    their own splitting equations.
    """
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for mixer",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))
    CONFIG.declare("outlet_list", ConfigValue(
        domain=list_of_strings,
        description="List of outlet names",
        doc="""A list containing names of outlets,
**default** - None.
**Valid values:** {
**None** - use num_outlets argument,
**list** - a list of names to use for outlets.}"""))
    CONFIG.declare("num_outlets", ConfigValue(
        domain=int,
        description="Number of outlets to unit",
        doc="""Argument indicating number (int) of outlets to construct, not
used if outlet_list arg is provided,
**default** - None.
**Valid values:** {
**None** - use outlet_list arg instead, or default to 2 if neither argument
provided,
**int** - number of outlets to create (will be named with sequential integers
from 1 to num_outlets).}"""))
    CONFIG.declare("split_basis", ConfigValue(
        default=SplittingType.totalFlow,
        domain=SplittingType,
        description="Basis for splitting stream",
        doc="""Argument indicating basis to use for splitting mixed stream,
**default** - SplittingType.totalFlow.
**Valid values:** {
**SplittingType.totalFlow** - split based on total flow (split
fraction indexed only by time and outlet),
**SplittingType.phaseFlow** - split based on phase flows (split fraction
indexed by time, outlet and phase),
**SplittingType.componentFlow** - split based on component flows (split
fraction indexed by time, outlet and components),
**SplittingType.phaseComponentFlow** - split based on phase-component flows (
split fraction indexed by both time, outlet, phase and components).}"""))
    CONFIG.declare("energy_split_basis", ConfigValue(
        default=EnergySplittingType.equal_temperature,
        domain=EnergySplittingType,
        description="Type of constraint to write for energy splitting",
        doc="""Argument indicating basis to use for splitting energy,
**default** - EnergySplittingType.equal_temperature.
**Valid values:** {
**EnergySplittingType.equal_temperature** - outlet temperatures equal inlet
**EnergySplittingType.equal_molar_enthalpy** - oulet molar enthalpies equal
inlet}"""))
    CONFIG.declare("ideal_separation", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Ideal splitting flag",
        doc="""Argument indicating whether ideal splitting should be used.
Ideal splitting assumes perfect spearation of material, and attempts to
avoid duplication of StateBlocks by directly partitioning outlet flows to
ports,
**default** - True.
**Valid values:** {
**True** - use ideal splitting methods,
**False** - use explicit splitting equations with split fractions.}"""))
    CONFIG.declare("ideal_split_map", ConfigValue(
        domain=dict,
        description="Ideal splitting partitioning map",
        doc="""Dictionary containing information on how extensive variables
should be partitioned when using ideal splitting (ideal_separation = True).
**default** - None.
**Valid values:** {
**dict** with keys of indexing set members and values indicating which outlet
this combination of keys should be partitioned to.
E.g. {("Vap", "H2"): "outlet_1"}}"""))
    CONFIG.declare("mixed_state_block", ConfigValue(
        domain=is_state_block,
        description="Existing StateBlock to use as mixed stream",
        doc="""An existing state block to use as the source stream from the
Separator block,
**default** - None.
**Valid values:** {
**None** - create a new StateBlock for the mixed stream,
**StateBlock** - a StateBock to use as the source for the mixed stream.}"""))
    CONFIG.declare("construct_ports", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Construct inlet and outlet Port objects",
        doc="""Argument indicating whether model should construct Port objects
linked the mixed state and all outlet states,
**default** - True.
**Valid values:** {
**True** - construct Ports for all states,
**False** - do not construct Ports."""))

    def build(self):
        """
        General build method for SeparatorData. This method calls a number
        of sub-methods which automate the construction of expected attributes
        of unit models.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        # Call super.build()
        super(SeparatorData, self).build()

        # Call setup methods from ControlVolumeBlockData
        self._get_property_package()
        self._get_indexing_sets()

        # Create list of inlet names
        outlet_list = self.create_outlet_list()

        if self.config.mixed_state_block is None:
            mixed_block = self.add_mixed_state_block()
        else:
            mixed_block = self.get_mixed_state_block()

        # Add inlet port
        self.add_inlet_port_objects(mixed_block)

        # Construct splitter based on ideal_separation argument
        if self.config.ideal_separation:
            # Use ideal partitioning method
            self.partition_outlet_flows(mixed_block, outlet_list)
        else:
            # Otherwise, Build StateBlocks for outlet
            outlet_blocks = self.add_outlet_state_blocks(outlet_list)

            # Add split fractions
            self.add_split_fractions(outlet_list)

            # Construct splitting equations
            self.add_material_splitting_constraints(mixed_block)
            self.add_energy_splitting_constraints(mixed_block)
            self.add_momentum_splitting_constraints(mixed_block)

            # Construct outlet port objects
            self.add_outlet_port_objects(outlet_list, outlet_blocks)

    def create_outlet_list(self):
        """
        Create list of outlet stream names based on config arguments.

        Returns:
            list of strings
        """
        if (self.config.outlet_list is not None and
                self.config.num_outlets is not None):
            # If both arguments provided and not consistent, raise Exception
            if len(self.config.outlet_list) != self.config.num_outlets:
                raise ConfigurationError(
                        "{} Separator provided with both outlet_list and "
                        "num_outlets arguments, which were not consistent ("
                        "length of outlet_list was not equal to num_outlets). "
                        "PLease check your arguments for consistency, and "
                        "note that it is only necessry to provide one of "
                        "these arguments.".format(self.name))
        elif (self.config.outlet_list is None and
              self.config.num_outlets is None):
            # If no arguments provided for outlets, default to num_outlets = 2
            self.config.num_outlets = 2

        # Create a list of names for outlet StateBlocks
        if self.config.outlet_list is not None:
            outlet_list = self.config.outlet_list
        else:
            outlet_list = ['outlet_' + str(n)
                           for n in range(1, self.config.num_outlets+1)]

        return outlet_list

    def add_outlet_state_blocks(self, outlet_list):
        """
        Construct StateBlocks for all outlet streams.

        Args:
            list of strings to use as StateBlock names

        Returns:
            list of StateBlocks
        """
        # Setup StateBlock argument dict
        tmp_dict = self.config.property_package_args
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        # Create empty list to hold StateBlocks for return
        outlet_blocks = []

        # Create an instance of StateBlock for all outlets
        for o in outlet_list:
            o_obj = self.config.property_package.state_block_class(
                        self.time_ref,
                        doc="Material properties at outlet",
                        default=tmp_dict)

            setattr(self, o+"_state", o_obj)

            outlet_blocks.append(getattr(self, o+"_state"))

        return outlet_blocks

    def add_mixed_state_block(self):
        """
        Constructs StateBlock to represent mixed stream.

        Returns:
            New StateBlock object
        """
        # Setup StateBlock argument dict
        tmp_dict = self.config.property_package_args
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True

        self.mixed_state = self.config.property_package.state_block_class(
                                self.time_ref,
                                doc="Material properties of mixed stream",
                                default=tmp_dict)

        return self.mixed_state

    def get_mixed_state_block(self):
        """
        Validates StateBlock provided in user arguments for mixed stream.

        Returns:
            The user-provided StateBlock or an Exception
        """
        # Sanity check to make sure method is not called when arg missing
        if self.config.mixed_state_block is None:
            raise BurntToast("{} get_mixed_state_block method called when "
                             "mixed_state_block argument is None. This should "
                             "not happen.".format(self.name))

        # Check that the user-provided StateBlock uses the same prop pack
        if (self.config.mixed_state_block[
                    self.time_ref.first()].config.parameters
                != self.config.property_package):
            raise ConfigurationError(
                    "{} StateBlock provided in mixed_state_block argument "
                    " does not come from the same property package as "
                    "provided in the property_package argument. All "
                    "StateBlocks within a Separator must use the same "
                    "property package.".format(self.name))

        return self.config.mixed_state_block

    def add_inlet_port_objects(self, mixed_block):
        """
        Adds inlet Port object if required.

        Args:
            a mixed state StateBlock object

        Returns:
            None
        """
        if self.config.construct_ports is True:
            self.add_port(name="inlet", block=mixed_block, doc="Inlet Port")

    def add_outlet_port_objects(self, outlet_list, outlet_blocks):
        """
        Adds outlet Port objects if required.

        Args:
            a list of outlet StateBlock objects

        Returns:
            None
        """
        if self.config.construct_ports is True:
            # Add ports
            for p in outlet_list:
                o_state = getattr(self, p+"_state")
                self.add_port(name=p, block=o_state, doc="Outlet Port")

    def add_split_fractions(self, outlet_list):
        """
        Creates outlet Port objects and tries to partiton mixed stream flows
        between these

        Args:
            StateBlock representing the mixed flow to be split
            a list of names for outlets

        Returns:
            None
        """
        self.outlet_idx = Set(initialize=outlet_list)

        if self.config.split_basis == SplittingType.totalFlow:
            sf_idx = [self.time_ref, self.outlet_idx]
            sf_sum_idx = [self.time_ref]
        elif self.config.split_basis == SplittingType.phaseFlow:
            sf_idx = [self.time_ref, self.outlet_idx, self.phase_list_ref]
            sf_sum_idx = [self.time_ref, self.phase_list_ref]
        elif self.config.split_basis == SplittingType.componentFlow:
            sf_idx = [self.time_ref, self.outlet_idx, self.component_list_ref]
            sf_sum_idx = [self.time_ref, self.component_list_ref]
        elif self.config.split_basis == SplittingType.phaseComponentFlow:
            sf_idx = [self.time_ref,
                      self.outlet_idx,
                      self.phase_list_ref,
                      self.component_list_ref]
            sf_sum_idx = [self.time_ref,
                          self.phase_list_ref,
                          self.component_list_ref]
        else:
            raise BurntToast("{} split_basis has unexpected value. This "
                             "should not happen.".format(self.name))

        # Create split fraction variable
        self.split_fraction = Var(*sf_idx,
                                  initialize=0.5,
                                  doc="Outlet split fractions")

        # Add constraint that split fractions sum to 1
        def sum_sf_rule(b, t, *args):
            return 1 == sum(b.split_fraction[t, o, args]
                            for o in self.outlet_idx)
        self.sum_split_frac = Constraint(*sf_sum_idx, rule = sum_sf_rule)

    def add_material_splitting_constraints(self, mixed_block):
        """
        Creates constraints for splitting the material flows
        """
        def sf(t, o, p, j):
            if self.config.split_basis == SplittingType.totalFlow:
                return self.split_fraction[t, o]
            elif self.config.split_basis == SplittingType.phaseFlow:
                return self.split_fraction[t, o, p]
            elif self.config.split_basis == SplittingType.componentFlow:
                return self.split_fraction[t, o, j]
            elif self.config.split_basis == SplittingType.phaseComponentFlow:
                return self.split_fraction[t, o, p, j]

        @self.Constraint(self.time_ref,
                         self.outlet_idx,
                         self.phase_list_ref,
                         self.component_list_ref,
                         doc="Material splitting equations")
        def material_splitting_eqn(b, t, o, p, j):
            o_block = getattr(self, o+"_state")
            return (sf(t, o, p, j) *
                    mixed_block[t].get_material_flow_terms(p, j) ==
                    o_block[t].get_material_flow_terms(p, j))

    def add_energy_splitting_constraints(self, mixed_block):
        """
        Creates constraints for splitting the energy flows - done by equating
        temperatures in outlets.
        """
        if self.config.energy_split_basis == \
            EnergySplittingType.equal_temperature:
            @self.Constraint(self.time_ref,
                             self.outlet_idx,
                             doc="Temperature equality constraint")
            def temperature_equality_eqn(b, t, o):
                o_block = getattr(self, o+"_state")
                return mixed_block[t].temperature == o_block[t].temperature
        elif self.config.energy_split_basis == \
            EnergySplittingType.equal_molar_enthalpy:
            @self.Constraint(self.time_ref,
                             self.outlet_idx,
                             doc="Molar enthalpy equality constraint")
            def temperature_equality_eqn(b, t, o):
                o_block = getattr(self, o+"_state")
                return mixed_block[t].enth_mol == o_block[t].enth_mol

    def add_momentum_splitting_constraints(self, mixed_block):
        """
        Creates constraints for splitting the momentum flows - done by equating
        pressures in outlets.
        """
        @self.Constraint(self.time_ref,
                         self.outlet_idx,
                         doc="Pressure equality constraint")
        def pressure_equality_eqn(b, t, o):
            o_block = getattr(self, o+"_state")
            return mixed_block[t].pressure == o_block[t].pressure

    def partition_outlet_flows(self, mb, outlet_list):
        """
        Creates outlet Port objects and tries to partiton mixed stream flows
        between these

        Args:
            StateBlock representing the mixed flow to be split
            a list of names for outlets

        Returns:
            None
        """
        # Check arguments
        if self.config.construct_ports is False:
            raise ConfigurationError("{} cannot have and ideal separator "
                                     "(ideal_separation = True) with "
                                     "construct_ports = False."
                                     .format(self.name))
        if self.config.split_basis == SplittingType.totalFlow:
            raise ConfigurationError("{} cannot do an ideal separation based "
                                     "on total flow.".format(self.name))
        if self.config.ideal_split_map is None:
            raise ConfigurationError("{} was not provided with an "
                                     "ideal_split_map argument which is "
                                     "necessary for doing an ideal_separation."
                                     .format(self.name))

        # Validate split map
        split_map = self.config.ideal_split_map
        idx_list = []
        if self.config.split_basis == SplittingType.phaseFlow:
            for p in self.phase_list_ref:
                idx_list.append((p))

            if len(idx_list) != len(split_map):
                raise ConfigurationError(
                        "{} ideal_split_map does not match with "
                        "split_basis chosen. ideal_split_map must"
                        " have a key for each combination of indices."
                        .format(self.name))
            for k in idx_list:
                if k not in split_map:
                    raise ConfigurationError(
                        "{} ideal_split_map does not match with "
                        "split_basis chosen. ideal_split_map must"
                        " have a key for each combination of indices."
                        .format(self.name))

        elif self.config.split_basis == SplittingType.componentFlow:
            for j in self.component_list_ref:
                idx_list.append((j))

            if len(idx_list) != len(split_map):
                raise ConfigurationError(
                        "{} ideal_split_map does not match with "
                        "split_basis chosen. ideal_split_map must"
                        " have a key for each component."
                        .format(self.name))
        elif self.config.split_basis == SplittingType.phaseComponentFlow:
            for p in self.phase_list_ref:
                for j in self.component_list_ref:
                    idx_list.append((p, j))

            if len(idx_list) != len(split_map):
                raise ConfigurationError(
                        "{} ideal_split_map does not match with "
                        "split_basis chosen. ideal_split_map must"
                        " have a key for each phase-component pair."
                        .format(self.name))

        # Check that no. outlets matches split_basis
        if len(outlet_list) != len(idx_list):
            raise ConfigurationError(
                        "{} Cannot perform ideal separation. Must have one "
                        "outlet for each possible combination of the "
                        "chosen split_basis."
                        .format(self.name))

        # Get list of port members
        s_vars = mb[self.time_ref.first()].define_port_members()

        # Add empty Port objects
        for o in outlet_list:
            p_obj = Port(noruleinit=True,
                         doc="Outlet Port")
            setattr(self, o, p_obj)

            # Iterate over members to create References or Expressions
            for s in s_vars:
                # Get local variable name of component
                l_name = s_vars[s].local_name

                if l_name == "pressure" or l_name == "temperature":
                    # Assume outlets same as mixed flow - make Reference
                    e_obj = Reference(mb[:].component(l_name))

                elif (l_name.startswith("mole_frac") or
                      l_name.startswith("mass_frac")):
                    # Mole and mass frac need special handling
                    if l_name.endswith("_phase"):
                        def e_rule(b, t, p, j):
                            if self.config.split_basis == \
                                        SplittingType.phaseFlow:
                                s_check = split_map[p]
                            elif self.config.split_basis == \
                                    SplittingType.componentFlow:
                                s_check = split_map[j]
                            elif self.config.split_basis == \
                                    SplittingType.phaseComponentFlow:
                                s_check = split_map[p, j]
                            else:
                                raise BurntToast(
                                        "{} This should not happen. Please "
                                        "report this bug to the IDAES "
                                        "developers.".format(self.name))

                            if s_check == o:
                                return mb[t].component(l_name)[p, j]
                            else:
                                return 0

                        e_obj = Expression(self.time_ref,
                                           self.phase_list_ref,
                                           self.component_list_ref,
                                           rule=e_rule)

                    else:
                        if self.config.split_basis == \
                                    SplittingType.componentFlow:
                            def e_rule(b, t, j):
                                if split_map[j] == o:
                                    return mb[t].component(l_name)[j]
                                # else:
                                return 0

                        else:
                            def e_rule(b, t, j):
                                try:
                                    mfp = mb[t].component(l_name+"_phase")
                                except AttributeError:
                                    raise AttributeError(
                                        "{} Cannot use ideal splitting with "
                                        "this property package. Package uses "
                                        "indexed port member {} which does not"
                                        " have the correct indexing sets, and "
                                        "an equivalent variable with correct "
                                        "indexing sets is not available."
                                        .format(self.name, s))

                                for p in self.phase_list_ref:
                                    if self.config.split_basis == \
                                            SplittingType.phaseFlow:
                                        s_check = split_map[p]
                                    elif self.config.split_basis == \
                                            SplittingType.phaseComponentFlow:
                                        s_check = split_map[p, j]
                                    else:
                                        raise BurntToast(
                                            "{} This should not happen. Please"
                                            " report this bug to the IDAES "
                                            "developers.".format(self.name))

                                    if s_check == o:
                                        return mfp[p, j]
                                # else:
                                return 0

                        e_obj = Expression(self.time_ref,
                                           self.component_list_ref,
                                           rule=e_rule)

                elif l_name.endswith("_phase_comp"):
                    def e_rule(b, t, p, j):
                        if self.config.split_basis == \
                                SplittingType.phaseFlow:
                            s_check = split_map[p]
                        elif self.config.split_basis == \
                                SplittingType.componentFlow:
                            s_check = split_map[j]
                        elif self.config.split_basis == \
                                SplittingType.phaseComponentFlow:
                            s_check = split_map[p, j]
                        else:
                            raise BurntToast(
                                    "{} This should not happen. Please"
                                    " report this bug to the IDAES "
                                    "developers.".format(self.name))

                        if s_check == o:
                            return mb[t].component(l_name)[p, j]
                        else:
                            return 0

                    e_obj = Expression(self.time_ref,
                                       self.phase_list_ref,
                                       self.component_list_ref,
                                       rule=e_rule)

                elif l_name.endswith("_phase"):
                    if self.config.split_basis == \
                                    SplittingType.phaseFlow:
                        def e_rule(b, t, p):
                            if split_map[p] == o:
                                return mb[t].component(l_name)[p]
                            else:
                                return 0

                    else:
                        def e_rule(b, t, p):
                            try:
                                mfp = mb[t].component(l_name+"_comp")
                            except AttributeError:
                                raise AttributeError(
                                    "{} Cannot use ideal splitting with this "
                                    "property package. Package uses indexed "
                                    "port member {} which does not have the "
                                    "correct indexing sets, and an equivalent "
                                    "variable with correct indexing sets is "
                                    "not available."
                                    .format(self.name, s))

                            for j in self.component_list_ref:
                                if self.config.split_basis == \
                                        SplittingType.componentFlow:
                                    s_check = split_map[j]
                                elif self.config.split_basis == \
                                        SplittingType.phaseComponentFlow:
                                    s_check = split_map[p, j]
                                else:
                                    raise BurntToast(
                                        "{} This should not happen. Please"
                                        " report this bug to the IDAES "
                                        "developers.".format(self.name))

                                if s_check == o:
                                    return mfp[p, j]
                            # else:
                            return 0

                    e_obj = Expression(self.time_ref,
                                       self.phase_list_ref,
                                       rule=e_rule)

                elif l_name.endswith("_comp"):
                    if self.config.split_basis == \
                            SplittingType.componentFlow:
                        def e_rule(b, t, j):
                            if split_map[j] == o:
                                return mb[t].component(l_name)[j]
                            else:
                                return 0

                    elif self.config.split_basis == \
                            SplittingType.phaseFlow:

                        def e_rule(b, t, j):
                            try:
                                mfp = mb[t].component(
                                        "{0}_phase{1}"
                                        .format(l_name[:-5], s[-5:]))
                            except AttributeError:
                                raise AttributeError(
                                    "{} Cannot use ideal splitting with this "
                                    "property package. Package uses indexed "
                                    "port member {} which does not have the "
                                    "correct indexing sets, and an equivalent "
                                    "variable with correct indexing sets is "
                                    "not available."
                                    .format(self.name, s))

                            for p in self.phase_list_ref:
                                if self.config.split_basis == \
                                        SplittingType.phaseFlow:
                                    s_check = split_map[p]
                                elif self.config.split_basis == \
                                        SplittingType.phaseComponentFlow:
                                    s_check = split_map[p, j]
                                else:
                                    raise BurntToast(
                                        "{} This should not happen. Please"
                                        " report this bug to the IDAES "
                                        "developers.".format(self.name))

                                if s_check == o:
                                    return mfp[p, j]
                            # else:
                            return 0

                    e_obj = Expression(self.time_ref,
                                       self.component_list_ref,
                                       rule=e_rule)

                else:
                    # Not a recognised state, check for indexing sets
                    if mb[self.time_ref.first()].component(
                            l_name).is_indexed():
                        # Is indexed, assume indexes match and partition

                        def e_rule(b, t, k):
                            if split_map[k] == o:
                                try:
                                    return mb[t].component(l_name)[k]
                                except KeyError:
                                    raise KeyError(
                                        "{} Cannot use ideal splitting with"
                                        " this property package. Package uses "
                                        "indexed port member {} which does not"
                                        " have suitable indexing set(s)."
                                        .format(self.name, s))
                            else:
                                return 0

                        # TODO : Reusing indexing set from first port member.
                        # TODO : Not sure how good of an idea this is.
                        e_obj = Expression(
                                    self.time_ref,
                                    mb[self.time_ref.first()]
                                        .component(l_name).index_set(),
                                    rule=e_rule)

                    else:
                        # Is not indexed, look for indexed equivalent
                        try:
                            if self.config.split_basis == \
                                    SplittingType.phaseFlow:
                                def e_rule(b, t):
                                    for p in self.phase_list_ref:
                                        if split_map[p] == o:
                                            return mb[t].component(
                                                    l_name+"_phase")[p]
                                    # else
                                    return 0

                            elif self.config.split_basis == \
                                    SplittingType.componentFlow:
                                def e_rule(b, t):
                                    for j in self.component_list_ref:
                                        if split_map[j] == o:
                                            return mb[t].component(
                                                    l_name+"_comp")[j]
                                    # else
                                    return 0

                            elif self.config.split_basis == \
                                    SplittingType.phaseComponentFlow:
                                def e_rule(b, t):
                                    for p in self.phase_list_ref:
                                        for j in self.component_list_ref:
                                            if split_map[p, j] == o:
                                                return (mb[t].component(
                                                        l_name+"_phase_comp")
                                                        [p, j])
                                    # else
                                    return 0

                        except AttributeError:
                            raise AttributeError(
                                "{} Cannot use ideal splitting with this "
                                "property package. Package uses unindexed "
                                "port member {} which does not have an "
                                "equivalent indexed form."
                                .format(self.name, s))

                    e_obj = Expression(self.time_ref,
                                       rule=e_rule)

                # Add Reference/Expression object to Separator model object
                setattr(self, "_"+o+"_"+l_name+"_ref", e_obj)

                # Add member to Port object
                p_obj.add(e_obj, s)

    def model_check(blk):
        """
        This method executes the model_check methods on the associated state
        blocks (if they exist). This method is generally called by a unit model
        as part of the unit's model_check method.

        Args:
            None

        Returns:
            None
        """
        # Try property block model check
        for t in blk.time_ref:
            try:
                if blk.config.mixed_state_block is None:
                    blk.mixed_state[t].model_check()
                else:
                    blk.config.mixed_state_block.model_check()
            except AttributeError:
                _log.warning('{} Separator inlet state block has no '
                             'model check. To correct this, add a '
                             'model_check method to the associated '
                             'StateBlock class.'.format(blk.name))

            try:
                outlet_list = blk.create_outlet_list()
                for o in outlet_list:
                    o_block = getattr(blk, o+"_state")
                    o_block[t].model_check()
            except AttributeError:
                _log.warning('{} Separator outlet state block has no '
                             'model checks. To correct this, add a model_check'
                             ' method to the associated StateBlock class.'
                             .format(blk.name))

    def initialize(blk, outlvl=0, optarg=None,
                   solver='ipopt', hold_state=True):
        '''
        Initialisation routine for separator (default solver ipopt)

        Keyword Arguments:
            outlvl : sets output level of initialisation routine. **Valid
                     values:** **0** - no output (default), **1** - return
                     solver state for each step in routine, **2** - include
                     solver output infomation (tee=True)
            optarg : solver options dictionary object (default=None)
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')
            hold_state : flag indicating whether the initialization routine
                     should unfix any state variables fixed during
                     initialization, **default** - True. **Valid values:**
                     **True** - states variables are not unfixed, and a dict of
                     returned containing flags for which states were fixed
                     during initialization, **False** - state variables are
                     unfixed after initialization by calling the release_state
                     method.

        Returns:
            If hold_states is True, returns a dict containing flags for which
            states were fixed during initialization.
        '''
        # Set solver options
        if outlvl > 1:
            stee = True
        else:
            stee = False

        opt = SolverFactory(solver)
        opt.options = optarg

        # Initialize mixed state block
        if blk.config.mixed_state_block is not None:
            mblock = blk.config.mixed_state_block
        else:
            mblock = blk.mixed_state
        flags = mblock.initialize(outlvl=outlvl-1,
                                  optarg=optarg,
                                  solver=solver,
                                  hold_state=hold_state)

        if blk.config.ideal_separation:
            # If using ideal splitting, initialisation should be complete
            return flags

        # Initialize outlet StateBlocks
        outlet_list = blk.create_outlet_list()

        for o in outlet_list:
            # Get corresponding outlet StateBlock
            o_block = getattr(blk, o+"_state")

            for t in blk.time_ref:

                # Calculate values for state variables
                s_vars = o_block[t].define_state_vars()

                for v in s_vars:
                    m_var = getattr(mblock[t], s_vars[v].local_name)

                    if "flow" in v:
                        # If a "flow" variable, is extensive
                        # Apply split fraction
                        try:
                            for k in s_vars[v]:
                                if k is None:
                                    s_vars[v][k].value = value(
                                    m_var[k]*blk.split_fraction[(t, o)])
                                else:
                                    s_vars[v][k].value = value(
                                            m_var[k]*blk.split_fraction[
                                                    (t, o) + k])
                        except KeyError:
                            raise KeyError(
                                    "{} state variable and split fraction "
                                    "indexing sets do not match. The in-built"
                                    " initialization routine for Separators "
                                    "relies on the split fraction and state "
                                    "variable indexing sets matching to "
                                    "calculate initial guesses for extensive "
                                    "variables. In other cases users will "
                                    "need to provide their own initial "
                                    "guesses".format(blk.name))
                    else:
                        # Otherwise intensive, equate to mixed stream
                        for k in s_vars[v]:
                            s_vars[v][k].value = m_var[k].value

                # Call initialization routine for outlet StateBlock
                o_block.initialize(outlvl=outlvl-1,
                                   optarg=optarg,
                                   solver=solver,
                                   hold_state=False)

        if blk.config.mixed_state_block is None:
            results = opt.solve(blk, tee=stee)

            if outlvl > 0:
                if results.solver.termination_condition == \
                        TerminationCondition.optimal:
                    _log.info('{} Initialisation Complete.'.format(blk.name))
                else:
                    _log.warning('{} Initialisation Failed.'.format(blk.name))
        else:
            _log.info('{} Initialisation Complete.'.format(blk.name))

        return flags

    def release_state(blk, flags, outlvl=0):
        '''
        Method to release state variables fixed during initialisation.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state = True.
            outlvl : sets output level of logging

        Returns:
            None
        '''
        if blk.config.mixed_state_block is None:
            mblock = blk.mixed_state
        else:
            mblock = blk.config.mixed_state_block

        mblock.release_state(flags, outlvl=outlvl-1)
