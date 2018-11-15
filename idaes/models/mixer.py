##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
General purpose mixer block for IDAES models
"""
from __future__ import absolute_import  # disable implicit relative imports
from __future__ import division, print_function

import logging

from pyomo.environ import Constraint, Param, PositiveReals, Reals, Set, Var
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyutilib.enum import Enum

from idaes.core import (declare_process_block_class,
                        UnitBlockData,
                        useDefault)
from idaes.core.util.config import (is_physical_parameter_block,
                                    is_state_block,
                                    list_of_strings)
from idaes.core.util.exceptions import (BurntToast,
                                        ConfigurationError,
                                        PropertyNotSupportedError)
from idaes.core.util.math import smooth_min

__author__ = "Andrew Lee"


# Set up logger
_log = logging.getLogger(__name__)


# Enumerate options for material balances
MomentumMixingType = Enum(
    'minimize',
    'equality')


@declare_process_block_class("MixerBlock")
class MixerBlockData(UnitBlockData):
    """
    This is a general purpose model for a Mixer block with the IDAES modeling
    framework. This block can be used either as a stand-alone Mixer unit
    operation, or as a sub-model within another unit operation.

    This model creates a number of StateBlocks to represent the incoming
    streams, then writes a set of phase-component material balances, an
    overall enthalpy balance and a momentum balance (2 options) linked to a
    mixed-state StateBlock. The mixed-state StateBlock can either be specified
    by the user (allowing use as a sub-model), or created by the MixerBlock.

    When being used as a sub-model, MixerBlock should only be used when a set
    of new StateBlocks are required for the streams to be mixed. It should not
    be used to mix streams from mutiple ControlVolumes in a single unit model -
    in these cases the unit model developer should write their own mixing
    equations.
    """
    CONFIG = UnitBlockData.CONFIG()
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for mixer",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault. **Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
 and used when constructing these, **default** - None. **Valid values:** {
see property package for documentation.}"""))
    CONFIG.declare("inlet_list", ConfigValue(
        domain=list_of_strings,
        description="List of inlet names",
        doc="""A list containing names of inlets, **default** - None.
**Valid values:** {**None** - use num_inlets argument, **list** - a list of
names to use for inlets.}"""))
    CONFIG.declare("num_inlets", ConfigValue(
        domain=int,
        description="Number of inlets to unit",
        doc="""Argument indication number (int) of inlets to construct, not
used if inlet_list arg is provided, **default** - None. **Valid values:** {
**None** - use inlet_list arg instead, or default to 2 if neither argument
provided, **int** - number of inlets to creat (will be named with sequential
integers from 1 to num_inlets).}"""))
    CONFIG.declare("calculate_phase_equilibrium", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Calculate phase equilibrium in mixed stream",
        doc="""Argument indicating whether phase equilibrium should be
calculated for the resulting mixed stream, **default** - False. **Valid values:
** {**True** - calculate phase equilibrium in mixed stream, **False** -
do not calculate equilibrium in mixed stream.}"""))
    CONFIG.declare("momentum_mixing_type", ConfigValue(
        default=MomentumMixingType.minimize,
        domain=MomentumMixingType,
        description="Method to use when mxing momentum/pressure",
        doc="""Argument indicating what method to use when mixing momentum/
pressure of incoming streams, **default** - MomentumMixingType.minimize.
**Valid values:** {**MomentumMixingType.minimize** - mixed stream has
pressure equal to the minimimum pressure of the incoming streams (uses
smoothMin operator), **MomentumMixingType.equality** - enforces equality of
pressure in mixed and all incoming streams.}"""))
    CONFIG.declare("mixed_state_block", ConfigValue(
        domain=is_state_block,
        description="Existing StateBlock to use as mixed stream",
        doc="""An existing state block to use as the outlet stream from the
Mixer block, **default** - None. **Valid values:** {**None** - create a new
StateBlock for the mixed stream, **StateBlock** - a StateBock to use as the
destination for the mixed stream.}"""))
    CONFIG.declare("construct_ports", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Construct inlet and outlet Port objects",
        doc="""Argument indicating whether model should construct Port objects
linked to all inlet states and the mixed state, **default** - True.
**Valid values:** {**True** - construct Ports for all states, **False** - do
not construct Ports."""))

    def build(self):
        """
        General build method for MixerBlockData. This method calls a number
        of sub-methods which automate the construction of expected attributes
        of unit models.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        # Call super.build()
        super(MixerBlockData, self).build()

        # Call setup methods from ControlVolumeBase
        self._get_property_package()
        self._get_indexing_sets()

        # Create list of inlet names
        inlet_list = self.create_inlet_list()

        # Build StateBlocks
        inlet_blocks = self.add_inlet_state_blocks(inlet_list)

        if self.config.mixed_state_block is None:
            mixed_block = self.add_mixed_state_block()
        else:
            mixed_block = self.get_mixed_state_block()

        self.add_material_mixing_equations(inlet_blocks=inlet_blocks,
                                           mixed_block=mixed_block)
#        self.add_energy_mixing_equations(inlet_blocks=inlet_blocks,
#                                         mixed_block=mixed_block)
#
#        if self.config.momentum_mixing_type == MomentumMixingType.minimize:
#            self.add_pressure_minimization_equations(inlet_blocks=inlet_blocks,
#                                                     mixed_block=mixed_block)
#        elif self.config.momentum_mixing_type == MomentumMixingType.equality:
#            self.add_pressure_equality_equations(inlet_blocks=inlet_blocks,
#                                                 mixed_block=mixed_block)
#        else:
#            raise ConfigurationError("{} recieved unrecognised value for "
#                                     "momentum_mixing_type argument. This "
#                                     "should not occur, so please contact "
#                                     "the IDAES developers with this bug."
#                                     .format(self.name))
#
        self.add_port_objects(inlet_blocks, mixed_block)

    def create_inlet_list(self):
        """
        Create list on inlet stream names based on config arguments.

        Returns:
            list of strings
        """
        if (self.config.inlet_list is not None and
                self.config.num_inlets is not None):
            # If both arguments provided and not consistent, raise Exception
            if len(self.config.inlet_list) != self.config.num_inlets:
                raise ConfigurationError(
                        "{} MixerBlock provided with both inlet_list and "
                        "num_inlets arguments, which were not consistent ("
                        "length of inlet_list was not equal to num_inlets). "
                        "PLease check your arguments for consistency, and "
                        "note that it is only necessry to provide one of "
                        "these arguments.".format(self.name))
        elif self.config.inlet_list is None and self.config.num_inlets is None:
            # If no arguments provided for inlets, default to num_inlets = 2
            self.config.num_inlets = 2

        # Create a list of names for inlet StateBlocks
        if self.config.inlet_list is not None:
            inlet_list = self.config.inlet_list
        else:
            inlet_list = ['inlet_' + str(n)
                          for n in range(1, self.config.num_inlets+1)]

        return inlet_list

    def add_inlet_state_blocks(self, inlet_list):
        """
        Construct StateBlocks for all inlet streams.

        Args:
            list of strings to use as StateBlock names

        Returns:
            list of StateBlocks
        """
        # Setup StateBlock argument dict
        tmp_dict = self.config.property_package_args
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True

        # Create empty list to hold StateBlocks for return
        inlet_blocks = []

        # Create an instance of StateBlock for all inlets
        for i in inlet_list:
            i_obj = self._property_module.StateBlock(
                        self.time,
                        doc="Material properties at inlet",
                        default=tmp_dict)

            setattr(self, i+"_state", i_obj)

            inlet_blocks.append(getattr(self, i+"_state"))

        return inlet_blocks

    def add_mixed_state_block(self):
        """
        Constructs StateBlock to represent mixed stream.

        Returns:
            New StateBlock object
        """
        # Setup StateBlock argument dict
        tmp_dict = self.config.property_package_args
        tmp_dict["has_phase_equilibrium"] = \
            self.config.calculate_phase_equilibrium
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.mixed_state = self._property_module.StateBlock(
                                self.time,
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
        if (self.config.mixed_state_block[self.time.first()].config.parameters
                != self.config.property_package):
            raise ConfigurationError(
                    "{} StateBlock provided in mixed_state_block argument "
                    " does not come from the same property package as "
                    "provided in the property_package argument. All "
                    "StateBlocks within a MixerBlock must use the same "
                    "property package.".format(self.name))

        return self.config.mixed_state_block

    def add_material_mixing_equations(self, inlet_blocks, mixed_block):
        """
        Add material mixing equations (phase-component balances).
        """
        # Create equilibrium generation term and constraints if required
        if self.config.calculate_phase_equilibrium is True:
            # Get units from property package
            units = {}
            for u in ['holdup', 'time']:
                try:
                    units[u] = (self.config.property_package
                                .get_metadata().default_units[u])
                except KeyError:
                    units[u] = '-'

            try:
                # TODO : replace with Reference
                object.__setattr__(
                    self,
                    "phase_equilibrium_idx",
                    self.config.property_package.phase_equilibrium_idx)
            except AttributeError:
                raise PropertyNotSupportedError(
                    "{} Property package does not contain a list of phase "
                    "equilibrium reactions (phase_equilibrium_idx), thus does "
                    "not support phase equilibrium.".format(self.name))
            self.phase_equilibrium_generation = Var(
                        self.time,
                        self.phase_equilibrium_idx,
                        domain=Reals,
                        doc="Amount of generation in unit by phase "
                            "equilibria [{}/{}]"
                            .format(units['holdup'], units['time']))

        # Define terms to use in mixing equation
        def phase_equilibrium_term(b, t, p, j):
            if self.config.calculate_phase_equilibrium:
                sd = {}
                sblock = mixed_block[t]
                for r in b.phase_equilibrium_idx:
                    if sblock.phase_equilibrium_list[r][0] == j:
                        if sblock.phase_equilibrium_list[r][1][0] == p:
                            sd[r] = 1
                        elif sblock.phase_equilibrium_list[r][1][1] == p:
                            sd[r] = -1
                        else:
                            sd[r] = 0
                    else:
                        sd[r] = 0

                return sum(b.phase_equilibrium_generation[t, r]*sd[r]
                           for r in b.phase_equilibrium_idx)
            else:
                return 0

        # Get phase component list(s)
        phase_component_list = self._get_phase_comp_list()

        # Write phase-component balances
        @self.Constraint(self.time,
                         self.phase_list,
                         self.component_list,
                         doc="Material mixing equations")
        def material_mixing_equations(b, t, p, j):
            if j in phase_component_list[p]:
                return 0 == (
                        sum(inlet_blocks[i][t].get_material_flow_terms(p, j)
                            for i in range(len(inlet_blocks))) -
                        mixed_block[t].get_material_flow_terms(p, j) +
                        phase_equilibrium_term(b, t, p, j))
            else:
                return Constraint.Skip

    def add_energy_mixing_equations(self, inlet_blocks, mixed_block):
        """
        Add energy mixing equations (total enthalpy balance).
        """
        @self.Constraint(self.time, doc="Energy balances")
        def enthalpy_mixing_equations(b, t):
            return 0 == (sum(sum(inlet_blocks[i][t].get_enthalpy_flow_terms(p)
                                 for p in b.phase_list)
                             for i in range(len(inlet_blocks))) -
                         sum(mixed_block[t].get_enthalpy_flow_terms(p)
                             for p in b.phase_list))

    def add_pressure_minimization_equations(self, inlet_blocks, mixed_block):
        """
        Add pressure minimization equations. This is done by sequential
        comparisons of each inlet to the minimum pressure so far, using
       the IDAES smooth minimum fuction.
        """
        # Add variables
        self.inlet_idx = Set(initialize=range(1, len(inlet_blocks)+1),
                             ordered=True)

        self.minimum_pressure = Var(self.time,
                                    self.inlet_idx,
                                    doc='Variable for calculating '
                                        'minimum inlet pressure')

        self.eps_pressure = Param(mutable=True,
                                  initialize=1e-3,
                                  domain=PositiveReals,
                                  doc='Smoothing term for '
                                      'minimum inlet pressure')

        # Calculate minimum inlet pressure
        @self.Constraint(self.time,
                         self.inlet_idx,
                         doc='Calculation for minimum inlet pressure')
        def minimum_pressure_constraint(b, t, i):
            if i == self.inlet_idx.first():
                return self.minimum_pressure[t, i] == (
                           inlet_blocks[i-1][t].pressure)
            else:
                return self.minimum_pressure[t, i] == (
                        smooth_min(self.minimum_pressure[t, i-1],
                                   inlet_blocks[i-1][t].pressure,
                                   self.eps_pressure))

        # Set inlet pressure to minimum pressure
        @self.Constraint(self.time, doc='link pressure to holdup')
        def mixture_pressure(b, t):
            return mixed_block[t].pressure == (
                       self.minimum_pressure[t,
                                             self.inlet_idx.last()])

    def add_port_objects(self, inlet_list, inlet_blocks, mixed_block):
        """
        Adds Port objects if required.

        Args:
            a list of inlet StateBlock objects
            a mixed state StateBlock object

        Returns:
            None
        """
        if self.config.construct_ports is True:
            # Add ports
            for p in inlet_list:
                i_state = getattr(self, p+"_state")
                self.add_port(name=p, block=i_state, doc="Inlet Port")
            self.add_port(name="outlet", block=mixed_block, doc="Outlet Port")
