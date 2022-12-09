#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
General purpose mixer block for IDAES models
"""
from pyomo.environ import Param, PositiveReals, value
from pyomo.common.config import ConfigBlock, ConfigValue, In, ListOf

from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError

from idaes.core.util.math import smooth_min
from idaes.models.unit_models import MomentumMixingType
from idaes.core.util import from_json, to_json, StoreSpec

import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import idaes.logger as idaeslog


__author__ = "John Eslick"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("HelmMixer")
class HelmMixerData(UnitModelBlockData):
    """
    This is a Helmholtz EOS specific mixed unit model.
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
**default** = False. Mixer blocks are always steady-state.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Mixer blocks do not contain holdup, thus this must be
False.""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for mixer",
            doc="""Property parameter object used to define property
calculations, **default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property
block(s) and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "inlet_list",
        ConfigValue(
            domain=ListOf(str),
            description="List of inlet names",
            doc="""A list containing names of inlets,
**default** - None.
**Valid values:** {
**None** - use num_inlets argument,
**list** - a list of names to use for inlets.}""",
        ),
    )
    CONFIG.declare(
        "num_inlets",
        ConfigValue(
            domain=int,
            description="Number of inlets to unit",
            doc="""Argument indicating number (int) of inlets to construct, not
used if inlet_list arg is provided,
**default** - None.
**Valid values:** {
**None** - use inlet_list arg instead, or default to 2 if neither argument
provided,
**int** - number of inlets to create (will be named with sequential integers
from 1 to num_inlets).}""",
        ),
    )
    CONFIG.declare(
        "momentum_mixing_type",
        ConfigValue(
            default=MomentumMixingType.minimize,
            domain=MomentumMixingType,
            description="Method to use when mixing momentum/pressure",
            doc="""Argument indicating what method to use when mixing momentum/
pressure of incoming streams,
**default** - MomentumMixingType.minimize.
**Valid values:** {
**MomentumMixingType.none** - do not include momentum mixing equations,
**MomentumMixingType.minimize** - mixed stream has pressure equal to the
minimimum pressure of the incoming streams (uses smoothMin operator),
**MomentumMixingType.equality** - enforces equality of pressure in mixed and
all incoming streams.,
**MomentumMixingType.minimize_and_equality** - add constraints for pressure
equal to the minimum pressure of the inlets and constraints for equality of
pressure in mixed and all incoming streams. When the model is initially built,
the equality constraints are deactivated.  This option is useful for switching
between flow and pressure driven simulations.}""",
        ),
    )

    def build(self):
        """
        General build method for MixerData. This method calls a number
        of sub-methods which automate the construction of expected attributes
        of unit models.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        # Call super.build()
        super().build()

        self._get_property_package()

        # Create list of inlet names
        # PYLINT-TODO: assigning the result of self.create_inlet_list() to unused local variable inlet_list
        # causes pylint error assignment-from-no-return; check if removing assignment is OK
        self.create_inlet_list()

        # Build StateBlocks
        self.add_inlet_state_blocks()
        self.add_mixed_state_block()

        @self.Constraint(self.flowsheet().time)
        def mass_balance(b, t):
            return self.mixed_state[t].flow_mol == sum(
                self.inlet_blocks[i][t].flow_mol for i in self.inlet_list
            )

        @self.Constraint(self.flowsheet().time)
        def energy_balance(b, t):
            return self.mixed_state[t].enth_mol * self.mixed_state[t].flow_mol == sum(
                self.inlet_blocks[i][t].enth_mol * self.inlet_blocks[i][t].flow_mol
                for i in self.inlet_list
            )

        mmx_type = self.config.momentum_mixing_type
        if mmx_type == MomentumMixingType.minimize:
            self.add_pressure_minimization_equations()
        elif mmx_type == MomentumMixingType.equality:
            self.add_pressure_equality_equations()
        elif mmx_type == MomentumMixingType.minimize_and_equality:
            self.add_pressure_minimization_equations()
            self.add_pressure_equality_equations()
            self.use_minimum_inlet_pressure_constraint()

        self.add_port_objects()

    def create_inlet_list(self):
        """
        Create list of inlet stream names based on config arguments.

        Returns:
            list of strings
        """
        if self.config.inlet_list is not None and self.config.num_inlets is not None:
            # If both arguments provided and not consistent, raise Exception
            if len(self.config.inlet_list) != self.config.num_inlets:
                raise ConfigurationError(
                    "{} Mixer provided with both inlet_list and "
                    "num_inlets arguments, which were not consistent ("
                    "length of inlet_list was not equal to num_inlets). "
                    "PLease check your arguments for consistency, and "
                    "note that it is only necessary to provide one of "
                    "these arguments.".format(self.name)
                )
        elif self.config.inlet_list is None and self.config.num_inlets is None:
            # If no arguments provided for inlets, default to num_inlets = 2
            self.config.num_inlets = 2

        # Create a list of names for inlet StateBlocks
        if self.config.inlet_list is not None:
            inlet_list = self.config.inlet_list
        else:
            inlet_list = [
                "inlet_{}".format(n) for n in range(1, self.config.num_inlets + 1)
            ]
        self.inlet_list = inlet_list

    def add_inlet_state_blocks(self):
        """
        Construct StateBlocks for all inlet streams.

        Args:
            list of strings to use as StateBlock names

        Returns:
            list of StateBlocks
        """
        # Setup StateBlock argument dict
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["defined_state"] = True

        # Create empty list to hold StateBlocks for return
        self.inlet_blocks = {}

        # Create an instance of StateBlock for all inlets
        for i in self.inlet_list:
            i_obj = self.config.property_package.build_state_block(
                self.flowsheet().time, doc="Material properties at inlet", **tmp_dict
            )
            setattr(self, "{}_state".format(i), i_obj)
            self.inlet_blocks[i] = i_obj

    def add_mixed_state_block(self):
        """
        Constructs StateBlock to represent mixed stream.

        Returns:
            New StateBlock object
        """
        # Setup StateBlock argument dict
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["defined_state"] = False

        self.mixed_state = self.config.property_package.build_state_block(
            self.flowsheet().time, doc="Material properties of mixed stream", **tmp_dict
        )
        return self.mixed_state

    def add_pressure_minimization_equations(self):
        """
        Add pressure minimization equations. This is done by sequential
        comparisons of each inlet to the minimum pressure so far, using
        the IDAES smooth minimum function.
        """
        units_meta = self.config.property_package.get_metadata()
        self.eps_pressure = Param(
            mutable=True,
            initialize=1e-3,
            domain=PositiveReals,
            doc="Smoothing term for minimum inlet pressure",
            units=units_meta.get_derived_units("pressure"),
        )

        # Calculate minimum inlet pressure
        @self.Expression(
            self.flowsheet().time,
            self.inlet_list,
            doc="Calculation for minimum inlet pressure",
        )
        def minimum_pressure(b, t, i):
            if i == self.inlet_list[0]:
                return self.inlet_blocks[i][t].pressure
            else:
                pi = self.inlet_list[self.inlet_list.index(i) - 1]
                prev_p = self.minimum_pressure[t, pi]
                this_p = self.inlet_blocks[i][t].pressure
                return smooth_min(this_p, prev_p, self.eps_pressure)

        # Set inlet pressure to minimum pressure
        @self.Constraint(self.flowsheet().time, doc="Link pressure to control volume")
        def minimum_pressure_constraint(b, t):
            return self.mixed_state[t].pressure == (
                self.minimum_pressure[t, self.inlet_list[-1]]
            )

    def add_pressure_equality_equations(self):
        """
        Add pressure equality equations. Note that this writes a number of
        constraints equal to the number of inlets, enforcing equality between
        all inlets and the mixed stream.
        """
        # Create equality constraints
        @self.Constraint(
            self.flowsheet().time,
            self.inlet_list,
            doc="Calculation for minimum inlet pressure",
        )
        def pressure_equality_constraints(b, t, i):
            return self.mixed_state[t].pressure == self.inlet_blocks[i][t].pressure

    def add_port_objects(self):
        """
        Adds Port objects if required.

        Args:
            a list of inlet StateBlock objects
            a mixed state StateBlock object

        Returns:
            None
        """
        for p in self.inlet_list:
            self.add_port(name=p, block=self.inlet_blocks[p], doc="Inlet Port")
        self.add_port(name="outlet", block=self.mixed_state, doc="Outlet Port")

    def use_minimum_inlet_pressure_constraint(self):
        """Activate the mixer pressure = mimimum inlet pressure constraint and
        deactivate the mixer pressure and all inlet pressures are equal
        constraints. This should only be used when momentum_mixing_type ==
        MomentumMixingType.minimize_and_equality.
        """
        if self.config.momentum_mixing_type != MomentumMixingType.minimize_and_equality:
            _log.warning(
                """use_minimum_inlet_pressure_constraint() can only be used
                when momentum_mixing_type ==
                MomentumMixingType.minimize_and_equality"""
            )
            return
        self.minimum_pressure_constraint.activate()
        self.pressure_equality_constraints.deactivate()

    def use_equal_pressure_constraint(self):
        """Deactivate the mixer pressure = mimimum inlet pressure constraint
        and activate the mixer pressure and all inlet pressures are equal
        constraints. This should only be used when momentum_mixing_type ==
        MomentumMixingType.minimize_and_equality.
        """
        if self.config.momentum_mixing_type != MomentumMixingType.minimize_and_equality:
            _log.warning(
                """use_equal_pressure_constraint() can only be used when
                momentum_mixing_type ==
                MomentumMixingType.minimize_and_equality"""
            )
            return
        self.minimum_pressure_constraint.deactivate()
        self.pressure_equality_constraints.activate()

    def initialize_build(self, outlvl=idaeslog.NOTSET, optarg=None, solver=None):
        """
        Initialization routine for mixer.

        Keyword Arguments:
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # This shouldn't require too much initializtion, just fixing inlets
        # and solving should always work.

        # sp is what to save to make sure state after init is same as the start
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        for b in self.inlet_blocks.values():
            for bdat in b.values():
                bdat.pressure.fix()
                bdat.enth_mol.fix()
                bdat.flow_mol.fix()

        for t, v in self.outlet.pressure.items():
            if not v.fixed:
                v.value = min(
                    [value(self.inlet_blocks[i][t].pressure) for i in self.inlet_blocks]
                )
        self.outlet.unfix()

        if (
            hasattr(self, "pressure_equality_constraints")
            and self.pressure_equality_constraints.active
        ):
            # If using the equal pressure constraint fix the outlet and free
            # the inlet pressures, this is typical for pressure driven flow
            for i, b in self.inlet_blocks.items():
                for bdat in b.values():
                    bdat.pressure.unfix()
            self.outlet.pressure.fix()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))
        from_json(self, sd=istate, wts=sp)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        for t, c in self.mass_balance.items():
            s = iscale.get_scaling_factor(self.mixed_state[t].flow_mol)
            iscale.constraint_scaling_transform(c, s, overwrite=False)
        for t, c in self.energy_balance.items():
            s = iscale.get_scaling_factor(self.mixed_state[t].enth_mol)
            s *= iscale.get_scaling_factor(self.mixed_state[t].flow_mol)
            iscale.constraint_scaling_transform(c, s, overwrite=False)
        if hasattr(self, "minimum_pressure_constraint"):
            for t, c in self.minimum_pressure_constraint.items():
                s = iscale.get_scaling_factor(self.mixed_state[t].pressure)
                iscale.constraint_scaling_transform(c, s, overwrite=False)
        if hasattr(self, "pressure_equality_constraints"):
            for (t, i), c in self.pressure_equality_constraints.items():
                s = iscale.get_scaling_factor(self.mixed_state[t].pressure)
                iscale.constraint_scaling_transform(c, s, overwrite=False)
