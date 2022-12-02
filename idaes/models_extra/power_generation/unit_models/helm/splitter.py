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
Flow splitter depending on a Helmholtz EOS property package.  This just
multiplies the flow by the split fractions for the outlets and has a constraint
(sum_split) that requires the split fractions to sum to 1.  The usual way to
specify this unit is to have fix or have constraints that set the inlet, and
have n_outlets - 1 specified split fractions or outlet flows.

This model is psuedo-steady-state when used in dynamic mode.
"""

from pyomo.environ import Var, value
from pyomo.common.config import ConfigBlock, ConfigValue, In, ListOf

from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale


__author__ = "John Eslick"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("HelmSplitter")
class HelmSplitterData(UnitModelBlockData):
    """
    This is a basic stream splitter which splits flow into outlet streams based
    on split fractions. This does not do phase seperation, and assumes that you
    are using a Helmholtz EOS propery package with P-H state variables. In
    dynamic mode this uses a pseudo-steady-state model.

    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for mixer",
            doc="""Property parameter object used to define property
calculations,
**default** - useDefault.
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
        "outlet_list",
        ConfigValue(
            domain=ListOf(str),
            description="List of outlet names",
            doc="""A list containing names of outlets,
**default** - None.
**Valid values:** {
**None** - use num_outlets argument,
**list** - a list of names to use for outlets.}""",
        ),
    )
    CONFIG.declare(
        "num_outlets",
        ConfigValue(
            domain=int,
            description="Number of outlets to unit",
            doc="""Argument indicating number (int) of outlets to construct,
not used if outlet_list arg is provided,
**default** - None.
**Valid values:** {
**None** - use outlet_list arg instead, or default to 2 if neither argument
provided,
**int** - number of outlets to create (will be named with sequential integers
from 1 to num_outlets).}""",
        ),
    )

    def build(self):
        """
        Build a splitter.

        Args:
            None

        Returns:
            None
        """
        time = self.flowsheet().time
        super().build()

        self._get_property_package()

        self.create_outlet_list()
        self.add_inlet_state_and_port()
        self.add_outlet_state_blocks()
        self.add_outlet_port_objects()

        self.split_fraction = Var(
            time,
            self.outlet_list,
            initialize=1.0 / len(self.outlet_list),
            doc="Split fractions for outlet streams",
        )

        @self.Constraint(time, doc="Splt constraint")
        def sum_split(b, t):
            return 1 == sum(self.split_fraction[t, o] for o in self.outlet_list)

        @self.Constraint(time, self.outlet_list, doc="Pressure constraint")
        def pressure_eqn(b, t, o):
            o_block = getattr(self, "{}_state".format(o))
            return self.mixed_state[t].pressure == o_block[t].pressure

        @self.Constraint(time, self.outlet_list, doc="Enthalpy constraint")
        def enthalpy_eqn(b, t, o):
            o_block = getattr(self, "{}_state".format(o))
            return self.mixed_state[t].enth_mol == o_block[t].enth_mol

        @self.Constraint(time, self.outlet_list, doc="Flow constraint")
        def flow_eqn(b, t, o):
            o_block = getattr(self, "{}_state".format(o))
            sf = self.split_fraction[t, o]
            return self.mixed_state[t].flow_mol * sf == o_block[t].flow_mol

    def add_inlet_state_and_port(self):
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["defined_state"] = True
        self.mixed_state = self.config.property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties of mixed (inlet) stream",
            **tmp_dict
        )
        self.add_port(name="inlet", block=self.mixed_state, doc="Inlet Port")

    def create_outlet_list(self):
        """
        Create list of outlet stream names based on config arguments.

        Returns:
            list of strings
        """
        config = self.config
        if config.outlet_list is not None and config.num_outlets is not None:
            # If both arguments provided and not consistent, raise Exception
            if len(config.outlet_list) != config.num_outlets:
                raise ConfigurationError(
                    "{} Splitter provided with both outlet_list and "
                    "num_outlets arguments, which were not consistent ("
                    "length of outlet_list was not equal to num_outlets). "
                    "Please check your arguments for consistency, and "
                    "note that it is only necessry to provide one of "
                    "these arguments.".format(self.name)
                )
        elif config.outlet_list is None and config.num_outlets is None:
            # If no arguments provided for outlets, default to num_outlets = 2
            config.num_outlets = 2

        # Create a list of names for outlet StateBlocks
        if config.outlet_list is not None:
            outlet_list = self.config.outlet_list
        else:
            outlet_list = [
                "outlet_{}".format(n) for n in range(1, config.num_outlets + 1)
            ]
        self.outlet_list = outlet_list

    def add_outlet_state_blocks(self):
        """
        Construct StateBlocks for all outlet streams.

        Args:
            None

        Returns:
            list of StateBlocks
        """
        # Setup StateBlock argument dict
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = False

        # Create empty list to hold StateBlocks for return
        self.outlet_blocks = {}

        # Create an instance of StateBlock for all outlets
        for o in self.outlet_list:
            o_obj = self.config.property_package.build_state_block(
                self.flowsheet().time, doc="Material properties at outlet", **tmp_dict
            )
            setattr(self, o + "_state", o_obj)
            self.outlet_blocks[o] = o_obj

    def add_outlet_port_objects(self):
        """
        Adds outlet Port objects if required.

        Args:
            None

        Returns:
            None
        """
        self.outlet_ports = {}
        for p in self.outlet_list:
            self.add_port(name=p, block=self.outlet_blocks[p], doc="Outlet")
            self.outlet_ports[p] = getattr(self, p)

    def initialize_build(self, outlvl=idaeslog.NOTSET, optarg=None, solver=None):
        """
        Initialization routine for splitter

        Keyword Arguments:
            outlvl: sets output level of initialization routine
            optarg: solver options dictionary object (default=None, use
                    default solver options)
            solver: str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            If hold_states is True, returns a dict containing flags for which
            states were fixed during initialization.
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # sp is what to save to make sure state after init is same as the start
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        # check for fixed outlet flows and use them to calculate fixed split
        # fractions
        for t in self.flowsheet().time:
            for o in self.outlet_list:
                if self.outlet_blocks[o][t].flow_mol.fixed:
                    self.split_fraction[t, o].fix(
                        value(
                            self.mixed_state[t].flow_mol
                            / self.outlet_blocks[o][t].flow_mol
                        )
                    )

        # fix or unfix split fractions so n - 1 are fixed
        for t in self.flowsheet().time:
            # see how many split fractions are fixed
            n = sum(1 for o in self.outlet_list if self.split_fraction[t, o].fixed)
            # if number of outlets - 1 we're good
            if n == len(self.outlet_list) - 1:
                continue
            # if too mant are fixed un fix the first, generally assume that is
            # the main flow, and is the calculated split fraction
            if n == len(self.outlet_list):
                self.split_fraction[t, self.outlet_list[0]].unfix()
            # if not enough fixed, start fixing from the back until there are
            # are enough
            for o in reversed(self.outlet_list):
                if not self.split_fraction[t, o].fixed:
                    self.split_fraction[t, o].fix()
                    n += 1
                if n == len(self.outlet_list) - 1:
                    break

        # This model is really simple so it should easily solve without much
        # effort to initialize
        self.inlet.fix()
        for o, p in self.outlet_ports.items():
            p.unfix()
        assert degrees_of_freedom(self) == 0
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        from_json(self, sd=istate, wts=sp)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        for (t, i), c in self.pressure_eqn.items():
            o_block = getattr(self, "{}_state".format(i))
            s = iscale.get_scaling_factor(o_block[t].pressure)
            iscale.constraint_scaling_transform(c, s, overwrite=False)
        for (t, i), c in self.enthalpy_eqn.items():
            o_block = getattr(self, "{}_state".format(i))
            s = iscale.get_scaling_factor(o_block[t].enth_mol)
            iscale.constraint_scaling_transform(c, s, overwrite=False)
        for (t, i), c in self.flow_eqn.items():
            o_block = getattr(self, "{}_state".format(i))
            s = iscale.get_scaling_factor(o_block[t].flow_mol)
            iscale.constraint_scaling_transform(c, s, overwrite=False)
