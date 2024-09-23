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
Unit model to adjust size of streams to represent, for example, a stream being split across several identical units,
which are then all modeled as a single IDAES unit
"""
from enum import Enum
from functools import partial

from pyomo.environ import (
    Block,
    check_optimal_termination,
    Param,
    PositiveReals,
    Reals,
    RangeSet,
    units as pyunits,
    Var,
)
from pyomo.network import Port
from pyomo.common.config import ConfigBlock, ConfigValue, In, ListOf, Bool

from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
    MaterialBalanceType,
    MaterialFlowBasis,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_state_block,
)
from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    PropertyNotSupportedError,
    InitializationError,
)
from idaes.core.base.var_like_expression import VarLikeExpression
from idaes.core.util.math import smooth_min
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

__author__ = "Douglas Allan"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("StreamScaler")
class StreamScalerData(UnitModelBlockData):
    """
    TODO
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
                **default** = False. Scaler blocks are always steady-state.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="Scaler blocks do not contain holdup, thus this must be False.",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for StreamScaler",
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

    def build(self):
        """
        General build method for StreamScalerData. This method calls a number
        of sub-methods which automate the construction of expected attributes
        of unit models.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        # Call super.build()
        super(StreamScalerData, self).build()

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = True

        # Call setup methods from ControlVolumeBlockData
        self._get_property_package()
        self._get_indexing_sets()

        self.inlet_block = self.config.property_package.build_state_block(
                self.flowsheet().time, doc="Material properties at inlet", **tmp_dict
        )
        self.outlet_block = Block()
        self.multiplier = Var(
            initialize=1,
            domain=PositiveReals,
            units=pyunits.dimensionless,
            doc="Factor by which to scale dimensionless streams"
        )
        self.add_inlet_port(name="inlet", block=self.inlet_block)
        self.outlet = Port(doc="Outlet port")

        def rule_scale_var(b, *args, var=None):
            return self.multiplier * var[args]

        def rule_no_scale_var(b, *args, var=None):
            return var[args]

        for var_name in self.inlet.vars.keys():
            var = getattr(self.inlet, var_name)
            if "flow" in var_name:
                rule=partial(rule_scale_var, var=var)
            else:
                rule=partial(rule_no_scale_var, var=var)
            self.outlet_block.add_component(
                var_name,
                VarLikeExpression(
                    var.index_set(),
                    rule=rule
                )
            )
            expr = getattr(self.outlet_block, var_name)
            self.outlet.add(expr, var_name)

    def initialize_build(
        blk, outlvl=idaeslog.NOTSET, optarg=None, solver=None, hold_state=False
    ):
        """
        Initialization routine for StreamScaler.

        Keyword Arguments:
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)
            hold_state : flag indicating whether the initialization routine
                     should unfix any state variables fixed during
                     initialization, **default** - False. **Valid values:**
                     **True** - states variables are not unfixed, and a dict of
                     returned containing flags for which states were fixed
                     during initialization, **False** - state variables are
                     unfixed after initialization by calling the release_state
                     method.

        Returns:
            If hold_states is True, returns a dict containing flags for which
            states were fixed during initialization.
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver

        # Initialize inlet state blocks
        flags = blk.inlet_block.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        if hold_state is True:
            return flags
        else:
            blk.release_state(flags, outlvl=outlvl)

    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        """
        Method to release state variables fixed during initialization.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state = True.
            outlvl : sets output level of logging

        Returns:
            None
        """
        blk.inlet_block.release_state(flags, outlvl=outlvl)

    def _get_stream_table_contents(self, time_point=0):
        io_dict = {
            "Inlet": self.inlet,
            "Outlet": self.outlet,
        }
        return create_stream_table_dataframe(io_dict, time_point=time_point)

    def calculate_scaling_factors(self):
        # Scaling factors for the property block are calculated automatically
        super().calculate_scaling_factors()

        # Need to pass on scaling factors from the property block to the outlet
        # VarLikeExpressions so arcs get scaled right
        scale = 1/self.multiplier.value
        for var_name in self.inlet.vars.keys():
            var = getattr(self.inlet, var_name)
            outlet_expr = getattr(self.outlet, var_name)
            for key, subvar in var.items():
                sf = iscale.get_scaling_factor(subvar, default=1, warning=True)
                iscale.set_scaling_factor(outlet_expr[key],scale*sf)
        