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
Generic template for a translator block.
"""
# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
from idaes.core import declare_process_block_class, UnitModelBlockData
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

__author__ = "Andrew Lee"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Translator")
class TranslatorData(UnitModelBlockData):
    """
    Standard Translator Block Class
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Translator blocks are always steady-state.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Translator blocks do not contain holdup.""",
        ),
    )
    CONFIG.declare(
        "outlet_state_defined",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Indicated whether outlet state will be fully defined",
            doc="""Indicates whether unit model will fully define outlet state.
If False, the outlet property package will enforce constraints such as sum
of mole fractions and phase equilibrium.
**default** - True.
**Valid values:** {
**True** - outlet state will be fully defined,
**False** - outlet property package should enforce sumation and equilibrium
constraints.}""",
        ),
    )
    CONFIG.declare(
        "has_phase_equilibrium",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Indicates whether outlet is in phase equilibrium",
            doc="""Indicates whether outlet property package should enforce
phase equilibrium constraints.
**default** - False.
**Valid values:** {
**True** - outlet property package should calculate phase equilibrium,
**False** - outlet property package should notcalculate phase equilibrium.}
""",
        ),
    )
    CONFIG.declare(
        "inlet_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for incoming stream",
            doc="""Property parameter object used to define property
calculations for the incoming stream,
**default** - None.
**Valid values:** {
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "inlet_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property package "
            "of the incoming stream",
            doc="""A ConfigBlock with arguments to be passed to the property
block associated with the incoming stream,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "outlet_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for outgoing stream",
            doc="""Property parameter object used to define property
calculations for the outgoing stream,
**default** - None.
**Valid values:** {
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "outlet_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property package "
            "of the outgoing stream",
            doc="""A ConfigBlock with arguments to be passed to the property
block associated with the outgoing stream,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )

    def build(self):
        """
        Begin building model.

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(TranslatorData, self).build()

        # Check construction argumnet consistency
        if self.config.outlet_state_defined and self.config.has_phase_equilibrium:
            raise ConfigurationError(
                "{} cannot calcuate phase equilibrium (has_phase_equilibrium "
                "= True) when outlet state is set to be fully defined ("
                "outlet_state_defined = True).".format(self.name)
            )

        # Add State Blocks
        self.properties_in = self.config.inlet_property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties in incoming stream",
            defined_state=True,
            has_phase_equilibrium=False,
            **self.config.inlet_property_package_args
        )

        self.properties_out = self.config.outlet_property_package.build_state_block(
            self.flowsheet().time,
            doc="Material properties in outgoing stream",
            defined_state=self.config.outlet_state_defined,
            has_phase_equilibrium=self.config.has_phase_equilibrium,
            **self.config.outlet_property_package_args
        )

        # Add outlet port
        self.add_port(name="inlet", block=self.properties_in, doc="Inlet Port")
        self.add_port(name="outlet", block=self.properties_out, doc="Outlet Port")

    def initialize_build(
        blk,
        state_args_in=None,
        state_args_out=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        This method calls the initialization method of the state blocks.

        Keyword Arguments:
            state_args_in : a dict of arguments to be passed to the inlet
                            property package (to provide an initial state for
                            initialization (see documentation of the specific
                            property package) (default = None).
            state_args_out : a dict of arguments to be passed to the outlet
                             property package (to provide an initial state for
                             initialization (see documentation of the specific
                             property package) (default = None).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize state block
        flags = blk.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_in,
            hold_state=True,
        )

        blk.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )

        if degrees_of_freedom(blk) == 0:
            with idaeslog.solver_log(init_log, idaeslog.DEBUG) as slc:
                res = opt.solve(blk, tee=slc.tee)

            init_log.info("Initialization Complete {}.".format(idaeslog.condition(res)))
        else:
            init_log.warning(
                "Initialization incomplete. Degrees of freedom "
                "were not zero. Please provide sufficient number "
                "of constraints linking the state variables "
                "between the two state blocks."
            )

        blk.properties_in.release_state(flags=flags, outlvl=outlvl)
