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
The basic functionality is to allow a user to add any model of choice that
can interface with other IDAES unit models that leverage the ControlVolume
blocks. This is done by letting the user create ports and populate with
members that match with state variables defined in an IDAES unit model.
Also, a default initialize method is provided with a flexibility to pass a
callback from the user.
"""

from pyomo.network import Port
from pyomo.common.config import ConfigValue, In

from idaes.core.base.process_base import declare_process_block_class, ProcessBlockData
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

__author__ = "Jaffer Ghouse"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("SkeletonUnitModel")
class SkeletonUnitModelData(ProcessBlockData):
    """
    This is the class for a skeleteon unit model. This will provide a user
    with a generic skeleton to add a custom or a surrogate unit model
    when controlvolumes are not required.
    """

    # This is a staticmethod that will be the default callable set for the
    # initializer flag in the config block.
    @staticmethod
    def _default_initializer(
        model, opt=None, init_log=None, solve_log=None, initial_guess=None
    ):
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(model, tee=slc.tee)
        init_log.info(
            "Initialization completed using default method {}.".format(
                idaeslog.condition(res)
            )
        )

    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Dynamic model flag",
            doc="""Model does not support dynamics or holdup automatically.
        The variables declared in this unit will need to be indexed with the
        flowsheet time domain to facilitate connection with other
        unit models.""",
        ),
    )
    CONFIG.declare(
        "initializer",
        ConfigValue(
            # Note: staticmethod is unbound, extract with __func__
            default=_default_initializer.__func__,
            description="Initializer to be used",
            doc="""Flag to set the callback from user for initialization.
        Default is set to use the default_initialize method.""",
        ),
    )

    def build(self):
        """
        General build method for SkeletonUnitModelBlockData.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        super().build()

    def add_ports(self, name, member_dict, doc=None):
        """
        This is a method to build Port objects in the SkeletonUnitModel and
        populate them with appropriate port members as specified. User can add
        as many inlet and outlet ports as required.

        Keyword Args:
            name : name to use for Port object.
            dict : dictionary containing variables to be added to the port
            doc : doc string for Port object

        Returns:
            A Pyomo Port object and associated components.
        """
        # Validate that member_dict is a dict
        if not isinstance(member_dict, dict):
            raise ConfigurationError(
                "member_dict should be a dictionary "
                "with the keys being the name assigned(strings) and "
                "values being the variable objects declared in the model."
            )

        # Create empty Port
        p = Port(noruleinit=True, doc=doc)

        # Add port object to model
        setattr(self, name, p)

        # Populate port and map names to actual variables as defined
        for k in member_dict.keys():
            p.add(member_dict[k], name=k)

    # TODO : Work out how to make this work with new UnitModel initialization
    def initialize(
        self, outlvl=idaeslog.NOTSET, solver=None, optarg=None, initial_guess=None
    ):
        """Initialize method for the SkeletonUnitModel. If a custom function
        is provided via the initializer argument in the config block,
        then, this method will use it. If no custom function is
        provided, the _default_initializer method is used. Irrespective
        of which method is being used (default or custom), the expectation is
        that the degrees of freedom is zero at the start of initialization.

        Args:
            outlvl (idaes logger, optional):
                Set idaes logger level. Defaults to idaeslog.NOTSET.
            solver (string, optional):
                solver to be used for solve. Defaults to None but if None,
                ipopt from IDAES is used.
            optarg ({dict}, optional): solver arguments. Defaults to None
                but if None, default args for ipopt from IDAES is used.
            initial_guess ({dict}, optional): initial guess that is passed
            to custom initialize. Not used for the default initialize method.

        Raises:
            ConfigurationError: If degrees of freedom is not zero at the start
            of initialization.
        """

        # Set solver options
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        opt = get_solver(solver=solver, options=optarg)

        if degrees_of_freedom(self) == 0:
            self.config.initializer(
                self,
                opt=opt,
                init_log=init_log,
                solve_log=solve_log,
                initial_guess=initial_guess,
            )
        else:
            raise ConfigurationError(
                "Degrees of freedom is not zero during start of "
                "initialization. Fix/unfix appropriate number of variables "
                "to result in zero degrees of freedom for initialization."
            )
