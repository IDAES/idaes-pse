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
Generic base template for a custom model. The basic functionality is to allow
interfacing with other IDAES unit models that leverage the ControlVolume
blocks. This is done by letting the user create ports and populate with
members that can connect with an IDAES unit model. Also, a default initialize
method is provided with a flexibility to pass a callback from the user.
"""

from pyomo.environ import Reference
from pyomo.network import Port
from pyomo.common.config import ConfigValue, In

from idaes.core.process_base import declare_process_block_class, \
    ProcessBlockData
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util import get_solver
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

    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare("dynamic", ConfigValue(
        default=False,
        domain=In([False]),
        description="Dynamic model flag",
        doc="""Model assumed to be steady-state."""))
    CONFIG.declare("custom_initializer", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Custom initializer",
        doc="""Custom initiliazer provided by the user"""))

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

        if self.config.custom_initializer:
            self.custom_initializer = None

        # TODO: Build checks for time domain and a warning if flowhseet is
        # dynamic, let the user know that the time index has to match
        # at ports

    def add_ports(self, name, member_dict, doc=None):
        """
        This is a method to build Port objects in a surrogate model and
        populate this with appropriate port members as specified. User can add
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
                "values being the variable objects declared in the model.")

        # Create empty Port
        p = Port(noruleinit=True, doc=doc)

        # Add port object to model
        setattr(self, name, p)

        # Populate port and map names to actual variables as defined
        for k in member_dict.keys():
            p.add(member_dict[k], name=k)

    def _custom_initialize_method(self):
        if self.custom_initializer is None:
            raise ConfigurationError(
                "_custom_initializer attribute was not set to a method. "
                "Set a custom method before calling the initialize method.")
        else:
            self.custom_initializer()

    def initialize(self, outlvl=idaeslog.NOTSET,
                   solver=None, optarg=None):
        """Initialize method for the custom model. If a custom initializer is
        provided then, this method will use that before a final solve
        is executed after checking that degrees of freedom is zero. If no
        custom initializer is provided, then the final solve is the only
        solve executed after checking for zero degrees of freedom.

        Args:
            outlvl ([idaes logger], optional):
                [Set idaes logger level]. Defaults to idaeslog.NOTSET.
            solver ([string], optional):
                [solver to be used for solve]. Defaults to None but if None,
                ipopt from IDAES is used.
            optarg ([dict], optional): [solver arguments]. Defaults to None
                but if None, default args for ipopt from IDAES is used.

        Raises:
            ConfigurationError: [Checks for degrees of freedom before final
            solve and raises exception when not zero.]
        """
        # Set solver options
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        opt = get_solver(solver=solver, options=optarg)

        if self.config.custom_initializer:
            self._custom_initialize_method()

        if degrees_of_freedom(self) == 0:
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(self, tee=slc.tee)
            init_log.info_high("Initialization completed {}.".
                               format(idaeslog.condition(res)))
        else:
            raise ConfigurationError(
                "Degrees of freedom is not zero during initialization. "
                "Fix/unfix appropriate number of variables to result "
                "in zero degrees of freedom for initialization.")

