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
Generic template for a surrogate unit model.
"""

from pyomo.environ import Reference, SolverFactory
from pyomo.network import Port
from pyomo.common.config import ConfigValue, In

from idaes.core.process_base import (declare_process_block_class,
                           ProcessBlockData,
                           useDefault)
from idaes.core.util.exceptions import (BurntToast,
                                        ConfigurationError,
                                        PropertyPackageError,
                                        BalanceTypeNotSupportedError)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.logger as idaeslog

__author__ = "Jaffer Ghouse"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("SurrogateModel")
class SurrogateModelData(ProcessBlockData):
    """
    This is the class for a surrogate unit model. This will provide a user
    with a generic skeleton to add a surrogate unit model.
    """

    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare("dynamic", ConfigValue(
        default=False,
        domain=In([False]),
        description="Dynamic model flag",
        doc="""Surrogate model assumed to be steady-state."""))

    def build(self):
        """
        General build method for SurrogateModelBlockData.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        super(SurrogateModelData, self).build()

    def add_ports(self, name=None, member_list=None, doc=None):
        """
        This is a method to build Port objects in a surrogate model and
        populate this with appropriate port members as specified.

        Keyword Args:
            name : name to use for Port object.
            dict : dictionary containing variables to be added to the port
            doc : doc string for Port object

        Returns:
            A Pyomo Port object and associated components.
        """
        # Validate that members is a dict
        if not isinstance(member_list, dict):
            raise ConfigurationError("{} block object provided to add_port "
                                     "method is not an instance of a "
                                     "StateBlock object. IDAES port objects "
                                     "should only be associated with "
                                     "StateBlocks.".format(self.name))

        # Create empty Port
        p = Port(noruleinit=True, doc=doc)
        setattr(self, name, p)

        # Create References for port members
        for k in member_list.keys():
            if not member_list[k].is_indexed():
                local_name = member_list[k]
            else:
                local_name = member_list[k]

            # Add Reference to Port
            p.add(Reference(local_name), k)

    def _get_stream_table_contents(self, time_point=0):
        """
        Assume unit has standard configuration of 1 inlet and 1 outlet.

        Developers should overload this as appropriate.
        """
        try:
            return create_stream_table_dataframe({"Inlet": self.inlet,
                                                  "Outlet": self.outlet},
                                                 time_point=time_point)
        except AttributeError:
            raise ConfigurationError(
                "Unit model {self.name} does not have the standard Port "
                "names (inet and outlet). Please contact the unit model "
                "developer to develop a unit specific stream table.")

    def initialize(self, custom_initialize = None, outlvl=idaeslog.NOTSET,
                   solver='ipopt', optarg={'tol': 1e-6}):
        '''
        This is a simple initialization routine for surrogate
        models. If the user, does not provide a custom callback function, this
        method will check for degrees of freedom and if zero, will attempt
        to solve the model as is.

        Keyword Arguments:
            custom_initialize : custom call back for intialization
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating which solver to use during
                     initialization (default = 'ipopt')

        Returns:
            None
        '''
        # Set solver options
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        opt = SolverFactory(solver)
        opt.options = optarg

        if custom_initialize is not None:
            custom_initialize
        else:
            if degrees_of_freedom(self) == 0:
                with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                    res = opt.solve(self, tee=slc.tee)
                init_log.info_high("Initialization completed {}.".
                    format(idaeslog.condition(res)))
            else:
                raise ConfigurationError(
                    "Degrees of freedom is not 0 during initialization. "
                    "Fix/unfix appropriate number of variables to result "
                    "in 0 degrees of freedom or provide a custom callback "
                    "function for initialization.")

