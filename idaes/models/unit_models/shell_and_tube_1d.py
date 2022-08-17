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
1D Single pass shell and tube HX model with 0D wall conduction model.

This model derives from the HeatExchanger1D unit model.
"""
# Import Pyomo libraries
from pyomo.environ import (
    Var,
    check_optimal_termination,
    Constraint,
    value,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
from idaes.core import declare_process_block_class
from idaes.models.unit_models.heat_exchanger_1D import (
    HeatExchanger1DData,
    HeatExchangerFlowPattern,
)

from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants as c
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

import idaes.logger as idaeslog


__author__ = "Jaffer Ghouse"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ShellAndTube")
class ShellAndTubeData(HeatExchanger1DData):
    """1D Shell and Tube HX Unit Model Class."""

    CONFIG = HeatExchanger1DData.CONFIG(implicit=True)

    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super().build()

    # def _make_performance(self):
    #     """
    #     Constraints for unit model.
    #
    #     Args:
    #         None
    #
    #     Returns:
    #         None
    #     """
    #     hot_side_units = (
    #         self.config.hot_side.property_package.get_metadata().get_derived_units
    #     )
    #     cold_side_units = (
    #         self.config.cold_side.property_package.get_metadata().get_derived_units
    #     )
    #
    #     # Unit model variables
    #     # HX dimensions
    #     self.d_hot_side = Var(
    #         initialize=1, doc="Diameter of hot side", units=hot_side_units("length")
    #     )
    #     self.d_cold_side_outer = Var(
    #         initialize=0.011,
    #         doc="Outer diameter of cold side",
    #         units=hot_side_units("length"),
    #     )
    #     self.d_cold_side_inner = Var(
    #         initialize=0.010,
    #         doc="Inner diameter of cold side",
    #         units=hot_side_units("length"),
    #     )
    #     self.N_tubes = Var(
    #         initialize=1, doc="Number of tubes", units=pyunits.dimensionless
    #     )
    #
    #     # Note: In addition to the above variables, "hot_side_length" and
    #     # "cold_side_length" need to be fixed at the flowsheet level
    #
    #     # Performance variables
    #     self.hot_side_heat_transfer_coefficient = Var(
    #         self.flowsheet().time,
    #         self.hot_side.length_domain,
    #         initialize=50,
    #         doc="Heat transfer coefficient",
    #         units=hot_side_units("heat_transfer_coefficient"),
    #     )
    #     self.cold_side_heat_transfer_coefficient = Var(
    #         self.flowsheet().time,
    #         self.cold_side.length_domain,
    #         initialize=50,
    #         doc="Heat transfer coefficient",
    #         units=cold_side_units("heat_transfer_coefficient"),
    #     )
    #
    #     # Wall 0D model (Q_hot_side = Q_cold_side*N_tubes)
    #     if self.config.has_wall_conduction == WallConductionType.zero_dimensional:
    #         self.temperature_wall = Var(
    #             self.flowsheet().time,
    #             self.cold_side.length_domain,
    #             initialize=298.15,
    #             units=hot_side_units("temperature"),
    #         )
    #
    #         # Performance equations
    #         # Energy transfer between hot_side and cold_side wall
    #
    #         @self.Constraint(
    #             self.flowsheet().time,
    #             self.hot_side.length_domain,
    #             doc="Heat transfer between hot_side and cold_side",
    #         )
    #         def hot_side_heat_transfer_eq(self, t, x):
    #             return self.hot_side.heat[t, x] == -self.N_tubes * (
    #                 self.hot_side_heat_transfer_coefficient[t, x]
    #                 * c.pi
    #                 * self.d_cold_side_outer
    #                 * (
    #                     self.hot_side.properties[t, x].temperature
    #                     - self.temperature_wall[t, x]
    #                 )
    #             )
    #
    #         # Energy transfer between cold_side wall and cold_side
    #         @self.Constraint(
    #             self.flowsheet().time,
    #             self.cold_side.length_domain,
    #             doc="Convective heat transfer",
    #         )
    #         def cold_side_heat_transfer_eq(self, t, x):
    #             return self.cold_side.heat[
    #                 t, x
    #             ] == self.cold_side_heat_transfer_coefficient[
    #                 t, x
    #             ] * c.pi * pyunits.convert(
    #                 self.d_cold_side_inner, to_units=cold_side_units("length")
    #             ) * (
    #                 pyunits.convert(
    #                     self.temperature_wall[t, x],
    #                     to_units=cold_side_units("temperature"),
    #                 )
    #                 - self.cold_side.properties[t, x].temperature
    #             )
    #
    #         if hot_side_units("length") is None:
    #             # Backwards compatability check
    #             q_units = None
    #         else:
    #             q_units = hot_side_units("power") / hot_side_units("length")
    #         # Wall 0D model
    #         @self.Constraint(
    #             self.flowsheet().time,
    #             self.hot_side.length_domain,
    #             doc="wall 0D model",
    #         )
    #         def wall_0D_model(self, t, x):
    #             return pyunits.convert(
    #                 self.cold_side.heat[t, x], to_units=q_units
    #             ) == -(self.hot_side.heat[t, x] / self.N_tubes)
    #
    #     else:
    #         raise NotImplementedError(
    #             "{} HeatExchanger1D has not yet implemented support for "
    #             "wall conduction models."
    #         )
    #
    #     # Define cold_side area in terms of tube diameter
    #     self.area_calc_cold_side = Constraint(
    #         expr=4 * self.cold_side_area
    #         == c.pi
    #         * pyunits.convert(
    #             self.d_cold_side_inner, to_units=cold_side_units("length")
    #         )
    #         ** 2
    #     )
    #
    #     # Define hot_side area in terms of hot_side and tube diameter
    #     self.area_calc_hot_side = Constraint(
    #         expr=4 * self.hot_side_area
    #         == c.pi
    #         * (self.d_hot_side**2 - self.N_tubes * self.d_cold_side_outer**2)
    #     )

    # def initialize_build(
    #     self,
    #     hot_side_state_args=None,
    #     cold_side_state_args=None,
    #     outlvl=idaeslog.NOTSET,
    #     solver=None,
    #     optarg=None,
    # ):
    #     """
    #     Initialization routine for the unit.
    #
    #     Keyword Arguments:
    #         state_args : a dict of arguments to be passed to the property
    #                      package(s) to provide an initial state for
    #                      initialization (see documentation of the specific
    #                      property package) (default = {}).
    #         outlvl : sets output level of initialization routine
    #         optarg : solver options dictionary object (default=None, use
    #                  default solver options)
    #         solver : str indicating which solver to use during
    #                  initialization (default = None, use default solver)
    #
    #     Returns:
    #         None
    #     """
    #     init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
    #     solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
    #
    #     # Create solver
    #     opt = get_solver(solver, optarg)
    #
    #     # ---------------------------------------------------------------------
    #     # Initialize hot_side block
    #     flags_hot_side = self.hot_side.initialize(
    #         outlvl=outlvl,
    #         optarg=optarg,
    #         solver=solver,
    #         state_args=hot_side_state_args,
    #     )
    #
    #     flags_cold_side = self.cold_side.initialize(
    #         outlvl=outlvl,
    #         optarg=optarg,
    #         solver=solver,
    #         state_args=cold_side_state_args,
    #     )
    #
    #     init_log.info_high("Initialization Step 1 Complete.")
    #
    #     # ---------------------------------------------------------------------
    #     # Solve unit
    #     # Wall 0D
    #     if self.config.has_wall_conduction == WallConductionType.zero_dimensional:
    #         hot_side_units = (
    #             self.config.hot_side.property_package.get_metadata().get_derived_units
    #         )
    #         for t in self.flowsheet().time:
    #             for z in self.hot_side.length_domain:
    #                 self.temperature_wall[t, z].fix(
    #                     value(
    #                         0.5
    #                         * (
    #                             self.hot_side.properties[t, 0].temperature
    #                             + pyunits.convert(
    #                                 self.cold_side.properties[t, 0].temperature,
    #                                 to_units=hot_side_units("temperature"),
    #                             )
    #                         )
    #                     )
    #                 )
    #
    #         self.cold_side.deactivate()
    #         self.cold_side_heat_transfer_eq.deactivate()
    #         self.wall_0D_model.deactivate()
    #
    #         with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
    #             res = opt.solve(self, tee=slc.tee)
    #         init_log.info_high(
    #             "Initialization Step 2 {}.".format(idaeslog.condition(res))
    #         )
    #
    #         self.cold_side.activate()
    #         self.cold_side_heat_transfer_eq.activate()
    #
    #         with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
    #             res = opt.solve(self, tee=slc.tee)
    #         init_log.info_high(
    #             "Initialization Step 3 {}.".format(idaeslog.condition(res))
    #         )
    #
    #         self.wall_0D_model.activate()
    #         self.temperature_wall.unfix()
    #
    #         with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
    #             res = opt.solve(self, tee=slc.tee)
    #         init_log.info_high(
    #             "Initialization Step 4 {}.".format(idaeslog.condition(res))
    #         )
    #     else:
    #         res = None
    #
    #     self.hot_side.release_state(flags_hot_side)
    #     self.cold_side.release_state(flags_cold_side)
    #
    #     if res is not None and not check_optimal_termination(res):
    #         raise InitializationError(
    #             f"{self.name} failed to initialize successfully. Please check "
    #             f"the output logs for more information."
    #         )
    #
    #     init_log.info("Initialization Complete.")

    # def _get_performance_contents(self, time_point=0):
    #     # TODO: Set this up to use user names if available
    #     var_dict = {}
    #     var_dict["Hot Side Area"] = self.hot_side.area
    #     var_dict["Hot Side Diameter"] = self.d_hot_side
    #     var_dict["Hot Side Length"] = self.hot_side.length
    #     var_dict["Cold Side Area"] = self.cold_side.area
    #     var_dict["Cold Side Outer Diameter"] = self.d_cold_side_outer
    #     var_dict["Cold Side Inner Diameter"] = self.d_cold_side_inner
    #     var_dict["Cold Side Length"] = self.cold_side.length
    #     var_dict["Number of Tubes"] = self.N_tubes
    #
    #     return {"vars": var_dict}
    #
    # def _get_stream_table_contents(self, time_point=0):
    #     # TODO : Set this up to use user provided names if available
    #     return create_stream_table_dataframe(
    #         {
    #             "Hot Side Inlet": self.hot_side_inlet,
    #             "Hot Side Outlet": self.hot_side_outlet,
    #             "Cold Side Inlet": self.cold_side_inlet,
    #             "Cold Side Outlet": self.cold_side_outlet,
    #         },
    #         time_point=time_point,
    #     )
    #
    # def calculate_scaling_factors(self):
    #     super().calculate_scaling_factors()
    #
    #     for i, c in self.hot_side_heat_transfer_eq.items():
    #         iscale.constraint_scaling_transform(
    #             c,
    #             iscale.get_scaling_factor(
    #                 self.hot_side.heat[i], default=1, warning=True
    #             ),
    #             overwrite=False,
    #         )
    #
    #     for i, c in self.cold_side_heat_transfer_eq.items():
    #         iscale.constraint_scaling_transform(
    #             c,
    #             iscale.get_scaling_factor(
    #                 self.cold_side.heat[i], default=1, warning=True
    #             ),
    #             overwrite=False,
    #         )
