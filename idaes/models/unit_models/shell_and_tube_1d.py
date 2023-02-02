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
from pyomo.common.config import ConfigValue, Bool

# Import IDAES cores
from idaes.core import declare_process_block_class, UnitModelBlockData
from idaes.models.unit_models.heat_exchanger_1D import (
    HeatExchanger1DData,
)
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.constants import Constants as c
from idaes.core.util import scaling as iscale
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.solvers import get_solver

import idaes.logger as idaeslog


__author__ = "Jaffer Ghouse"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ShellAndTube1D")
class ShellAndTube1DData(HeatExchanger1DData):
    """1D Shell and Tube HX Unit Model Class."""

    CONFIG = HeatExchanger1DData.CONFIG(implicit=True)
    CONFIG.declare(
        "shell_is_hot",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Shell side contains hot fluid",
            doc="""Boolean flag indicating whether shell side contains hot fluid (default=True).
    If True, shell side will be the hot_side, if False shell side will be cold_side.""",
        ),
    )

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

    def _process_config(self):
        super()._process_config()

        # Check for custom names, and if not present assign defaults
        if self.config.hot_side_name is None:
            if self.config.shell_is_hot:
                self.config.hot_side_name = "Shell"
            else:
                self.config.hot_side_name = "Tube"

        if self.config.cold_side_name is None:
            if self.config.shell_is_hot:
                self.config.cold_side_name = "Tube"
            else:
                self.config.cold_side_name = "Shell"

    def _make_geometry(self):
        # Add reference to control volume geometry
        add_object_reference(self, "hot_side_area", self.hot_side.area)
        add_object_reference(self, "length", self.hot_side.length)

        # Equate hot and cold side geometries
        hot_side_units = (
            self.config.hot_side.property_package.get_metadata().get_derived_units
        )
        cold_side_units = (
            self.config.cold_side.property_package.get_metadata().get_derived_units
        )

        @self.Constraint(
            self.flowsheet().time, doc="Equating hot and cold side lengths"
        )
        def length_equality(self, t):
            return (
                pyunits.convert(
                    self.cold_side.length, to_units=hot_side_units("length")
                )
                == self.hot_side.length
            )

        # Get hot and cold sides
        if self.config.shell_is_hot:
            shell = self.hot_side
            tube = self.cold_side
        else:
            shell = self.cold_side
            tube = self.hot_side

        # New Unit model variables and constraints
        self.shell_diameter = Var(
            initialize=0.011,
            doc="Diameter of shell",
            units=hot_side_units("length"),
        )
        self.tube_outer_diameter = Var(
            initialize=1,
            doc="Outer diameter of tubes",
            units=hot_side_units("length"),
        )
        self.tube_inner_diameter = Var(
            initialize=0.010,
            doc="Inner diameter of tubes",
            units=hot_side_units("length"),
        )
        self.number_of_tubes = Var(
            initialize=1, doc="Number of tubes", units=pyunits.dimensionless
        )

        # Calculate cross-sectional area of control volumes
        self.tube_side_xsec_area_calc = Constraint(
            expr=4 * tube.area
            == self.number_of_tubes
            * c.pi
            * pyunits.convert(
                self.tube_inner_diameter, to_units=cold_side_units("length")
            )
            ** 2,
            doc="Tube side cross-sectional area",
        )
        # Need to account for area occupied by tubes
        self.shell_side_xsec_area_calc = Constraint(
            expr=4 * shell.area
            == c.pi
            * (
                self.shell_diameter**2
                - self.number_of_tubes * self.tube_outer_diameter**2
            ),
            doc="Shell side cross-sectional area",
        )

    def _make_performance(self):
        hot_side_units = (
            self.config.hot_side.property_package.get_metadata().get_derived_units
        )

        # Performance variables
        self.hot_side_heat_transfer_coefficient = Var(
            self.flowsheet().time,
            self.hot_side.length_domain,
            initialize=50,
            doc="Hot side heat transfer coefficient",
            units=hot_side_units("heat_transfer_coefficient"),
        )
        self.cold_side_heat_transfer_coefficient = Var(
            self.flowsheet().time,
            self.cold_side.length_domain,
            initialize=50,
            doc="Cold side heat transfer coefficient",
            units=hot_side_units("heat_transfer_coefficient"),
        )
        self.temperature_wall = Var(
            self.flowsheet().time,
            self.hot_side.length_domain,
            initialize=298.15,
            units=hot_side_units("temperature"),
        )

        @self.Constraint(
            self.flowsheet().time,
            self.hot_side.length_domain,
            doc="Heat transfer between hot_side and wall",
        )
        def hot_side_heat_transfer_eq(self, t, x):
            return self.hot_side.heat[t, x] == -(
                self.hot_side_heat_transfer_coefficient[t, x]
                * self.number_of_tubes
                * c.pi
                * self.tube_outer_diameter
                * (
                    self.hot_side.properties[t, x].temperature
                    - self.temperature_wall[t, x]
                )
            )

        @self.Constraint(
            self.flowsheet().time,
            self.cold_side.length_domain,
            doc="Heat transfer between cold_side and wall",
        )
        def cold_side_heat_transfer_eq(self, t, x):
            return pyunits.convert(
                self.cold_side.heat[t, x],
                to_units=hot_side_units("power") / hot_side_units("length"),
            ) == (
                self.cold_side_heat_transfer_coefficient[t, x]
                * self.number_of_tubes
                * c.pi
                * self.tube_inner_diameter
                * (
                    self.temperature_wall[t, x]
                    - pyunits.convert(
                        self.cold_side.properties[t, x].temperature,
                        to_units=hot_side_units("temperature"),
                    )
                )
            )

    def initialize_build(
        self,
        hot_side_state_args=None,
        cold_side_state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Initialization routine for the unit.

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
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

        # ---------------------------------------------------------------------
        # Initialize hot_side block
        flags_hot_side = self.hot_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=hot_side_state_args,
        )

        flags_cold_side = self.cold_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=cold_side_state_args,
        )

        init_log.info_high("Initialization Step 1 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        hot_side_units = (
            self.config.hot_side.property_package.get_metadata().get_derived_units
        )
        for t in self.flowsheet().time:
            for z in self.hot_side.length_domain:
                self.temperature_wall[t, z].fix(
                    value(
                        0.5
                        * (
                            self.hot_side.properties[t, 0].temperature
                            + pyunits.convert(
                                self.cold_side.properties[t, 0].temperature,
                                to_units=hot_side_units("temperature"),
                            )
                        )
                    )
                )

        self.cold_side.deactivate()
        self.cold_side_heat_transfer_eq.deactivate()
        self.heat_conservation.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        self.cold_side.activate()
        self.cold_side_heat_transfer_eq.activate()
        self.heat_conservation.activate()
        self.temperature_wall.unfix()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        self.hot_side.release_state(flags_hot_side)
        self.cold_side.release_state(flags_cold_side)

        if res is not None and not check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        init_log.info("Initialization Complete.")

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Shell Diameter"] = self.shell_diameter
        var_dict["Length"] = self.length
        var_dict["Tube Outer Diameter"] = self.tube_outer_diameter
        var_dict["Tube Inner Diameter"] = self.tube_inner_diameter
        var_dict["Number of Tubes"] = self.number_of_tubes

        return {"vars": var_dict}

    def _get_stream_table_contents(self, time_point=0):
        # Get names for hot and cold sides
        hot_name = self.config.hot_side_name
        cold_name = self.config.cold_side_name
        return create_stream_table_dataframe(
            {
                f"{hot_name} Inlet": self.hot_side_inlet,
                f"{hot_name} Outlet": self.hot_side_outlet,
                f"{cold_name} Inlet": self.cold_side_inlet,
                f"{cold_name} Outlet": self.cold_side_outlet,
            },
            time_point=time_point,
        )

    def calculate_scaling_factors(self):
        super(UnitModelBlockData, self).calculate_scaling_factors()

        for i, c in self.hot_side_heat_transfer_eq.items():
            print(
                iscale.get_scaling_factor(
                    self.hot_side.heat[i], default=1, warning=False
                )
            )
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.hot_side.heat[i], default=1, warning=True
                ),
                overwrite=False,
            )

        for i, c in self.cold_side_heat_transfer_eq.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.hot_side.heat[i], default=1, warning=True
                ),
                overwrite=False,
            )

        for i, c in self.heat_conservation.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.hot_side.heat[i], default=1, warning=True
                ),
                overwrite=False,
            )
