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
NTU Condenser model, for use with Helmholtz EOS property packages.
"""

__author__ = "Jinliang Ma"

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.environ import units as pyunits
from pyomo.common.config import ConfigBlock, ConfigValue

# Import IDAES cores
from idaes.core import declare_process_block_class, UnitModelBlockData
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.core.solvers import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.models.unit_models.heater import (
    _make_heater_config_block,
    _make_heater_control_volume,
)
from idaes.models.unit_models.heat_exchanger import hx_process_config, add_hx_references
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog


def _make_heat_exchanger_config(config):
    """
    Declare configuration options for HeatExchangerData block.
    """
    config.declare(
        "hot_side_name",
        ConfigValue(
            default="shell",
            domain=str,
            doc="Hot side name, sets control volume and inlet and outlet names",
        ),
    )
    config.declare(
        "cold_side_name",
        ConfigValue(
            default="tube",
            domain=str,
            doc="Cold side name, sets control volume and inlet and outlet names",
        ),
    )
    config.declare(
        "hot_side",
        ConfigBlock(
            implicit=True,
            description="Config block for hot side",
            doc="""A config block used to construct the hot side control volume.
This config can be given by the hot side name instead of hot_side.""",
        ),
    )
    config.declare(
        "cold_side",
        ConfigBlock(
            implicit=True,
            description="Config block for cold side",
            doc="""A config block used to construct the cold side control volume.
This config can be given by the cold side name instead of cold_side.""",
        ),
    )
    _make_heater_config_block(config.hot_side)
    _make_heater_config_block(config.cold_side)


@declare_process_block_class("HelmNtuCondenser", doc="Simple 0D condenser model.")
class HelmNtuCondenserData(UnitModelBlockData):
    """
    Simple NTU condenser unit model.  This model assumes the property pacakages
    specified are Helmholtz EOS type.
    """

    CONFIG = UnitModelBlockData.CONFIG(implicit=True)
    _make_heat_exchanger_config(CONFIG)

    def build(self):
        """
        Building model

        Args:
            None
        Returns:
            None
        """
        ########################################################################
        #  Call UnitModel.build to setup dynamics and configure                #
        ########################################################################
        super().build()
        hx_process_config(self)
        config = self.config
        time = self.flowsheet().time

        ########################################################################
        # Add control volumes                                                  #
        ########################################################################
        hot_side = _make_heater_control_volume(
            self,
            "hot_side",
            config.hot_side,
            dynamic=config.dynamic,
            has_holdup=config.has_holdup,
        )
        cold_side = _make_heater_control_volume(
            self,
            "cold_side",
            config.cold_side,
            dynamic=config.dynamic,
            has_holdup=config.has_holdup,
        )
        ########################################################################
        # Add variables                                                        #
        ########################################################################
        # Use hot side units as basis
        s1_metadata = config.hot_side.property_package.get_metadata()

        f_units = s1_metadata.get_derived_units("flow_mole")
        cp_units = s1_metadata.get_derived_units("heat_capacity_mole")
        q_units = s1_metadata.get_derived_units("power")
        u_units = s1_metadata.get_derived_units("heat_transfer_coefficient")
        a_units = s1_metadata.get_derived_units("area")
        temp_units = s1_metadata.get_derived_units("temperature")

        self.overall_heat_transfer_coefficient = pyo.Var(
            time,
            domain=pyo.PositiveReals,
            initialize=100.0,
            doc="Overall heat transfer coefficient",
            units=u_units,
        )
        self.area = pyo.Var(
            domain=pyo.PositiveReals,
            initialize=1000.0,
            doc="Heat exchange area",
            units=a_units,
        )
        self.heat_duty = pyo.Reference(cold_side.heat)
        ########################################################################
        # Add ports                                                            #
        ########################################################################
        self.add_inlet_port(name="hot_side_inlet", block=hot_side, doc="Hot side inlet")
        self.add_inlet_port(
            name="cold_side_inlet",
            block=cold_side,
            doc="Cold side inlet",
        )
        self.add_outlet_port(
            name="hot_side_outlet", block=hot_side, doc="Hot side outlet"
        )
        self.add_outlet_port(
            name="cold_side_outlet",
            block=cold_side,
            doc="Cold side outlet",
        )
        ########################################################################
        # Add aliases                                                          #
        ########################################################################
        add_hx_references(self)

        ########################################################################
        # Add a unit level energy balance                                      #
        ########################################################################
        @self.Constraint(time, doc="Heat balance equation")
        def unit_heat_balance(b, t):
            return 0 == (
                hot_side.heat[t] + pyunits.convert(cold_side.heat[t], to_units=q_units)
            )

        ########################################################################
        # Add some useful expressions for condenser performance                #
        ########################################################################

        @self.Expression(time, doc="Inlet temperature difference")
        def delta_temperature_in(b, t):
            return hot_side.properties_in[t].temperature - pyunits.convert(
                cold_side.properties_in[t].temperature, temp_units
            )

        @self.Expression(time, doc="Outlet temperature difference")
        def delta_temperature_out(b, t):
            return hot_side.properties_out[t].temperature - pyunits.convert(
                cold_side.properties_out[t].temperature, temp_units
            )

        @self.Expression(time, doc="NTU Based temperature difference")
        def delta_temperature_ntu(b, t):
            return hot_side.properties_in[t].temperature_sat - pyunits.convert(
                cold_side.properties_in[t].temperature, temp_units
            )

        @self.Expression(
            time,
            doc="Minimum product of flow rate and heat "
            "capacity (always on tube side since shell side has phase change)",
        )
        def mcp_min(b, t):
            return pyunits.convert(
                cold_side.properties_in[t].flow_mol
                * cold_side.properties_in[t].cp_mol_phase["Liq"],
                f_units * cp_units,
            )

        @self.Expression(time, doc="Number of transfer units (NTU)")
        def ntu(b, t):
            return b.overall_heat_transfer_coefficient[t] * b.area / b.mcp_min[t]

        @self.Expression(time, doc="Condenser effectiveness factor")
        def effectiveness(b, t):
            return 1 - pyo.exp(-self.ntu[t])

        @self.Expression(time, doc="Heat treansfer")
        def heat_transfer(b, t):
            return b.effectiveness[t] * b.mcp_min[t] * b.delta_temperature_ntu[t]

        ########################################################################
        # Add Equations to calculate heat duty based on NTU method             #
        ########################################################################
        @self.Constraint(time, doc="Heat transfer rate equation based on NTU method")
        def heat_transfer_equation(b, t):
            return pyunits.convert(cold_side.heat[t], q_units) == self.heat_transfer[t]

        @self.Constraint(
            time, doc="Shell side outlet enthalpy is saturated water enthalpy"
        )
        def saturation_eqn(b, t):
            return (
                hot_side.properties_out[t].enth_mol
                == hot_side.properties_in[t].enth_mol_sat_phase["Liq"]
            )

    def set_initial_condition(self):
        if self.config.dynamic is True:
            self.hot_side.material_accumulation[:, :, :].value = 0
            self.hot_side.energy_accumulation[:, :].value = 0
            self.hot_side.material_accumulation[0, :, :].fix(0)
            self.hot_side.energy_accumulation[0, :].fix(0)
            self.cold_side.material_accumulation[:, :, :].value = 0
            self.cold_side.energy_accumulation[:, :].value = 0
            self.cold_side.material_accumulation[0, :, :].fix(0)
            self.cold_side.energy_accumulation[0, :].fix(0)

    def initialize_build(
        self,
        state_args_1=None,
        state_args_2=None,
        unfix="hot_flow",
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Condenser initialization method. The initialization routine assumes
        fixed area and heat transfer coefficient and adjusts the cooling water
        flow to condense steam to saturated water at shell side pressure.

        Args:
            state_args_1 : a dict of arguments to be passed to the property
                initialization for hot side (see documentation of the specific
                property package) (default = None).
            state_args_2 : a dict of arguments to be passed to the property
                initialization for cold side (see documentation of the specific
                property package) (default = None).
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None

        """
        if unfix not in {"hot_flow", "cold_flow", "pressure"}:
            raise Exception(
                "Condenser free variable must be in 'hot_flow', "
                "'cold_flow', or 'pressure'"
            )
        # Set solver options
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        hot_side = getattr(self, self.config.hot_side_name)
        cold_side = getattr(self, self.config.cold_side_name)

        # Store initial model specs, restored at the end of initializtion, so
        # the problem is not altered.  This can restore fixed/free vars,
        # active/inactive constraints, and fixed variable values.
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        # Create solver
        opt = get_solver(solver, optarg)

        flags1 = hot_side.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_1
        )
        flags2 = cold_side.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_2
        )
        init_log.info_high("Initialization Step 1 Complete.")

        # Solve with all constraints activated
        self.saturation_eqn.activate()
        if unfix == "pressure":
            hot_side.properties_in[:].pressure.unfix()
        elif unfix == "hot_flow":
            hot_side.properties_in[:].flow_mol.unfix()
        elif unfix == "cold_flow":
            cold_side.properties_in[:].flow_mol.unfix()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 4 {}.".format(idaeslog.condition(res)))

        # Release Inlet state
        hot_side.release_state(flags1, outlvl)
        cold_side.release_state(flags2, outlvl)
        from_json(self, sd=istate, wts=sp)

    def _get_performance_contents(self, time_point=0):
        var_dict = {
            "HX Coefficient": self.overall_heat_transfer_coefficient[time_point]
        }
        var_dict["HX Area"] = self.area
        var_dict["Heat Duty"] = self.heat_duty[time_point]
        expr_dict = {}
        expr_dict["Delta T In"] = self.delta_temperature_in[time_point]
        expr_dict["Delta T Out"] = self.delta_temperature_out[time_point]

        return {"vars": var_dict, "exprs": expr_dict}

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Hot Inlet": self.hot_inlet,
                "Hot Outlet": self.hot_outlet,
                "Cold Inlet": self.cold_inlet,
                "Cold Outlet": self.cold_outlet,
            },
            time_point=time_point,
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        area_sf_default = 1e-2
        overall_heat_transfer_coefficient_sf_default = 1e-2

        # Function to set defaults so I don't need to reproduce the same code
        def _fill_miss_with_default(name, s):
            try:
                c = getattr(self, name)
            except AttributeError:
                return  # it's okay if the attribute doesn't exist, spell careful
            if iscale.get_scaling_factor(c) is None:
                for ci in c.values():
                    if iscale.get_scaling_factor(ci) is None:
                        iscale.set_scaling_factor(ci, s)

        # Set defaults where scale factors are missing
        _fill_miss_with_default("area", area_sf_default)
        _fill_miss_with_default(
            "overall_heat_transfer_coefficient",
            overall_heat_transfer_coefficient_sf_default,
        )

        for t, c in self.heat_transfer_equation.items():
            sf = iscale.get_scaling_factor(self.cold_side.heat[t])
            iscale.constraint_scaling_transform(c, sf, overwrite=False)

        for t, c in self.unit_heat_balance.items():
            sf = iscale.get_scaling_factor(self.cold_side.heat[t])
            iscale.constraint_scaling_transform(c, sf, overwrite=False)

        for t, c in self.saturation_eqn.items():
            sf = iscale.get_scaling_factor(self.hot_side.properties_out[t].enth_mol)
            iscale.constraint_scaling_transform(c, sf, overwrite=False)
