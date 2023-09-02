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
Reboiler model for distillation.

While the reboiler model, is fairly simple, a major
portion of this code has gone into making this generic and be able to handle
different state variables and the associated splits.
"""

__author__ = "Jaffer Ghouse"

from pandas import DataFrame

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool
from pyomo.network import Port
from pyomo.environ import (
    check_optimal_termination,
    Reference,
    Var,
    Constraint,
    value,
    Set,
)

# Import IDAES cores
import idaes.logger as idaeslog
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    EnergyBalanceType,
    MomentumBalanceType,
    MaterialBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import (
    PropertyNotSupportedError,
    ConfigurationError,
    InitializationError,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models_extra.column_models.util import make_phase_split


_log = idaeslog.getIdaesLogger(__name__)


@declare_process_block_class("Reboiler")
class ReboilerData(UnitModelBlockData):
    """
    Reboiler unit for distillation model.
    Unit model to reboil the liquid from the bottom tray of
    the distillation column.
    """

    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        "has_boilup_ratio",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Boilup ratio term construction flag",
            doc="""Indicates whether terms for boilup ratio should be
constructed,
**default** - False.
**Valid values:** {
**True** - include construction of boilup ratio constraint,
**False** - exclude construction of boilup ratio constraint}""",
        ),
    )
    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.componentPhase.
**Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
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
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )

    def build(self):
        """Build the model.
        Args:
            None
        Returns:
            None
        """
        # Setup model build logger
        model_log = idaeslog.getModelLogger(self.name, tag="unit")

        # Call UnitModel.build to setup dynamics
        super(ReboilerData, self).build()

        # Add Control Volume for the Reboiler
        self.control_volume = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.control_volume.add_state_blocks(has_phase_equilibrium=True)

        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type, has_phase_equilibrium=True
        )

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type, has_heat_transfer=True
        )

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # Get liquid and vapor phase objects from the property package
        # to be used below. Avoids repetition.
        _liquid_list = []
        _vapor_list = []
        for p in self.control_volume.config.property_package.phase_list:
            pobj = self.control_volume.config.property_package.get_phase(p)
            if pobj.is_vapor_phase():
                _vapor_list.append(p)
            elif pobj.is_liquid_phase():
                _liquid_list.append(p)
            else:
                _liquid_list.append(p)
                model_log.warning(
                    "A non-liquid/non-vapor phase was detected but will "
                    "be treated as a liquid."
                )

        # Create a pyomo set for indexing purposes. This set is appended to
        # model otherwise results in an abstract set.
        self._liquid_set = Set(initialize=_liquid_list)
        self._vapor_set = Set(initialize=_vapor_list)

        if self.config.has_boilup_ratio is True:

            self.boilup_ratio = Var(initialize=0.5, doc="Boilup ratio for reboiler")

            def rule_boilup_ratio(self, t):
                if hasattr(self.control_volume.properties_out[t], "flow_mol_phase"):
                    return self.boilup_ratio * sum(
                        self.control_volume.properties_out[t].flow_mol_phase[p]
                        for p in self._liquid_set
                    ) == sum(
                        self.control_volume.properties_out[t].flow_mol_phase["Vap"]
                        for p in self._vapor_set
                    )
                elif hasattr(
                    self.control_volume.properties_out[t], "flow_mol_phase_comp"
                ):
                    return self.boilup_ratio * sum(
                        self.control_volume.properties_out[t].flow_mol_phase_comp[p, i]
                        for p in self._liquid_set
                        for i in self.control_volume.properties_out[
                            t
                        ].params.component_list
                    ) == sum(
                        self.control_volume.properties_out[t].flow_mol_phase_comp[p, i]
                        for p in self._vapor_set
                        for i in self.control_volume.properties_out[
                            t
                        ].params.component_list
                    )
                else:
                    raise PropertyNotSupportedError(
                        "Unrecognized names for flow variables encountered "
                        "while building the constraint for reboiler."
                    )

            self.eq_boilup_ratio = Constraint(
                self.flowsheet().time, rule=rule_boilup_ratio
            )

        self.add_inlet_port()

        # Outlet ports that always exist irrespective of reboiler type
        self.bottoms = Port(noruleinit=True, doc="Bottoms stream.")

        self.vapor_reboil = Port(
            noruleinit=True,
            doc="Vapor outlet stream that is returned to " "to the bottom tray.",
        )

        make_phase_split(
            self.control_volume,
            port=self.bottoms,
            phase=self._liquid_set,
            side_sf=1,
            equipmentType="Reboiler",
        )

        make_phase_split(
            self.control_volume,
            port=self.vapor_reboil,
            phase=self._vapor_set,
            side_sf=1,
            equipmentType="Reboiler",
        )

        # Add object reference to variables of the control volume
        # Reference to the heat duty
        self.heat_duty = Reference(self.control_volume.heat[:])

        # Reference to the pressure drop (if set to True)
        if self.config.has_pressure_change:
            self.deltaP = Reference(self.control_volume.deltaP[:])

    def _make_ports(self):

        # Add Ports for the reboiler
        # Inlet port (the vapor from the top tray)
        self.add_inlet_port()

        # Outlet ports that always exist irrespective of reboiler type
        self.bottoms = Port(noruleinit=True, doc="Bottoms stream.")

        self.vapor_reboil = Port(
            noruleinit=True,
            doc="Vapor outlet stream that is returned to " "to the bottom tray.",
        )

    def initialize(
        self, state_args=None, solver=None, optarg=None, outlvl=idaeslog.NOTSET
    ):

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        solverobj = get_solver(solver, optarg)

        # Initialize the inlet and outlet state blocks. Calling the state
        # blocks initialize methods directly so that custom set of state args
        # can be passed to the inlet and outlet state blocks as control_volume
        # initialize method initializes the state blocks with the same
        # state conditions.
        flags = self.control_volume.properties_in.initialize(
            state_args=state_args,
            solver=solver,
            optarg=optarg,
            outlvl=outlvl,
            hold_state=True,
        )

        # Initialize outlet state block at same conditions of inlet except
        # the temperature. Set the temperature to a temperature guess based
        # on the desired boilup_ratio.

        # Get index for bubble point temperature and and assume it
        # will have only a single phase equilibrium pair. This is to
        # support the generic property framework where the T_bubble
        # is indexed by the phases_in_equilibrium. In distillation,
        # the assumption is that there will only be a single pair
        # i.e. vap-liq.
        idx = next(iter(self.control_volume.properties_in[0].temperature_bubble))
        temp_guess = 0.5 * (
            value(self.control_volume.properties_in[0].temperature_dew[idx])
            - value(self.control_volume.properties_in[0].temperature_bubble[idx])
        ) + value(self.control_volume.properties_in[0].temperature_bubble[idx])

        state_args_outlet = {}
        state_dict_outlet = self.control_volume.properties_in[
            self.flowsheet().time.first()
        ].define_port_members()

        for k in state_dict_outlet.keys():
            if state_dict_outlet[k].is_indexed():
                state_args_outlet[k] = {}
                for m in state_dict_outlet[k].keys():
                    state_args_outlet[k][m] = value(state_dict_outlet[k][m])
            else:
                if k != "temperature":
                    state_args_outlet[k] = value(state_dict_outlet[k])
                else:
                    state_args_outlet[k] = temp_guess

        self.control_volume.properties_out.initialize(
            state_args=state_args_outlet,
            solver=solver,
            optarg=optarg,
            outlvl=outlvl,
            hold_state=False,
        )

        if degrees_of_freedom(self) == 0:
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solverobj.solve(self, tee=slc.tee)
            init_log.info(
                "Initialization Complete, {}.".format(idaeslog.condition(res))
            )
        else:
            raise ConfigurationError(
                "State vars fixed but degrees of freedom "
                "for reboiler is not zero during "
                "initialization. Please ensure that the boilup_ratio "
                "or the outlet temperature is fixed."
            )

        if not check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        self.control_volume.properties_in.release_state(flags=flags, outlvl=outlvl)

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        if hasattr(self, "heat_duty"):
            var_dict["Heat Duty"] = self.heat_duty[time_point]

        return {"vars": var_dict}

    def _get_stream_table_contents(self, time_point=0):
        stream_attributes = {}

        stream_dict = {
            "Inlet": "inlet",
            "Vapor Reboil": "vapor_reboil",
            "Bottoms": "bottoms",
        }

        for n, v in stream_dict.items():
            port_obj = getattr(self, v)

            stream_attributes[n] = {}

            for k in port_obj.vars:
                for i in port_obj.vars[k].keys():
                    if isinstance(i, float):
                        stream_attributes[n][k] = value(port_obj.vars[k][time_point])
                    else:
                        if len(i) == 2:
                            kname = str(i[1])
                        else:
                            kname = str(i[1:])
                        stream_attributes[n][k + " " + kname] = value(
                            port_obj.vars[k][time_point, i[1:]]
                        )

        return DataFrame.from_dict(stream_attributes, orient="columns")
