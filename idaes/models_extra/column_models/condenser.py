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
Condenser model for distillation.

While the condenser model (both total and partial), is fairly simple, a major
portion of this code has gone into making this generic and be able to handle
different state variables and the associated splits.
"""

__author__ = "Jaffer Ghouse"

from enum import Enum
from pandas import DataFrame

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.network import Port
from pyomo.environ import (
    Reference,
    Var,
    Constraint,
    value,
    Set,
    check_optimal_termination,
)

# Import IDAES cores
import idaes.logger as idaeslog
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    EnergyBalanceType,
    MaterialBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import (
    ConfigurationError,
    InitializationError,
)
from idaes.core.solvers import get_solver
from idaes.models_extra.column_models.util import make_phase_split

_log = idaeslog.getLogger(__name__)


class CondenserType(Enum):
    """
    Enum for supported condenser types.
    """

    totalCondenser = 0
    partialCondenser = 1


class TemperatureSpec(Enum):
    """
    Enum for temperature specifications.
    """

    atBubblePoint = 0
    customTemperature = 1


@declare_process_block_class("Condenser")
class CondenserData(UnitModelBlockData):
    """
    Condenser unit for distillation model.
    Unit model to condense (total/partial) the vapor from the top tray of
    the distillation column.
    """

    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        "condenser_type",
        ConfigValue(
            default=CondenserType.totalCondenser,
            domain=In(CondenserType),
            description="Type of condenser flag",
            doc="""Indicates what type of condenser should be constructed,
**default** - CondenserType.totalCondenser.
**Valid values:** {
**CondenserType.totalCondenser** - Incoming vapor from top tray is condensed
to all liquid,
**CondenserType.partialCondenser** - Incoming vapor from top tray is
partially condensed to a vapor and liquid stream.}""",
        ),
    )
    CONFIG.declare(
        "temperature_spec",
        ConfigValue(
            default=None,
            domain=In(TemperatureSpec),
            description="Temperature spec for the condenser",
            doc="""Temperature specification for the condenser,
**default** - TemperatureSpec.none
**Valid values:** {
**TemperatureSpec.none** - No spec is selected,
**TemperatureSpec.atBubblePoint** - Condenser temperature set at
bubble point i.e. total condenser,
**TemperatureSpec.customTemperature** - Condenser temperature at
user specified temperature.}""",
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
        super(CondenserData, self).build()

        # Check config arguments
        if self.config.temperature_spec is None:
            raise ConfigurationError(
                "temperature_spec config argument "
                "has not been specified. Please select "
                "a valid option."
            )
        if (self.config.condenser_type == CondenserType.partialCondenser) and (
            self.config.temperature_spec == TemperatureSpec.atBubblePoint
        ):
            raise ConfigurationError(
                "condenser_type set to partial but "
                "temperature_spec set to atBubblePoint. "
                "Select customTemperature and specify "
                "outlet temperature."
            )

        # Add Control Volume for the condenser
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

        # Note: No momentum balance added for the condenser as the condenser
        # outlet pressure is a spec set by the user.

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

        self._make_ports()

        @self.Expression(doc="Split fraction calculation calculations")
        def reflux_split_fraction(self):
            return self.reflux_ratio * (1 + self.reflux_ratio) ** -1

        if self.config.condenser_type == CondenserType.totalCondenser:

            make_phase_split(
                self.control_volume,
                port=self.reflux,
                phase=self._liquid_set,
                side_sf=self.reflux_split_fraction,
                equipmentType=CondenserType.totalCondenser,
            )

            make_phase_split(
                self.control_volume,
                port=self.distillate,
                phase=self._liquid_set,
                side_sf=1 - self.reflux_split_fraction,
                equipmentType=CondenserType.totalCondenser,
            )

            if self.config.temperature_spec == TemperatureSpec.atBubblePoint:
                # Option 1: if true, condition for total condenser
                # (T_cond = T_bubble)
                # Option 2: if this is false, then user has selected
                # custom temperature spec and needs to fix an outlet
                # temperature.

                # Get index for bubble point temperature and and assume it
                # will have only a single phase equilibrium pair. This is to
                # support the generic property framework where the T_bubble
                # is indexed by the phases_in_equilibrium. In distillation,
                # the assumption is that there will only be a single pair
                # i.e. vap-liq.
                idx = next(
                    iter(
                        self.control_volume.properties_out[
                            self.flowsheet().time.first()
                        ].temperature_bubble
                    )
                )

                def rule_total_cond(self, t):
                    return (
                        self.control_volume.properties_out[t].temperature
                        == self.control_volume.properties_out[t].temperature_bubble[idx]
                    )

                self.eq_total_cond_spec = Constraint(
                    self.flowsheet().time, rule=rule_total_cond
                )

        else:

            make_phase_split(
                self.control_volume,
                port=self.reflux,
                phase=self._liquid_set,
                side_sf=self.reflux_split_fraction,
                equipmentType=CondenserType.partialCondenser,
            )

            make_phase_split(
                self.control_volume,
                port=self.distillate,
                phase=self._liquid_set,
                side_sf=1 - self.reflux_split_fraction,
                equipmentType=CondenserType.partialCondenser,
            )
            make_phase_split(
                self.control_volume,
                port=self.vapor_outlet,
                phase=self._vapor_set,
                # Split fraction set to 1 as all vapor from condenser
                # returns to column
                side_sf=1,
                equipmentType=CondenserType.partialCondenser,
            )

        # Add object reference to variables of the control volume
        # Reference to the heat duty
        self.heat_duty = Reference(self.control_volume.heat[:])

        self.condenser_pressure = Reference(
            self.control_volume.properties_out[:].pressure
        )

    def _make_ports(self):

        # Add Ports for the condenser
        # Inlet port (the vapor from the top tray)
        self.add_inlet_port()

        # Outlet ports that always exist irrespective of condenser type
        self.reflux = Port(
            noruleinit=True, doc="Reflux stream that is" " returned to the top tray."
        )
        self.distillate = Port(
            noruleinit=True, doc="Distillate stream that is" " the top product."
        )

        if self.config.condenser_type == CondenserType.partialCondenser:
            self.vapor_outlet = Port(
                noruleinit=True, doc="Vapor outlet port from a " "partial condenser"
            )
        # Add codnenser specific variables
        self.reflux_ratio = Var(initialize=0.5, doc="Reflux ratio for the condenser")

    def initialize(
        self, state_args=None, solver=None, optarg=None, outlvl=idaeslog.NOTSET
    ):

        # TODO: Fix the inlets to the condenser to the vapor flow from
        # the top tray or take it as an argument to this method.

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        if self.config.temperature_spec == TemperatureSpec.customTemperature:
            if degrees_of_freedom(self) != 0:
                raise ConfigurationError(
                    "Degrees of freedom is not 0 during initialization. "
                    "Check if outlet temperature has been fixed in addition "
                    "to the other inputs required as customTemperature was "
                    "selected for temperature_spec config argument."
                )

        solverobj = get_solver(solver, optarg)

        if state_args is None:
            state_args = {}
            state_dict = self.control_volume.properties_in[
                self.flowsheet().time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = value(state_dict[k][m])
                else:
                    state_args[k] = value(state_dict[k])

        if self.config.condenser_type == CondenserType.totalCondenser:
            self.eq_total_cond_spec.deactivate()

        # Initialize the inlet and outlet state blocks
        flags = self.control_volume.initialize(
            state_args=state_args,
            solver=solver,
            optarg=optarg,
            outlvl=outlvl,
            hold_state=True,
        )

        # Activate the total condenser spec
        if self.config.condenser_type == CondenserType.totalCondenser:
            self.eq_total_cond_spec.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solverobj.solve(self, tee=slc.tee)
        init_log.info("Initialization Complete, {}.".format(idaeslog.condition(res)))
        if not check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        self.control_volume.release_state(flags=flags)

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        if hasattr(self, "heat_duty"):
            var_dict["Heat Duty"] = self.heat_duty[time_point]
        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]

        return {"vars": var_dict}

    def _get_stream_table_contents(self, time_point=0):
        stream_attributes = {}

        if self.config.condenser_type == CondenserType.totalCondenser:
            stream_dict = {
                "Inlet": "inlet",
                "Reflux": "reflux",
                "Distillate": "distillate",
            }
        else:
            stream_dict = {
                "Inlet": "inlet",
                "Vapor Outlet": "vapor_outlet",
                "Reflux": "reflux",
                "Distillate": "distillate",
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
