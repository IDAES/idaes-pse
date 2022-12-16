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
Condenser model for distillation.

While the condenser model (both total and partial), is fairly simple, a major
portion of this code has gone into making this generic and be able to handle
different state variables and the associated splits.
"""

__author__ = "Jaffer Ghouse"

from pandas import DataFrame
from enum import Enum

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.network import Port
from pyomo.environ import (
    Reference,
    Expression,
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
    PropertyPackageError,
    ConfigurationError,
    PropertyNotSupportedError,
    InitializationError,
)
from idaes.core.solvers import get_solver

_log = idaeslog.getLogger(__name__)


class CondenserType(Enum):
    totalCondenser = 0
    partialCondenser = 1


class TemperatureSpec(Enum):
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
        # to be used below. Avoids repition.
        _liquid_list = []
        _vapor_list = []
        for p in self.config.property_package.phase_list:
            pobj = self.config.property_package.get_phase(p)
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

        if self.config.condenser_type == CondenserType.totalCondenser:

            self._make_splits_total_condenser()

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
            self._make_splits_partial_condenser()

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
        self.reflux_ratio = Var(initialize=1, doc="Reflux ratio for the condenser")

    def _make_splits_total_condenser(self):

        # Get dict of Port members and names
        member_list = self.control_volume.properties_out[0].define_port_members()

        # Create references and populate the reflux, distillate ports
        for k in member_list:

            local_name = member_list[k].local_name
            # Create references and populate the intensive variables
            if "flow" not in local_name:
                if not member_list[k].is_indexed():
                    var = self.control_volume.properties_out[:].component(local_name)
                else:
                    var = self.control_volume.properties_out[:].component(local_name)[
                        ...
                    ]

                # add the reference and variable name to the reflux port
                self.reflux.add(Reference(var), k)

                # add the reference and variable name to the distillate port
                self.distillate.add(Reference(var), k)

            elif "flow" in local_name:
                # Create references and populate the extensive variables
                # This is for vars that are not indexed
                if not member_list[k].is_indexed():
                    # Expression for reflux flow and relation to the
                    # reflux_ratio variable

                    def rule_reflux_flow(self, t):
                        return self.control_volume.properties_out[t].component(
                            local_name
                        ) * (self.reflux_ratio / (1 + self.reflux_ratio))

                    self.e_reflux_flow = Expression(
                        self.flowsheet().time, rule=rule_reflux_flow
                    )
                    self.reflux.add(self.e_reflux_flow, k)

                    # Expression for distillate flow and relation to the
                    # reflux_ratio variable
                    def rule_distillate_flow(self, t):
                        return self.control_volume.properties_out[t].component(
                            local_name
                        ) / (1 + self.reflux_ratio)

                    self.e_distillate_flow = Expression(
                        self.flowsheet().time, rule=rule_distillate_flow
                    )
                    self.distillate.add(self.e_distillate_flow, k)
                else:
                    # Create references and populate the extensive variables
                    # This is for vars that are indexed by phase, comp or both.
                    index_set = member_list[k].index_set()

                    def rule_reflux_flow(self, t, *args):
                        return self.control_volume.properties_out[t].component(
                            local_name
                        )[args] * (self.reflux_ratio / (1 + self.reflux_ratio))

                    self.e_reflux_flow = Expression(
                        self.flowsheet().time, index_set, rule=rule_reflux_flow
                    )
                    self.reflux.add(self.e_reflux_flow, k)

                    def rule_distillate_flow(self, t, *args):
                        return self.control_volume.properties_out[t].component(
                            local_name
                        )[args] / (1 + self.reflux_ratio)

                    self.e_distillate_flow = Expression(
                        self.flowsheet().time, index_set, rule=rule_distillate_flow
                    )
                    self.distillate.add(self.e_distillate_flow, k)
            else:
                raise PropertyNotSupportedError(
                    "Unrecognized names for flow variables encountered while "
                    "building the condenser ports."
                )

    def _make_splits_partial_condenser(self):
        # Get dict of Port members and names
        member_list = self.control_volume.properties_out[0].define_port_members()

        # Create references and populate the reflux, distillate ports
        for k in member_list:

            local_name = member_list[k].local_name

            # Create references and populate the intensive variables
            if (
                "flow" not in local_name
                and "frac" not in local_name
                and "enth" not in local_name
            ):
                if not member_list[k].is_indexed():
                    var = self.control_volume.properties_out[:].component(local_name)
                else:
                    var = self.control_volume.properties_out[:].component(local_name)[
                        ...
                    ]

                # add the reference and variable name to the reflux port
                self.reflux.add(Reference(var), k)

                # add the reference and variable name to the distillate port
                self.distillate.add(Reference(var), k)

                # add the reference and variable name to the
                # vapor outlet port
                self.vapor_outlet.add(Reference(var), k)

            elif "frac" in local_name:

                # Mole/mass frac is typically indexed
                index_set = member_list[k].index_set()

                # if state var is not mole/mass frac by phase
                if "phase" not in local_name:
                    if "mole" in local_name:  # check mole basis/mass basis

                        # The following conditionals are required when a
                        # mole frac or mass frac is a state var i.e. will be
                        # a port member. This gets a bit tricky when handling
                        # non-conventional systems when you have more than one
                        # liquid or vapor phase. Hence, the logic here is that
                        # the mole frac that should be present in the liquid or
                        # vapor port should be computed by accounting for
                        # multiple liquid or vapor phases if present. For the
                        # classical VLE system, this holds too.
                        if hasattr(
                            self.control_volume.properties_out[0],
                            "mole_frac_phase_comp",
                        ) and hasattr(
                            self.control_volume.properties_out[0], "flow_mol_phase"
                        ):
                            flow_phase_comp = False
                            local_name_frac = "mole_frac_phase_comp"
                            local_name_flow = "flow_mol_phase"
                        elif hasattr(
                            self.control_volume.properties_out[0], "flow_mol_phase_comp"
                        ):
                            flow_phase_comp = True
                            local_name_flow = "flow_mol_phase_comp"
                        else:
                            raise PropertyNotSupportedError(
                                "No mole_frac_phase_comp or flow_mol_phase or"
                                " flow_mol_phase_comp variables encountered "
                                "while building ports for the condenser. "
                            )
                    elif "mass" in local_name:
                        if hasattr(
                            self.control_volume.properties_out[0],
                            "mass_frac_phase_comp",
                        ) and hasattr(
                            self.control_volume.properties_out[0], "flow_mass_phase"
                        ):
                            flow_phase_comp = False
                            local_name_frac = "mass_frac_phase_comp"
                            local_name_flow = "flow_mass_phase"
                        elif hasattr(
                            self.control_volume.properties_out[0],
                            "flow_mass_phase_comp",
                        ):
                            flow_phase_comp = True
                            local_name_flow = "flow_mass_phase_comp"
                        else:
                            raise PropertyNotSupportedError(
                                "No mass_frac_phase_comp or flow_mass_phase or"
                                " flow_mass_phase_comp variables encountered "
                                "while building ports for the condenser."
                            )
                    else:
                        raise PropertyNotSupportedError(
                            "No mass frac or mole frac variables encountered "
                            " while building ports for the condenser. "
                            "phase_frac as a state variable is not "
                            "supported with distillation unit models."
                        )

                    # Rule for liquid phase mole fraction
                    def rule_liq_frac(self, t, i):
                        if not flow_phase_comp:
                            sum_flow_comp = sum(
                                self.control_volume.properties_out[t].component(
                                    local_name_frac
                                )[p, i]
                                * self.control_volume.properties_out[t].component(
                                    local_name_flow
                                )[p]
                                for p in self._liquid_set
                            )

                            return sum_flow_comp / sum(
                                self.control_volume.properties_out[t].component(
                                    local_name_flow
                                )[p]
                                for p in self._liquid_set
                            )
                        else:
                            sum_flow_comp = sum(
                                self.control_volume.properties_out[t].component(
                                    local_name_flow
                                )[p, i]
                                for p in self._liquid_set
                            )

                            return sum_flow_comp / sum(
                                self.control_volume.properties_out[t].component(
                                    local_name_flow
                                )[p, i]
                                for p in self._liquid_set
                                for i in self.config.property_package.component_list
                            )

                    self.e_liq_frac = Expression(
                        self.flowsheet().time, index_set, rule=rule_liq_frac
                    )

                    # Rule for vapor phase mass/mole fraction
                    def rule_vap_frac(self, t, i):
                        if not flow_phase_comp:
                            sum_flow_comp = sum(
                                self.control_volume.properties_out[t].component(
                                    local_name_frac
                                )[p, i]
                                * self.control_volume.properties_out[t].component(
                                    local_name_flow
                                )[p]
                                for p in self._vapor_set
                            )
                            return sum_flow_comp / sum(
                                self.control_volume.properties_out[t].component(
                                    local_name_flow
                                )[p]
                                for p in self._vapor_set
                            )
                        else:
                            sum_flow_comp = sum(
                                self.control_volume.properties_out[t].component(
                                    local_name_flow
                                )[p, i]
                                for p in self._vapor_set
                            )

                            return sum_flow_comp / sum(
                                self.control_volume.properties_out[t].component(
                                    local_name_flow
                                )[p, i]
                                for p in self._vapor_set
                                for i in self.config.property_package.component_list
                            )

                    self.e_vap_frac = Expression(
                        self.flowsheet().time, index_set, rule=rule_vap_frac
                    )

                    # add the reference and variable name to the reflux port
                    self.reflux.add(self.e_liq_frac, k)

                    # add the reference and variable name to the
                    # distillate port
                    self.distillate.add(self.e_liq_frac, k)

                    # add the reference and variable name to the
                    # vapor port
                    self.vapor_outlet.add(self.e_vap_frac, k)
                else:

                    # Assumes mole_frac_phase or mass_frac_phase exist as
                    # state vars in the port and therefore access directly
                    # from the state block.
                    var = self.control_volume.properties_out[:].component(local_name)[
                        ...
                    ]

                    # add the reference and variable name to the reflux port
                    self.reflux.add(Reference(var), k)

                    # add the reference and variable name to the distillate port
                    self.distillate.add(Reference(var), k)
            elif "flow" in local_name:
                if "phase" not in local_name:

                    # Assumes that here the var is total flow or component
                    # flow. However, need to extract the flow by phase from
                    # the state block. Expects to find the var
                    # flow_mol_phase or flow_mass_phase in the state block.

                    # Check if it is not indexed by component list and this
                    # is total flow
                    if not member_list[k].is_indexed():
                        # if state var is not flow_mol/flow_mass by phase
                        local_name_flow = local_name + "_phase"

                        # Rule for vap phase flow
                        def rule_vap_flow(self, t):
                            return sum(
                                self.control_volume.properties_out[t].component(
                                    local_name_flow
                                )[p]
                                for p in self._vapor_set
                            )

                        self.e_vap_flow = Expression(
                            self.flowsheet().time, rule=rule_vap_flow
                        )

                        # Rule to link the liq phase flow to the reflux
                        def rule_reflux_flow(self, t):
                            return sum(
                                self.control_volume.properties_out[t].component(
                                    local_name_flow
                                )[p]
                                for p in self._liquid_set
                            ) * (self.reflux_ratio / (1 + self.reflux_ratio))

                        self.e_reflux_flow = Expression(
                            self.flowsheet().time, rule=rule_reflux_flow
                        )

                        # Rule to link the liq flow to the distillate
                        def rule_distillate_flow(self, t):
                            return sum(
                                self.control_volume.properties_out[t].component(
                                    local_name_flow
                                )[p]
                                for p in self._liquid_set
                            ) / (1 + self.reflux_ratio)

                        self.e_distillate_flow = Expression(
                            self.flowsheet().time, rule=rule_distillate_flow
                        )

                    else:
                        # when it is flow comp indexed by component list
                        str_split = local_name.split("_")
                        if len(str_split) == 3 and str_split[-1] == "comp":
                            local_name_flow = (
                                str_split[0] + "_" + str_split[1] + "_phase_" + "comp"
                            )

                        # Get the indexing set i.e. component list
                        index_set = member_list[k].index_set()

                        # Rule for vap phase flow to the vapor outlet
                        def rule_vap_flow(self, t, i):
                            return sum(
                                self.control_volume.properties_out[t].component(
                                    local_name_flow
                                )[p, i]
                                for p in self._vapor_set
                            )

                        self.e_vap_flow = Expression(
                            self.flowsheet().time, index_set, rule=rule_vap_flow
                        )

                        # Rule to link the liq flow to the reflux
                        def rule_reflux_flow(self, t, i):
                            return sum(
                                self.control_volume.properties_out[t].component(
                                    local_name_flow
                                )[p, i]
                                for p in self._liquid_set
                            ) * (self.reflux_ratio / (1 + self.reflux_ratio))

                        self.e_reflux_flow = Expression(
                            self.flowsheet().time, index_set, rule=rule_reflux_flow
                        )

                        # Rule to link the liq flow to the distillate
                        def rule_distillate_flow(self, t, i):
                            return sum(
                                self.control_volume.properties_out[t].component(
                                    local_name_flow
                                )[p, i]
                                for p in self._liquid_set
                            ) / (1 + self.reflux_ratio)

                        self.e_distillate_flow = Expression(
                            self.flowsheet().time, index_set, rule=rule_distillate_flow
                        )

                    # add the reference and variable name to the reflux port
                    self.reflux.add(self.e_reflux_flow, k)

                    # add the reference and variable name to the
                    # distillate port
                    self.distillate.add(self.e_distillate_flow, k)

                    # add the reference and variable name to the
                    # distillate port
                    self.vapor_outlet.add(self.e_vap_flow, k)
            elif "enth" in local_name:
                if "phase" not in local_name:
                    # assumes total mixture enthalpy (enth_mol or enth_mass)
                    # and hence should not be indexed by phase
                    if not member_list[k].is_indexed():
                        # if state var is not enth_mol/enth_mass
                        # by phase, add _phase string to extract the right
                        # value from the state block
                        local_name_enth = local_name + "_phase"
                    else:
                        raise PropertyPackageError(
                            "Enthalpy is indexed but the variable "
                            "name does not reflect the presence of an index. "
                            "Please follow the naming convention outlined "
                            "in the documentation for state variables."
                        )

                    # NOTE:pass phase index when generating expression only
                    # when multiple liquid or vapor phases detected
                    # else ensure consistency with state vars and do not
                    # add phase index to the port members. Hence, the check
                    # for length of local liq and vap phase sets.

                    # Rule for vap enthalpy. Setting the enthalpy to the
                    # enth_mol_phase['Vap'] value from the state block
                    def rule_vap_enth(self, t):
                        return sum(
                            self.control_volume.properties_out[t].component(
                                local_name_enth
                            )[p]
                            for p in self._vapor_set
                        )

                    self.e_vap_enth = Expression(
                        self.flowsheet().time, rule=rule_vap_enth
                    )

                    # Rule to link the liq enthalpy to the reflux.
                    # Setting the enthalpy to the
                    # enth_mol_phase['Liq'] value from the state block
                    def rule_reflux_enth(self, t):
                        return sum(
                            self.control_volume.properties_out[t].component(
                                local_name_enth
                            )[p]
                            for p in self._liquid_set
                        )

                    self.e_reflux_enth = Expression(
                        self.flowsheet().time, rule=rule_reflux_enth
                    )

                    # Rule to link the liq flow to the distillate.
                    # Setting the enthalpy to the
                    # enth_mol_phase['Liq'] value from the state block
                    def rule_distillate_enth(self, t):
                        return sum(
                            self.control_volume.properties_out[t].component(
                                local_name_enth
                            )[p]
                            for p in self._liquid_set
                        )

                    self.e_distillate_enth = Expression(
                        self.flowsheet().time, rule=rule_distillate_enth
                    )

                    # add the reference and variable name to the reflux port
                    self.reflux.add(self.e_reflux_enth, k)

                    # add the reference and variable name to the
                    # distillate port
                    self.distillate.add(self.e_distillate_enth, k)

                    # add the reference and variable name to the
                    # distillate port
                    self.vapor_outlet.add(self.e_vap_enth, k)
                elif "phase" in local_name:
                    # assumes enth_mol_phase or enth_mass_phase.
                    # This is an intensive property, you create a direct
                    # reference irrespective of the reflux, distillate and
                    # vap_outlet

                    # Rule for vap flow
                    if not member_list[k].is_indexed():
                        var = self.control_volume.properties_out[:].component(
                            local_name
                        )
                    else:
                        var = self.control_volume.properties_out[:].component(
                            local_name
                        )[...]

                    # add the reference and variable name to the reflux port
                    self.reflux.add(Reference(var), k)

                    # add the reference and variable name to the distillate port
                    self.distillate.add(Reference(var), k)

                    # add the reference and variable name to the
                    # vapor outlet port
                    self.vapor_outlet.add(Reference(var), k)
                else:
                    raise PropertyNotSupportedError(
                        "Unrecognized enthalpy state variable encountered "
                        "while building ports for the condenser. Only total "
                        "mixture enthalpy or enthalpy by phase are supported."
                    )

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
