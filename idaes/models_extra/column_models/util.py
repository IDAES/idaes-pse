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
Utility functions for column models
"""

__author__ = "Jaffer Ghouse"

import idaes.logger as idaeslog

# Import Pyomo libraries
from pyomo.network import Port
from pyomo.environ import (
    Reference,
    Expression,
)

# Import IDAES cores
from idaes.core.util.exceptions import (
    PropertyPackageError,
    PropertyNotSupportedError,
)

_log = idaeslog.getLogger(__name__)


def make_phase_split(
    self,
    port=None,
    phase=None,
    has_liquid_side_draw=False,
    has_vapor_side_draw=False,
    side_sf=None,
    equipmentType=None,
):
    """Method to split and populate the outlet ports with corresponding
    phase values from the mixed stream outlet block."""

    member_list = self.properties_out[0].define_port_members()

    for k in member_list:

        local_name = member_list[k].local_name

        # Create references and populate the intensive variables
        if (
            "flow" not in local_name
            and "frac" not in local_name
            and "enth" not in local_name
        ):
            if not member_list[k].is_indexed():
                var = self.properties_out[:].component(local_name)
            else:
                var = self.properties_out[:].component(local_name)[...]

            # add the reference and variable name to the port
            port.add(Reference(var), k)

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
                        self.properties_out[0], "mole_frac_phase_comp"
                    ) and hasattr(self.properties_out[0], "flow_mol_phase"):
                        flow_phase_comp = False
                        local_name_frac = "mole_frac_phase_comp"
                        local_name_flow = "flow_mol_phase"
                    elif hasattr(self.properties_out[0], "flow_mol_phase_comp"):
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
                        self.properties_out[0], "mass_frac_phase_comp"
                    ) and hasattr(self.properties_out[0], "flow_mass_phase"):
                        flow_phase_comp = False
                        local_name_frac = "mass_frac_phase_comp"
                        local_name_flow = "flow_mass_phase"
                    elif hasattr(self.properties_out[0], "flow_mass_phase_comp"):
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

                # Rule for mole fraction
                def rule_mole_frac(self, t, i):
                    if not flow_phase_comp:
                        sum_flow_comp = sum(
                            self.properties_out[t].component(local_name_frac)[p, i]
                            * self.properties_out[t].component(local_name_flow)[p]
                            for p in phase
                        )

                        return sum_flow_comp / sum(
                            self.properties_out[t].component(local_name_flow)[p]
                            for p in phase
                        )

                    else:
                        sum_flow_comp = sum(
                            self.properties_out[t].component(local_name_flow)[p, i]
                            for p in phase
                        )

                        return sum_flow_comp / sum(
                            self.properties_out[t].component(local_name_flow)[p, i]
                            for p in phase
                            for i in self.config.property_package.component_list
                        )

                # add the reference and variable name to the port
                expr = Expression(self.flowsheet().time, index_set, rule=rule_mole_frac)
                self.add_component("e_mole_frac_" + port.local_name, expr)
                port.add(expr, k)
            else:

                # Assumes mole_frac_phase or mass_frac_phase exist as
                # state vars in the port and therefore access directly
                # from the state block.
                var = self.properties_out[:].component(local_name)[...]

                # add the reference and variable name to the port
                ref = Reference(var)
                setattr(self, "_" + k + "_" + port.local_name + "_ref", ref)
                port.add(ref, k)
        elif "flow" in local_name:
            if "phase" not in local_name:

                # Assumes that here the var is total flow or component
                # flow. However, need to extract the flow by phase from
                # the state block. Expects to find the var
                # flow_mol_phase or flow_mass_phase in the state block.

                # Check if it is not indexed by component list and this
                # is total flow
                if "total" in str(equipmentType):
                    if not member_list[k].is_indexed():

                        def rule_flow(self, t):
                            return (self.properties_out[t].component(local_name)) * (
                                side_sf
                            )

                        # add the reference and variable name to the port
                        expr = Expression(self.flowsheet().time, rule=rule_flow)
                        self.add_component("e_flow_" + port.local_name, expr)
                        port.add(expr, k)
                    else:
                        # Create references and populate the extensive variables
                        # This is for vars that are indexed by phase, comp or both.
                        index_set = member_list[k].index_set()

                        # Rule to link the flow to the port
                        def rule_flow(self, t, i):
                            return (self.properties_out[t].component(local_name)[i]) * (
                                side_sf
                            )

                        # add the reference and variable name to the port
                        expr = Expression(
                            self.flowsheet().time, index_set, rule=rule_flow
                        )
                        self.add_component("e_flow_" + port.local_name, expr)
                        port.add(expr, k)

                else:
                    if not member_list[k].is_indexed():
                        # if state var is not flow_mol/flow_mass
                        # by phase
                        local_name_flow = local_name + "_phase"

                        # Rule to link the flow to the port
                        def rule_flow(self, t):
                            return sum(
                                self.properties_out[t].component(local_name_flow)[p]
                                for p in phase
                            ) * (side_sf)

                        # add the reference and variable name to the port
                        expr = Expression(self.flowsheet().time, rule=rule_flow)
                        self.add_component("e_flow_" + port.local_name, expr)
                        port.add(expr, k)
                    else:
                        # when it is flow comp indexed by component list
                        str_split = local_name.split("_")
                        if len(str_split) == 3 and str_split[-1] == "comp":
                            local_name_flow = (
                                str_split[0] + "_" + str_split[1] + "_phase_" + "comp"
                            )

                        # Get the indexing set i.e. component list
                        index_set = member_list[k].index_set()

                        # Rule to link the flow to the port
                        def rule_flow(self, t, i):
                            return sum(
                                self.properties_out[t].component(local_name_flow)[p, i]
                                for p in phase
                            ) * (side_sf)

                        expr = Expression(
                            self.flowsheet().time, index_set, rule=rule_flow
                        )
                        self.add_component("e_flow_" + port.local_name, expr)
                        port.add(expr, local_name)
            elif "phase" in local_name:
                # flow is indexed by phase and comp
                # Get the indexing sets i.e. component list and phase list
                component_set = self.config.property_package.component_list

                phase_set = self.config.property_package.phase_list

                # If statement to skip in case equimpment is not a Tray
                if equipmentType == None:

                    def rule_flow(self, t, p, i):
                        if (phase is self._liquid_set and p in self._liquid_set) or (
                            phase is self._vapor_set and p in self._vapor_set
                        ):
                            # pass appropriate phase flow values to port
                            return (
                                self.properties_out[t].component(local_name)[p, i]
                            ) * (side_sf)

                        else:
                            # return small number for phase that should not
                            # be in the appropriate port. For example,
                            # the state vars will be flow_mol_phase_comp
                            # which will include all phases. The liq port
                            # should have the correct references to the liq
                            # phase flow but the vapor phase flow should be 0.
                            return 1e-8

                else:

                    def rule_flow(self, t, p, i):
                        return (self.properties_out[t].component(local_name)[p, i]) * (
                            side_sf
                        )

                expr = Expression(
                    self.flowsheet().time, phase_set, component_set, rule=rule_flow
                )

                self.add_component("e_" + local_name + port.local_name, expr)
                port.add(expr, k)
            else:
                raise PropertyPackageError(
                    "Unrecognized flow state variable encountered "
                    "while building ports for the tray. Please follow "
                    "the naming convention outlined in the documentation "
                    "for state variables."
                )
        elif "enth" in local_name:
            if "phase" not in local_name:
                # assumes total mixture enthalpy (enth_mol or enth_mass)
                if not member_list[k].is_indexed():
                    # if state var is not enth_mol/enth_mass
                    # by phase, add _phase string to extract the right
                    # value from the state block
                    local_name_phase = local_name + "_phase"
                else:
                    raise PropertyPackageError(
                        "Enthalpy is indexed but the variable "
                        "name does not reflect the presence of an index. "
                        "Please follow the naming convention outlined "
                        "in the documentation for state variables."
                    )

                # Rule to link the phase enthalpy to the port.
                def rule_enth(self, t):
                    return sum(
                        self.properties_out[t].component(local_name_phase)[p]
                        for p in phase
                    )

                expr = Expression(self.flowsheet().time, rule=rule_enth)
                self.add_component("e_enth_" + port.local_name, expr)
                # add the reference and variable name to the port
                port.add(expr, k)

            elif "phase" in local_name:
                # assumes enth_mol_phase or enth_mass_phase.
                # This is an intensive property, you create a direct
                # reference irrespective of the reflux, distillate and
                # vap_outlet

                if not member_list[k].is_indexed():
                    var = self.properties_out[:].component(local_name)
                else:
                    var = self.properties_out[:].component(local_name)[...]

                # add the reference and variable name to the port
                ref = Reference(var)
                setattr(self, "_" + k + "_" + port.local_name + "_ref", ref)
                port.add(ref, k)
            else:
                raise PropertyNotSupportedError(
                    "Unrecognized enthalpy state variable encountered "
                    "while building ports for the tray. Only total "
                    "mixture enthalpy or enthalpy by phase are supported."
                )
