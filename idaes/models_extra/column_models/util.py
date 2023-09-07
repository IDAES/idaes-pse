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
Utility functions for column models
"""

__author__ = "Jaffer Ghouse, Alejandro Garcia-Diego"

# TODO: look into protected access issues - probably need to refactor
# pylint: disable=protected-access

from functools import partial

# Import Pyomo libraries
from pyomo.environ import (
    Reference,
    Expression,
)

# Import IDAES cores
from idaes.core.util.exceptions import (
    PropertyPackageError,
    PropertyNotSupportedError,
)
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def make_phase_split(
    model,
    port=None,
    phase=None,
    side_sf=None,
    equipmentType=None,
):
    """Method to split and populate the outlet ports with corresponding
    phase values from the mixed stream outlet block."""
    t0 = model.flowsheet().time.first()

    member_list = model.properties_out[t0].define_port_members()

    for k in member_list:

        local_name = member_list[k].local_name

        # Create references and populate the intensive variables
        if (
            "flow" not in local_name
            and "frac" not in local_name
            and "enth" not in local_name
        ):
            if not member_list[k].is_indexed():
                var = model.properties_out[:].component(local_name)
            else:
                var = model.properties_out[:].component(local_name)[...]

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
                        model.properties_out[t0], "mole_frac_phase_comp"
                    ) and hasattr(model.properties_out[t0], "flow_mol_phase"):
                        flow_phase_comp = False
                        local_name_frac = "mole_frac_phase_comp"
                        local_name_flow = "flow_mol_phase"
                    elif hasattr(model.properties_out[t0], "flow_mol_phase_comp"):
                        flow_phase_comp = True
                        local_name_flow = "flow_mol_phase_comp"
                    else:
                        raise PropertyNotSupportedError(
                            "No mole_frac_phase_comp or flow_mol_phase or"
                            " flow_mol_phase_comp variables encountered "
                            "while building ports. "
                        )
                elif "mass" in local_name:
                    if hasattr(
                        model.properties_out[t0], "mass_frac_phase_comp"
                    ) and hasattr(model.properties_out[t0], "flow_mass_phase"):
                        flow_phase_comp = False
                        local_name_frac = "mass_frac_phase_comp"
                        local_name_flow = "flow_mass_phase"
                    elif hasattr(model.properties_out[t0], "flow_mass_phase_comp"):
                        flow_phase_comp = True
                        local_name_flow = "flow_mass_phase_comp"
                    else:
                        raise PropertyNotSupportedError(
                            "No mass_frac_phase_comp or flow_mass_phase or"
                            " flow_mass_phase_comp variables encountered "
                            "while building ports."
                        )
                else:
                    raise PropertyNotSupportedError(
                        "No mass frac or mole frac variables encountered "
                        " while building ports. "
                        "phase_frac as a state variable is not "
                        "supported with distillation unit models."
                    )

                if not flow_phase_comp:
                    tmp_rule = partial(
                        _rule_mole_frac_0,
                        local_name_frac=local_name_frac,
                        local_name_flow=local_name_flow,
                        phase=phase,
                    )
                else:
                    tmp_rule = partial(
                        _rule_mole_frac_1,
                        local_name_flow=local_name_flow,
                        phase=phase,
                    )
                # add the reference and variable name to the port
                expr = Expression(model.flowsheet().time, index_set, rule=tmp_rule)
                model.add_component("e_mole_frac_" + port.local_name, expr)
                port.add(expr, k)
            else:

                # Assumes mole_frac_phase or mass_frac_phase exist as
                # state vars in the port and therefore access directly
                # from the state block.
                var = model.properties_out[:].component(local_name)[...]

                # add the reference and variable name to the port
                ref = Reference(var)
                setattr(model, "_" + k + "_" + port.local_name + "_ref", ref)
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

                        # add the reference and variable name to the port
                        expr = Expression(
                            model.flowsheet().time,
                            rule=partial(
                                _rule_flow_0, local_name=local_name, side_sf=side_sf
                            ),
                        )
                        model.add_component("e_flow_" + port.local_name, expr)
                        port.add(expr, k)
                    else:
                        # Create references and populate the extensive variables
                        # This is for vars that are indexed by phase, comp or both.
                        index_set = member_list[k].index_set()

                        # add the reference and variable name to the port
                        expr = Expression(
                            model.flowsheet().time,
                            index_set,
                            rule=partial(
                                _rule_flow_1, local_name=local_name, side_sf=side_sf
                            ),
                        )
                        model.add_component("e_flow_" + port.local_name, expr)
                        port.add(expr, k)

                else:
                    if not member_list[k].is_indexed():
                        # if state var is not flow_mol/flow_mass
                        # by phase
                        local_name_flow = local_name + "_phase"

                        # add the reference and variable name to the port
                        expr = Expression(
                            model.flowsheet().time,
                            rule=partial(
                                _rule_flow_2,
                                local_name_flow=local_name_flow,
                                side_sf=side_sf,
                                phase=phase,
                            ),
                        )
                        model.add_component("e_flow_" + port.local_name, expr)
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

                        expr = Expression(
                            model.flowsheet().time,
                            index_set,
                            rule=partial(
                                _rule_flow_3,
                                local_name_flow=local_name_flow,
                                side_sf=side_sf,
                                phase=phase,
                            ),
                        )
                        model.add_component("e_flow_" + port.local_name, expr)
                        port.add(expr, local_name)
            else:
                # flow is indexed by phase and comp
                # Get the indexing sets i.e. component list and phase list
                component_set = model.config.property_package.component_list

                phase_set = model.config.property_package.phase_list

                expr = Expression(
                    model.flowsheet().time,
                    phase_set,
                    component_set,
                    rule=partial(
                        _rule_flow_4,
                        local_name=local_name,
                        side_sf=side_sf,
                        phase=phase,
                        equipmentType=equipmentType,
                    ),
                )

                model.add_component("e_" + local_name + port.local_name, expr)
                port.add(expr, k)

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

                expr = Expression(
                    model.flowsheet().time,
                    rule=partial(
                        _rule_enth_0, local_name_phase=local_name_phase, phase=phase
                    ),
                )
                model.add_component("e_enth_" + port.local_name, expr)
                # add the reference and variable name to the port
                port.add(expr, k)

            else:
                # assumes enth_mol_phase or enth_mass_phase.
                # This is an intensive property, you create a direct
                # reference irrespective of the reflux, distillate and
                # vap_outlet

                if not member_list[k].is_indexed():
                    var = model.properties_out[:].component(local_name)
                else:
                    var = model.properties_out[:].component(local_name)[...]

                # add the reference and variable name to the port
                ref = Reference(var)
                setattr(model, "_" + k + "_" + port.local_name + "_ref", ref)
                port.add(ref, k)


# Rule for mole fraction
def _rule_mole_frac_0(model, t, i, local_name_frac, local_name_flow, phase):
    sum_flow_comp = sum(
        model.properties_out[t].component(local_name_frac)[p, i]
        * model.properties_out[t].component(local_name_flow)[p]
        for p in phase
    )

    return sum_flow_comp / sum(
        model.properties_out[t].component(local_name_flow)[p] for p in phase
    )


# Rule for mole fraction
def _rule_mole_frac_1(model, t, i, local_name_flow, phase):
    sum_flow_comp = sum(
        model.properties_out[t].component(local_name_flow)[p, i] for p in phase
    )

    return sum_flow_comp / sum(
        model.properties_out[t].component(local_name_flow)[p, i]
        for p in phase
        for i in model.config.property_package.component_list
    )


def _rule_flow_0(model, t, local_name, side_sf):
    return (model.properties_out[t].component(local_name)) * (side_sf)


# Rule to link the flow to the port
def _rule_flow_1(model, t, i, local_name, side_sf):
    return (model.properties_out[t].component(local_name)[i]) * (side_sf)


# Rule to link the flow to the port
def _rule_flow_2(model, t, local_name_flow, side_sf, phase):
    return sum(model.properties_out[t].component(local_name_flow)[p] for p in phase) * (
        side_sf
    )


# Rule to link the flow to the port
def _rule_flow_3(model, t, i, local_name_flow, side_sf, phase):
    return sum(
        model.properties_out[t].component(local_name_flow)[p, i] for p in phase
    ) * (side_sf)


def _rule_flow_4(model, t, p, i, local_name, side_sf, phase, equipmentType):
    # If statement to skip in case equipment is not a Tray
    if equipmentType is None:
        if (phase is model._liquid_set and p in model._liquid_set) or (
            phase is model._vapor_set and p in model._vapor_set
        ):
            # pass appropriate phase flow values to port
            return (model.properties_out[t].component(local_name)[p, i]) * (side_sf)

        else:
            # return small number for phase that should not
            # be in the appropriate port. For example,
            # the state vars will be flow_mol_phase_comp
            # which will include all phases. The liq port
            # should have the correct references to the liq
            # phase flow but the vapor phase flow should be 0.
            # TODO this should be a param and probably also have units
            return 1e-8
    else:
        return (model.properties_out[t].component(local_name)[p, i]) * (side_sf)


# Rule to link the phase enthalpy to the port.
def _rule_enth_0(model, t, local_name_phase, phase):
    return sum(model.properties_out[t].component(local_name_phase)[p] for p in phase)
