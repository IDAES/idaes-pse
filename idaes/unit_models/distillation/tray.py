##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Tray model for distillation.
"""

__author__ = "Jaffer Ghouse"

import logging

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.network import Port
from pyomo.environ import Reference, Expression, Var, Constraint, \
    TerminationCondition, value

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        MaterialBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import PropertyPackageError

_log = logging.getLogger(__name__)


@declare_process_block_class("Tray")
class TrayData(UnitModelBlockData):
    """
    Tray unit for distillation model.
    """
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare("is_top_tray", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Flag to indicate top tray.",
        doc="""Indicates if this is a top tray and constructs
corresponding ports,
**default** - False.
**Valid values:** {
**True** - top tray,
**False** - conventional tray}"""))
    CONFIG.declare("is_bottom_tray", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Flag to indicate bottom tray.",
        doc="""Indicates if this is a bottom tray and constructs
corresponding ports,
**default** - False.
**Valid values:** {
**True** - bottom tray,
**False** - conventional tray}"""))
    CONFIG.declare("has_liquid_side_draw", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Liquid side draw construction flag.",
        doc="""Indicates if there is a liquid side draw from the tray,
**default** - False.
**Valid values:** {
**True** - include a liquid side draw from the tray,
**False** - exclude a liquid side draw from the tray.}"""))
    CONFIG.declare("has_vapor_side_draw", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Vapor side draw construction flag.",
        doc="""Indicates if there is a vapor side draw from the tray,
**default** - False.
**Valid values:** {
**True** - include a vapor side draw from the tray,
**False** - exclude a vapor side draw from the tray.}"""))
    CONFIG.declare("has_heat_transfer", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Heat transfer to/from tray construction flag.",
        doc="""Indicates if there is heat transfer to/from the tray,
**default** - False.
**Valid values:** {
**True** - include a heat transfer term,
**False** - exclude a heat transfer term.}"""))
    CONFIG.declare("has_pressure_change", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
    constructed,
    **default** - False.
    **Valid values:** {
    **True** - include pressure change terms,
    **False** - exclude pressure change terms.}"""))
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))

    def build(self):
        """Build the model.

        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(TrayData, self).build()

        # Create the inlets for the tray
        inlet_list = ["liq", "vap"]

        # Setup StateBlock argument dict
        state_block_args = dict(**self.config.property_package_args)
        state_block_args["has_phase_equilibrium"] = False
        state_block_args["parameters"] = self.config.property_package
        state_block_args["defined_state"] = True

        for i in inlet_list:
            state_obj = self.config.property_package.state_block_class(
                self.flowsheet().config.time,
                doc="State block for inlet to tray",
                default=state_block_args)

            setattr(self, "properties_in_" + i, state_obj)

        # add mixed outlet state blocks which is the feed to the tray

        # Setup StateBlock argument dict
        mixed_block_args = dict(**self.config.property_package_args)
        mixed_block_args["has_phase_equilibrium"] = True
        mixed_block_args["parameters"] = self.config.property_package
        mixed_block_args["defined_state"] = False

        self.properties_out = self.config.property_package.\
            state_block_class(self.flowsheet().config.time,
                              doc="State block for inlet to tray",
                              default=mixed_block_args)

        self._add_material_balance()
        self._add_energy_balance()

        self._add_pressure_balance()
        self._add_ports()

    def _add_material_balance(self):

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.component_list,
                         doc="Material mixing equations")
        def material_mixing_equations(b, t, j):
            return 0 == sum(
                self.properties_in_liq[t].get_material_flow_terms(p, j) +
                self.properties_in_vap[t].get_material_flow_terms(p, j) -
                self.properties_out[t].get_material_flow_terms(p, j)
                for p in b.config.property_package.phase_list)

    def _add_energy_balance(self):

        if self.config.has_heat_transfer:
            self.heat_duty = Var(initialize=0,
                                 doc="heat duty for the tray")

        @self.Constraint(self.flowsheet().config.time, doc="Energy balances")
        def enthalpy_mixing_equations(b, t):
            if self.config.has_heat_transfer:
                return 0 == (
                    sum(self.properties_in_liq[t].get_enthalpy_flow_terms(p)
                        for p in b.config.property_package.phase_list) +
                    sum(self.properties_in_vap[t].get_enthalpy_flow_terms(p)
                        for p in b.config.property_package.phase_list) -
                    sum(self.properties_out[t].get_enthalpy_flow_terms(p)
                        for p in b.config.property_package.phase_list)) + \
                    self.heat_duty
            else:
                return 0 == (
                    sum(self.properties_in_liq[t].get_enthalpy_flow_terms(p)
                        for p in b.config.property_package.phase_list) +
                    sum(self.properties_in_vap[t].get_enthalpy_flow_terms(p)
                        for p in b.config.property_package.phase_list) -
                    sum(self.properties_out[t].get_enthalpy_flow_terms(p)
                        for p in b.config.property_package.phase_list))

    def _add_pressure_balance(self):
        if self.config.has_pressure_change:
            self.deltaP = Var(initialize=0,
                              doc="pressure drop across tray")

        @self.Constraint(self.flowsheet().config.time,
                         doc="Pressure drop constraint for tray")
        def pressure_drop_equation(self, t):
            if self.config.has_pressure_change:
                return self.properties_out[t].pressure == \
                    self.properties_in_vap[t].pressure - self.deltaP
            else:
                return self.properties_out[t].pressure == \
                    self.properties_in_vap[t].pressure

    def _add_ports(self):

        if self.config.has_liquid_side_draw:
            self.liq_side_draw = Port(noruleinit=True, doc="liquid side draw.")
            self._make_liq_side_split()
        if self.config.has_vapor_side_draw:
            self.vap_side_draw = Port(noruleinit=True, doc="vapor side draw.")
            self._make_vap_side_split()
        if self.config.is_top_tray and not self.config.is_bottom_tray:
            self.vap_out = Port(noruleinit=True, doc="vapor from top tray.")
            self._make_vap_out_top_tray()
        elif self.config.is_bottom_tray and not self.config.is_top_tray:
            self.liq_out = Port(noruleinit=True, doc="liquid from bottom tray")
            self._make_liq_out_bottom_tray()
        elif self.config.is_top_tray and self.config.is_bottom_tray:
            raise Exception("Configuration Error")

    def _make_liq_side_split(self):

        self.liq_side_sf = Var(initialize=0.01,
                               doc="split fraction for the liquid side draw")

        member_list = self.properties_in_liq[0].define_port_members()

        for k in member_list:
            # Create references and populate the intensive variables
            if "flow" not in k and "frac" not in k and "enth" not in k:
                if not member_list[k].is_indexed():
                    var = self.properties_out[:].\
                        component(member_list[k].local_name)
                else:
                    var = self.properties_out[:].\
                        component(member_list[k].local_name)[...]

                # add the reference and variable name to the liq_side_draw port
                self.liq_side_draw.add(Reference(var), k)

            elif "frac" in k and ("mole" in k or "mass" in k):

                # Mole/mass frac is typically indexed
                index_set = member_list[k].index_set()

                # if state var is not mole/mass frac by phase
                if "phase" not in k:
                    # Assuming the state block has the var
                    # "mole_frac_phase_comp". Valid if VLE is supported
                    # Create a string "mole_frac_phase_comp" or
                    # "mass_frac_phase_comp". Cannot directly append phase
                    # to k as the naming convention is phase followed
                    # by comp
                    str_split = k.split('_')
                    local_name = '_'.join(str_split[0:2]) + \
                        "_phase" + "_" + str_split[2]

                    # Rule for liquid fraction
                    def rule_liq_frac(self, t, i):
                        return self.properties_out[t].\
                            component(local_name)["Liq", i]
                    self.e_liq_frac = Expression(
                        self.flowsheet().time, index_set,
                        rule=rule_liq_frac)

                    # add the reference and variable name to the liq_side_draw port
                    self.liq_side_draw.add(self.e_liq_frac, k)

                else:

                    # Assumes mole_frac_phase or mass_frac_phase exist as
                    # state vars in the port and therefore access directly
                    # from the state block.
                    var = self.properties_out[:].\
                        component(member_list[k].local_name)[...]

                    # add the reference and variable name to the
                    # liq_side_draw port
                    self.liq_side_draw.add(Reference(var), k)
            elif "flow" in k:
                if "phase" not in k:

                    # Assumes that here the var is total flow or component
                    # flow. However, need to extract the flow by phase from
                    # the state block. Expects to find the var
                    # flow_mol_phase or flow_mass_phase in the state block.

                    # Check if it is not indexed by component list and this
                    # is total flow
                    if not member_list[k].is_indexed():
                        # if state var is not flow_mol/flow_mass
                        # by phase
                        local_name = str(member_list[k].local_name) + \
                            "_phase"

                        # Rule to link the liq flow to the liq_side_draw
                        def rule_liq_flow(self, t):
                            return self.properties_out[t].\
                                component(local_name)["Liq"] * \
                                (self.liq_side_sf)
                        self.e_liq_flow = Expression(
                            self.flowsheet().time,
                            rule=rule_liq_flow)
                    else:
                        # when it is flow comp indexed by component list
                        str_split = \
                            str(member_list[k].local_name).split("_")
                        if len(str_split) == 3 and str_split[-1] == "comp":
                            local_name = str_split[0] + "_" + \
                                str_split[1] + "_phase_" + "comp"

                        # Get the indexing set i.e. component list
                        index_set = member_list[k].index_set()

                        # Rule to link the liq flow to the liq_side_draw
                        def rule_liq_flow(self, t, i):
                            return self.properties_out[t].\
                                component(local_name)["Liq", i] * \
                                (self.liq_side_sf)
                        self.e_liq_flow = Expression(
                            self.flowsheet().time, index_set,
                            rule=rule_liq_flow)

                    # add the reference and variable name to the
                    # liq_side_draw port
                    self.liq_side_draw.add(self.e_liq_flow, k)
            elif "enth" in k:
                if "phase" not in k:
                    # assumes total mixture enthalpy (enth_mol or enth_mass)
                    if not member_list[k].is_indexed():
                        # if state var is not enth_mol/enth_mass
                        # by phase, add _phase string to extract the right
                        # value from the state block
                        local_name = str(member_list[k].local_name) + \
                            "_phase"
                    else:
                        raise PropertyPackageError(
                            "Expected an unindexed variable.")

                    # Rule to link the liq enthalpy to the reflux.
                    # Setting the enthalpy to the
                    # enth_mol_phase['Liq'] value from the state block
                    def rule_liq_enth(self, t):
                        return self.properties_out[t].\
                            component(local_name)["Liq"]
                    self.e_liq_enth = Expression(
                        self.flowsheet().time,
                        rule=rule_liq_enth)

                    # add the reference and variable name to the reflux port
                    self.liq_side_draw.add(self.e_liq_enth, k)
                elif "phase" in k:
                    # assumes enth_mol_phase or enth_mass_phase.
                    # This is an intensive property, you create a direct
                    # reference irrespective of the reflux, distillate and
                    # vap_outlet

                    # Rule for liq flow
                    if not member_list[k].is_indexed():
                        var = self.properties_out[:].\
                            component(member_list[k].local_name)
                    else:
                        var = self.properties_out[:].\
                            component(member_list[k].local_name)[...]

                    # add the reference and variable name to the reflux port
                    self.liq_side_draw.add(Reference(var), k)
                else:
                    raise Exception(
                        "Unrecognized enthalpy state variable. "
                        "Only total mixture enthalpy or enthalpy by "
                        "phase are supported.")

    def _make_vap_side_split(self):
        self.vap_side_sf = Var(initialize=0.01,
                               doc="split fraction for the vapor side draw")

        member_list = self.properties_in_vap[0].define_port_members()

        for k in member_list:
            # Create references and populate the intensive variables
            if "flow" not in k and "frac" not in k and "enth" not in k:
                if not member_list[k].is_indexed():
                    var = self.properties_out[:].\
                        component(member_list[k].local_name)
                else:
                    var = self.properties_out[:].\
                        component(member_list[k].local_name)[...]

                # add the reference and variable name to the liq_side_draw port
                self.vap_side_draw.add(Reference(var), k)

            elif "frac" in k and ("mole" in k or "mass" in k):

                # Mole/mass frac is typically indexed
                index_set = member_list[k].index_set()

                # if state var is not mole/mass frac by phase
                if "phase" not in k:
                    # Assuming the state block has the var
                    # "mole_frac_phase_comp". Valid if VLE is supported
                    # Create a string "mole_frac_phase_comp" or
                    # "mass_frac_phase_comp". Cannot directly append phase
                    # to k as the naming convention is phase followed
                    # by comp
                    str_split = k.split('_')
                    local_name = '_'.join(str_split[0:2]) + \
                        "_phase" + "_" + str_split[2]

                    # Rule for vapor fraction
                    def rule_vap_frac(self, t, i):
                        return self.properties_out[t].\
                            component(local_name)["Vap", i]
                    self.e_vap_frac = Expression(
                        self.flowsheet().time, index_set,
                        rule=rule_vap_frac)

                    # add the reference and variable name to the liq_side_draw port
                    self.vap_side_draw.add(self.e_vap_frac, k)

                else:

                    # Assumes mole_frac_phase or mass_frac_phase exist as
                    # state vars in the port and therefore access directly
                    # from the state block.
                    var = self.properties_out[:].\
                        component(member_list[k].local_name)[...]

                    # add the reference and variable name to the
                    # liq_side_draw port
                    self.vap_side_draw.add(Reference(var), k)
            elif "flow" in k:
                if "phase" not in k:

                    # Assumes that here the var is total flow or component
                    # flow. However, need to extract the flow by phase from
                    # the state block. Expects to find the var
                    # flow_mol_phase or flow_mass_phase in the state block.

                    # Check if it is not indexed by component list and this
                    # is total flow
                    if not member_list[k].is_indexed():
                        # if state var is not flow_mol/flow_mass
                        # by phase
                        local_name = str(member_list[k].local_name) + \
                            "_phase"

                        # Rule to link the liq flow to the liq_side_draw
                        def rule_vap_flow(self, t):
                            return self.properties_out[t].\
                                component(local_name)["Vap"] * \
                                (self.vap_side_sf)
                        self.e_vap_flow = Expression(
                            self.flowsheet().time,
                            rule=rule_vap_flow)
                    else:
                        # when it is flow comp indexed by component list
                        str_split = \
                            str(member_list[k].local_name).split("_")
                        if len(str_split) == 3 and str_split[-1] == "comp":
                            local_name = str_split[0] + "_" + \
                                str_split[1] + "_phase_" + "comp"

                        # Get the indexing set i.e. component list
                        index_set = member_list[k].index_set()

                        # Rule to link the liq flow to the liq_side_draw
                        def rule_vap_flow(self, t, i):
                            return self.properties_out[t].\
                                component(local_name)["Vap", i] * \
                                (self.vap_side_sf)
                        self.e_vap_flow = Expression(
                            self.flowsheet().time, index_set,
                            rule=rule_vap_flow)

                    # add the reference and variable name to the
                    # liq_side_draw port
                    self.vap_side_draw.add(self.e_vap_flow, k)
            elif "enth" in k:
                if "phase" not in k:
                    # assumes total mixture enthalpy (enth_mol or enth_mass)
                    if not member_list[k].is_indexed():
                        # if state var is not enth_mol/enth_mass
                        # by phase, add _phase string to extract the right
                        # value from the state block
                        local_name = str(member_list[k].local_name) + \
                            "_phase"
                    else:
                        raise PropertyPackageError(
                            "Expected an unindexed variable.")

                    # Rule to link the liq enthalpy to the reflux.
                    # Setting the enthalpy to the
                    # enth_mol_phase['Liq'] value from the state block
                    def rule_vap_enth(self, t):
                        return self.properties_out[t].\
                            component(local_name)["Vap"]
                    self.e_vap_enth = Expression(
                        self.flowsheet().time,
                        rule=rule_vap_enth)

                    # add the reference and variable name to the reflux port
                    self.vap_side_draw.add(self.e_vap_enth, k)
                elif "phase" in k:
                    # assumes enth_mol_phase or enth_mass_phase.
                    # This is an intensive property, you create a direct
                    # reference irrespective of the reflux, distillate and
                    # vap_outlet

                    # Rule for liq flow
                    if not member_list[k].is_indexed():
                        var = self.properties_out[:].\
                            component(member_list[k].local_name)
                    else:
                        var = self.properties_out[:].\
                            component(member_list[k].local_name)[...]

                    # add the reference and variable name to the reflux port
                    self.vap_side_draw.add(Reference(var), k)
                else:
                    raise Exception(
                        "Unrecognized enthalpy state variable. "
                        "Only total mixture enthalpy or enthalpy by "
                        "phase are supported.")

    def _make_vap_out_top_tray(self):
        if not self.config.has_vapor_side_draw:
            self.vap_side_sf = 0

        member_list = self.properties_in_vap[0].define_port_members()

        for k in member_list:
            # Create references and populate the intensive variables
            if "flow" not in k and "frac" not in k and "enth" not in k:
                if not member_list[k].is_indexed():
                    var = self.properties_out[:].\
                        component(member_list[k].local_name)
                else:
                    var = self.properties_out[:].\
                        component(member_list[k].local_name)[...]

                # add the reference and variable name to the liq_side_draw port
                self.vap_out.add(Reference(var), k)

            elif "frac" in k and ("mole" in k or "mass" in k):

                # Mole/mass frac is typically indexed
                index_set = member_list[k].index_set()

                # if state var is not mole/mass frac by phase
                if "phase" not in k:
                    # Assuming the state block has the var
                    # "mole_frac_phase_comp". Valid if VLE is supported
                    # Create a string "mole_frac_phase_comp" or
                    # "mass_frac_phase_comp". Cannot directly append phase
                    # to k as the naming convention is phase followed
                    # by comp
                    str_split = k.split('_')
                    local_name = '_'.join(str_split[0:2]) + \
                        "_phase" + "_" + str_split[2]

                    # Rule for vapor fraction
                    def rule_vap_frac(self, t, i):
                        return self.properties_out[t].\
                            component(local_name)["Vap", i]
                    self.e_vap_frac_top = Expression(
                        self.flowsheet().time, index_set,
                        rule=rule_vap_frac)

                    # add the reference and variable name to the liq_side_draw port
                    self.vap_out.add(self.e_vap_frac_top, k)

                else:

                    # Assumes mole_frac_phase or mass_frac_phase exist as
                    # state vars in the port and therefore access directly
                    # from the state block.
                    var = self.properties_out[:].\
                        component(member_list[k].local_name)[...]

                    # add the reference and variable name to the
                    # liq_side_draw port
                    self.vap_out.add(Reference(var), k)
            elif "flow" in k:
                if "phase" not in k:

                    # Assumes that here the var is total flow or component
                    # flow. However, need to extract the flow by phase from
                    # the state block. Expects to find the var
                    # flow_mol_phase or flow_mass_phase in the state block.

                    # Check if it is not indexed by component list and this
                    # is total flow
                    if not member_list[k].is_indexed():
                        # if state var is not flow_mol/flow_mass
                        # by phase
                        local_name = str(member_list[k].local_name) + \
                            "_phase"

                        # Rule to link the liq flow to the liq_side_draw
                        def rule_vap_flow(self, t):
                            return self.properties_out[t].\
                                component(local_name)["Vap"] * \
                                (1 - self.vap_side_sf)
                        self.e_vap_flow_top = Expression(
                            self.flowsheet().time,
                            rule=rule_vap_flow)
                    else:
                        # when it is flow comp indexed by component list
                        str_split = \
                            str(member_list[k].local_name).split("_")
                        if len(str_split) == 3 and str_split[-1] == "comp":
                            local_name = str_split[0] + "_" + \
                                str_split[1] + "_phase_" + "comp"

                        # Get the indexing set i.e. component list
                        index_set = member_list[k].index_set()

                        # Rule to link the liq flow to the liq_side_draw
                        def rule_vap_flow(self, t, i):
                            return self.properties_out[t].\
                                component(local_name)["Vap", i] * \
                                (1 - self.vap_side_sf)
                        self.e_vap_flow_top = Expression(
                            self.flowsheet().time, index_set,
                            rule=rule_vap_flow)

                    # add the reference and variable name to the
                    # liq_side_draw port
                    self.vap_out.add(self.e_vap_flow_top, k)
            elif "enth" in k:
                if "phase" not in k:
                    # assumes total mixture enthalpy (enth_mol or enth_mass)
                    if not member_list[k].is_indexed():
                        # if state var is not enth_mol/enth_mass
                        # by phase, add _phase string to extract the right
                        # value from the state block
                        local_name = str(member_list[k].local_name) + \
                            "_phase"
                    else:
                        raise PropertyPackageError(
                            "Expected an unindexed variable.")

                    # Rule to link the liq enthalpy to the reflux.
                    # Setting the enthalpy to the
                    # enth_mol_phase['Liq'] value from the state block
                    def rule_vap_enth(self, t):
                        return self.properties_out[t].\
                            component(local_name)["Vap"]
                    self.e_vap_enth_top = Expression(
                        self.flowsheet().time,
                        rule=rule_vap_enth)

                    # add the reference and variable name to the reflux port
                    self.vap_out.add(self.e_vap_enth_top, k)
                elif "phase" in k:
                    # assumes enth_mol_phase or enth_mass_phase.
                    # This is an intensive property, you create a direct
                    # reference irrespective of the reflux, distillate and
                    # vap_outlet

                    # Rule for liq flow
                    if not member_list[k].is_indexed():
                        var = self.properties_out[:].\
                            component(member_list[k].local_name)
                    else:
                        var = self.properties_out[:].\
                            component(member_list[k].local_name)[...]

                    # add the reference and variable name to the reflux port
                    self.vap_out.add(Reference(var), k)
                else:
                    raise Exception(
                        "Unrecognized enthalpy state variable. "
                        "Only total mixture enthalpy or enthalpy by "
                        "phase are supported.")

    def initialize(self, solver=None, outlvl=None):

        # Initialize the inlet state blocks
        self.properties_in_liq.initialize(outlvl=outlvl)
        self.properties_in_vap.initialize(outlvl=outlvl)

        # Initialize the mixed outlet state block
        self.properties_out.initialize(outlvl=outlvl)

        # Deactivate energy balance
        self.enthalpy_mixing_equations.deactivate()
        average_temperature = \
            0.5 * (self.properties_in_liq[0].temperature.value
                   + self.properties_in_vap[0].temperature.value)

        self.properties_out[:].temperature.fix(average_temperature)

        # Deactivate pressure balance
        self.pressure_drop_equation.deactivate()
        self.properties_out[:].pressure.\
            fix(self.properties_in_vap[0].pressure.value)

        if solver is not None:
            if outlvl > 2:
                tee = True
            else:
                tee = False

        solver_output = solver.solve(self, tee=tee)

        if solver_output.solver.termination_condition == \
                TerminationCondition.optimal:
            _log.info('{} Mass balance solve successful.'
                      .format(self.name))
        else:
            _log.info('{} Mass balance solve failed.'
                      .format(self.name))

        # Activate energy balance
        self.enthalpy_mixing_equations.activate()
        self.properties_out[:].temperature.unfix()

        solver_output = solver.solve(self, tee=tee)

        if solver_output.solver.termination_condition == \
                TerminationCondition.optimal:
            _log.info('{} Mass/Energy balance solve successful.'
                      .format(self.name))
        else:
            _log.info('{} Mass/Energy balance solve failed.'
                      .format(self.name))

        # Activate pressure balance
        self.pressure_drop_equation.activate()
        self.properties_out[:].pressure.unfix()

        if solver_output.solver.termination_condition == \
                TerminationCondition.optimal:
            _log.info('{} Tray initialisation complete.'
                      .format(self.name))
        else:
            _log.info('{} Mass/Energy/Pressure balance solve failed.'
                      .format(self.name))
