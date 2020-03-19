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
IDAES Component objects

@author: alee
"""
from pyomo.environ import Set, Param, Var
from pyomo.common.config import ConfigBlock, ConfigValue

from .process_base import (declare_process_block_class,
                           ProcessBlockData)


@declare_process_block_class("Component")
class ComponentData(ProcessBlockData):
    CONFIG = ConfigBlock()

    CONFIG.declare("mw", ConfigValue(
        domain=float,
        description="Molecular weight of component"))
    CONFIG.declare("pressure_crit", ConfigValue(
        domain=float,
        description="Critical pressure of component"))
    CONFIG.declare("temperature_crit", ConfigValue(
        domain=float,
        description="Critical temperature of component"))

    CONFIG.declare("dens_mol_liq_comp", ConfigValue(
        description="Method to use to calculate liquid phase molar density"))
    CONFIG.declare("enth_mol_liq_comp", ConfigValue(
        description="Method to calculate liquid component molar enthalpies"))
    CONFIG.declare("enth_mol_ig_comp", ConfigValue(
        description="Method to calculate ideal gas component molar enthalpies"
        ))
    CONFIG.declare("entr_mol_liq_comp", ConfigValue(
        description="Method to calculate liquid component molar entropies"))
    CONFIG.declare("entr_mol_ig_comp", ConfigValue(
        description="Method to calculate ideal gas component molar entropies"))
    CONFIG.declare("pressure_sat_comp", ConfigValue(
        description="Method to use to calculate saturation pressure"))

    CONFIG.declare("phase_equilibrium_form", ConfigValue(
        domain=dict,
        description="Form of phase equilibrium constraints for component"))

    CONFIG.declare("_component_list_exists", ConfigValue(
            default=False,
            doc="Internal config argument indicating whether component_list "
            "needs to be populated."))

    def build(self):
        super(ComponentData, self).build()

        # If the component_list does not exist, add reference to new Component
        # The IF is mostly for backwards compatability, to allow for old-style
        # property packages where the component_list already exists but we
        # need to add new Component objects
        if not self.config._component_list_exists:
            self.__add_to_component_list()

        # Create Param for molecular weight if provided
        if self.config.mw is not None:
            self.mw = Param(initialize=self.config.mw)

        # Create Vars for common parameters
        for p in ["pressure_crit", "temperature_crit"]:
            if self.config[p] is not None:
                self.add_component(p, Var(
                    initialize=self.config[p],
                    doc=self.config.get(p)._description))

    def __add_to_component_list(self):
        """
        Method to add reference to new Component in component_list
        """
        parent = self.parent_block()
        try:
            comp_list = getattr(parent, "component_list")
            comp_list.add(self.local_name)
        except AttributeError:
            # Parent does not have a component_list yet, so create one
            parent.component_list = Set(initialize=[self.local_name],
                                        ordered=True)
