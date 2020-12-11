##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Methods for setting up FpcTP as the state variables in a generic property
package
"""
from pyomo.environ import Constraint, Expression, NonNegativeReals, Var, value

from idaes.core import (MaterialFlowBasis,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.generic_models.properties.core.generic.utility import \
    get_bounds_from_config
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)


def set_metadata(b):
    # The default metadata should be correct in this case, so no need to update
    pass


def define_state(b):
    # FpcTP contains full information on the phase equilibrium, so flash
    # calculations re not always needed
    b.always_flash = False

    # Check that only necessary state_bounds are defined
    expected_keys = ["flow_mol_phase_comp", "enth_mol",
                     "temperature", "pressure"]
    if (b.params.config.state_bounds is not None and
            any(b.params.config.state_bounds.keys()) not in expected_keys):
        for k in b.params.config.state_bounds.keys():
            if k not in expected_keys:
                raise ConfigurationError(
                    "{} - found unexpected state_bounds key {}. Please ensure "
                    "bounds are provided only for expected state variables "
                    "and that you have typed the variable names correctly."
                    .format(b.name, k))

    units = b.params.get_metadata().derived_units
    # Get bounds and initial values from config args
    f_bounds, f_init = get_bounds_from_config(
        b, "flow_mol_phase_comp", units["flow_mole"])
    t_bounds, t_init = get_bounds_from_config(
        b, "temperature", units["temperature"])
    p_bounds, p_init = get_bounds_from_config(
        b, "pressure", units["pressure"])

    # Add state variables
    b.flow_mol_phase_comp = Var(b.params._phase_component_set,
                                initialize=f_init,
                                domain=NonNegativeReals,
                                bounds=f_bounds,
                                doc='Phase-component molar flowrate',
                                units=units["flow_mole"])
    b.pressure = Var(initialize=p_init,
                     domain=NonNegativeReals,
                     bounds=p_bounds,
                     doc='State pressure',
                     units=units["pressure"])
    b.temperature = Var(initialize=t_init,
                        domain=NonNegativeReals,
                        bounds=t_bounds,
                        doc='State temperature',
                        units=units["temperature"])

    # Add supporting variables
    b.flow_mol = Expression(
        expr=sum(b.flow_mol_phase_comp[i]
                 for i in b.params._phase_component_set),
        doc="Total molar flowrate")

    def flow_mol_phase(b, p):
        return sum(b.flow_mol_phase_comp[p, j]
                   for j in b.params.component_list
                   if (p, j) in b.params._phase_component_set)
    b.flow_mol_phase = Expression(b.params.phase_list,
                                  rule=flow_mol_phase,
                                  doc='Phase molar flow rates')

    def rule_flow_mol_comp(b, j):
        return sum(b.flow_mol_phase_comp[p, j]
                   for p in b.params.phase_list
                   if (p, j) in b.params._phase_component_set)
    b.flow_mol_comp = Expression(b.params.component_list,
                                 rule=rule_flow_mol_comp,
                                 doc='Component molar flow rates')

    def mole_frac_comp(b, j):
        return (sum(b.flow_mol_phase_comp[p, j]
                    for p in b.params.phase_list
                    if (p, j) in b.params._phase_component_set) / b.flow_mol)
    b.mole_frac_comp = Expression(b.params.component_list,
                                  rule=mole_frac_comp,
                                  doc='Mixture mole fractions')

    b.mole_frac_phase_comp = Var(
            b.params._phase_component_set,
            initialize=1/len(b.params.component_list),
            doc='Phase mole fractions',
            units=None)

    def rule_mole_frac_phase_comp(b, p, j):
        return b.mole_frac_phase_comp[p, j] * b.flow_mol_phase[p] == \
            b.flow_mol_phase_comp[p, j]
    b.mole_frac_phase_comp_eq = Constraint(
        b.params._phase_component_set, rule=rule_mole_frac_phase_comp)

    def rule_phase_frac(b, p):
        if len(b.params.phase_list) == 1:
            return 1
        else:
            return b.flow_mol_phase[p] / b.flow_mol
    b.phase_frac = Expression(
        b.params.phase_list,
        rule=rule_phase_frac,
        doc='Phase fractions')

    # -------------------------------------------------------------------------
    # General Methods
    def get_material_flow_terms_FpcTP(p, j):
        """Create material flow terms for control volume."""
        return b.flow_mol_phase_comp[p, j]
    b.get_material_flow_terms = get_material_flow_terms_FpcTP

    def get_enthalpy_flow_terms_FpcTP(p):
        """Create enthalpy flow terms."""
        return b.flow_mol_phase[p] * b.enth_mol_phase[p]
    b.get_enthalpy_flow_terms = get_enthalpy_flow_terms_FpcTP

    def get_material_density_terms_FpcTP(p, j):
        """Create material density terms."""
        if j in b.params.component_list:
            return b.dens_mol_phase[p] * b.mole_frac_phase_comp[p, j]
        else:
            return 0
    b.get_material_density_terms = get_material_density_terms_FpcTP

    def get_energy_density_terms_FpcTP(p):
        """Create energy density terms."""
        return b.dens_mol_phase[p] * b.enth_mol_phase[p]
    b.get_energy_density_terms = get_energy_density_terms_FpcTP

    def default_material_balance_type_FpcTP():
        return MaterialBalanceType.componentTotal
    b.default_material_balance_type = default_material_balance_type_FpcTP

    def default_energy_balance_type_FpcTP():
        return EnergyBalanceType.enthalpyTotal
    b.default_energy_balance_type = default_energy_balance_type_FpcTP

    def get_material_flow_basis_FpcTP():
        return MaterialFlowBasis.molar
    b.get_material_flow_basis = get_material_flow_basis_FpcTP

    def define_state_vars_FpcTP():
        """Define state vars."""
        return {"flow_mol_phase_comp": b.flow_mol_phase_comp,
                "temperature": b.temperature,
                "pressure": b.pressure}
    b.define_state_vars = define_state_vars_FpcTP

    def define_display_vars_FpcTP():
        """Define display vars."""
        return {"Molar Flowrate": b.flow_mol_phase_comp,
                "Temperature": b.temperature,
                "Pressure": b.pressure}
    b.define_display_vars = define_display_vars_FpcTP


def state_initialization(b):
    for i in b.params._phase_component_set:
        b.mole_frac_phase_comp[i].value = value(
            b.flow_mol_phase_comp[i] / b.flow_mol_phase[i[0]])


do_not_initialize = []
