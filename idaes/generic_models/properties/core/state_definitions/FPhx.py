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
Methods for setting up FPhx as the state variables in a generic property
package
"""
from pyomo.environ import Constraint, NonNegativeReals, Var

from idaes.core import (MaterialFlowBasis,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.generic_models.properties.core.generic.utility import \
    get_bounds_from_config
from idaes.generic_models.properties.core.state_definitions.FTPx import (
    state_initialization)
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)


# TODO : Need a way to get a guess for T during initialization
def set_metadata(b):
    # Need to update metadata so that enth_mol is recorded as being part of the
    # state variables, and to ensure that getattr does not try to build it
    # using the default method.
    b.get_metadata().properties['enth_mol'] = {
        'method': None}


def define_state(b):
    # FTPx formulation always requires a flash, so set flag to True
    # TODO: should have some checking to make sure developers implement this properly
    b.always_flash = True

    # Check that only necessary state_bounds are defined
    expected_keys = ["flow_mol", "enth_mol", "temperature", "pressure"]
    if (b.params.config.state_bounds is not None and
            any(b.params.config.state_bounds.keys()) not in expected_keys):
        for k in b.params.config.state_bounds.keys():
            if "mole_frac" in k:
                _log.warning("{} - found state_bounds argument for {}."
                             " Mole fraction bounds are set automatically and "
                             "this argument will be ignored."
                             .format(b.name, k))
            elif k not in expected_keys:
                raise ConfigurationError(
                    "{} - found unexpected state_bounds key {}. Please ensure "
                    "bounds are provided only for expected state variables "
                    "and that you have typed the variable names correctly."
                    .format(b.name, k))

    units = b.params.get_metadata().derived_units
    # Get bounds and initial values from config args
    f_bounds, f_init = get_bounds_from_config(
        b, "flow_mol", units["flow_mole"])
    h_bounds, h_init = get_bounds_from_config(
        b, "enth_mol", units["energy_mole"])
    p_bounds, p_init = get_bounds_from_config(
        b, "pressure", units["pressure"])
    t_bounds, t_init = get_bounds_from_config(
        b, "temperature", units["temperature"])

    # Add state variables
    b.flow_mol = Var(initialize=f_init,
                     domain=NonNegativeReals,
                     bounds=f_bounds,
                     doc=' Total molar flowrate',
                     units=units["flow_mole"])
    b.mole_frac_comp = Var(b.params.component_list,
                           bounds=(0, None),
                           initialize=1 / len(b.params.component_list),
                           doc='Mixture mole fractions',
                           units=None)
    b.pressure = Var(initialize=p_init,
                     domain=NonNegativeReals,
                     bounds=p_bounds,
                     doc='State pressure',
                     units=units["pressure"])

    b.enth_mol = Var(initialize=h_init,
                     bounds=h_bounds,
                     doc='State molar enthalpy',
                     units=units["energy_mole"])

    # Add supporting variables
    if f_init is None:
        fp_init = None
    else:
        fp_init = f_init / len(b.params.phase_list)

    b.flow_mol_phase = Var(b.params.phase_list,
                           initialize=fp_init,
                           domain=NonNegativeReals,
                           bounds=f_bounds,
                           doc='Phase molar flow rates',
                           units=units["flow_mole"])

    b.mole_frac_phase_comp = Var(
        b.params._phase_component_set,
        initialize=1/len(b.params.component_list),
        bounds=(0, None),
        doc='Phase mole fractions',
        units=None)

    b.temperature = Var(initialize=t_init,
                        domain=NonNegativeReals,
                        bounds=t_bounds,
                        doc='Temperature',
                        units=units["temperature"])

    b.phase_frac = Var(
        b.params.phase_list,
        initialize=1/len(b.params.phase_list),
        bounds=(0, None),
        doc='Phase fractions',
        units=None)

    # Add supporting constraints
    if b.config.defined_state is False:
        # applied at outlet only
        b.sum_mole_frac_out = Constraint(
            expr=1e3 == 1e3*sum(b.mole_frac_comp[i]
                                for i in b.params.component_list))

    def rule_enth_mol(b):
        return b.enth_mol == sum(b.enth_mol_phase[p]*b.phase_frac[p]
                                 for p in b.params.phase_list)
    b.enth_mol_eq = Constraint(rule=rule_enth_mol,
                               doc="Total molar enthalpy mixing rule")

    if len(b.params.phase_list) == 1:
        def rule_total_mass_balance(b):
            return b.flow_mol_phase[b.params.phase_list[1]] == b.flow_mol
        b.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return 1e3*b.mole_frac_comp[i] == \
                1e3*b.mole_frac_phase_comp[b.params.phase_list[1], i]
        b.component_flow_balances = Constraint(b.params.component_list,
                                               rule=rule_comp_mass_balance)

        def rule_phase_frac(b, p):
            return b.phase_frac[p] == 1
        b.phase_fraction_constraint = Constraint(b.params.phase_list,
                                                 rule=rule_phase_frac)

    elif len(b.params.phase_list) == 2:
        # For two phase, use Rachford-Rice formulation
        def rule_total_mass_balance(b):
            return sum(b.flow_mol_phase[p] for p in b.params.phase_list) == \
                b.flow_mol
        b.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.flow_mol*b.mole_frac_comp[i] == sum(
                b.flow_mol_phase[p]*b.mole_frac_phase_comp[p, i]
                for p in b.params.phase_list
                if (p, i) in b.params._phase_component_set)
        b.component_flow_balances = Constraint(b.params.component_list,
                                               rule=rule_comp_mass_balance)

        def rule_mole_frac(b):
            return 1e3*sum(b.mole_frac_phase_comp[b.params.phase_list[1], i]
                           for i in b.params.component_list
                           if (b.params.phase_list[1], i)
                           in b.params._phase_component_set) -\
                1e3*sum(b.mole_frac_phase_comp[b.params.phase_list[2], i]
                        for i in b.params.component_list
                        if (b.params.phase_list[2], i)
                        in b.params._phase_component_set) == 0
        b.sum_mole_frac = Constraint(rule=rule_mole_frac)

        def rule_phase_frac(b, p):
            return b.phase_frac[p]*b.flow_mol == b.flow_mol_phase[p]
        b.phase_fraction_constraint = Constraint(b.params.phase_list,
                                                 rule=rule_phase_frac)

    else:
        # Otherwise use a general formulation
        def rule_comp_mass_balance(b, i):
            return b.flow_mol*b.mole_frac_comp[i] == sum(
                b.flow_mol_phase[p]*b.mole_frac_phase_comp[p, i]
                for p in b.params.phase_list
                if (p, i) in b.params._phase_component_set)
        b.component_flow_balances = Constraint(b.params.component_list,
                                               rule=rule_comp_mass_balance)

        def rule_mole_frac(b, p):
            return 1e3*sum(b.mole_frac_phase_comp[p, i]
                           for i in b.params.component_list
                           if (p, i) in b.params._phase_component_set) == 1e3
        b.sum_mole_frac = Constraint(b.params.phase_list,
                                     rule=rule_mole_frac)

        def rule_phase_frac(b, p):
            return b.phase_frac[p]*b.flow_mol == b.flow_mol_phase[p]
        b.phase_fraction_constraint = Constraint(b.params.phase_list,
                                                 rule=rule_phase_frac)

    # -------------------------------------------------------------------------
    # General Methods
    def get_material_flow_terms_FTPx(p, j):
        """Create material flow terms for control volume."""
        if j in b.params.component_list:
            return b.flow_mol_phase[p] * b.mole_frac_phase_comp[p, j]
        else:
            return 0
    b.get_material_flow_terms = get_material_flow_terms_FTPx

    def get_enthalpy_flow_terms_FTPx(p):
        """Create enthalpy flow terms."""
        return b.flow_mol_phase[p] * b.enth_mol_phase[p]
    b.get_enthalpy_flow_terms = get_enthalpy_flow_terms_FTPx

    def get_material_density_terms_FTPx(p, j):
        """Create material density terms."""
        if j in b.params.component_list:
            return b.dens_mol_phase[p] * b.mole_frac_phase_comp[p, j]
        else:
            return 0
    b.get_material_density_terms = get_material_density_terms_FTPx

    def get_energy_density_terms_FTPx(p):
        """Create energy density terms."""
        return b.dens_mol_phase[p] * b.enth_mol_phase[p]
    b.get_energy_density_terms = get_energy_density_terms_FTPx

    def default_material_balance_type_FTPx():
        return MaterialBalanceType.componentTotal
    b.default_material_balance_type = default_material_balance_type_FTPx

    def default_energy_balance_type_FTPx():
        return EnergyBalanceType.enthalpyTotal
    b.default_energy_balance_type = default_energy_balance_type_FTPx

    def get_material_flow_basis_FTPx():
        return MaterialFlowBasis.molar
    b.get_material_flow_basis = get_material_flow_basis_FTPx

    def define_state_vars_FPhx():
        """Define state vars."""
        return {"flow_mol": b.flow_mol,
                "mole_frac_comp": b.mole_frac_comp,
                "enth_mol": b.enth_mol,
                "pressure": b.pressure}
    b.define_state_vars = define_state_vars_FPhx

    def define_display_vars_FPhx():
        """Define display vars."""
        return {"Total Molar Flowrate": b.flow_mol,
                "Total Mole Fraction": b.mole_frac_comp,
                "Molar Enthalpy": b.enth_mol,
                "Pressure": b.pressure}
    b.define_display_vars = define_display_vars_FPhx


# Inherit state_initialization from FTPX form, as the process is the same


do_not_initialize = ["sum_mole_frac_out"]
