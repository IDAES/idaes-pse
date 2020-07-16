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
Methods for setting up FcTP as the state variables in a generic property
package
"""
from pyomo.environ import Constraint, Expression, NonNegativeReals, Var

from idaes.core import (MaterialFlowBasis,
                        MaterialBalanceType,
                        EnergyBalanceType)

from idaes.generic_models.properties.core.state_definitions.FTPx import (
    state_initialization)


def set_metadata(b):
    # The default metadata should be correct in this case, so no need to update
    pass


def define_state(b):
    # FcTP formulation always requires a flash, so set flag to True
    # TODO: should have some checking to make sure developers implement this properly
    b.always_flash = True

    # Get bounds if provided
    try:
        f_bounds = b.params.config.state_bounds["flow_mol_comp"]
    except (KeyError, TypeError):
        f_bounds = (None, None)

    try:
        t_bounds = b.params.config.state_bounds["temperature"]
    except (KeyError, TypeError):
        t_bounds = (None, None)

    try:
        p_bounds = b.params.config.state_bounds["pressure"]
    except (KeyError, TypeError):
        p_bounds = (None, None)

    # Set an initial value for each state var
    if f_bounds == (None, None):
        # No bounds, default to 1
        f_init = 1
    elif f_bounds[1] is None and f_bounds[0] is not None:
        # Only lower bound, use lower bound + 10
        f_init = f_bounds[0] + 10
    elif f_bounds[1] is not None and f_bounds[0] is None:
        # Only upper bound, use half upper bound
        f_init = f_bounds[1]/2
    else:
        # Both bounds, use mid point
        f_init = (f_bounds[0] + f_bounds[1])/2

    if t_bounds == (None, None):
        # No bounds, default to 298.15
        t_init = 298.15
    elif t_bounds[1] is None and t_bounds[0] is not None:
        # Only lower bound, use lower bound + 10
        t_init = t_bounds[0] + 10
    elif t_bounds[1] is not None and t_bounds[0] is None:
        # Only upper bound, use half upper bound
        t_init = t_bounds[1]/2
    else:
        # Both bounds, use mid point
        t_init = (t_bounds[0] + t_bounds[1])/2

    if p_bounds == (None, None):
        # No bounds, default to 101325
        p_init = 101325
    elif p_bounds[1] is None and p_bounds[0] is not None:
        # Only lower bound, use lower bound + 10
        p_init = p_bounds[0] + 10
    elif p_bounds[1] is not None and p_bounds[0] is None:
        # Only upper bound, use half upper bound
        p_init = p_bounds[1]/2
    else:
        # Both bounds, use mid point
        p_init = (p_bounds[0] + p_bounds[1])/2

    # Add state variables
    b.flow_mol_comp = Var(b.params.component_list,
                          initialize=f_init,
                          domain=NonNegativeReals,
                          bounds=f_bounds,
                          doc=' Component molar flowrate')
    b.pressure = Var(initialize=p_init,
                     domain=NonNegativeReals,
                     bounds=p_bounds,
                     doc='State pressure')
    b.temperature = Var(initialize=t_init,
                        domain=NonNegativeReals,
                        bounds=t_bounds,
                        doc='State temperature')

    # Add supporting variables
    b.flow_mol = Expression(
        expr=sum(b.flow_mol_comp[j] for j in b.params.component_list),
        doc="Total molar flowrate")

    b.flow_mol_phase = Var(b.params.phase_list,
                           initialize=f_init / len(b.params.phase_list),
                           domain=NonNegativeReals,
                           bounds=f_bounds,
                           doc='Phase molar flow rates')

    b.mole_frac_comp = Var(b.params.component_list,
                           bounds=(0, None),
                           initialize=1 / len(b.params.component_list),
                           doc='Mixture mole fractions')

    b.mole_frac_phase_comp = Var(
        b.params._phase_component_set,
        initialize=1/len(b.params.component_list),
        bounds=(0, None),
        doc='Phase mole fractions')

    b.phase_frac = Var(
        b.params.phase_list,
        initialize=1/len(b.params.phase_list),
        bounds=(0, None),
        doc='Phase fractions')

    # Add supporting constraints
    def rule_mole_frac_comp(b, j):
        if len(b.params.component_list) > 1:
            return b.flow_mol_comp[j] == b.mole_frac_comp[j]*sum(
                b.flow_mol_comp[k] for k in b.params.component_list)
        else:
            return b.mole_frac_comp[j] == 1
    b.mole_frac_comp_eq = Constraint(b.params.component_list,
                                     rule=rule_mole_frac_comp)

    if len(b.params.phase_list) == 1:
        def rule_total_mass_balance(b):
            return b.flow_mol_phase[b.params.phase_list[1]] == b.flow_mol
        b.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.mole_frac_comp[i]*1e3 == \
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
            return b.flow_mol_comp[i] == sum(
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
            return b.flow_mol_comp[i] == sum(
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

    def define_state_vars_FTPx():
        """Define state vars."""
        return {"flow_mol_comp": b.flow_mol_comp,
                "temperature": b.temperature,
                "pressure": b.pressure}
    b.define_state_vars = define_state_vars_FTPx

    def define_display_vars_FTPx():
        """Define display vars."""
        return {"Molar Flowrate": b.flow_mol_comp,
                "Temperature": b.temperature,
                "Pressure": b.pressure}
    b.define_display_vars = define_display_vars_FTPx


# Inherit state_initialization from FTPX form, as the process is the same


do_not_initialize = []
