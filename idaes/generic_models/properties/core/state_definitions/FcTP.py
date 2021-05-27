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
from pyomo.environ import \
    Constraint, Expression, NonNegativeReals, Var, units as pyunits

from idaes.core import (MaterialFlowBasis,
                        MaterialBalanceType,
                        EnergyBalanceType)

from idaes.generic_models.properties.core.state_definitions.FTPx import (
    state_initialization)
from idaes.generic_models.properties.core.generic.utility import \
    get_bounds_from_config
from .electrolyte_states import \
    define_electrolyte_state, calculate_electrolyte_scaling
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

# Set up logger
_log = idaeslog.getLogger(__name__)


def set_metadata(b):
    # The default metadata should be correct in this case, so no need to update
    pass


def define_state(b):
    # FcTP formulation always requires a flash, so set flag to True
    # TODO: should have some checking to make sure developers implement this properly
    b.always_flash = True

    # Check that only necessary state_bounds are defined
    expected_keys = ["flow_mol_comp", "temperature", "pressure"]
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
        b, "flow_mol_comp", units["flow_mole"])
    t_bounds, t_init = get_bounds_from_config(
        b, "temperature", units["temperature"])
    p_bounds, p_init = get_bounds_from_config(
        b, "pressure", units["pressure"])

    # Add state variables
    b.flow_mol_comp = Var(b.component_list,
                          initialize=f_init,
                          domain=NonNegativeReals,
                          bounds=f_bounds,
                          doc=' Component molar flowrate',
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
        expr=sum(b.flow_mol_comp[j] for j in b.component_list),
        doc="Total molar flowrate")

    if f_init is None:
        fp_init = None
    else:
        fp_init = f_init / len(b.phase_list)

    b.flow_mol_phase = Var(b.phase_list,
                           initialize=fp_init,
                           domain=NonNegativeReals,
                           bounds=f_bounds,
                           doc='Phase molar flow rates',
                           units=units["flow_mole"])

    b.mole_frac_comp = Var(b.component_list,
                           bounds=(0, None),
                           initialize=1 / len(b.component_list),
                           doc='Mixture mole fractions',
                           units=None)

    b.mole_frac_phase_comp = Var(
        b.phase_component_set,
        initialize=1/len(b.component_list),
        bounds=(0, None),
        doc='Phase mole fractions',
        units=None)

    def flow_mol_phase_comp_rule(b, p, j):
        return b.flow_mol_phase[p] * b.mole_frac_phase_comp[p, j]
    b.flow_mol_phase_comp = Expression(b.phase_component_set,
                                       rule=flow_mol_phase_comp_rule)

    b.phase_frac = Var(
        b.phase_list,
        initialize=1/len(b.phase_list),
        bounds=(0, None),
        doc='Phase fractions',
        units=None)

    # Add electrolye state vars if required
    if b.params._electrolyte:
        define_electrolyte_state(b)

    # Add supporting constraints
    def rule_mole_frac_comp(b, j):
        if len(b.component_list) > 1:
            return b.flow_mol_comp[j] == b.mole_frac_comp[j]*sum(
                b.flow_mol_comp[k] for k in b.component_list)
        else:
            return b.mole_frac_comp[j] == 1
    b.mole_frac_comp_eq = Constraint(b.component_list,
                                     rule=rule_mole_frac_comp)

    if len(b.phase_list) == 1:
        def rule_total_mass_balance(b):
            return b.flow_mol_phase[b.phase_list[1]] == b.flow_mol
        b.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.mole_frac_comp[i] == \
                b.mole_frac_phase_comp[b.phase_list[1], i]
        b.component_flow_balances = Constraint(b.component_list,
                                               rule=rule_comp_mass_balance)

        def rule_phase_frac(b, p):
            return b.phase_frac[p] == 1
        b.phase_fraction_constraint = Constraint(b.phase_list,
                                                 rule=rule_phase_frac)

    elif len(b.phase_list) == 2:
        # For two phase, use Rachford-Rice formulation
        def rule_total_mass_balance(b):
            return sum(b.flow_mol_phase[p] for p in b.phase_list) == \
                b.flow_mol
        b.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.flow_mol_comp[i] == sum(
                b.flow_mol_phase[p]*b.mole_frac_phase_comp[p, i]
                for p in b.phase_list
                if (p, i) in b.phase_component_set)
        b.component_flow_balances = Constraint(b.component_list,
                                               rule=rule_comp_mass_balance)

        def rule_mole_frac(b):
            return sum(b.mole_frac_phase_comp[b.phase_list[1], i]
                       for i in b.component_list
                       if (b.phase_list[1], i) in b.phase_component_set) -\
                sum(b.mole_frac_phase_comp[b.phase_list[2], i]
                    for i in b.component_list
                    if (b.phase_list[2], i) in b.phase_component_set) == 0
        b.sum_mole_frac = Constraint(rule=rule_mole_frac)

        def rule_phase_frac(b, p):
            return b.phase_frac[p]*b.flow_mol == b.flow_mol_phase[p]
        b.phase_fraction_constraint = Constraint(b.phase_list,
                                                 rule=rule_phase_frac)

    else:
        # Otherwise use a general formulation
        def rule_comp_mass_balance(b, i):
            return b.flow_mol_comp[i] == sum(
                b.flow_mol_phase[p]*b.mole_frac_phase_comp[p, i]
                for p in b.phase_list
                if (p, i) in b.phase_component_set)
        b.component_flow_balances = Constraint(b.component_list,
                                               rule=rule_comp_mass_balance)

        def rule_mole_frac(b, p):
            return sum(b.mole_frac_phase_comp[p, i]
                       for i in b.component_list
                       if (p, i) in b.phase_component_set) == 1
        b.sum_mole_frac = Constraint(b.phase_list,
                                     rule=rule_mole_frac)

        def rule_phase_frac(b, p):
            return b.phase_frac[p]*b.flow_mol == b.flow_mol_phase[p]
        b.phase_fraction_constraint = Constraint(b.phase_list,
                                                 rule=rule_phase_frac)

    # -------------------------------------------------------------------------
    # General Methods
    def get_material_flow_terms_FcTP(p, j):
        """Create material flow terms for control volume."""
        return b.flow_mol_phase_comp[p, j]
    b.get_material_flow_terms = get_material_flow_terms_FcTP

    def get_enthalpy_flow_terms_FcTP(p):
        """Create enthalpy flow terms."""
        # enth_mol_phase probably does not exist when this is created
        # Use try/except to build flow term if not present
        try:
            eflow = b._enthalpy_flow_term
        except AttributeError:
            def rule_eflow(b, p):
                return b.flow_mol_phase[p] * b.enth_mol_phase[p]
            eflow = b._enthalpy_flow_term = Expression(
                b.phase_list, rule=rule_eflow)
        return eflow[p]
    b.get_enthalpy_flow_terms = get_enthalpy_flow_terms_FcTP

    def get_material_density_terms_FcTP(p, j):
        """Create material density terms."""
        # dens_mol_phase probably does not exist when this is created
        # Use try/except to build term if not present
        try:
            mdens = b._material_density_term
        except AttributeError:
            def rule_mdens(b, p, j):
                return b.dens_mol_phase[p] * b.mole_frac_phase_comp[p, j]
            mdens = b._material_density_term = Expression(
                b.phase_component_set, rule=rule_mdens)
        return mdens[p, j]
    b.get_material_density_terms = get_material_density_terms_FcTP

    def get_energy_density_terms_FcTP(p):
        """Create energy density terms."""
        # Density and energy terms probably do not exist when this is created
        # Use try/except to build term if not present
        try:
            edens = b._energy_density_term
        except AttributeError:
            def rule_edens(b, p):
                return b.dens_mol_phase[p] * b.energy_internal_mol_phase[p]
            edens = b._energy_density_term = Expression(
                b.phase_list, rule=rule_edens)
        return edens[p]
    b.get_energy_density_terms = get_energy_density_terms_FcTP

    def default_material_balance_type_FcTP():
        return MaterialBalanceType.componentTotal
    b.default_material_balance_type = default_material_balance_type_FcTP

    def default_energy_balance_type_FcTP():
        return EnergyBalanceType.enthalpyTotal
    b.default_energy_balance_type = default_energy_balance_type_FcTP

    def get_material_flow_basis_FcTP():
        return MaterialFlowBasis.molar
    b.get_material_flow_basis = get_material_flow_basis_FcTP

    def define_state_vars_FcTP():
        """Define state vars."""
        return {"flow_mol_comp": b.flow_mol_comp,
                "temperature": b.temperature,
                "pressure": b.pressure}
    b.define_state_vars = define_state_vars_FcTP

    def define_display_vars_FcTP():
        """Define display vars."""
        return {"Molar Flowrate": b.flow_mol_comp,
                "Temperature": b.temperature,
                "Pressure": b.pressure}
    b.define_display_vars = define_display_vars_FcTP


def define_default_scaling_factors(b):
    """
    Method to set default scaling factors for the property package. Scaling
    factors are based on the default initial value for each variable provided
    in the state_bounds config argument.
    """
    # Get bounds and initial values from config args
    units = b.get_metadata().derived_units
    state_bounds = b.config.state_bounds

    if state_bounds is None:
        return

    try:
        f_bounds = state_bounds["flow_mol_comp"]
        if len(f_bounds) == 4:
            f_init = pyunits.convert_value(f_bounds[1],
                                           from_units=f_bounds[3],
                                           to_units=units["flow_mole"])
        else:
            f_init = f_bounds[1]
    except KeyError:
        f_init = 1

    try:
        p_bounds = state_bounds["pressure"]
        if len(p_bounds) == 4:
            p_init = pyunits.convert_value(p_bounds[1],
                                           from_units=p_bounds[3],
                                           to_units=units["pressure"])
        else:
            p_init = p_bounds[1]
    except KeyError:
        p_init = 1

    try:
        t_bounds = state_bounds["temperature"]
        if len(t_bounds) == 4:
            t_init = pyunits.convert_value(t_bounds[1],
                                           from_units=t_bounds[3],
                                           to_units=units["temperature"])
        else:
            t_init = t_bounds[1]
    except KeyError:
        t_init = 1

    # Set default scaling factors
    b.set_default_scaling("flow_mol", 1/f_init)
    b.set_default_scaling("flow_mol_phase", 1/f_init)
    b.set_default_scaling("flow_mol_comp", 1/f_init)
    b.set_default_scaling("flow_mol_phase_comp", 1/f_init)
    b.set_default_scaling("pressure", 1/p_init)
    b.set_default_scaling("temperature", 1/t_init)


def calculate_scaling_factors(b):
    sf_flow = iscale.get_scaling_factor(b.flow_mol, default=1, warning=True)
    sf_mf = {}
    for i, v in b.mole_frac_phase_comp.items():
        sf_mf[i] = iscale.get_scaling_factor(v, default=1e3, warning=True)

    for j in b.component_list:
        sf_j = iscale.get_scaling_factor(
            b.mole_frac_comp[j], default=1e3, warning=True)
        iscale.constraint_scaling_transform(
            b.mole_frac_comp_eq[j], sf_j, overwrite=False)

    if len(b.phase_list) == 1:
        iscale.constraint_scaling_transform(
            b.total_flow_balance, sf_flow, overwrite=False)

        for j in b.component_list:
            sf_j = iscale.get_scaling_factor(
                b.mole_frac_comp[j], default=1e3, warning=True)
            iscale.constraint_scaling_transform(
                b.component_flow_balances[j], sf_j, overwrite=False)

        # b.phase_fraction_constraint is well scaled

    elif len(b.phase_list) == 2:
        iscale.constraint_scaling_transform(
            b.total_flow_balance, sf_flow, overwrite=False)

        for j in b.component_list:
            sf_fc = iscale.get_scaling_factor(
                b.flow_mol_comp[j], default=1, warning=True)
            iscale.constraint_scaling_transform(
                b.component_flow_balances[j], sf_fc, overwrite=False)

        iscale.constraint_scaling_transform(
            b.sum_mole_frac, min(sf_mf.values()), overwrite=False)

        for p in b.phase_list:
            sf_p = iscale.get_scaling_factor(
                b.flow_mol_phase[p], default=1, warning=True)
            iscale.constraint_scaling_transform(
                b.phase_fraction_constraint[p], sf_p, overwrite=False)

    else:
        for j in b.component_list:
            sf_fc = iscale.get_scaling_factor(
                b.flow_mol_comp[j], default=1, warning=True)
            iscale.constraint_scaling_transform(
                b.component_flow_balances[j], sf_fc, overwrite=False)

        for p in b.phase_list:
            sf_p = iscale.get_scaling_factor(
                b.flow_mol_phase[p], default=1, warning=True)
            iscale.constraint_scaling_transform(
                b.sum_mole_frac[p], min(sf_mf[p, :].values()), overwrite=False)
            iscale.constraint_scaling_transform(
                b.phase_fraction_constraint[p], sf_p, overwrite=False)

    if b.params._electrolyte:
        calculate_electrolyte_scaling(b)

# Inherit state_initialization from FTPx form, as the process is the same


do_not_initialize = []


class FcTP(object):
    set_metadata = set_metadata
    define_state = define_state
    state_initialization = state_initialization
    do_not_initialize = do_not_initialize
    define_default_scaling_factors = define_default_scaling_factors
    calculate_scaling_factors = calculate_scaling_factors
