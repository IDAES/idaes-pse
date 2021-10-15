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
Methods for setting up FTPx as the state variables in a generic property
package
"""

from pyomo.environ import \
    Constraint, Expression, NonNegativeReals, Var, value, units as pyunits

from idaes.core import (MaterialFlowBasis,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.generic_models.properties.core.generic.utility import \
    get_bounds_from_config, get_method, GenericPropertyPackageError
from idaes.generic_models.properties.core.phase_equil.bubble_dew import \
    _valid_VL_component_list
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog
from .electrolyte_states import \
    define_electrolyte_state, calculate_electrolyte_scaling
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import BurntToast, UserModelError


# Set up logger
_log = idaeslog.getLogger(__name__)


def set_metadata(b):
    # This is the default assumption for state vars, so we don't need to change
    # the metadata dict
    pass


def define_state(b):
    # FTPx formulation always requires a flash, so set flag to True
    # TODO: should have some checking to make sure developers implement this properly
    b.always_flash = True

    units = b.params.get_metadata().derived_units

    # Check that only necessary state_bounds are defined
    expected_keys = ["flow_mol", "temperature", "pressure"]
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

    # Get bounds and initial values from config args
    f_bounds, f_init = get_bounds_from_config(
        b, "flow_mol", units["flow_mole"])
    t_bounds, t_init = get_bounds_from_config(
        b, "temperature", units["temperature"])
    p_bounds, p_init = get_bounds_from_config(
        b, "pressure", units["pressure"])

    # Add state variables
    b.flow_mol = Var(initialize=f_init,
                     domain=NonNegativeReals,
                     bounds=f_bounds,
                     doc=' Total molar flowrate',
                     units=units["flow_mole"])
    b.mole_frac_comp = Var(b.component_list,
                           bounds=(0, None),
                           initialize=1 / len(b.component_list),
                           doc='Mixture mole fractions',
                           units=None)
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

    b.mole_frac_phase_comp = Var(
        b.phase_component_set,
        initialize=1/len(b.component_list),
        bounds=(0, None),
        doc='Phase mole fractions',
        units=None)

    def Fpc_expr(b, p, j):
        return b.flow_mol_phase[p] * b.mole_frac_phase_comp[p, j]
    b.flow_mol_phase_comp = Expression(
        b.phase_component_set,
        rule=Fpc_expr,
        doc='Phase-component molar flowrates')

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
    if b.config.defined_state is False:
        # applied at outlet only
        b.sum_mole_frac_out = Constraint(
            expr=1 == sum(b.mole_frac_comp[i] for i in b.component_list))

    if len(b.phase_list) == 1:
        def rule_total_mass_balance(b):
            return b.flow_mol_phase[b.phase_list.first()] == b.flow_mol
        b.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.mole_frac_comp[i] == \
                b.mole_frac_phase_comp[b.phase_list.first(), i]
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
            return b.flow_mol*b.mole_frac_comp[i] == sum(
                b.flow_mol_phase[p]*b.mole_frac_phase_comp[p, i]
                for p in b.phase_list
                if (p, i) in b.phase_component_set)
        b.component_flow_balances = Constraint(b.component_list,
                                               rule=rule_comp_mass_balance)

        def rule_mole_frac(b):
            return sum(b.mole_frac_phase_comp[b.phase_list.first(), i]
                       for i in b.component_list
                       if (b.phase_list.first(), i) in b.phase_component_set) -\
                sum(b.mole_frac_phase_comp[b.phase_list.last(), i]
                    for i in b.component_list
                    if (b.phase_list.last(), i) in b.phase_component_set) == 0
        b.sum_mole_frac = Constraint(rule=rule_mole_frac)

        def rule_phase_frac(b, p):
            return b.phase_frac[p]*b.flow_mol == b.flow_mol_phase[p]
        b.phase_fraction_constraint = Constraint(b.phase_list,
                                                 rule=rule_phase_frac)

    else:
        # Otherwise use a general formulation
        def rule_comp_mass_balance(b, i):
            return b.flow_mol*b.mole_frac_comp[i] == sum(
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
    def get_material_flow_terms_FTPx(p, j):
        """Create material flow terms for control volume."""
        return b.flow_mol_phase_comp[p, j]
    b.get_material_flow_terms = get_material_flow_terms_FTPx

    def get_enthalpy_flow_terms_FTPx(p):
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
    b.get_enthalpy_flow_terms = get_enthalpy_flow_terms_FTPx

    def get_material_density_terms_FTPx(p, j):
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
    b.get_material_density_terms = get_material_density_terms_FTPx

    def get_energy_density_terms_FTPx(p):
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
        return {"flow_mol": b.flow_mol,
                "mole_frac_comp": b.mole_frac_comp,
                "temperature": b.temperature,
                "pressure": b.pressure}
    b.define_state_vars = define_state_vars_FTPx

    def define_display_vars_FTPx():
        """Define display vars."""
        return {"Total Molar Flowrate": b.flow_mol,
                "Total Mole Fraction": b.mole_frac_comp,
                "Temperature": b.temperature,
                "Pressure": b.pressure}
    b.define_display_vars = define_display_vars_FTPx


def state_initialization(b):
    for j in b.component_list:
        if value(b.mole_frac_comp[j]) <= -1E-6:
            raise ValueError(f"Component {j} has a negative "
                             f"mole fraction in block {b.getname()}. "
                             "Check your initialization.")
    
    if sum([value(b.mole_frac_comp[j]) for j in b.component_list]) >= 1.01:
        raise ValueError(f"Mole fractions in block {b.getname()} sum to value "
                         "greater than 1.")

    vl_comps = []
    henry_comps = []
    _pe_pairs = []
    if hasattr(b.params,'_pe_pairs'):
        _pe_pairs = b.params._pe_pairs
    init_VLE = False
    for pp in _pe_pairs:
        # Look for a VLE pair with this phase - should only be 1
        if ((b.params.get_phase(pp[0]).is_liquid_phase() and
             b.params.get_phase(pp[1]).is_vapor_phase()) or
            (b.params.get_phase(pp[1]).is_liquid_phase() and
             b.params.get_phase(pp[0]).is_vapor_phase())):
            init_VLE = True
            # Get bubble and dew points
            tbub = None
            tdew = None
            if hasattr(b, "eq_temperature_bubble"):
                try:
                    tbub = b.temperature_bubble[pp].value
                except KeyError:
                    pass
            if hasattr(b, "eq_temperature_dew"):
                try:
                    tdew = b.temperature_dew[pp].value
                except KeyError:
                    pass
            if len(vl_comps) > 0:
                # More than one VLE. Just use the default initialization for
                # now
                init_VLE = False
            (l_phase,
             v_phase,
             vl_comps,
             henry_comps,
             l_only_comps,
             v_only_comps) = _valid_VL_component_list(b, pp)
            
            H = {j: value(get_method(b, "henry_component", j, l_phase)(
                        b, l_phase, j, b.temperature))
                        for j in henry_comps}
            try:
                psat = {j: value(get_method(b, "pressure_sat_comp", j)(
                        b, b.params.get_component(j),b.temperature)) 
                        for j in vl_comps}
                raoult_init = True
            except GenericPropertyPackageError:
                # No method for calculating Psat, use default values
                raoult_init = False
            if raoult_init:
                for j in vl_comps:
                    if psat[j] <= -1E-6:
                        raise UserModelError(f"Component {j} has a negative "
                                         f"saturation pressure in block "
                                         f"{b.getname()}. Check "
                                         f"your implementation and parameters.")
                for j in henry_comps:
                    if H[j] <= -1E-6:
                        raise ConfigurationError(f"Component {j} has a negative "
                                         f"Henry's Law constant in block "
                                         f"{b.getname()}. Check "
                                         f"your implementation and parameters.")

                # Calculate bubble and dew pressures for an ideal mixture
                if len(l_only_comps) == 0:
                    p_dew_ideal = 1/(sum([value(b.mole_frac_comp[j])/psat[j]
                                        for j in vl_comps]) + sum(
                                        [value(b.mole_frac_comp[j])/H[j]
                                        for j in henry_comps]))
                else:
                    p_dew_ideal = 0
                if len(v_only_comps) == 0:
                    p_bubble_ideal = (sum([value(b.mole_frac_comp[j])*psat[j]
                                        for j in vl_comps]) + sum(
                                        [value(b.mole_frac_comp[j])*H[j]
                                        for j in henry_comps]))
                else:
                    p_bubble_ideal = float('inf')
                    
                # Sanity check for the block's pressure. The method to 
                # calculate the ideal vapor fraction converges may not
                # converge to a reasonable value if the system pressure does
                # not lie between the ideal dew and bubble pressures. If 
                # this is the case, we give up and assume that the 
                # is predominantly one phase or the other
                if value(b.pressure) <= p_dew_ideal:
                    vapFrac = 0.95
                elif value(b.pressure) >= p_bubble_ideal:
                    vapFrac = 0.05
                else:
                    # Equation defining vapor fraction for an ideal mixture.
                    # This equation has two roots: a trivial one at vapFrac = 1
                    # and a nontrivial one for some 0<vapFrac<1. B is a
                    # convex function
                    def B(vapFrac):
                        return (
                            sum([value(b.mole_frac_comp[j]*psat[j]
                                  /(vapFrac*psat[j]+(1-vapFrac)*b.pressure))
                                for j in vl_comps])
                            + sum([value(b.mole_frac_comp[j]/vapFrac)
                                   for j in v_only_comps]) 
                            + sum([value(b.mole_frac_comp[j]*H[j]
                                  /(vapFrac*H[j]+(1-vapFrac)*b.pressure))
                                for j in henry_comps]) - 1
                            )
                    def dB_dvapFrac(vapFrac):
                        return (
                            sum([value(b.mole_frac_comp[j]*psat[j]
                                       *(b.pressure - psat[j])
                                  /(vapFrac*psat[j]+(1-vapFrac)*b.pressure)**2)
                                for j in vl_comps])
                            + sum([value(b.mole_frac_comp[j]*H[j]
                                       *(b.pressure - H[j])
                                  /(vapFrac*H[j]+(1-vapFrac)*b.pressure)**2)
                                for j in henry_comps])
                            - sum([value(b.mole_frac_comp[j]/vapFrac**2)
                                   for j in v_only_comps])
                            )
                    # Newton's method will always undershoot the root of a 
                    # convex equation if one starts from a value from which
                    # the function is positive. However, there is a singularity
                    # at vapFrac=0 if there are noncondensable components.
                    # Therefore, use a binary search to find a value where
                    # B(vapFrac) is positive
                    k = 0
                    vapFrac = 0.5
                    while B(vapFrac) < 0:
                        vapFrac /= 2
                        k += 1
                        if k > 40:
                            raise_runtime_error(b)
                    # Now use Newton's method to calculate vapFrac
                    k = 0
                    while abs(B(vapFrac)) > 1E-6:
                        vapFrac -= B(vapFrac)/dB_dvapFrac(vapFrac)
                        k += 1
                        if k > 40 or vapFrac<0 or vapFrac>1:
                            raise_runtime_error(b)

                
    for p in b.phase_list:
        # Start with phase mole fractions equal to total mole fractions
        for j in b.component_list:
            if (p, j) in b.phase_component_set:
                b.mole_frac_phase_comp[p, j].value = b.mole_frac_comp[j].value

        b.flow_mol_phase[p].value = value(
                        b.flow_mol / len(b.phase_list))

        # Try to refine guesses - Check phase type
        pobj = b.params.get_phase(p)

        if not hasattr(b.params, "_pe_pairs"):
            return #TODO: Are we sure we want to end initialization instead of
                    # breaking?
        
        if pobj.is_liquid_phase() and init_VLE:
            if tbub is None and tdew is None:
                # No VLE pair found, or no bubble and dew point
                # Do nothing
                pass
            elif tdew is not None and b.temperature.value > tdew:
                # Pure vapour
                b.flow_mol_phase[p].value = value(1e-5*b.flow_mol)
                b.phase_frac[p].value = 1e-5

                for j in b.component_list:
                    if (p, j) in b.phase_component_set:
                        b.mole_frac_phase_comp[p, j].value = \
                            b._mole_frac_tdew[pp, j].value
            elif tbub is not None and b.temperature.value < tbub:
                # Pure liquid
                b.flow_mol_phase[p].value = value(b.flow_mol)
                b.phase_frac[p].value = 1

                for j in b.component_list:
                    if (p, j) in b.phase_component_set:
                        b.mole_frac_phase_comp[p, j].value = \
                            b.mole_frac_comp[j].value
            elif raoult_init:
                if tbub is not None and tdew is not None:
                    # Two-phase with bounds two-phase region
                    # Thanks to Rahul Gandhi for the method
                    vapFrac = value((b.temperature-tbub) / (tdew-tbub))
                else:
                    # Two-phase, but with non-vaporizables and/or non-condensables
                    # Cannot take difference of tdew and tbub, so use vapor
                    # fraction calculated with Raoult's law.
                    pass
                b.flow_mol_phase[p].value = value((1-vapFrac)*b.flow_mol)
                for j in vl_comps:
                    kfact = value(psat[j] / b.pressure)
                    b.mole_frac_phase_comp[p, j].value = value(
                        b.mole_frac_comp[j]/(1+vapFrac*(kfact-1)))
                for j in henry_comps:
                    kfact = value(H[j] / b.pressure)
                    b.mole_frac_phase_comp[p, j].value = value(
                        b.mole_frac_comp[j]/(1+vapFrac*(kfact-1)))
                for j in l_only_comps:
                    b.mole_frac_phase_comp[p, j].value = value(
                        b.mole_frac_comp[j]/(1-vapFrac))
                else:
                    # Cannot do anything without a method to calculate Psat
                    pass

        elif pobj.is_vapor_phase() and init_VLE:            
            if tbub is None and tdew is None:
                # No VLE pair found, or no bubble and dew point
                # Do nothing
                pass
            elif tdew is not None and b.temperature.value > tdew:
                # Pure vapour
                b.flow_mol_phase[p].value = value(b.flow_mol)
                b.phase_frac[p].value = 1

                for j in b.component_list:
                    if (p, j) in b.phase_component_set:
                        b.mole_frac_phase_comp[p, j].value = \
                            b.mole_frac_comp[j].value
            elif tbub is not None and b.temperature.value < tbub:
                # Pure liquid
                b.flow_mol_phase[p].value = value(1e-5*b.flow_mol)
                b.phase_frac[p].value = 1e-5

                for j in b.component_list:
                    if (p, j) in b.phase_component_set:
                        b.mole_frac_phase_comp[p, j].value = \
                            b._mole_frac_tbub[pp, j].value
            elif raoult_init:
                if tbub is not None and tdew is not None:
                    # Two-phase with bounds two-phase region
                    # Thanks to Rahul Gandhi for the method
                    vapFrac = value((b.temperature-tbub) / (tdew-tbub))
                else:
                    # Two-phase, but with non-vaporizables and/or non-condensables
                    # Cannot take difference of tdew and tbub, so use vapor
                    # fraction calculated with Raoult's law.
                    pass
                b.flow_mol_phase[p].value = value(vapFrac*b.flow_mol)
                for j in vl_comps:
                    kfact = value(psat[j] / b.pressure)
                    b.mole_frac_phase_comp[p, j].value = value(
                        b.mole_frac_comp[j]*kfact/(kfact+1-vapFrac))
                for j in henry_comps:
                    kfact = value(H[j] / b.pressure)
                    b.mole_frac_phase_comp[p, j].value = value(
                        b.mole_frac_comp[j]*kfact/(kfact+1-vapFrac))
                for j in v_only_comps:
                    b.mole_frac_phase_comp[p, j].value = value(
                        b.mole_frac_comp[j]/vapFrac)
            else:
                # Cannot do anything without a method to calculate Psat
                pass


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
        f_bounds = state_bounds["flow_mol"]
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
    sf_flow = iscale.get_scaling_factor(
        b.flow_mol, default=1, warning=True)
    sf_mf = {}
    for i, v in b.mole_frac_phase_comp.items():
        sf_mf[i] = iscale.get_scaling_factor(v, default=1e3, warning=True)

    if b.config.defined_state is False:
        iscale.constraint_scaling_transform(
            b.sum_mole_frac_out, min(sf_mf.values()), overwrite=False)

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
            sf_j = iscale.get_scaling_factor(
                b.mole_frac_comp[j], default=1e3, warning=True)
            iscale.constraint_scaling_transform(
                b.component_flow_balances[j], sf_j*sf_flow, overwrite=False)

        iscale.constraint_scaling_transform(
            b.sum_mole_frac, min(sf_mf.values()), overwrite=False)

        for p in b.phase_list:
            iscale.constraint_scaling_transform(
                b.phase_fraction_constraint[p], sf_flow, overwrite=False)

    else:
        iscale.constraint_scaling_transform(
            b.total_flow_balance, sf_flow, overwrite=False)

        for j in b.component_list:
            sf_j = iscale.get_scaling_factor(
                b.mole_frac_comp[j], default=1e3, warning=True)
            iscale.constraint_scaling_transform(
                b.component_flow_balances[j], sf_j*sf_flow, overwrite=False)

        for p in b.phase_list:
            iscale.constraint_scaling_transform(
                b.sum_mole_frac[p], min(sf_mf[p, :].values()), overwrite=False)
            iscale.constraint_scaling_transform(
                b.phase_fraction_constraint[p], sf_flow, overwrite=False)

    if b.params._electrolyte:
        calculate_electrolyte_scaling(b)


do_not_initialize = ["sum_mole_frac_out"]


class FTPx(object):
    set_metadata = set_metadata
    define_state = define_state
    state_initialization = state_initialization
    do_not_initialize = do_not_initialize
    define_default_scaling_factors = define_default_scaling_factors
    calculate_scaling_factors = calculate_scaling_factors

# Mathematically, it should be impossible for the ideal vapor fraction
# calculations to not converge, given nonnegative saturation pressures, Henry's
# law constants, and mole fractions (that also add up to 1). However, the
# methods that check those properties relax the bounds slightly to avoid
# raising exceptions over roundoff errors, and, of course, floating point
# calculations do not perfectly capture mathematical reality. Therefore, 
# the error message is put in a separate function so we can test formatting.
def raise_runtime_error(b):
    raise RuntimeError("Could not calculate ideal vapor fraction "
                        f"in initialization of block {b.getname()}. "
                        "Check that mole fractions and saturation "
                        "pressures have reasonable values.")