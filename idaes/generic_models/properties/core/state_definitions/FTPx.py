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
Methods for setting up FTPx as the state variables in a generic property
package
"""
from pyomo.environ import Constraint, NonNegativeReals, Var, value

from idaes.core import (MaterialFlowBasis,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.generic_models.properties.core.generic.utility import \
    get_bounds_from_config, get_method, GenericPropertyPackageError
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog

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
    b.temperature = Var(initialize=t_init,
                        domain=NonNegativeReals,
                        bounds=t_bounds,
                        doc='State temperature',
                        units=units["temperature"])

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
    for p in b.params.phase_list:
        # Start with phase mole fractions equal to total mole fractions
        for j in b.components_in_phase(p):
            b.mole_frac_phase_comp[p, j].value = b.mole_frac_comp[j].value

        b.flow_mol_phase[p].value = value(
                        b.flow_mol / len(b.params.phase_list))

        # Try to refine guesses - Check phase type
        pobj = b.params.get_phase(p)

        if not hasattr(b.params, "_pe_pairs"):
            return

        if pobj.is_liquid_phase():
            tbub = None
            tdew = None
            for pp in b.params._pe_pairs:
                # Look for a VLE pair with this phase - should only be 1
                if ((pp[0] == p and
                     b.params.get_phase(pp[1]).is_vapor_phase()) or
                    (pp[1] == p and
                     b.params.get_phase(pp[0]).is_vapor_phase())):
                    # Get bubble and dew points
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
                    break

            if tbub is None and tdew is None:
                # No VLE pair found, or no bubble and dew point
                # Do nothing
                pass
            elif tdew is not None and b.temperature.value > tdew:
                # Pure vapour
                b.flow_mol_phase[p].value = value(1e-5*b.flow_mol)
                b.phase_frac[p].value = 1e-5

                for j in b.params.component_list:
                    if (p, j) in b.params._phase_component_set:
                        b.mole_frac_phase_comp[p, j].value = \
                            b._mole_frac_tdew[pp, j].value
            elif tbub is not None and b.temperature.value < tbub:
                # Pure liquid
                b.flow_mol_phase[p].value = value(b.flow_mol)
                b.phase_frac[p].value = 1

                for j in b.params.component_list:
                    if (p, j) in b.params._phase_component_set:
                        b.mole_frac_phase_comp[p, j].value = \
                            b.mole_frac_comp[j].value
            elif tbub is not None and tdew is not None:
                # Two-phase with bounds two-phase region
                # Thanks to Rahul Gandhi for the method
                vapRatio = value((b.temperature-tbub) / (tdew-tbub))

                b.flow_mol_phase["Liq"].value = value((1-vapRatio)*b.flow_mol)

                try:
                    for p2, j in b.params._phase_component_set:
                        if p2 == p:
                            psat_j = value(get_method(
                                b, "pressure_sat_comp", j)(
                                    b,
                                    b.params.get_component(j),
                                    b.temperature))
                            kfact = value(psat_j / b.pressure)

                            b.mole_frac_phase_comp["Liq", j].value = value(
                                b.mole_frac_comp[j]/(1+vapRatio*(kfact-1)))
                except GenericPropertyPackageError:
                    # No method for calculating Psat, use default values
                    pass
            else:
                # Two-phase, but with non-vaporizables and/or non-condensables
                pass

        elif pobj.is_vapor_phase():
            # Look for a VLE pair with this phase - will go with 1st found
            tbub = None
            tdew = None
            for pp in b.params._pe_pairs:
                if ((pp[0] == p and
                     b.params.get_phase(pp[1]).is_liquid_phase()) or
                    (pp[1] == p and
                     b.params.get_phase(pp[0]).is_liquid_phase())):
                    # Get bubble and dew points
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
                    break

            if tbub is None and tdew is None:
                # No VLE pair found, or no bubble and dew point
                # Do nothing
                pass
            elif tdew is not None and b.temperature.value > tdew:
                # Pure vapour
                b.flow_mol_phase[p].value = value(b.flow_mol)
                b.phase_frac[p].value = 1

                for j in b.params.component_list:
                    if (p, j) in b.params._phase_component_set:
                        b.mole_frac_phase_comp[p, j].value = \
                            b.mole_frac_comp[j].value
            elif tbub is not None and b.temperature.value < tbub:
                # Pure liquid
                b.flow_mol_phase[p].value = value(1e-5*b.flow_mol)
                b.phase_frac[p].value = 1e-5

                for j in b.params.component_list:
                    if (p, j) in b.params._phase_component_set:
                        b.mole_frac_phase_comp[p, j].value = \
                            b._mole_frac_tbub[pp, j].value
            elif tbub is not None and tdew is not None:
                # Two-phase with bounds two-phase region
                # Thanks to Rahul Gandhi for the method
                vapRatio = value((b.temperature-tbub) / (tdew-tbub))

                b.flow_mol_phase["Liq"].value = value((1-vapRatio)*b.flow_mol)

                try:
                    for p2, j in b.params._phase_component_set:
                        if p2 == p:
                            psat_j = value(get_method(
                                b, "pressure_sat_comp", j)(
                                    b,
                                    b.params.get_component(j),
                                    b.temperature))
                            kfact = value(psat_j / b.pressure)

                            b.mole_frac_phase_comp["Vap", j].value = value(
                                b.mole_frac_comp[j] /
                                (1+vapRatio*(kfact-1))*kfact)
                except GenericPropertyPackageError:
                    # No method for calculating Psat, use default values
                    pass
            else:
                # Two-phase, but with non-vaporizables and/or non-condensables
                pass


do_not_initialize = ["sum_mole_frac_out"]
