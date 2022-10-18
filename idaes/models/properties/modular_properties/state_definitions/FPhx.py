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
Methods for setting up FPhx as the state variables in a generic property
package
"""
from pyomo.environ import (
    Constraint,
    Expression,
    NonNegativeReals,
    Var,
    units as pyunits,
)

from idaes.core import MaterialFlowBasis, MaterialBalanceType, EnergyBalanceType
from idaes.models.properties.modular_properties.base.utility import (
    get_bounds_from_config,
)
from idaes.models.properties.modular_properties.state_definitions.FTPx import (
    state_initialization,
)
from .electrolyte_states import define_electrolyte_state, calculate_electrolyte_scaling
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

# Set up logger
_log = idaeslog.getLogger(__name__)


# TODO : Need a way to get a guess for T during initialization
def set_metadata(b):
    # Need to update metadata so that enth_mol is recorded as being part of the
    # state variables, and to ensure that getattr does not try to build it
    # using the default method.
    b.get_metadata().properties["enth_mol"] = {"method": None}


def define_state(b):
    # FTPx formulation always requires a flash, so set flag to True
    # TODO: should have some checking to make sure developers implement this properly
    b.always_flash = True

    # Check that only necessary state_bounds are defined
    expected_keys = ["flow_mol", "enth_mol", "temperature", "pressure"]
    if (
        b.params.config.state_bounds is not None
        and any(b.params.config.state_bounds.keys()) not in expected_keys
    ):
        for k in b.params.config.state_bounds.keys():
            if "mole_frac" in k:
                _log.warning(
                    "{} - found state_bounds argument for {}."
                    " Mole fraction bounds are set automatically and "
                    "this argument will be ignored.".format(b.name, k)
                )
            elif k not in expected_keys:
                raise ConfigurationError(
                    "{} - found unexpected state_bounds key {}. Please ensure "
                    "bounds are provided only for expected state variables "
                    "and that you have typed the variable names correctly.".format(
                        b.name, k
                    )
                )

    units = b.params.get_metadata().derived_units
    # Get bounds and initial values from config args
    f_bounds, f_init = get_bounds_from_config(b, "flow_mol", units.FLOW_MOLE)
    h_bounds, h_init = get_bounds_from_config(b, "enth_mol", units.ENERGY_MOLE)
    p_bounds, p_init = get_bounds_from_config(b, "pressure", units.PRESSURE)
    t_bounds, t_init = get_bounds_from_config(b, "temperature", units.TEMPERATURE)

    # Add state variables
    b.flow_mol = Var(
        initialize=f_init,
        domain=NonNegativeReals,
        bounds=f_bounds,
        doc=" Total molar flowrate",
        units=units.FLOW_MOLE,
    )
    b.mole_frac_comp = Var(
        b.component_list,
        bounds=(1e-20, 1.001),
        initialize=1 / len(b.component_list),
        doc="Mixture mole fractions",
        units=pyunits.dimensionless,
    )
    b.pressure = Var(
        initialize=p_init,
        domain=NonNegativeReals,
        bounds=p_bounds,
        doc="State pressure",
        units=units.PRESSURE,
    )

    b.enth_mol = Var(
        initialize=h_init,
        bounds=h_bounds,
        doc="State molar enthalpy",
        units=units.ENERGY_MOLE,
    )

    # Add supporting variables
    if f_init is None:
        fp_init = None
    else:
        fp_init = f_init / len(b.phase_list)

    b.flow_mol_phase = Var(
        b.phase_list,
        initialize=fp_init,
        domain=NonNegativeReals,
        bounds=f_bounds,
        doc="Phase molar flow rates",
        units=units.FLOW_MOLE,
    )

    b.mole_frac_phase_comp = Var(
        b.phase_component_set,
        initialize=1 / len(b.component_list),
        bounds=(1e-20, 1.001),
        doc="Phase mole fractions",
        units=pyunits.dimensionless,
    )

    def flow_mol_phase_comp_rule(b, p, j):
        return b.flow_mol_phase[p] * b.mole_frac_phase_comp[p, j]

    b.flow_mol_phase_comp = Expression(
        b.phase_component_set, rule=flow_mol_phase_comp_rule
    )

    b.temperature = Var(
        initialize=t_init,
        domain=NonNegativeReals,
        bounds=t_bounds,
        doc="Temperature",
        units=units.TEMPERATURE,
    )

    b.phase_frac = Var(
        b.phase_list,
        initialize=1 / len(b.phase_list),
        bounds=(0, None),
        doc="Phase fractions",
        units=pyunits.dimensionless,
    )

    # Add electrolye state vars if required
    # This must occur before adding the enthalpy constraint, as it needs true
    # species mole fractions
    if b.params._electrolyte:
        define_electrolyte_state(b)

    # Add supporting constraints
    if b.config.defined_state is False:
        # applied at outlet only
        b.sum_mole_frac_out = Constraint(
            expr=1 == sum(b.mole_frac_comp[i] for i in b.component_list)
        )

    def rule_enth_mol(b):
        return b.enth_mol == sum(
            b.enth_mol_phase[p] * b.phase_frac[p] for p in b.phase_list
        )

    b.enth_mol_eq = Constraint(
        rule=rule_enth_mol, doc="Total molar enthalpy mixing rule"
    )

    if len(b.phase_list) == 1:

        def rule_total_mass_balance(b):
            return b.flow_mol_phase[b.phase_list.first()] == b.flow_mol

        b.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return (
                b.mole_frac_comp[i] == b.mole_frac_phase_comp[b.phase_list.first(), i]
            )

        b.component_flow_balances = Constraint(
            b.component_list, rule=rule_comp_mass_balance
        )

        def rule_phase_frac(b, p):
            return b.phase_frac[p] == 1

        b.phase_fraction_constraint = Constraint(b.phase_list, rule=rule_phase_frac)

    elif len(b.phase_list) == 2:
        # For two phase, use Rachford-Rice formulation
        def rule_total_mass_balance(b):
            return sum(b.flow_mol_phase[p] for p in b.phase_list) == b.flow_mol

        b.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.flow_mol * b.mole_frac_comp[i] == sum(
                b.flow_mol_phase[p] * b.mole_frac_phase_comp[p, i]
                for p in b.phase_list
                if (p, i) in b.phase_component_set
            )

        b.component_flow_balances = Constraint(
            b.component_list, rule=rule_comp_mass_balance
        )

        def rule_mole_frac(b):
            return (
                sum(
                    b.mole_frac_phase_comp[b.phase_list.first(), i]
                    for i in b.component_list
                    if (b.phase_list.first(), i) in b.phase_component_set
                )
                - sum(
                    b.mole_frac_phase_comp[b.phase_list.last(), i]
                    for i in b.component_list
                    if (b.phase_list.last(), i) in b.phase_component_set
                )
                == 0
            )

        b.sum_mole_frac = Constraint(rule=rule_mole_frac)

        def rule_phase_frac(b, p):
            return b.phase_frac[p] * b.flow_mol == b.flow_mol_phase[p]

        b.phase_fraction_constraint = Constraint(b.phase_list, rule=rule_phase_frac)

    else:
        # Otherwise use a general formulation
        def rule_comp_mass_balance(b, i):
            return b.flow_mol * b.mole_frac_comp[i] == sum(
                b.flow_mol_phase[p] * b.mole_frac_phase_comp[p, i]
                for p in b.phase_list
                if (p, i) in b.phase_component_set
            )

        b.component_flow_balances = Constraint(
            b.component_list, rule=rule_comp_mass_balance
        )

        def rule_mole_frac(b, p):
            return (
                sum(
                    b.mole_frac_phase_comp[p, i]
                    for i in b.component_list
                    if (p, i) in b.phase_component_set
                )
                == 1
            )

        b.sum_mole_frac = Constraint(b.phase_list, rule=rule_mole_frac)

        def rule_phase_frac(b, p):
            return b.phase_frac[p] * b.flow_mol == b.flow_mol_phase[p]

        b.phase_fraction_constraint = Constraint(b.phase_list, rule=rule_phase_frac)

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

            eflow = b._enthalpy_flow_term = Expression(b.phase_list, rule=rule_eflow)
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
                b.phase_component_set, rule=rule_mdens
            )
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

            edens = b._energy_density_term = Expression(b.phase_list, rule=rule_edens)
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

    def define_state_vars_FPhx():
        """Define state vars."""
        return {
            "flow_mol": b.flow_mol,
            "mole_frac_comp": b.mole_frac_comp,
            "enth_mol": b.enth_mol,
            "pressure": b.pressure,
        }

    b.define_state_vars = define_state_vars_FPhx

    def define_display_vars_FPhx():
        """Define display vars."""
        return {
            "Total Molar Flowrate": b.flow_mol,
            "Total Mole Fraction": b.mole_frac_comp,
            "Molar Enthalpy": b.enth_mol,
            "Pressure": b.pressure,
        }

    b.define_display_vars = define_display_vars_FPhx


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
            f_init = pyunits.convert_value(
                f_bounds[1], from_units=f_bounds[3], to_units=units.FLOW_MOLE
            )
        else:
            f_init = f_bounds[1]
    except KeyError:
        f_init = 1

    try:
        p_bounds = state_bounds["pressure"]
        if len(p_bounds) == 4:
            p_init = pyunits.convert_value(
                p_bounds[1], from_units=p_bounds[3], to_units=units.PRESSURE
            )
        else:
            p_init = p_bounds[1]
    except KeyError:
        p_init = 1

    try:
        h_bounds = state_bounds["enth_mol"]
        if len(h_bounds) == 4:
            h_init = pyunits.convert_value(
                h_bounds[1], from_units=h_bounds[3], to_units=units.ENERGY_MOLE
            )
        else:
            h_init = h_bounds[1]
    except KeyError:
        h_init = 1

    try:
        t_bounds = state_bounds["temperature"]
        if len(t_bounds) == 4:
            t_init = pyunits.convert_value(
                t_bounds[1], from_units=t_bounds[3], to_units=units.TEMPERATURE
            )
        else:
            t_init = t_bounds[1]
    except KeyError:
        t_init = 1

    # Set default scaling factors
    b.set_default_scaling("flow_mol", 1 / f_init)
    b.set_default_scaling("flow_mol_phase", 1 / f_init)
    b.set_default_scaling("flow_mol_comp", 1 / f_init)
    b.set_default_scaling("flow_mol_phase_comp", 1 / f_init)
    b.set_default_scaling("pressure", 1 / p_init)
    b.set_default_scaling("temperature", 1 / t_init)
    b.set_default_scaling("enth_mol", 1 / h_init)


def calculate_scaling_factors(b):
    sf_flow = iscale.get_scaling_factor(b.flow_mol, default=1, warning=True)
    sf_mf = {}
    for i, v in b.mole_frac_phase_comp.items():
        sf_mf[i] = iscale.get_scaling_factor(v, default=1e3, warning=True)
    sf_h = iscale.get_scaling_factor(b.enth_mol, default=1e-4, warning=True)

    if b.config.defined_state is False:
        iscale.constraint_scaling_transform(
            b.sum_mole_frac_out, min(sf_mf.values()), overwrite=False
        )

    iscale.constraint_scaling_transform(b.enth_mol_eq, sf_h, overwrite=False)

    if len(b.phase_list) == 1:
        iscale.constraint_scaling_transform(
            b.total_flow_balance, sf_flow, overwrite=False
        )

        for j in b.component_list:
            sf_j = iscale.get_scaling_factor(
                b.mole_frac_comp[j], default=1e3, warning=True
            )
            iscale.constraint_scaling_transform(
                b.component_flow_balances[j], sf_j, overwrite=False
            )

        # b.phase_fraction_constraint is well scaled

    elif len(b.phase_list) == 2:
        iscale.constraint_scaling_transform(
            b.total_flow_balance, sf_flow, overwrite=False
        )

        for j in b.component_list:
            sf_j = iscale.get_scaling_factor(
                b.mole_frac_comp[j], default=1e3, warning=True
            )
            iscale.constraint_scaling_transform(
                b.component_flow_balances[j], sf_j * sf_flow, overwrite=False
            )

        iscale.constraint_scaling_transform(
            b.sum_mole_frac, min(sf_mf.values()), overwrite=False
        )

        for p in b.phase_list:
            iscale.constraint_scaling_transform(
                b.phase_fraction_constraint[p], sf_flow, overwrite=False
            )

    else:
        iscale.constraint_scaling_transform(
            b.total_flow_balance, sf_flow, overwrite=False
        )

        for j in b.component_list:
            sf_j = iscale.get_scaling_factor(
                b.mole_frac_comp[j], default=1e3, warning=True
            )
            iscale.constraint_scaling_transform(
                b.component_flow_balances[j], sf_j * sf_flow, overwrite=False
            )

        for p in b.phase_list:
            iscale.constraint_scaling_transform(
                b.sum_mole_frac[p], min(sf_mf[p, :].values()), overwrite=False
            )
            iscale.constraint_scaling_transform(
                b.phase_fraction_constraint[p], sf_flow, overwrite=False
            )

    if b.params._electrolyte:
        calculate_electrolyte_scaling(b)


# Inherit state_initialization from FTPX form, as the process is the same


do_not_initialize = ["sum_mole_frac_out"]


class FPhx(object):
    set_metadata = set_metadata
    define_state = define_state
    state_initialization = state_initialization
    do_not_initialize = do_not_initialize
    define_default_scaling_factors = define_default_scaling_factors
    calculate_scaling_factors = calculate_scaling_factors
