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
Methods for setting up FpcTP as the state variables in a generic property
package
"""
from pyomo.environ import (
    Constraint,
    Expression,
    NonNegativeReals,
    Var,
    value,
    units as pyunits,
)

from idaes.core import MaterialFlowBasis, MaterialBalanceType, EnergyBalanceType
from idaes.models.properties.modular_properties.base.utility import (
    get_bounds_from_config,
)
from .electrolyte_states import define_electrolyte_state, calculate_electrolyte_scaling
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

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
    expected_keys = ["flow_mol_phase_comp", "enth_mol", "temperature", "pressure"]
    if (
        b.params.config.state_bounds is not None
        and any(b.params.config.state_bounds.keys()) not in expected_keys
    ):
        for k in b.params.config.state_bounds.keys():
            if k not in expected_keys:
                raise ConfigurationError(
                    "{} - found unexpected state_bounds key {}. Please ensure "
                    "bounds are provided only for expected state variables "
                    "and that you have typed the variable names correctly.".format(
                        b.name, k
                    )
                )

    units = b.params.get_metadata().derived_units
    # Get bounds and initial values from config args
    f_bounds, f_init = get_bounds_from_config(b, "flow_mol_phase_comp", units.FLOW_MOLE)
    t_bounds, t_init = get_bounds_from_config(b, "temperature", units.TEMPERATURE)
    p_bounds, p_init = get_bounds_from_config(b, "pressure", units.PRESSURE)

    # Add state variables
    b.flow_mol_phase_comp = Var(
        b.phase_component_set,
        initialize=f_init,
        domain=NonNegativeReals,
        bounds=f_bounds,
        doc="Phase-component molar flowrate",
        units=units.FLOW_MOLE,
    )
    b.pressure = Var(
        initialize=p_init,
        domain=NonNegativeReals,
        bounds=p_bounds,
        doc="State pressure",
        units=units.PRESSURE,
    )
    b.temperature = Var(
        initialize=t_init,
        domain=NonNegativeReals,
        bounds=t_bounds,
        doc="State temperature",
        units=units.TEMPERATURE,
    )

    # Add supporting variables
    b.flow_mol = Expression(
        expr=sum(b.flow_mol_phase_comp[i] for i in b.phase_component_set),
        doc="Total molar flowrate",
    )

    def flow_mol_phase(b, p):
        return sum(
            b.flow_mol_phase_comp[p, j]
            for j in b.component_list
            if (p, j) in b.phase_component_set
        )

    b.flow_mol_phase = Expression(
        b.phase_list, rule=flow_mol_phase, doc="Phase molar flow rates"
    )

    def rule_flow_mol_comp(b, j):
        return sum(
            b.flow_mol_phase_comp[p, j]
            for p in b.phase_list
            if (p, j) in b.phase_component_set
        )

    b.flow_mol_comp = Expression(
        b.component_list, rule=rule_flow_mol_comp, doc="Component molar flow rates"
    )

    def mole_frac_comp(b, j):
        return (
            sum(
                b.flow_mol_phase_comp[p, j]
                for p in b.phase_list
                if (p, j) in b.phase_component_set
            )
            / b.flow_mol
        )

    b.mole_frac_comp = Expression(
        b.component_list, rule=mole_frac_comp, doc="Mixture mole fractions"
    )

    b.mole_frac_phase_comp = Var(
        b.phase_component_set,
        bounds=(1e-20, 1.001),
        initialize=1 / len(b.component_list),
        doc="Phase mole fractions",
        units=pyunits.dimensionless,
    )

    def rule_mole_frac_phase_comp(b, p, j):
        # Calcualting mole frac phase comp is degenerate if there is only one
        # component in phase.
        # Count components
        comp_count = 0
        for p1, j1 in b.phase_component_set:
            if p1 == p:
                comp_count += 1

        if comp_count > 1:
            return (
                b.mole_frac_phase_comp[p, j] * b.flow_mol_phase[p]
                == b.flow_mol_phase_comp[p, j]
            )
        else:
            return b.mole_frac_phase_comp[p, j] == 1

    b.mole_frac_phase_comp_eq = Constraint(
        b.phase_component_set, rule=rule_mole_frac_phase_comp
    )

    def rule_phase_frac(b, p):
        if len(b.phase_list) == 1:
            return 1
        else:
            return b.flow_mol_phase[p] / b.flow_mol

    b.phase_frac = Expression(b.phase_list, rule=rule_phase_frac, doc="Phase fractions")

    # Add electrolye state vars if required
    if b.params._electrolyte:
        define_electrolyte_state(b)

    # -------------------------------------------------------------------------
    # General Methods
    def get_material_flow_terms_FpcTP(p, j):
        """Create material flow terms for control volume."""
        return b.flow_mol_phase_comp[p, j]

    b.get_material_flow_terms = get_material_flow_terms_FpcTP

    def get_enthalpy_flow_terms_FpcTP(p):
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

    b.get_enthalpy_flow_terms = get_enthalpy_flow_terms_FpcTP

    def get_material_density_terms_FpcTP(p, j):
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

    b.get_material_density_terms = get_material_density_terms_FpcTP

    def get_energy_density_terms_FpcTP(p):
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
        return {
            "flow_mol_phase_comp": b.flow_mol_phase_comp,
            "temperature": b.temperature,
            "pressure": b.pressure,
        }

    b.define_state_vars = define_state_vars_FpcTP

    def define_display_vars_FpcTP():
        """Define display vars."""
        return {
            "Molar Flowrate": b.flow_mol_phase_comp,
            "Temperature": b.temperature,
            "Pressure": b.pressure,
        }

    b.define_display_vars = define_display_vars_FpcTP


def state_initialization(b):
    for i in b.phase_component_set:
        b.mole_frac_phase_comp[i].value = value(
            b.flow_mol_phase_comp[i] / b.flow_mol_phase[i[0]]
        )


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
        f_bounds = state_bounds["flow_mol_phase_comp"]
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


def calculate_scaling_factors(b):
    for p, j in b.phase_component_set:
        sf = iscale.get_scaling_factor(
            b.flow_mol_phase_comp[p, j], default=1, warning=True
        )
        iscale.constraint_scaling_transform(
            b.mole_frac_phase_comp_eq[p, j], sf, overwrite=False
        )

    if b.params._electrolyte:
        calculate_electrolyte_scaling(b)


do_not_initialize = []


class FpcTP(object):
    set_metadata = set_metadata
    define_state = define_state
    state_initialization = state_initialization
    do_not_initialize = do_not_initialize
    define_default_scaling_factors = define_default_scaling_factors
    calculate_scaling_factors = calculate_scaling_factors
