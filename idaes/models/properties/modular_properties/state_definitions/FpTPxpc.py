#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Methods for setting up FpTPxpc as the state variables in a generic property
package
"""
# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

# TODO: Look into protected access issues
# pylint: disable=protected-access

from types import MethodType
from pyomo.environ import (
    Constraint,
    Expression,
    NonNegativeReals,
    Var,
    value,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import MaterialFlowBasis, MaterialBalanceType, EnergyBalanceType
from idaes.models.properties.modular_properties.base.utility import (
    get_bounds_from_config,
    get_method,
    GenericPropertyPackageError,
)
from idaes.models.properties.modular_properties.base.utility import (
    identify_VL_component_list,
)
from idaes.models.properties.modular_properties.phase_equil.henry import (
    HenryType,
    henry_equilibrium_ratio,
)
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from .electrolyte_states import define_electrolyte_state, calculate_electrolyte_scaling


# Set up logger
_log = idaeslog.getLogger(__name__)


def set_metadata(blk):
    # This is the default assumption for state vars, so we don't need to change
    # the metadata dict
    pass


def define_state(blk):
    blk.always_flash = False

    units = blk.params.get_metadata().derived_units

    # Check that only necessary state_bounds are defined
    expected_keys = ["flow_mol_phase", "temperature", "pressure"]
    if (
        blk.params.config.state_bounds is not None
        and any(blk.params.config.state_bounds.keys()) not in expected_keys
    ):
        for k in blk.params.config.state_bounds.keys():
            if "mole_frac" in k:
                _log.warning(
                    "{} - found state_bounds argument for {}."
                    " Mole fraction bounds are set automatically and "
                    "this argument will be ignored.".format(blk.name, k)
                )
            elif k not in expected_keys:
                raise ConfigurationError(
                    "{} - found unexpected state_bounds key {}. Please ensure "
                    "bounds are provided only for expected state variables "
                    "and that you have typed the variable names correctly.".format(
                        blk.name, k
                    )
                )

    # Get bounds and initial values from config args

    f_bounds, f_init = get_bounds_from_config(blk, "flow_mol_phase", units.FLOW_MOLE)
    t_bounds, t_init = get_bounds_from_config(blk, "temperature", units.TEMPERATURE)
    p_bounds, p_init = get_bounds_from_config(blk, "pressure", units.PRESSURE)

    # Add state variables
    blk.flow_mol_phase = Var(
        blk.phase_list,
        initialize=f_init,
        domain=NonNegativeReals,
        bounds=f_bounds,
        doc="Total molar flowrate",
        units=units.FLOW_MOLE,
    )
    blk.mole_frac_phase_comp = Var(
        blk.phase_component_set,
        bounds=(1e-20, 1.001),
        initialize=1 / len(blk.component_list),
        doc="Mixture mole fractions",
        units=pyunits.dimensionless,
    )
    blk.pressure = Var(
        initialize=p_init,
        domain=NonNegativeReals,
        bounds=p_bounds,
        doc="State pressure",
        units=units.PRESSURE,
    )
    blk.temperature = Var(
        initialize=t_init,
        domain=NonNegativeReals,
        bounds=t_bounds,
        doc="State temperature",
        units=units.TEMPERATURE,
    )

    # Add supporting variables and expressions

    @blk.Expression(doc="Total molar flow")
    def flow_mol(b):
        return sum(b.flow_mol_phase[p] for p in b.phase_list)

    blk.mole_frac_comp = Var(
        blk.component_list,
        initialize=1 / len(blk.component_list),
        bounds=(1e-20, 1.001),
        doc="Mole fractions",
        units=pyunits.dimensionless,
    )

    @blk.Expression(blk.phase_component_set, doc="Phase-component molar flowrates")
    def flow_mol_phase_comp(b, p, j):
        return b.flow_mol_phase[p] * b.mole_frac_phase_comp[p, j]

    blk.phase_frac = Var(
        blk.phase_list,
        initialize=1.0 / len(blk.phase_list),
        bounds=(0, None),
        doc="Phase fractions",
        units=pyunits.dimensionless,
    )

    # Add electrolyte state vars if required
    if blk.params._electrolyte:
        # TODO do we need special handling here?
        # raise NotImplementedError()
        define_electrolyte_state(blk)

    # Add supporting constraints
    if blk.config.defined_state is False:
        # applied at outlet only
        @blk.Constraint(blk.phase_list)
        def sum_mole_frac_out(b, p):
            return 1.0 == sum(blk.mole_frac_phase_comp[p, j]
                              for j in b.component_list
                              if (p, j) in b.phase_component_set)

    @blk.Constraint(blk.component_list, doc="Defines mole_frac_comp")
    def mole_frac_comp_eq(b, j):
        return b.mole_frac_comp[j] * b.flow_mol == sum(
            b.mole_frac_phase_comp[p, j] * b.flow_mol_phase[p]
            for p in b.phase_list
            if (p, j) in b.phase_component_set
        )

    if len(blk.phase_list) == 1:
        @blk.Constraint(blk.phase_list, doc="Defines phase_frac")
        def phase_frac_eqn(b, p):
            return b.phase_frac[p] == 1.0

    else:
        def rule_phase_frac(b, p):
            return b.phase_frac[p] * b.flow_mol == b.flow_mol_phase[p]

        blk.phase_frac_eqn = Constraint(blk.phase_list, rule=rule_phase_frac)

    # -------------------------------------------------------------------------
    # General Methods
    def get_material_flow_terms_FpTPxpc(b, p, j):
        """Create material flow terms for control volume."""
        return b.flow_mol_phase_comp[p, j]

    blk.get_material_flow_terms = MethodType(get_material_flow_terms_FpTPxpc, blk)

    def get_enthalpy_flow_terms_FpTPxpc(b, p):
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

    blk.get_enthalpy_flow_terms = MethodType(get_enthalpy_flow_terms_FpTPxpc, blk)

    def get_material_density_terms_FpTPxpc(b, p, j):
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

    blk.get_material_density_terms = MethodType(get_material_density_terms_FpTPxpc, blk)

    def get_energy_density_terms_FpTPxpc(b, p):
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

    blk.get_energy_density_terms = MethodType(get_energy_density_terms_FpTPxpc, blk)

    def default_material_balance_type_FpTPxpc():
        return MaterialBalanceType.componentTotal

    blk.default_material_balance_type = default_material_balance_type_FpTPxpc

    def default_energy_balance_type_FpTPxpc():
        return EnergyBalanceType.enthalpyTotal

    blk.default_energy_balance_type = default_energy_balance_type_FpTPxpc

    def get_material_flow_basis_FpTPxpc():
        return MaterialFlowBasis.molar

    blk.get_material_flow_basis = get_material_flow_basis_FpTPxpc

    def define_state_vars_FpTPxpc(b):
        """Define state vars."""
        return {
            "flow_mol_phase": b.flow_mol_phase,
            "mole_frac_phase_comp": b.mole_frac_phase_comp,
            "temperature": b.temperature,
            "pressure": b.pressure,
        }

    blk.define_state_vars = MethodType(define_state_vars_FpTPxpc, blk)

    def define_display_vars_FpTPxpc(b):
        """Define display vars."""
        return {
            "Total Molar Flowrate": b.flow_mol,
            "Total Mole Fraction": b.mole_frac_comp,
            "Temperature": b.temperature,
            "Pressure": b.pressure,
        }

    blk.define_display_vars = MethodType(define_display_vars_FpTPxpc, blk)

# TODO flash initialization---needs to be a separate method because
# we need to build a ControlVolume0D to enforce material balances
def state_initialization(b):
    """
    Initialize state variables that need to be initialize irrespective
    of whether or not there's a flash.
    """
    for j in b.component_list:
        # Linear in mole_frac_comp when all other variables
        # are fixed---this should converge in one iteration
        calculate_variable_from_constraint(
            b.mole_frac_comp[j],
            b.mole_frac_comp_eq[j]
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
        f_bounds = state_bounds["flow_mol_phase"]
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


def calculate_scaling_factors(blk):
    sf_flow_phase = {}
    for p, v in blk.flow_mol_phase.items():
        sf_flow_phase[p] = iscale.get_scaling_factor(blk.flow_mol_phase[p], default=1, warning=True)
    sf_mf = {}
    for i, v in blk.mole_frac_phase_comp.items():
        sf_mf[i] = iscale.get_scaling_factor(v, default=1e3, warning=True)

    for i in blk.mole_frac_comp:
        sf = iscale.get_scaling_factor(blk.mole_frac_comp[i])
        if sf is None:
            sf = min(sf_mf[p, j] for j, p in blk.phase_component_set if i==j)
            iscale.set_scaling_factor(sf, blk.mole_frac_comp[i])
        iscale.constraint_scaling_transform(blk.mole_frac_comp_eq[i], sf, overwrite=False)

    if blk.config.defined_state is False:
        for p in blk.phase_list:
            sf_eqn = min(sf_mf[p, j] for j in blk.components_in_phase(p))
            iscale.constraint_scaling_transform(
                blk.sum_mole_frac_out[p], sf_eqn, overwrite=False
            )

    if len(blk.phase_list) == 1:
        for p in blk.phase_list:
            iscale.constraint_scaling_transform(
                blk.phase_frac_eqn[p],
                1,
                overwrite=False
            )
    else:
        for p in blk.phase_list:
            iscale.constraint_scaling_transform(
                blk.phase_frac_eqn[p],
                sf_flow_phase[p],
                overwrite=False
            )

    if blk.params._electrolyte:
        # raise NotImplementedError()
        calculate_electrolyte_scaling(blk)


do_not_initialize = ["sum_mole_frac_out"]


class FpTPxpc(object):
    """Total flow, temperature, pressure, mole fraction state."""

    set_metadata = set_metadata
    define_state = define_state
    state_initialization = state_initialization
    do_not_initialize = do_not_initialize
    define_default_scaling_factors = define_default_scaling_factors
    calculate_scaling_factors = calculate_scaling_factors
