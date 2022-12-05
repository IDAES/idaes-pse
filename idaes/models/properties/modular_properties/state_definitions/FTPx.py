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
    get_method,
    GenericPropertyPackageError,
)
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    _valid_VL_component_list,
)
from idaes.models.properties.modular_properties.phase_equil.henry import (
    HenryType,
    henry_equilibrium_ratio,
)
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog
from .electrolyte_states import define_electrolyte_state, calculate_electrolyte_scaling
import idaes.core.util.scaling as iscale


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

    # Get bounds and initial values from config args
    f_bounds, f_init = get_bounds_from_config(b, "flow_mol", units.FLOW_MOLE)
    t_bounds, t_init = get_bounds_from_config(b, "temperature", units.TEMPERATURE)
    p_bounds, p_init = get_bounds_from_config(b, "pressure", units.PRESSURE)

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
    b.temperature = Var(
        initialize=t_init,
        domain=NonNegativeReals,
        bounds=t_bounds,
        doc="State temperature",
        units=units.TEMPERATURE,
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

    def Fpc_expr(b, p, j):
        return b.flow_mol_phase[p] * b.mole_frac_phase_comp[p, j]

    b.flow_mol_phase_comp = Expression(
        b.phase_component_set, rule=Fpc_expr, doc="Phase-component molar flowrates"
    )

    b.phase_frac = Var(
        b.phase_list,
        initialize=1 / len(b.phase_list),
        bounds=(0, None),
        doc="Phase fractions",
        units=pyunits.dimensionless,
    )

    # Add electrolye state vars if required
    if b.params._electrolyte:
        define_electrolyte_state(b)

    # Add supporting constraints
    if b.config.defined_state is False:
        # applied at outlet only
        b.sum_mole_frac_out = Constraint(
            expr=1 == sum(b.mole_frac_comp[i] for i in b.component_list)
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

    def define_state_vars_FTPx():
        """Define state vars."""
        return {
            "flow_mol": b.flow_mol,
            "mole_frac_comp": b.mole_frac_comp,
            "temperature": b.temperature,
            "pressure": b.pressure,
        }

    b.define_state_vars = define_state_vars_FTPx

    def define_display_vars_FTPx():
        """Define display vars."""
        return {
            "Total Molar Flowrate": b.flow_mol,
            "Total Mole Fraction": b.mole_frac_comp,
            "Temperature": b.temperature,
            "Pressure": b.pressure,
        }

    b.define_display_vars = define_display_vars_FTPx


def state_initialization(b):
    # Need to do sanity checking of mole fractions for ideal phase fraction
    # calculations
    for j in b.component_list:
        if value(b.mole_frac_comp[j]) < 0:
            raise ValueError(
                f"Component {j} has a negative "
                f"mole fraction in block {b.name}. "
                "Check your initialization."
            )

    for p in b.phase_list:
        # Start with phase mole fractions equal to total mole fractions
        for j in b.component_list:
            if (p, j) in b.phase_component_set:
                b.mole_frac_phase_comp[p, j].value = b.mole_frac_comp[j].value

        b.flow_mol_phase[p].value = value(b.phase_frac[p] * b.flow_mol)

    if not hasattr(b.params, "_pe_pairs"):
        return
    else:
        _pe_pairs = b.params._pe_pairs

    vl_comps = []
    henry_comps = []
    init_VLE = False
    for pp in _pe_pairs:
        # Look for a VLE pair with this phase - should only be 1
        if (
            b.params.get_phase(pp[0]).is_liquid_phase()
            and b.params.get_phase(pp[1]).is_vapor_phase()
        ) or (
            b.params.get_phase(pp[1]).is_liquid_phase()
            and b.params.get_phase(pp[0]).is_vapor_phase()
        ):
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
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                v_only_comps,
            ) = _valid_VL_component_list(b, pp)
            pp_VLE = pp

    if init_VLE:
        henry_mole_frac = []
        henry_conc = []
        henry_other = []
        K = {}

        for j in henry_comps:
            henry_type = b.params.get_component(j).config.henry_component[l_phase][
                "type"
            ]
            if henry_type in {HenryType.Hxp, HenryType.Kpx}:
                # Quick and dirty way of getting a Henry constant in the right
                # form/units is by temporarily setting its mole fraction to
                # one, then calling the henry_pressure method. Reset mole
                # fraction to original value to "do no harm"
                henry_mole_frac.append(j)
                K[j] = value(henry_equilibrium_ratio(b, l_phase, j))

            elif henry_type in {HenryType.Hcp, HenryType.Kpc}:
                # Can't evaluate concentration without density, can't evaluate
                # density without liquid phase concentration. Treat as
                # noncondensible and correct value during a second pass
                henry_conc.append(j)
            else:
                _log.warning(
                    "{} - Components with unsupported Henry's Law "
                    "type {} are present, and will be treated as "
                    "noncondensible by the vapor-liquid equilibrium "
                    "initialization for {} state variables.".format(
                        b.name, str(henry_type), str(b.params.config.state_definition)
                    )
                )
                henry_other.append(j)
        for j in vl_comps:
            try:
                K[j] = value(
                    get_method(b, "pressure_sat_comp", j)(
                        b, b.params.get_component(j), b.temperature
                    )
                    / b.pressure
                )
            except GenericPropertyPackageError:
                # No method for calculating Psat, use default values
                K = None
                break

    if init_VLE:
        raoult_init = False
        if tdew is not None and b.temperature.value > tdew:
            # Pure vapour
            vap_frac = 1 - 1e-5
        elif tbub is not None and b.temperature.value < tbub:
            # Pure liquid
            vap_frac = 1e-5
        elif tbub is not None and tdew is not None:
            # Two-phase with bounds two-phase region
            # Thanks to Rahul Gandhi for the method
            vap_frac = value(b.temperature - tbub) / (tdew - tbub)
        elif K is not None:
            raoult_init = True
            vap_frac = _modified_rachford_rice(
                b,
                K,
                vl_comps + henry_mole_frac,
                l_only_comps,
                v_only_comps + henry_conc + henry_other,
            )
        else:
            # No way to estimate phase fraction
            vap_frac = None

    if vap_frac is not None:
        b.phase_frac[v_phase] = vap_frac
        b.phase_frac[l_phase] = 1 - vap_frac

        for p in pp_VLE:
            b.flow_mol_phase[p].value = value(b.phase_frac[p] * b.flow_mol)

        # If dew/bubble calculations or MRR gives a nearly single-phase vapor
        # fraction, try to initialize mole fraction using dew/bubble comps
        if tbub is not None and vap_frac < 1e-4:
            # Pure liquid
            for j in b.component_list:
                if (l_phase, j) in b.phase_component_set:
                    b.mole_frac_phase_comp[l_phase, j].value = b.mole_frac_comp[j].value
                if (v_phase, j) in b.phase_component_set:
                    b.mole_frac_phase_comp[v_phase, j].value = b._mole_frac_tbub[
                        pp_VLE, j
                    ].value
        elif tdew is not None and vap_frac > 1 - 1e-4:
            # Pure Vapor
            for j in b.component_list:
                if (v_phase, j) in b.phase_component_set:
                    b.mole_frac_phase_comp[v_phase, j].value = b.mole_frac_comp[j].value
                if (l_phase, j) in b.phase_component_set:
                    b.mole_frac_phase_comp[l_phase, j].value = b._mole_frac_tdew[
                        pp_VLE, j
                    ].value
        elif K is not None:
            _set_mole_fractions_vle(
                b,
                K,
                vap_frac,
                l_phase,
                v_phase,
                vl_comps + henry_mole_frac,
                l_only_comps,
                v_only_comps + henry_conc + henry_other,
            )
            if len(henry_conc) > 0:
                # With initial guesses for the mole fractions, there is a
                # decent value for density to calculate the concentration with
                for j in henry_conc:
                    K[j] = value(henry_equilibrium_ratio(b, l_phase, j))

                # Only recompute vapor fraction if we didn't get it from
                # dew and bubble point calculations
                if raoult_init:
                    vap_frac = _modified_rachford_rice(
                        b,
                        K,
                        vl_comps + henry_mole_frac + henry_conc,
                        l_only_comps,
                        v_only_comps + henry_other,
                    )
                    # If MRR failed with K[j] set for j in henry_conc, just
                    # leave the preliminary values assigned
                    if vap_frac is not None:
                        b.phase_frac[v_phase] = vap_frac
                        b.phase_frac[l_phase] = 1 - vap_frac

                        for p in pp_VLE:
                            b.flow_mol_phase[p].value = value(
                                b.phase_frac[p] * b.flow_mol
                            )
                        _set_mole_fractions_vle(
                            b,
                            K,
                            vap_frac,
                            l_phase,
                            v_phase,
                            vl_comps + henry_mole_frac + henry_conc,
                            l_only_comps,
                            v_only_comps + henry_other,
                        )
                else:
                    # Assign mole fractions for vap_frac from Rahul Gahndi's
                    # method
                    _set_mole_fractions_vle(
                        b,
                        K,
                        vap_frac,
                        l_phase,
                        v_phase,
                        vl_comps + henry_mole_frac + henry_conc,
                        l_only_comps,
                        v_only_comps + henry_other,
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
    sf_flow = iscale.get_scaling_factor(b.flow_mol, default=1, warning=True)
    sf_mf = {}
    for i, v in b.mole_frac_phase_comp.items():
        sf_mf[i] = iscale.get_scaling_factor(v, default=1e3, warning=True)

    if b.config.defined_state is False:
        iscale.constraint_scaling_transform(
            b.sum_mole_frac_out, min(sf_mf.values()), overwrite=False
        )

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


do_not_initialize = ["sum_mole_frac_out"]


class FTPx(object):
    set_metadata = set_metadata
    define_state = define_state
    state_initialization = state_initialization
    do_not_initialize = do_not_initialize
    define_default_scaling_factors = define_default_scaling_factors
    calculate_scaling_factors = calculate_scaling_factors


def _set_mole_fractions_vle(
    b, K, vap_frac, l_phase, v_phase, vl_comps, l_only_comps, v_only_comps
):
    # The only reason this method exists is so we can use it twice because
    # concentration based Henry law coefficients make everything difficult
    for j in l_only_comps:
        b.mole_frac_phase_comp[l_phase, j].value = (
            1 / (1 - vap_frac) * b.mole_frac_comp[j]
        )

    for j in v_only_comps:
        b.mole_frac_phase_comp[v_phase, j].value = 1 / vap_frac * b.mole_frac_comp[j]

    for j in vl_comps:
        b.mole_frac_phase_comp[l_phase, j].value = value(
            b.mole_frac_comp[j] / (1 + vap_frac * (K[j] - 1))
        )
        b.mole_frac_phase_comp[v_phase, j].value = value(
            b.mole_frac_comp[j] * K[j] / (K[j] + 1 - vap_frac)
        )


def _modified_rachford_rice(b, K, vl_comps, l_only_comps, v_only_comps, eps=1e-5):
    """
    Uses a modified version of the Rachford Rice method to compute the vapor
    fraction for a two-phase vapor liquid equilibrium.

    Arguments:
        b: Property block for which this calculation is taking place
        K: Dictionary of equilibrium constants in the relationship
            y_j = K_j x_j for all j in vl_comps.
        vl_comps: List of components present in both phases
        l_only_comps: List of components present in only the liquid phase
        v_only_comps: List of components present in only the vapor phase
        eps (Optional): Offset from 0 or 1 in the event of a single phase
            solution

    Returns:
        phase_frac: Either a float between 0 and 1 if a root is found or
            RR predicts a single phase solution, or None if no root is found
            or some K_j < 0. Roots outside the range of 0 and 1 are clipped to
            be within that range.

    Notes:
        Mathematically, one can show that this method is guaranteed to give a
        meaningful answer given that:
            1. All K > 0 (any component with K_j = 0 should
                                            be designated an l_only_comp)
            2. All mole_frac_comp[j] > 0 in block b
            3. sum(mole_frac_comp[j] for j in b.component_list) == 1
        If the harmonic mean of K[j] is greater than one, the stream is pure
        vapor. If the arithmatic mean of K[j] is less than one, the stream is
        pure liquid. If neither of those conditions hold, a root exists between
        zero and one and the root-finding method will converge to it.

        However, validating user input to ensure that all these conditions
        hold turned out to be more trouble than it was worth (frequently users
        end up initializing with a sum of mole fractions less than one in
        reasonable use cases). Additionally, ill-conditioned problems open
        a gap between mathematical certainty and the results of floating point
        arithmetic.
    """
    # Validate split ratios K are nonnegative
    for j in vl_comps:
        if K[j] < 0:
            _log.warning(
                "While initializing block {}, the vapor/liquid split "
                "ratio of Component {} was calculated to be negative. "
                "Check the implementation of the saturation pressure, "
                "Henry's law method, or liquid density.".format(b.name, j)
            )
            return None

    # Calculate harmonic and arithmatic means of the split ratios KS
    if len(l_only_comps) == 0:
        tmp = sum([value(b.mole_frac_comp[j]) / K[j] for j in vl_comps])
        if tmp < 1e-14:
            K_bar_H = float("inf")
        else:
            K_bar_H = 1 / tmp
    else:
        K_bar_H = 0
    if len(v_only_comps) == 0:
        K_bar_A = sum([value(b.mole_frac_comp[j]) * K[j] for j in vl_comps])
    else:
        K_bar_A = float("inf")

    # If neither of these conditions hold, RR gives a single-phase solution
    if K_bar_H >= 1:
        # Vapor fraction is nearly one
        return 1 - eps
    elif K_bar_A <= 1:
        # Vapor fraction is nearly zero
        return eps

    # I discovered this method is a varient on solving the
    # Rachford-Rice equation, which has been known to ChemEs
    # since the 50s. I'm not sure if the link between Newton's
    # method and convexity was explicitly stated, though. For good
    # reason, the classic RR equation does not appear to necessarily
    # be convex

    # There are varients where the flash calculation is robust
    # to the calculated vapor fraction being greater than 1 or
    # less than zero. Future reader, check out
    # https://doi.org/10.1016/j.fluid.2011.12.005 and
    # https://doi.org/10.1016/S0378-3812(97)00179-9 if it
    # becomes necessary to increase robustness to very nonideal
    # systems.

    # Equation defining vapor fraction for an ideal mixture.
    # This equation has two roots: a trivial one at vap_frac = 1
    # and a nontrivial one for some 0<vap_frac<1. B is a
    # convex function
    def B(vap_frac):
        return (
            sum(
                [
                    value(b.mole_frac_comp[j] * K[j] / (1 + vap_frac * (K[j] - 1)))
                    for j in vl_comps
                ]
            )
            + sum([value(b.mole_frac_comp[j] / vap_frac) for j in v_only_comps])
            - 1
        )

    def dB_dvap_frac(vap_frac):
        return sum(
            [
                value(
                    -b.mole_frac_comp[j]
                    * K[j]
                    * (K[j] - 1)
                    / ((1 + vap_frac * (K[j] - 1)) ** 2)
                )
                for j in vl_comps
            ]
        ) - sum([value(b.mole_frac_comp[j] / vap_frac**2) for j in v_only_comps])

    # Newton's method will always undershoot the root of a
    # convex equation if one starts from a value from which
    # the function is positive. However, there is a singularity
    # at vap_frac=0 if there are noncondensable components.
    # Therefore, use a binary search to find a value where
    # B(vap_frac) is positive
    k = 0
    vap_frac = 0.5
    mRR_failed = False
    while B(vap_frac) < 0:
        vap_frac /= 2
        k += 1
        if k > 40:
            mRR_failed = True
            break
    # Now use Newton's method to calculate vap_frac
    k = 0
    if not mRR_failed:
        while abs(B(vap_frac)) > 1e-6:
            vap_frac -= B(vap_frac) / dB_dvap_frac(vap_frac)
            k += 1
            if k > 40:
                mRR_failed = True
                break
    if mRR_failed:
        # mRR did not even converge to a root. Do not attempt to assign
        # a vapor fraction
        vap_frac = None
    elif vap_frac < 0:
        # mRR converged to a negative root. Initialize as mostly liquid
        vap_frac = eps
        mRR_failed = True
    elif vap_frac > 1:
        # mRR converged to a root > 1. Initialize as mostly vapor
        vap_frac = 1 - eps
        mRR_failed = True
    # Regardless of the failure mode, log a warning
    if mRR_failed:
        _log.warning(
            "Block {} - phase faction initialization using "
            "modified Rachford-Rice failed. This could be "
            "because a component is essentially "
            "nonvolatile or noncondensible, or "
            "because mole fractions sum to more than "
            "one.".format(b.name)
        )
    return vap_frac
