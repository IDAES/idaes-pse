#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""Generic Helmholtz EOS StateBlock Class
"""

__author__ = "John Eslick"

import copy

import pyomo.environ as pyo
from pyomo.common.collections import ComponentSet
from pyomo.environ import units as pyunits

from idaes.core.util.math import smooth_max
from idaes.core import declare_process_block_class
from idaes.core import (
    StateBlock,
    StateBlockData,
    MaterialBalanceType,
    EnergyBalanceType,
    MaterialFlowBasis,
)
from idaes.models.properties.general_helmholtz.helmholtz_functions import (
    add_helmholtz_external_functions,
    HelmholtzThermoExpressions,
    AmountBasis,
    PhaseType,
    StateVars,
    _data_dir,
)
from idaes.models.properties.general_helmholtz.components import (
    viscosity_available,
    thermal_conductivity_available,
    surface_tension_available,
)
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.initialization import fix_state_vars
from idaes.core.initialization.initializer_base import (
    InitializerBase,
    InitializationStatus,
)

_log = idaeslog.getLogger(__name__)


class HelmholtzEoSInitializer(InitializerBase):
    """
    Initializer object for Helmholtz EoS packages using external functions.

    Due to the use of external functions, Helmholtz EoS StateBlocks have no constraints,
    thus there is nothing to initialize. This Initializer replaces the general initialize
    method with a no-op.

    """

    CONFIG = InitializerBase.CONFIG()

    def initialize(
        self,
        model: pyo.Block,
        output_level=None,
    ):
        """
        Initialize method for Helmholtz EoS state blocks. This is a no-op.

        Args:
            model: model to be initialized
            output_level: (optional) output level to use during initialization run (overrides global setting).

        Returns:
            InitializationStatus.Ok
        """
        self._update_summary(model, "status", InitializationStatus.Ok)
        return self.summary[model]["status"]


class _StateBlock(StateBlock):
    """
    This class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    # Set default initializer
    default_initializer = HelmholtzEoSInitializer

    @staticmethod
    def _set_fixed(v, f):
        if f:
            v.fix()
        else:
            v.unfix()

    @staticmethod
    def _set_not_fixed(v, state, key, hold):
        if state is not None:
            if not v.fixed:
                try:
                    v.value = state[key]
                except KeyError:
                    pass
        if hold:
            v.fix()

    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)

    def initialize(self, *args, **kwargs):
        flags = {}
        hold_state = kwargs.pop("hold_state", False)
        state_args = kwargs.pop("state_args", None)
        for i, v in self.items():
            pp = self[i].config.parameters.config.phase_presentation
            sv = self[i].state_vars
            ab = self[i].amount_basis
            if sv == StateVars.PH and ab == AmountBasis.MOLE:
                flags[i] = (v.flow_mol.fixed, v.enth_mol.fixed, v.pressure.fixed)
                self._set_not_fixed(v.flow_mol, state_args, "flow_mol", hold_state)
                self._set_not_fixed(v.enth_mol, state_args, "enth_mol", hold_state)
                self._set_not_fixed(v.pressure, state_args, "pressure", hold_state)
            elif sv == StateVars.PH and ab == AmountBasis.MASS:
                flags[i] = (v.flow_mass.fixed, v.enth_mass.fixed, v.pressure.fixed)
                self._set_not_fixed(v.flow_mass, state_args, "flow_mass", hold_state)
                self._set_not_fixed(v.enth_mass, state_args, "enth_mass", hold_state)
                self._set_not_fixed(v.pressure, state_args, "pressure", hold_state)
            elif sv == StateVars.PS and ab == AmountBasis.MOLE:
                flags[i] = (v.flow_mol.fixed, v.entr_mol.fixed, v.pressure.fixed)
                self._set_not_fixed(v.flow_mol, state_args, "flow_mol", hold_state)
                self._set_not_fixed(v.entr_mol, state_args, "enth_mol", hold_state)
                self._set_not_fixed(v.pressure, state_args, "pressure", hold_state)
            elif sv == StateVars.PS and ab == AmountBasis.MASS:
                flags[i] = (v.flow_mass.fixed, v.entr_mass.fixed, v.pressure.fixed)
                self._set_not_fixed(v.flow_mass, state_args, "flow_mass", hold_state)
                self._set_not_fixed(v.entr_mass, state_args, "enth_mass", hold_state)
                self._set_not_fixed(v.pressure, state_args, "pressure", hold_state)
            elif sv == StateVars.PU and ab == AmountBasis.MOLE:
                flags[i] = (
                    v.flow_mol.fixed,
                    v.energy_internal_mol.fixed,
                    v.pressure.fixed,
                )
                self._set_not_fixed(
                    v.energy_internal_mol, state_args, "energy_internal_mol", hold_state
                )
                self._set_not_fixed(
                    v.energy_internal_mol, state_args, "energy_internal_mol", hold_state
                )
                self._set_not_fixed(
                    v.energy_internal_mol, state_args, "energy_internal_mol", hold_state
                )
            elif sv == StateVars.PU and ab == AmountBasis.MASS:
                flags[i] = (
                    v.flow_mass.fixed,
                    v.energy_internal_mass.fixed,
                    v.pressure.fixed,
                )
                self._set_not_fixed(
                    v.energy_internal_mass,
                    state_args,
                    "energy_internal_mass",
                    hold_state,
                )
                self._set_not_fixed(
                    v.energy_internal_mass,
                    state_args,
                    "energy_internal_mass",
                    hold_state,
                )
                self._set_not_fixed(v.pressure, state_args, "pressure", hold_state)
            elif sv == StateVars.TPX and ab == AmountBasis.MOLE:
                # Hold the T-P-x vars
                if pp in (PhaseType.MIX, PhaseType.LG):
                    flags[i] = (
                        v.flow_mol.fixed,
                        v.temperature.fixed,
                        v.pressure.fixed,
                        v.vapor_frac.fixed,
                    )
                    self._set_not_fixed(v.flow_mol, state_args, "flow_mol", hold_state)
                    self._set_not_fixed(
                        v.temperature, state_args, "temperature", hold_state
                    )
                    self._set_not_fixed(v.pressure, state_args, "pressure", hold_state)
                    self._set_not_fixed(
                        v.vapor_frac, state_args, "vapor_frac", hold_state
                    )
                else:
                    flags[i] = (v.flow_mol.fixed, v.temperature.fixed, v.pressure.fixed)
                    self._set_not_fixed(v.flow_mol, state_args, "flow_mol", hold_state)
                    self._set_not_fixed(
                        v.temperature, state_args, "temperature", hold_state
                    )
                    self._set_not_fixed(v.pressure, state_args, "pressure", hold_state)
            elif sv == StateVars.TPX and ab == AmountBasis.MASS:
                # Hold the T-P-x vars
                if pp in (PhaseType.MIX, PhaseType.LG):
                    flags[i] = (
                        v.flow_mass.fixed,
                        v.temperature.fixed,
                        v.pressure.fixed,
                        v.vapor_frac.fixed,
                    )
                    self._set_not_fixed(
                        v.flow_mass, state_args, "flow_mass", hold_state
                    )
                    self._set_not_fixed(
                        v.temperature, state_args, "temperature", hold_state
                    )
                    self._set_not_fixed(v.pressure, state_args, "pressure", hold_state)
                    self._set_not_fixed(
                        v.vapor_frac, state_args, "vapor_frac", hold_state
                    )
                else:
                    flags[i] = (
                        v.flow_mass.fixed,
                        v.temperature.fixed,
                        v.pressure.fixed,
                    )
                    self._set_not_fixed(
                        v.flow_mass, state_args, "flow_mass", hold_state
                    )
                    self._set_not_fixed(
                        v.temperature, state_args, "temperature", hold_state
                    )
                    self._set_not_fixed(v.pressure, state_args, "pressure", hold_state)

        return flags

    def release_state(self, flags, **kwargs):
        """Set the variables back to there original fixed/free state
        after initialization.

        Args:
            flags (dict): Original variable states
        """
        for i, f in flags.items():
            pp = self[i].config.parameters.config.phase_presentation
            sv = self[i].state_vars
            ab = self[i].amount_basis
            if sv == StateVars.PH and ab == AmountBasis.MOLE:
                self._set_fixed(self[i].flow_mol, f[0])
                self._set_fixed(self[i].enth_mol, f[1])
                self._set_fixed(self[i].pressure, f[2])
            elif sv == StateVars.PH and ab == AmountBasis.MASS:
                self._set_fixed(self[i].flow_mass, f[0])
                self._set_fixed(self[i].enth_mass, f[1])
                self._set_fixed(self[i].pressure, f[2])
            elif sv == StateVars.PS and ab == AmountBasis.MOLE:
                self._set_fixed(self[i].flow_mol, f[0])
                self._set_fixed(self[i].entr_mol, f[1])
                self._set_fixed(self[i].pressure, f[2])
            elif sv == StateVars.PS and ab == AmountBasis.MASS:
                self._set_fixed(self[i].flow_mass, f[0])
                self._set_fixed(self[i].entr_mass, f[1])
                self._set_fixed(self[i].pressure, f[2])
            elif sv == StateVars.PU and ab == AmountBasis.MOLE:
                self._set_fixed(self[i].flow_mol, f[0])
                self._set_fixed(self[i].energy_internal_mol, f[1])
                self._set_fixed(self[i].pressure, f[2])
            elif sv == StateVars.PU and ab == AmountBasis.MASS:
                self._set_fixed(self[i].flow_mass, f[0])
                self._set_fixed(self[i].energy_internal_mass, f[1])
                self._set_fixed(self[i].pressure, f[2])
            elif sv == StateVars.TPX and ab == AmountBasis.MOLE:
                self._set_fixed(self[i].flow_mol, f[0])
                self._set_fixed(self[i].temperature, f[1])
                self._set_fixed(self[i].pressure, f[2])
                if pp in (PhaseType.MIX, PhaseType.LG):
                    self._set_fixed(self[i].vapor_frac, f[3])
            elif sv == StateVars.TPX and ab == AmountBasis.MASS:
                self._set_fixed(self[i].flow_mass, f[0])
                self._set_fixed(self[i].temperature, f[1])
                self._set_fixed(self[i].pressure, f[2])
                if pp in (PhaseType.MIX, PhaseType.LG):
                    self._set_fixed(self[i].vapor_frac, f[3])


@declare_process_block_class("HelmholtzStateBlock", block_class=_StateBlock)
class HelmholtzStateBlockData(StateBlockData):
    """
    This is a base class for Helmholtz equations of state using IDAES standard
    Helmholtz EOS external functions written in C++.
    """

    def _state_vars(self):
        """Create the state variables"""
        params = self.config.parameters
        cmp = params.pure_component
        phase_set = params.config.phase_presentation

        # Check for inconsistent phase presentation and phase equilibrium
        if phase_set == PhaseType.MIX and self.config.has_phase_equilibrium:
            _log.warning(
                "Helmholtz EoS packages using Mixed phase representation ignore the "
                "'has_phase_equilibrium' configuration argument. However, setting "
                "this to True can result in errors when constructing material balances "
                "due to only having a single phase (thus phase transfer terms cannot "
                "be constructed)."
            )

        # Add flow variable
        if self.amount_basis == AmountBasis.MOLE:
            self.flow_mol = pyo.Var(
                initialize=1, doc="Total mole flow", units=pyunits.mol / pyunits.s
            )
            self.flow_mass = pyo.Expression(
                expr=self.mw * self.flow_mol, doc="Total mass flow"
            )
        else:
            self.flow_mass = pyo.Var(
                initialize=1, doc="Total mass flow", units=pyunits.kg / pyunits.s
            )
            self.flow_mol = pyo.Expression(
                expr=self.flow_mass / self.mw, doc="Total mole flow"
            )
        # All supported state variable sets include pressure
        self.pressure = pyo.Var(
            domain=pyo.PositiveReals,
            initialize=params.default_pressure_value,
            doc="Pressure",
            bounds=params.default_pressure_bounds,
            units=pyunits.Pa,
        )
        self.p_kPa = pyo.Expression(expr=self.pressure * params.uc["Pa to kPa"])
        # If it's single phase use provide fixed expressions for vapor frac
        if phase_set == PhaseType.L:
            self.vapor_frac = pyo.Expression(
                expr=0.0, doc="Vapor mole fraction (mol vapor/mol total)"
            )
        elif phase_set == PhaseType.G:
            self.vapor_frac = pyo.Expression(
                expr=1.0,
                doc="Vapor mole fraction (mol vapor/mol total)",
            )
        # Add other state var
        if self.state_vars == StateVars.PH and self.amount_basis == AmountBasis.MOLE:
            self.enth_mol = pyo.Var(
                initialize=params.default_enthalpy_mol_value,
                doc="Total molar enthalpy",
                bounds=params.default_enthalpy_mol_bounds,
                units=pyunits.J / pyunits.mol,
            )
            self.h_kJ_per_kg = pyo.Expression(
                expr=self.enth_mol * params.uc["J/mol to kJ/kg"]
            )
            self.temperature = pyo.Expression(
                expr=self.t_hp_func(cmp, self.h_kJ_per_kg, self.p_kPa, _data_dir),
                doc="Temperature",
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = pyo.Expression(
                    expr=self.vf_hp_func(cmp, self.h_kJ_per_kg, self.p_kPa, _data_dir),
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "enth_mol": self.enth_mol,
                "pressure": self.pressure,
            }
            self.extensive_set = ComponentSet((self.flow_mol,))
            self.intensive_set = ComponentSet((self.enth_mol, self.pressure))
        elif self.state_vars == StateVars.PH and self.amount_basis == AmountBasis.MASS:
            self.enth_mass = pyo.Var(
                initialize=params.default_enthalpy_mass_value,
                doc="Total enthalpy per mass",
                bounds=params.default_enthalpy_mass_bounds,
                units=pyunits.J / pyunits.kg,
            )
            self.h_kJ_per_kg = pyo.Expression(
                expr=self.enth_mass * params.uc["J/kg to kJ/kg"]
            )
            self.temperature = pyo.Expression(
                expr=self.t_hp_func(cmp, self.h_kJ_per_kg, self.p_kPa, _data_dir),
                doc="Temperature",
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = pyo.Expression(
                    expr=self.vf_hp_func(cmp, self.h_kJ_per_kg, self.p_kPa, _data_dir),
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )
            self._state_vars_dict = {
                "flow_mass": self.flow_mass,
                "enth_mass": self.enth_mass,
                "pressure": self.pressure,
            }
            self.extensive_set = ComponentSet((self.flow_mass,))
            self.intensive_set = ComponentSet((self.enth_mass, self.pressure))
        elif self.state_vars == StateVars.PS and self.amount_basis == AmountBasis.MOLE:
            self.entr_mol = pyo.Var(
                initialize=params.default_entropy_mol_value,
                doc="Total molar entropy",
                bounds=params.default_entropy_mol_bounds,
                units=pyunits.J / pyunits.mol / pyunits.K,
            )
            self.s_kJ_per_kgK = pyo.Expression(
                expr=self.entr_mol * params.uc["J/mol/K to kJ/kg/K"]
            )
            self.temperature = pyo.Expression(
                expr=self.t_sp_func(cmp, self.s_kJ_per_kgK, self.p_kPa, _data_dir),
                doc="Temperature",
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = pyo.Expression(
                    expr=self.vf_sp_func(cmp, self.s_kJ_per_kgK, self.p_kPa, _data_dir),
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "entr_mol": self.entr_mol,
                "pressure": self.pressure,
            }
            self.extensive_set = ComponentSet((self.flow_mol,))
            self.intensive_set = ComponentSet((self.entr_mol, self.pressure))
        elif self.state_vars == StateVars.PS and self.amount_basis == AmountBasis.MASS:
            self.entr_mass = pyo.Var(
                initialize=params.default_entropy_mass_value,
                doc="Total entropy per mass",
                bounds=params.default_entropy_mass_bounds,
                units=pyunits.J / pyunits.kg / pyunits.K,
            )
            self.s_kJ_per_kgK = pyo.Expression(
                expr=self.entr_mass * params.uc["J/kg/K to kJ/kg/K"]
            )
            self.temperature = pyo.Expression(
                expr=self.t_sp_func(cmp, self.s_kJ_per_kgK, self.p_kPa, _data_dir),
                doc="Temperature",
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = pyo.Expression(
                    expr=self.vf_sp_func(cmp, self.s_kJ_per_kgK, self.p_kPa, _data_dir),
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )
            self._state_vars_dict = {
                "flow_mass": self.flow_mass,
                "entr_mass": self.entr_mass,
                "pressure": self.pressure,
            }
            self.extensive_set = ComponentSet((self.flow_mass,))
            self.intensive_set = ComponentSet((self.entr_mass, self.pressure))
        elif self.state_vars == StateVars.PU and self.amount_basis == AmountBasis.MOLE:
            self.energy_internal_mol = pyo.Var(
                initialize=params.default_energy_internal_mol_value,
                doc="Total molar internal energy",
                bounds=params.default_energy_internal_mol_bounds,
                units=pyunits.J / pyunits.mol,
            )
            self.u_kJ_per_kg = pyo.Expression(
                expr=self.energy_internal_mol * params.uc["J/mol to kJ/kg"]
            )
            self.temperature = pyo.Expression(
                expr=self.t_up_func(cmp, self.u_kJ_per_kg, self.p_kPa, _data_dir),
                doc="Temperature",
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = pyo.Expression(
                    expr=self.vf_up_func(cmp, self.u_kJ_per_kg, self.p_kPa, _data_dir),
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )
            self._state_vars_dict = {
                "flow_mol": self.flow_mol,
                "energy_internal_mol": self.energy_internal_mol,
                "pressure": self.pressure,
            }
            self.extensive_set = ComponentSet((self.flow_mol,))
            self.intensive_set = ComponentSet((self.energy_internal_mol, self.pressure))
        elif self.state_vars == StateVars.PU and self.amount_basis == AmountBasis.MASS:
            self.energy_internal_mass = pyo.Var(
                initialize=params.default_energy_internal_mass_value,
                doc="Total internal energy per mass",
                bounds=params.default_energy_internal_mass_bounds,
                units=pyunits.J / pyunits.kg,
            )
            self.u_kJ_per_kg = pyo.Expression(
                expr=self.energy_internal_mass * params.uc["J/kg to kJ/kg"]
            )
            self.temperature = pyo.Expression(
                expr=self.t_up_func(cmp, self.u_kJ_per_kg, self.p_kPa, _data_dir),
                doc="Temperature",
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = pyo.Expression(
                    expr=self.vf_up_func(cmp, self.u_kJ_per_kg, self.p_kPa, _data_dir),
                    doc="Vapor mole fraction (mol vapor/mol total)",
                )
            self._state_vars_dict = {
                "flow_mass": self.flow_mass,
                "energy_internal_mass": self.energy_internal_mass,
                "pressure": self.pressure,
            }
            self.extensive_set = ComponentSet((self.flow_mass,))
            self.intensive_set = ComponentSet(
                (self.energy_internal_mass, self.pressure)
            )
        if self.state_vars == StateVars.TPX:
            self.temperature = pyo.Var(
                domain=pyo.PositiveReals,
                initialize=params.default_temperature_value,
                doc="Temperature",
                bounds=params.default_temperature_bounds,
                units=pyunits.K,
            )
            if phase_set == PhaseType.MIX or phase_set == PhaseType.LG:
                self.vapor_frac = pyo.Var(
                    initialize=0.0,
                    units=pyunits.dimensionless,
                    doc="Vapor fraction"
                    # No bounds here, since it is often (usually) on it's bound
                    # and that's not the best for IPOPT
                )
                self.intensive_set = ComponentSet(
                    (self.temperature, self.pressure, self.vapor_frac)
                )
                self._state_vars_dict = {
                    "temperature": self.temperature,
                    "pressure": self.pressure,
                    "vapor_frac": self.vapor_frac,
                }
            else:
                self.intensive_set = ComponentSet((self.temperature, self.pressure))
                self._state_vars_dict = {
                    "temperature": self.temperature,
                    "pressure": self.pressure,
                }
            if self.amount_basis == AmountBasis.MOLE:
                self.extensive_set = ComponentSet((self.flow_mol,))
                self._state_vars_dict["flow_mol"] = self.flow_mol
            else:
                self.extensive_set = ComponentSet((self.flow_mass,))
                self._state_vars_dict["flow_mass"] = self.flow_mass

    def _tpx_phase_eq(self):
        params = self.config.parameters
        eps_pu = params.smoothing_pressure_under
        eps_po = params.smoothing_pressure_over
        priv_plist = params.private_phase_list
        # Convert pressures to kPa for external functions and nicer scaling in
        # the complementarity-type constraints.
        vf = self.vapor_frac
        # Terms for determining if you are above, below, or at the Psat
        self.pressure_under_sat = pyo.Expression(
            expr=smooth_max(0, self.pressure_sat - self.pressure, eps_pu),
            doc="pressure above Psat, 0 if liquid exists",
        )
        self.pressure_over_sat = pyo.Expression(
            expr=smooth_max(0, self.pressure - self.pressure_sat, eps_po),
            doc="pressure below Psat, 0 if vapor exists",
        )

        if not self.config.defined_state:
            self.eq_complementarity = pyo.Constraint(
                expr=0
                == (
                    vf * self.pressure_over_sat / 1000.0
                    - (1 - vf) * self.pressure_under_sat / 1000.0
                )
            )

        # Calculate liquid and vapor density.  If the phase doesn't exist,
        # density will be calculated at the saturation or critical pressure
        def rule_pressure_phase(b, p):
            if p == "Liq":
                return self.pressure + self.pressure_under_sat
            else:
                return self.pressure - self.pressure_over_sat

        self.pressure_phase = pyo.Expression(priv_plist, rule=rule_pressure_phase)

        # Constraint that can be activated to enforce that you are in the two-phase region
        self.eq_sat = pyo.Constraint(
            expr=self.pressure * 1e-6 == self.pressure_sat * 1e-6
        )
        self.eq_sat.deactivate()

    def build(self, *args):
        """
        Callable method for Block construction
        """
        # Call base class build
        super().build(*args)
        # Short path to the parameter block
        params = self.config.parameters
        # Check if the library is available, and add external functions.
        add_helmholtz_external_functions(self)
        cmp = params.pure_component
        # Which state vars to use
        self.state_vars = params.state_vars
        # Mass or Mole basis
        self.amount_basis = params.config.amount_basis
        # Private phase list
        phlist = params.private_phase_list
        # Public phase list
        pub_phlist = params.phase_list
        component_list = params.component_list
        self.phase_equilibrium_list = params.phase_equilibrium_list

        # Add component mole fraction for standardization
        def mole_frac_comp_rule(b, i):
            return 1.0

        self.mole_frac_comp = pyo.Expression(component_list, rule=mole_frac_comp_rule)
        # Expressions that link to some parameters in the param block, which
        # are commonly needed, this lets you get the parameters with scale
        # factors directly from the state block
        self.temperature_crit = pyo.Expression(
            expr=params.temperature_crit, doc="critical temperature"
        )
        self.temperature_star = pyo.Expression(
            expr=params.temperature_star, doc="temperature for tau calculation"
        )
        self.pressure_crit = pyo.Expression(
            expr=params.pressure_crit, doc="critical pressure"
        )
        self.dens_mass_crit = pyo.Expression(
            expr=params.dens_mass_crit, doc="critical mass density"
        )
        self.dens_mass_star = pyo.Expression(
            expr=params.dens_mass_star, doc="mass density for delta calculation"
        )
        self.dens_mol_crit = pyo.Expression(
            expr=params.dens_mass_crit / params.mw, doc="critical mole density"
        )
        self.dens_mol_star = pyo.Expression(
            expr=params.dens_mass_star / params.mw,
            doc="mole density for delta calculation",
        )
        self.mw = pyo.Expression(expr=params.mw, doc="molecular weight")
        # create the appropriate state variables and expressions for anything
        # that could be a state variable (H, S, U, P, T, X) on a mass and mole
        # basis if using TPx as state variables it also adds the complementarity
        # constraints. Beyond this everything else is common.
        self._state_vars()

        # P is pressure in kPa for external function calls
        T = self.temperature
        vf = self.vapor_frac
        # Saturation temperature expression
        self.temperature_sat = pyo.Expression(
            expr=self.t_sat_func(cmp, self.p_kPa, _data_dir),
            doc="Saturation temperature",
        )
        # Saturation pressure
        self.pressure_sat = pyo.Expression(
            expr=self.p_sat_t_func(cmp, T, _data_dir) * params.uc["kPa to Pa"],
            doc="Saturation pressure",
        )
        # Add the complementarity constraint for phase equilibrium with TPx.
        if self.state_vars == StateVars.TPX and len(phlist) > 1:
            self._tpx_phase_eq()

        # to make writing the remaining expressions simpler create a state var dict
        if self.state_vars == StateVars.PH:
            sv_dict = {"h": self.h_kJ_per_kg, "p": self.p_kPa}
        elif self.state_vars == StateVars.PS:
            sv_dict = {"s": self.s_kJ_per_kgK, "p": self.p_kPa}
        elif self.state_vars == StateVars.PU:
            sv_dict = {"u": self.u_kJ_per_kg, "p": self.p_kPa}
        elif self.state_vars == StateVars.TPX:
            sv_dict = {
                "T": self.temperature,
                "p": self.p_kPa,
                "x": self.vapor_frac,
            }
        sv_dict_liq = copy.copy(sv_dict)
        sv_dict_vap = copy.copy(sv_dict)
        if self.state_vars == StateVars.TPX and len(phlist) > 1:
            self.p_kPa_liq = pyo.Expression(
                expr=self.pressure_phase["Liq"] * params.uc["Pa to kPa"]
            )
            self.p_kPa_vap = pyo.Expression(
                expr=self.pressure_phase["Vap"] * params.uc["Pa to kPa"]
            )
            sv_dict_liq["p"] = self.p_kPa_liq
            sv_dict_vap["p"] = self.p_kPa_vap

        self.expression_writer = HelmholtzThermoExpressions(self, params)

        # Saturated Enthalpy molar
        def rule_enth_mol_sat_phase(b, p):
            if p == "Liq":
                return self.expression_writer.h_liq_sat(
                    p=b.p_kPa, result_basis=AmountBasis.MOLE, convert_args=False
                )
            else:
                return self.expression_writer.h_vap_sat(
                    p=b.p_kPa, result_basis=AmountBasis.MOLE, convert_args=False
                )

        self.enth_mol_sat_phase = pyo.Expression(
            phlist,
            rule=rule_enth_mol_sat_phase,
            doc="Saturated enthalpy of the phases at pressure",
        )

        # Saturated Enthalpy mass
        def rule_enth_mass_sat_phase(b, p):
            if p == "Liq":
                return self.expression_writer.h_liq_sat(
                    p=b.p_kPa, result_basis=AmountBasis.MASS, convert_args=False
                )
            else:
                return self.expression_writer.h_vap_sat(
                    p=b.p_kPa, result_basis=AmountBasis.MASS, convert_args=False
                )

        self.enth_mass_sat_phase = pyo.Expression(
            phlist,
            rule=rule_enth_mass_sat_phase,
            doc="Saturated enthalpy of the phases at pressure",
        )

        # Saturated Entropy molar
        def rule_entr_mol_sat_phase(b, p):
            if p == "Liq":
                return self.expression_writer.s_liq_sat(
                    p=b.p_kPa, result_basis=AmountBasis.MOLE, convert_args=False
                )
            else:
                return self.expression_writer.s_vap_sat(
                    p=b.p_kPa, result_basis=AmountBasis.MOLE, convert_args=False
                )

        self.entr_mol_sat_phase = pyo.Expression(
            phlist,
            rule=rule_entr_mol_sat_phase,
            doc="Saturated entropy of the phases at pressure",
        )

        # Saturated Entropy mass
        def rule_entr_mass_sat_phase(b, p):
            if p == "Liq":
                return self.expression_writer.s_liq_sat(
                    p=b.p_kPa, result_basis=AmountBasis.MASS, convert_args=False
                )
            else:
                return self.expression_writer.s_vap_sat(
                    p=b.p_kPa, result_basis=AmountBasis.MASS, convert_args=False
                )

        self.entr_mass_sat_phase = pyo.Expression(
            phlist,
            rule=rule_entr_mass_sat_phase,
            doc="Saturated entropy of the phases at pressure",
        )

        # Saturated internal energy molar
        def rule_energy_internal_mol_sat_phase(b, p):
            if p == "Liq":
                return self.expression_writer.u_liq_sat(
                    p=b.p_kPa, result_basis=AmountBasis.MOLE, convert_args=False
                )
            else:
                return self.expression_writer.u_vap_sat(
                    p=b.p_kPa, result_basis=AmountBasis.MOLE, convert_args=False
                )

        self.energy_internal_mol_sat_phase = pyo.Expression(
            phlist,
            rule=rule_energy_internal_mol_sat_phase,
            doc="Saturated internal energy of the phases at pressure",
        )

        # Saturated internal energy mass
        def rule_energy_internal_mass_sat_phase(b, p):
            if p == "Liq":
                return self.expression_writer.u_liq_sat(
                    p=b.p_kPa, result_basis=AmountBasis.MASS, convert_args=False
                )
            else:
                return self.expression_writer.u_vap_sat(
                    p=b.p_kPa, result_basis=AmountBasis.MASS, convert_args=False
                )

        self.energy_internal_mass_sat_phase = pyo.Expression(
            phlist,
            rule=rule_energy_internal_mass_sat_phase,
            doc="Saturated internal energy of the phases at pressure",
        )

        # Saturated molar volume
        def rule_volume_mol_sat_phase(b, p):
            if p == "Liq":
                return self.expression_writer.v_liq_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MOLE, convert_args=False
                )
            else:
                return self.expression_writer.v_vap_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MOLE, convert_args=False
                )

        self.volume_mol_sat_phase = pyo.Expression(
            phlist,
            rule=rule_volume_mol_sat_phase,
            doc="Saturated molar volume of the phases at pressure",
        )

        # Saturated specific volume
        def rule_volume_mass_sat_phase(b, p):
            if p == "Liq":
                return self.expression_writer.v_liq_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MASS, convert_args=False
                )
            else:
                return self.expression_writer.v_vap_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MASS, convert_args=False
                )

        self.volume_mass_sat_phase = pyo.Expression(
            phlist,
            rule=rule_volume_mass_sat_phase,
            doc="Saturated specific volume of the phases at pressure",
        )

        # delta h vap at P
        self.dh_vap_mol = pyo.Expression(
            expr=(
                self.expression_writer.h_vap_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MOLE, convert_args=False
                )
                - self.expression_writer.h_liq_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MOLE, convert_args=False
                )
            ),
            doc="Enthalpy of vaporization at pressure and saturation temperature",
        )

        # delta s vap at P
        self.ds_vap_mol = pyo.Expression(
            expr=(
                self.expression_writer.s_vap_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MOLE, convert_args=False
                )
                - self.expression_writer.s_liq_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MOLE, convert_args=False
                )
            ),
            doc="Entropy of vaporization at pressure and saturation temperature",
        )

        # delta u vap at P
        self.du_vap_mol = pyo.Expression(
            expr=(
                self.expression_writer.u_vap_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MOLE, convert_args=False
                )
                - self.expression_writer.u_liq_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MOLE, convert_args=False
                )
            ),
            doc="Internal energy of vaporization at pressure and saturation temperature",
        )

        # delta h vap at P
        self.dh_vap_mass = pyo.Expression(
            expr=(
                self.expression_writer.h_vap_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MASS, convert_args=False
                )
                - self.expression_writer.h_liq_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MASS, convert_args=False
                )
            ),
            doc="Enthalpy of vaporization at pressure and saturation temperature",
        )

        # delta s vap at P
        self.ds_vap_mass = pyo.Expression(
            expr=(
                self.expression_writer.s_vap_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MASS, convert_args=False
                )
                - self.expression_writer.s_liq_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MASS, convert_args=False
                )
            ),
            doc="Entropy of vaporization at pressure and saturation temperature",
        )

        # delta u vap at P
        self.du_vap_mass = pyo.Expression(
            expr=(
                self.expression_writer.u_vap_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MASS, convert_args=False
                )
                - self.expression_writer.u_liq_sat(
                    p=self.p_kPa, result_basis=AmountBasis.MASS, convert_args=False
                )
            ),
            doc="Internal energy of vaporization at pressure and saturation temperature",
        )

        # Phase fraction
        def rule_phase_frac(b, p):
            if p == "Vap":
                return vf
            elif p == "Liq":
                return 1.0 - vf

        self.phase_frac = pyo.Expression(
            phlist, rule=rule_phase_frac, doc="Phase fraction"
        )

        def rule_energy_internal_mol_phase(b, p):
            if p == "Liq":
                return self.expression_writer.u_liq(
                    **sv_dict_liq, result_basis=AmountBasis.MOLE, convert_args=False
                )
            else:
                return self.expression_writer.u_vap(
                    **sv_dict_vap, result_basis=AmountBasis.MOLE, convert_args=False
                )

        self.energy_internal_mol_phase = pyo.Expression(
            phlist,
            rule=rule_energy_internal_mol_phase,
            doc="Phase internal energy",
        )

        def rule_energy_internal_mass_phase(b, p):
            if p == "Liq":
                return self.expression_writer.u_liq(
                    **sv_dict_liq, result_basis=AmountBasis.MASS, convert_args=False
                )
            else:
                return self.expression_writer.u_vap(
                    **sv_dict_vap, result_basis=AmountBasis.MASS, convert_args=False
                )

        self.energy_internal_mass_phase = pyo.Expression(
            phlist,
            rule=rule_energy_internal_mass_phase,
            doc="Phase internal energy",
        )

        # Phase Enthalpy
        def rule_enth_mol_phase(b, p):
            if p == "Liq":
                return self.expression_writer.h_liq(
                    **sv_dict_liq, result_basis=AmountBasis.MOLE, convert_args=False
                )
            else:
                return self.expression_writer.h_vap(
                    **sv_dict_vap, result_basis=AmountBasis.MOLE, convert_args=False
                )

        self.enth_mol_phase = pyo.Expression(
            phlist,
            rule=rule_enth_mol_phase,
            doc="Phase enthalpy",
        )

        def rule_enth_mass_phase(b, p):
            if p == "Liq":
                return self.expression_writer.h_liq(
                    **sv_dict_liq, result_basis=AmountBasis.MASS, convert_args=False
                )
            else:
                return self.expression_writer.h_vap(
                    **sv_dict_vap, result_basis=AmountBasis.MASS, convert_args=False
                )

        self.enth_mass_phase = pyo.Expression(
            phlist,
            rule=rule_enth_mass_phase,
            doc="Phase enthalpy",
        )

        # Phase Entropy
        def rule_entr_mol_phase(b, p):
            if p == "Liq":
                return self.expression_writer.s_liq(
                    **sv_dict_liq, result_basis=AmountBasis.MOLE, convert_args=False
                )
            else:
                return self.expression_writer.s_vap(
                    **sv_dict_vap, result_basis=AmountBasis.MOLE, convert_args=False
                )

        self.entr_mol_phase = pyo.Expression(
            phlist,
            rule=rule_entr_mol_phase,
            doc="Phase entropy",
        )

        def rule_entr_mass_phase(b, p):
            if p == "Liq":
                return self.expression_writer.s_liq(
                    **sv_dict_liq, result_basis=AmountBasis.MASS, convert_args=False
                )
            else:
                return self.expression_writer.s_vap(
                    **sv_dict_vap, result_basis=AmountBasis.MASS, convert_args=False
                )

        self.entr_mass_phase = pyo.Expression(
            phlist,
            rule=rule_entr_mass_phase,
            doc="Phase entropy",
        )

        # Phase constant pressure heat capacity, cp
        def rule_cp_mol_phase(b, p):
            if p == "Liq":
                return self.expression_writer.cp_liq(
                    **sv_dict_liq, result_basis=AmountBasis.MOLE, convert_args=False
                )
            else:
                return self.expression_writer.cp_vap(
                    **sv_dict_vap, result_basis=AmountBasis.MOLE, convert_args=False
                )

        self.cp_mol_phase = pyo.Expression(
            phlist,
            rule=rule_cp_mol_phase,
            doc="Phase isobaric heat capacity",
        )

        def rule_cp_mass_phase(b, p):
            if p == "Liq":
                return self.expression_writer.cp_liq(
                    **sv_dict_liq, result_basis=AmountBasis.MASS, convert_args=False
                )
            else:
                return self.expression_writer.cp_vap(
                    **sv_dict_vap, result_basis=AmountBasis.MASS, convert_args=False
                )

        self.cp_mass_phase = pyo.Expression(
            phlist,
            rule=rule_cp_mass_phase,
            doc="Phase isobaric heat capacity",
        )

        # Phase constant volume heat capacity, cv
        def rule_cv_mol_phase(b, p):
            if p == "Liq":
                return self.expression_writer.cv_liq(
                    **sv_dict_liq, result_basis=AmountBasis.MOLE, convert_args=False
                )
            else:
                return self.expression_writer.cv_vap(
                    **sv_dict_vap, result_basis=AmountBasis.MOLE, convert_args=False
                )

        self.cv_mol_phase = pyo.Expression(
            phlist,
            rule=rule_cv_mol_phase,
            doc="Phase isochoric heat capacity",
        )

        def rule_cv_mass_phase(b, p):
            if p == "Liq":
                return self.expression_writer.cv_liq(
                    **sv_dict_liq, result_basis=AmountBasis.MASS, convert_args=False
                )
            else:
                return self.expression_writer.cv_vap(
                    **sv_dict_vap, result_basis=AmountBasis.MASS, convert_args=False
                )

        self.cv_mass_phase = pyo.Expression(
            phlist,
            rule=rule_cv_mass_phase,
            doc="Phase isochoric heat capacity",
        )

        # Phase speed of sound
        def rule_speed_sound_phase(b, p):
            if p == "Liq":
                return self.expression_writer.w_liq(**sv_dict_liq, convert_args=False)
            else:
                return self.expression_writer.w_vap(**sv_dict_liq, convert_args=False)

        self.speed_sound_phase = pyo.Expression(
            phlist,
            rule=rule_speed_sound_phase,
            doc="Phase speed of sound or saturated if phase doesn't exist",
        )

        # molar volume
        def rule_vol_mol_phase(b, p):
            if p == "Liq":
                return self.expression_writer.v_liq(
                    **sv_dict_liq, result_basis=AmountBasis.MOLE, convert_args=False
                )
            else:
                return self.expression_writer.v_vap(
                    **sv_dict_vap, result_basis=AmountBasis.MOLE, convert_args=False
                )

        self.vol_mol_phase = pyo.Expression(
            phlist,
            rule=rule_vol_mol_phase,
            doc="Molar volume of phase",
        )

        # specific volume
        def rule_vol_mass_phase(b, p):
            if p == "Liq":
                return self.expression_writer.v_liq(
                    **sv_dict_liq, result_basis=AmountBasis.MASS, convert_args=False
                )
            else:
                return self.expression_writer.v_vap(
                    **sv_dict_vap, result_basis=AmountBasis.MASS, convert_args=False
                )

        self.vol_mass_phase = pyo.Expression(
            phlist,
            rule=rule_vol_mass_phase,
            doc="Specific volume of phase",
        )

        # molar density
        def rule_dens_mol_phase(b, p):
            if p == "Liq":
                return 1.0 / self.expression_writer.v_liq(
                    **sv_dict_liq, result_basis=AmountBasis.MOLE, convert_args=False
                )
            else:
                return 1.0 / self.expression_writer.v_vap(
                    **sv_dict_vap, result_basis=AmountBasis.MOLE, convert_args=False
                )

        self.dens_mol_phase = pyo.Expression(
            phlist,
            rule=rule_dens_mol_phase,
            doc="Mole density of phase",
        )

        # specific density
        def rule_dens_mass_phase(b, p):
            if p == "Liq":
                return 1.0 / self.expression_writer.v_liq(
                    **sv_dict_liq, result_basis=AmountBasis.MASS, convert_args=False
                )
            else:
                return 1.0 / self.expression_writer.v_vap(
                    **sv_dict_vap, result_basis=AmountBasis.MASS, convert_args=False
                )

        self.dens_mass_phase = pyo.Expression(
            phlist,
            rule=rule_dens_mass_phase,
            doc="Mass density of phase",
        )

        # Component flow (for units that need it)
        def component_flow_mol(b, i):
            return self.flow_mol

        self.flow_mol_comp = pyo.Expression(
            component_list,
            rule=component_flow_mol,
            doc="Total mole flow (both phases) of component",
        )

        def component_flow_mass(b, i):
            return self.flow_mass

        self.flow_mass_comp = pyo.Expression(
            component_list,
            rule=component_flow_mass,
            doc="Total mass flow (both phases) of component",
        )

        #
        # Total (mixed phase) properties
        #

        # Enthalpy and pressure state variables two-phase properties
        if self.state_vars == StateVars.PH:
            if self.amount_basis == AmountBasis.MOLE:
                self.enth_mass = pyo.Expression(
                    expr=self.enth_mol * params.uc["J/mol to J/kg"]
                )
            else:
                self.enth_mol = pyo.Expression(
                    expr=self.enth_mass * params.uc["J/kg to J/mol"]
                )
            self.entr_mass = pyo.Expression(
                expr=self.s_hp_func(cmp, self.h_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/kg/K"]
            )
            self.entr_mol = pyo.Expression(
                expr=self.s_hp_func(cmp, self.h_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/mol/K"]
            )
            self.energy_internal_mass = pyo.Expression(
                expr=self.u_hp_func(cmp, self.h_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg to J/kg"]
            )
            self.energy_internal_mol = pyo.Expression(
                expr=self.u_hp_func(cmp, self.h_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg to J/mol"]
            )
            self.cp_mass = pyo.Expression(
                expr=self.cp_hp_func(cmp, self.h_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/kg/K"]
            )
            self.cp_mol = pyo.Expression(
                expr=self.cp_hp_func(cmp, self.h_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/mol/K"]
            )
            self.cv_mass = pyo.Expression(
                expr=self.cv_hp_func(cmp, self.h_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/kg/K"]
            )
            self.cv_mol = pyo.Expression(
                expr=self.cv_hp_func(cmp, self.h_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/mol/K"]
            )
            self.dens_mass = pyo.Expression(
                expr=1.0 / self.v_hp_func(cmp, self.h_kJ_per_kg, self.p_kPa, _data_dir)
            )
            self.dens_mol = pyo.Expression(
                expr=params.uc["kg/m3 to mol/m3"]
                / self.v_hp_func(cmp, self.h_kJ_per_kg, self.p_kPa, _data_dir)
            )
        # Entropy is a state variable
        elif self.state_vars == StateVars.PS:
            if self.amount_basis == AmountBasis.MOLE:
                self.entr_mass = pyo.Expression(
                    expr=self.entr_mol * params.uc["J/mol/K to J/kg/K"]
                )
            else:
                self.enth_mol = pyo.Expression(
                    expr=self.entr_mass * params.uc["J/kg/K to J/mol/K"]
                )
            self.enth_mass = pyo.Expression(
                expr=self.h_sp_func(cmp, self.s_kJ_per_kgK, self.p_kPa, _data_dir)
                * params.uc["kJ/kg to J/kg"]
            )
            self.enth_mol = pyo.Expression(
                expr=self.h_sp_func(cmp, self.s_kJ_per_kgK, self.p_kPa, _data_dir)
                * params.uc["kJ/kg to J/mol"]
            )
            self.energy_internal_mass = pyo.Expression(
                expr=self.u_sp_func(cmp, self.s_kJ_per_kgK, self.p_kPa, _data_dir)
                * params.uc["kJ/kg to J/kg"]
            )
            self.energy_internal_mol = pyo.Expression(
                expr=self.u_sp_func(cmp, self.s_kJ_per_kgK, self.p_kPa, _data_dir)
                * params.uc["kJ/kg to J/mol"]
            )
            self.cp_mass = pyo.Expression(
                expr=self.cp_sp_func(cmp, self.s_kJ_per_kgK, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/kg/K"]
            )
            self.cp_mol = pyo.Expression(
                expr=self.cp_sp_func(cmp, self.s_kJ_per_kgK, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/mol/K"]
            )
            self.cv_mass = pyo.Expression(
                expr=self.cv_sp_func(cmp, self.s_kJ_per_kgK, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/kg/K"]
            )
            self.cv_mol = pyo.Expression(
                expr=self.cv_sp_func(cmp, self.s_kJ_per_kgK, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/mol/K"]
            )
            self.dens_mass = pyo.Expression(
                expr=1.0 / self.v_sp_func(cmp, self.s_kJ_per_kgK, self.p_kPa, _data_dir)
            )
            self.dens_mol = pyo.Expression(
                expr=params.uc["kg/m3 to mol/m3"]
                / self.v_sp_func(cmp, self.s_kJ_per_kgK, self.p_kPa, _data_dir)
            )
        # Internal energy is a state variable
        elif self.state_vars == StateVars.PU:
            if self.amount_basis == AmountBasis.MOLE:
                self.energy_internal_mass = pyo.Expression(
                    expr=self.energy_internal_mol * params.uc["J/mol to J/kg"]
                )
            else:
                self.energy_internal_mol = pyo.Expression(
                    expr=self.energy_internal_mass * params.uc["J/kg to J/mol"]
                )
            self.enth_mass = pyo.Expression(
                expr=self.h_up_func(cmp, self.u_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg to J/kg"]
            )
            self.enth_mol = pyo.Expression(
                expr=self.h_up_func(cmp, self.u_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg to J/mol"]
            )
            self.entr_mass = pyo.Expression(
                expr=self.s_up_func(cmp, self.u_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/kg/K"]
            )
            self.entr_mol = pyo.Expression(
                expr=self.s_up_func(cmp, self.u_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/mol/K"]
            )
            self.cp_mass = pyo.Expression(
                expr=self.cp_up_func(cmp, self.u_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/kg/K"]
            )
            self.cp_mol = pyo.Expression(
                expr=self.cp_up_func(cmp, self.u_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/mol/K"]
            )
            self.cv_mass = pyo.Expression(
                expr=self.cv_up_func(cmp, self.u_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/kg/K"]
            )
            self.cv_mol = pyo.Expression(
                expr=self.cv_up_func(cmp, self.u_kJ_per_kg, self.p_kPa, _data_dir)
                * params.uc["kJ/kg/K to J/mol/K"]
            )
            self.dens_mass = pyo.Expression(
                expr=1.0 / self.v_up_func(cmp, self.u_kJ_per_kg, self.p_kPa, _data_dir)
            )
            self.dens_mol = pyo.Expression(
                expr=params.uc["kg/m3 to mol/m3"]
                / self.v_up_func(cmp, self.u_kJ_per_kg, self.p_kPa, _data_dir)
            )
        else:  # T, P, x
            # enthalpy
            self.enth_mol = pyo.Expression(
                expr=sum(self.phase_frac[p] * self.enth_mol_phase[p] for p in phlist)
            )

            self.enth_mass = pyo.Expression(
                expr=sum(self.phase_frac[p] * self.enth_mass_phase[p] for p in phlist)
            )

            # Entropy
            self.entr_mol = pyo.Expression(
                expr=sum(self.phase_frac[p] * self.entr_mol_phase[p] for p in phlist)
            )
            self.entr_mass = pyo.Expression(
                expr=sum(self.phase_frac[p] * self.entr_mass_phase[p] for p in phlist)
            )
            # Internal Energy
            self.energy_internal_mol = pyo.Expression(
                expr=sum(
                    self.phase_frac[p] * self.energy_internal_mol_phase[p]
                    for p in phlist
                )
            )
            self.energy_internal_mass = pyo.Expression(
                expr=sum(
                    self.phase_frac[p] * self.energy_internal_mass_phase[p]
                    for p in phlist
                )
            )
            # cp
            self.cp_mol = pyo.Expression(
                expr=sum(self.phase_frac[p] * self.cp_mol_phase[p] for p in phlist)
            )
            self.cp_mass = pyo.Expression(
                expr=sum(self.phase_frac[p] * self.cp_mass_phase[p] for p in phlist)
            )
            # cv
            self.cv_mol = pyo.Expression(
                expr=sum(self.phase_frac[p] * self.cv_mol_phase[p] for p in phlist)
            )
            self.cv_mass = pyo.Expression(
                expr=sum(self.phase_frac[p] * self.cv_mass_phase[p] for p in phlist)
            )
            # mass density
            self.dens_mass = pyo.Expression(
                expr=1.0
                / sum(
                    self.phase_frac[p] * 1.0 / self.dens_mass_phase[p] for p in phlist
                )
            )
            # mole density
            self.dens_mol = pyo.Expression(
                expr=1.0
                / sum(self.phase_frac[p] * 1.0 / self.dens_mol_phase[p] for p in phlist)
            )

        # heat capacity ratio
        self.heat_capacity_ratio = pyo.Expression(expr=self.cp_mol / self.cv_mol)
        # Volumetric flow
        self.flow_vol = pyo.Expression(
            expr=self.flow_mol / self.dens_mol,
            doc="Total liquid + vapor volumetric flow",
        )

        def rule_mole_frac_phase_comp(b, p, j):
            if p == "Mix":
                return 1.0
            else:
                return self.phase_frac[p]

        self.mole_frac_phase_comp = pyo.Expression(
            pub_phlist, component_list, rule=rule_mole_frac_phase_comp
        )
        self.mass_frac_phase_comp = pyo.Expression(  # since pure comp same as mole
            pub_phlist, component_list, rule=rule_mole_frac_phase_comp
        )

        if viscosity_available(cmp):

            def rule_visc_d_phase(b, p):
                if p == "Liq":
                    return self.expression_writer.viscosity_liq(
                        **sv_dict_liq, convert_args=False
                    )
                else:
                    return self.expression_writer.viscosity_vap(
                        **sv_dict_vap, convert_args=False
                    )

            self.visc_d_phase = pyo.Expression(
                phlist,
                rule=rule_visc_d_phase,
                doc="(Dynamic) viscosity of phase",
            )

            def rule_visc_k_phase(b, p):
                if p == "Liq":
                    return (
                        self.expression_writer.viscosity_liq(
                            **sv_dict_liq, convert_args=False
                        )
                        / self.dens_mass_phase[p]
                    )
                else:
                    return (
                        self.expression_writer.viscosity_vap(
                            **sv_dict_vap, convert_args=False
                        )
                        / self.dens_mass_phase[p]
                    )

            self.visc_k_phase = pyo.Expression(
                phlist,
                rule=rule_visc_k_phase,
                doc="Kinematic viscosity of phase",
            )

        if thermal_conductivity_available(cmp):

            def rule_therm_cond_phase(b, p):
                if p == "Liq":
                    return self.expression_writer.thermal_conductivity_liq(
                        **sv_dict_liq, convert_args=False
                    )
                else:
                    return self.expression_writer.thermal_conductivity_vap(
                        **sv_dict_vap, convert_args=False
                    )

            self.therm_cond_phase = pyo.Expression(
                phlist,
                rule=rule_therm_cond_phase,
                doc="Thermal conductivity of phase",
            )

        if surface_tension_available(cmp):
            self.surf_tens = pyo.Expression(
                phlist,
                expr=self.expression_writer.surface_tension(
                    **sv_dict_liq, convert_args=False
                ),
                doc="Thermal conductivity of phase",
            )

        # Define some expressions for the balance terms returned by functions
        # This is just to allow assigning scale factors to the expressions
        # returned
        #
        # Material flow term expressions
        if self.amount_basis == AmountBasis.MOLE:

            def rule_material_flow_terms(b, p):
                if p == "Mix":
                    return self.flow_mol
                else:
                    return self.flow_mol * self.phase_frac[p]

        else:

            def rule_material_flow_terms(b, p):
                if p == "Mix":
                    return self.flow_mass
                else:
                    return self.flow_mass * self.phase_frac[p]

        self.material_flow_terms = pyo.Expression(
            pub_phlist, rule=rule_material_flow_terms
        )

        # Enthalpy flow term expressions
        if self.amount_basis == AmountBasis.MOLE:

            def rule_enthalpy_flow_terms(b, p):
                if p == "Mix":
                    return self.enth_mol * self.flow_mol
                else:
                    return self.enth_mol_phase[p] * self.phase_frac[p] * self.flow_mol

        else:

            def rule_enthalpy_flow_terms(b, p):
                if p == "Mix":
                    return self.enth_mass * self.flow_mass
                else:
                    return self.enth_mass_phase[p] * self.phase_frac[p] * self.flow_mass

        self.enthalpy_flow_terms = pyo.Expression(
            pub_phlist, rule=rule_enthalpy_flow_terms
        )

        # Energy density term expressions
        if self.amount_basis == AmountBasis.MOLE:

            def rule_energy_density_terms(b, p):
                if p == "Mix":
                    return self.dens_mol * self.energy_internal_mol
                else:
                    return self.dens_mol_phase[p] * self.energy_internal_mol_phase[p]

        else:

            def rule_energy_density_terms(b, p):
                if p == "Mix":
                    return self.dens_mass * self.energy_internal_mass
                else:
                    return self.dens_mass_phase[p] * self.energy_internal_mass_phase[p]

        self.energy_density_terms = pyo.Expression(
            pub_phlist, rule=rule_energy_density_terms
        )

    def get_material_flow_terms(self, p, j):
        """Get material flow terms for phase

        Args:
            p (str): phase

        Returns:
            Expression
        """
        return self.material_flow_terms[p]

    def get_enthalpy_flow_terms(self, p):
        """Get enthalpy flow terms for phase

        Args:
            p (str): phase

        Returns:
            Expression
        """
        return self.enthalpy_flow_terms[p]

    def get_material_density_terms(self, p, j):
        """Get material density terms for phase

        Args:
            p (str): phase

        Returns:
            Expression
        """
        if self.amount_basis == AmountBasis.MOLE:
            if p == "Mix":
                return self.dens_mol
            else:
                return self.dens_mol_phase[p]
        else:
            if p == "Mix":
                return self.dens_mass
            else:
                return self.dens_mass_phase[p]

    def get_energy_density_terms(self, p):
        """Get energy density terms for phase

        Args:
            p (str): phase

        Returns:
            Expression
        """
        return self.energy_density_terms[p]

    def default_material_balance_type(self):
        """Get default material balance type suggestion"""
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        """Get default energy balance type suggestion"""
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_basis(b):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return self._state_vars_dict

    def define_display_vars(self):
        if self.amount_basis == AmountBasis.MOLE:
            return {
                "Molar Flow": self.flow_mol,
                "Mass Flow": self.flow_mass,
                "T": self.temperature,
                "P": self.pressure,
                "Vapor Fraction": self.vapor_frac,
                "Molar Enthalpy": self.enth_mol,
            }
        else:
            return {
                "Molar Flow": self.flow_mol,
                "Mass Flow": self.flow_mass,
                "T": self.temperature,
                "P": self.pressure,
                "Vapor Fraction": self.vapor_frac,
                "Mass Enthalpy": self.enth_mass,
            }

    def extensive_state_vars(self):
        """Return the set of extensive variables"""
        return self.extensive_set

    def intensive_state_vars(self):
        """Return the set of intensive variables"""
        return self.intensive_set

    def model_check(self):
        """Currently doesn't do anything, here for compatibility"""

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if self.params.config.amount_basis == AmountBasis.MOLE:
            sf_flow = iscale.get_scaling_factor(self.flow_mol, default=1)
            sf_enth = iscale.get_scaling_factor(self.enth_mol, default=1)
            sf_inte = iscale.get_scaling_factor(self.energy_internal_mol, default=1)
            sf_dens = iscale.get_scaling_factor(self.dens_mol, default=1)
        else:
            sf_flow = iscale.get_scaling_factor(self.flow_mass, default=1)
            sf_enth = iscale.get_scaling_factor(self.enth_mass, default=1)
            sf_inte = iscale.get_scaling_factor(self.energy_internal_mass, default=1)
            sf_dens = iscale.get_scaling_factor(self.dens_mass, default=1)
        sf_pres = iscale.get_scaling_factor(self.pressure, default=1)

        for v in self.material_flow_terms.values():
            iscale.set_scaling_factor(v, sf_flow)
        for v in self.enthalpy_flow_terms.values():
            iscale.set_scaling_factor(v, sf_enth * sf_flow)
        for k, v in self.energy_density_terms.items():
            if k == "Mix":
                iscale.set_scaling_factor(v, sf_inte * sf_dens)
            else:
                if self.params.config.amount_basis == AmountBasis.MOLE:
                    sf_inte_p = iscale.get_scaling_factor(
                        self.energy_internal_mol_phase[k], default=1
                    )
                    sf_dens_p = iscale.get_scaling_factor(
                        self.dens_mol_phase[k], default=1
                    )
                else:
                    sf_inte_p = iscale.get_scaling_factor(
                        self.energy_internal_mass_phase[k], default=1
                    )
                    sf_dens_p = iscale.get_scaling_factor(
                        self.dens_mass_phase[k], default=1
                    )
                iscale.set_scaling_factor(v, sf_inte_p * sf_dens_p)
        try:
            iscale.set_scaling_factor(self.eq_sat, sf_pres / 1000.0)
        except AttributeError:
            pass  # may not have eq_sat, and that's ok
        try:
            iscale.set_scaling_factor(self.eq_complementarity, sf_pres / 10)
        except AttributeError:
            pass  # may not have eq_complementarity which is fine
