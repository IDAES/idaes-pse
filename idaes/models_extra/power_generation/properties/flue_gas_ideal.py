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
Basic property package for flue gas.

Main assumptions:
    - ideal gas
    - components in flue gas: O2, N2, NO, CO2, H2O, SO2
"""
# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    Param,
    PositiveReals,
    Reals,
    value,
    log,
    sqrt,
    Var,
    Expression,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.common.config import ConfigValue

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    Component,
    VaporPhase,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_activated_constraints,
)
from idaes.core import MaterialBalanceType, EnergyBalanceType, MaterialFlowBasis
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
from idaes.core.util import constants
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import ConfigurationError, InitializationError

# Import Python libraries
import idaes.logger as idaeslog


# Some more inforation about this module
__author__ = "Boiler Subsystem Team  J. Ma, M. Zamarripa, T. Burgard"
__version__ = "3"

# Set up logger
_log = idaeslog.getLogger("idaes.unit_model.properties")


@declare_process_block_class("FlueGasParameterBlock")
class FlueGasParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class

    Contains parameters and indexing sets associated with properties for flue
    gas. The ideal gas assumption is applied.

    """

    CONFIG = PhysicalParameterBlock.CONFIG()
    CONFIG.declare(
        "components",
        ConfigValue(
            default=["N2", "O2", "NO", "CO2", "H2O", "SO2"],
            domain=list,
            description="Components to include",
        ),
    )

    def build(self):
        """Add contents to the block."""
        super().build()
        self._state_block_class = FlueGasStateBlock
        _valid_comps = ["N2", "O2", "NO", "CO2", "H2O", "SO2"]

        for j in self.config.components:
            if j not in _valid_comps:
                raise ConfigurationError(f"Component '{j}' is not supported")
            self.add_component(j, Component())

        # Create Phase object
        self.Vap = VaporPhase()

        # Molecular weight
        self.mw_comp = Param(
            self.component_list,
            initialize={
                k: v
                for k, v in {
                    "O2": 0.0319988,
                    "N2": 0.0280134,
                    "NO": 0.0300061,
                    "CO2": 0.0440095,
                    "H2O": 0.0180153,
                    "SO2": 0.064064,
                }.items()
                if k in self.component_list
            },
            doc="Molecular Weight [kg/mol]",
            units=pyunits.kg / pyunits.mol,
        )

        # Thermodynamic reference state
        self.pressure_ref = Param(
            within=PositiveReals,
            default=1.01325e5,
            doc="Reference pressure [Pa]",
            units=pyunits.Pa,
        )

        self.temperature_ref = Param(
            within=PositiveReals,
            default=298.15,
            doc="Reference temperature [K]",
            units=pyunits.K,
        )

        # Critical Properties
        self.pressure_crit = Param(
            self.component_list,
            within=PositiveReals,
            initialize={
                k: v
                for k, v in {
                    "O2": 50.45985e5,
                    "N2": 33.943875e5,
                    "NO": 64.85e5,
                    "CO2": 73.8e5,
                    "H2O": 220.64e5,
                    "SO2": 7.883e6,
                }.items()
                if k in self.component_list
            },
            doc="Critical pressure [Pa]",
            units=pyunits.Pa,
        )

        self.temperature_crit = Param(
            self.component_list,
            within=PositiveReals,
            initialize={
                k: v
                for k, v in {
                    "O2": 154.58,
                    "N2": 126.19,
                    "NO": 180.0,
                    "CO2": 304.18,
                    "H2O": 647,
                    "SO2": 430.8,
                }.items()
                if k in self.component_list
            },
            doc="Critical temperature [K]",
            units=pyunits.K,
        )

        # Constants for specific heat capacity, enthalpy, and entropy
        # calculations for ideal gas (from NIST 01/08/2020
        # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Units=SI&Mask=1#Thermo-Gas)
        cp_mol_ig_comp_coeff_parameter_A = {
            k: v
            for k, v in {
                "N2": 19.50583,
                "O2": 30.03235,
                "CO2": 24.99735,
                "H2O": 30.092,
                "NO": 23.83491,
                "SO2": 21.43049,
            }.items()
            if k in self.component_list
        }
        cp_mol_ig_comp_coeff_parameter_B = {
            k: v
            for k, v in {
                "N2": 19.88705,
                "O2": 8.772972,
                "CO2": 55.18696,
                "H2O": 6.832514,
                "NO": 12.58878,
                "SO2": 74.35094,
            }.items()
            if k in self.component_list
        }
        cp_mol_ig_comp_coeff_parameter_C = {
            k: v
            for k, v in {
                "N2": -8.598535,
                "O2": -3.98813,
                "CO2": -33.69137,
                "H2O": 6.793435,
                "NO": -1.139011,
                "SO2": -57.75217,
            }.items()
            if k in self.component_list
        }
        cp_mol_ig_comp_coeff_parameter_D = {
            k: v
            for k, v in {
                "N2": 1.369784,
                "O2": 0.788313,
                "CO2": 7.948387,
                "H2O": -2.53448,
                "NO": -1.497459,
                "SO2": 16.35534,
            }.items()
            if k in self.component_list
        }
        cp_mol_ig_comp_coeff_parameter_E = {
            k: v
            for k, v in {
                "N2": 0.527601,
                "O2": -0.7416,
                "CO2": -0.136638,
                "H2O": 0.082139,
                "NO": 0.214194,
                "SO2": 0.086731,
            }.items()
            if k in self.component_list
        }
        cp_mol_ig_comp_coeff_parameter_F = {
            k: v
            for k, v in {
                "N2": -4.935202,
                "O2": -11.3247,
                "CO2": -403.6075,
                "H2O": -250.881,
                "NO": 83.35783,
                "SO2": -305.7688,
            }.items()
            if k in self.component_list
        }
        cp_mol_ig_comp_coeff_parameter_G = {
            k: v
            for k, v in {
                "N2": 212.39,
                "O2": 236.1663,
                "CO2": 228.2431,
                "H2O": 223.3967,
                "NO": 237.1219,
                "SO2": 254.8872,
            }.items()
            if k in self.component_list
        }
        cp_mol_ig_comp_coeff_parameter_H = {
            k: v
            for k, v in {
                "N2": 0,
                "O2": 0,
                "CO2": -393.5224,
                "H2O": -241.8264,
                "NO": 90.29114,
                "SO2": -296.8422,
            }.items()
            if k in self.component_list
        }

        self.cp_mol_ig_comp_coeff_A = Param(
            self.component_list,
            initialize=cp_mol_ig_comp_coeff_parameter_A,
            doc="Constants for spec. heat capacity for ideal gas",
            units=pyunits.J / pyunits.mol / pyunits.K,
        )
        self.cp_mol_ig_comp_coeff_B = Param(
            self.component_list,
            initialize=cp_mol_ig_comp_coeff_parameter_B,
            doc="Constants for spec. heat capacity for ideal gas",
            units=pyunits.J / pyunits.mol / pyunits.K / pyunits.kK,
        )
        self.cp_mol_ig_comp_coeff_C = Param(
            self.component_list,
            initialize=cp_mol_ig_comp_coeff_parameter_C,
            doc="Constants for spec. heat capacity for ideal gas",
            units=pyunits.J / pyunits.mol / pyunits.K / pyunits.kK**2,
        )
        self.cp_mol_ig_comp_coeff_D = Param(
            self.component_list,
            initialize=cp_mol_ig_comp_coeff_parameter_D,
            doc="Constants for spec. heat capacity for ideal gas",
            units=pyunits.J / pyunits.mol / pyunits.K / pyunits.kK**3,
        )
        self.cp_mol_ig_comp_coeff_E = Param(
            self.component_list,
            initialize=cp_mol_ig_comp_coeff_parameter_E,
            doc="Constants for spec. heat capacity for ideal gas",
            units=pyunits.J / pyunits.mol / pyunits.K * pyunits.kK**2,
        )
        self.cp_mol_ig_comp_coeff_F = Param(
            self.component_list,
            initialize=cp_mol_ig_comp_coeff_parameter_F,
            doc="Constants for spec. heat capacity for ideal gas",
            units=pyunits.kJ / pyunits.mol,
        )
        self.cp_mol_ig_comp_coeff_G = Param(
            self.component_list,
            initialize=cp_mol_ig_comp_coeff_parameter_G,
            doc="Constants for spec. heat capacity for ideal gas",
            units=pyunits.J / pyunits.mol / pyunits.K,
        )
        self.cp_mol_ig_comp_coeff_H = Param(
            self.component_list,
            initialize=cp_mol_ig_comp_coeff_parameter_H,
            doc="Constants for spec. heat capacity for ideal gas",
            units=pyunits.kJ / pyunits.mol,
        )

        # Viscosity and thermal conductivity parameters
        self.ce_param = Param(
            initialize=2.6693e-5,
            units=(
                pyunits.g**0.5
                * pyunits.mol**0.5
                * pyunits.angstrom**2
                * pyunits.K**-0.5
                * pyunits.cm**-1
                * pyunits.s**-1
            ),
            doc="Parameter for the Chapman-Enskog viscosity correlation",
        )

        self.sigma = Param(
            self.component_list,
            initialize={
                k: v
                for k, v in {
                    "O2": 3.458,
                    "N2": 3.621,
                    "NO": 3.47,
                    "CO2": 3.763,
                    "H2O": 2.605,
                    "SO2": 4.29,
                }.items()
                if k in self.component_list
            },
            doc="collision diameter",
            units=pyunits.angstrom,
        )
        self.ep_Kappa = Param(
            self.component_list,
            initialize={
                k: v
                for k, v in {
                    "O2": 107.4,
                    "N2": 97.53,
                    "NO": 119.0,
                    "CO2": 244.0,
                    "H2O": 572.4,
                    "SO2": 252.0,
                }.items()
                if k in self.component_list
            },
            doc="Boltzmann constant divided by characteristic Lennard-Jones energy",
            units=pyunits.K,
        )

        self.set_default_scaling("flow_mol", 1e-4)
        self.set_default_scaling("flow_mass", 1e-3)
        self.set_default_scaling("flow_vol", 1e-3)
        # anything not explicitly listed
        self.set_default_scaling("mole_frac_comp", 1)
        self.set_default_scaling("mole_frac_comp", 1e3, index="NO")
        self.set_default_scaling("mole_frac_comp", 1e3, index="SO2")
        self.set_default_scaling("mole_frac_comp", 1e2, index="H2O")
        self.set_default_scaling("mole_frac_comp", 1e2, index="CO2")
        self.set_default_scaling("flow_vol", 1)

        # For flow_mol_comp, will calculate from flow_mol and mole_frac_comp
        # user should set a scale for both, and for each compoent of
        # mole_frac_comp
        self.set_default_scaling("pressure", 1e-5)
        self.set_default_scaling("temperature", 1e-1)
        self.set_default_scaling("pressure_red", 1e-3)
        self.set_default_scaling("temperature_red", 1)
        self.set_default_scaling("enth_mol_phase", 1e-3)
        self.set_default_scaling("enth_mol", 1e-3)
        self.set_default_scaling("entr_mol", 1e-2)
        self.set_default_scaling("entr_mol_phase", 1e-2)
        self.set_default_scaling("cp_mol", 0.1)
        self.set_default_scaling("cp_mol_phase", 0.1)
        self.set_default_scaling("compress_fact", 1)
        self.set_default_scaling("dens_mol_phase", 1)
        self.set_default_scaling("pressure_sat", 1e-4)
        self.set_default_scaling("visc_d_comp", 1e4)
        self.set_default_scaling("therm_cond_comp", 1e2)
        self.set_default_scaling("visc_d", 1e4)
        self.set_default_scaling("therm_cond", 1e2)
        self.set_default_scaling("mw", 1)
        self.set_default_scaling("mw_comp", 1)

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_mol_comp": {"method": None, "units": "mol/s"},
                "pressure": {"method": None, "units": "Pa"},
                "temperature": {"method": None, "units": "K"},
                "pressure_crit": {"method": None, "units": "Pa"},
                "temperature_crit": {"method": None, "units": "K"},
                "pressure_red": {"method": None, "units": None},
                "temperature_red": {"method": None, "units": None},
                "enth_mol_phase": {"method": "_enthalpy_calc", "units": "J/mol"},
                "entr_mol_phase": {"method": "_entropy_calc", "units": "J/mol/K"},
                "enth_mol": {"method": "_enthalpy_calc", "units": "J/mol"},
                "entr_mol": {"method": "_entropy_calc", "units": "J/mol.K"},
                "cp_mol": {"method": "_heat_cap_calc", "units": "J/mol.K"},
                "cp_mol_phase": {"method": "_heat_cap_calc", "units": "J/mol/K"},
                "compress_fact": {"method": "_compress_fact", "units": None},
                "dens_mol_phase": {"method": "_dens_mol_phase", "units": "mol/m^3"},
                "pressure_sat": {"method": "_vapor_pressure", "units": "Pa"},
                "flow_vol": {"method": "_flow_volume", "units": "m^3/s"},
                "visc_d": {"method": "_therm_cond", "units": "kg/m-s"},
                "therm_cond": {"method": "_therm_cond", "units": "W/m-K"},
                "mw_comp": {"method": None, "units": "kg/mol"},
                "mw": {"method": None, "units": "kg/mol"},
            }
        )

        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )


class _FlueGasStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(
        self,
        state_args=None,
        hold_state=False,
        state_vars_fixed=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """Initialisation routine for property package.

        Key values for the state_args dict:
            flow_mol_comp : value at which to initialize component flows
                (default=27.5e3 mol/s)
            pressure : value at which to initialize pressure
                (default=2.97e7 Pa)
            temperature : value at which to initialize temperature
                (default=866.5 K)

        Args:
            outlvl: sets logging level
            state_vars_fixed: Flag to denote state vars have been fixed.
                - True - states have been fixed by the control volume 1D.
                         Control volume 0D does not fix the state vars, so will
                         be False if this state block is used with 0D blocks.
                - False - states have not been fixed. The state block will deal
                          with fixing/unfixing.
            optarg: solver options dictionary object (default=None, use
                    default solver options)
            solver: str indicating which solver to use during
                     initialization (default = None, use default solver)
            hold_state: flag indicating whether the initialization routine
                should unfix any state variables fixed during initialization
                (default=False).
                - True - states varaibles are not unfixed, and a dict of
                         returned containing flags for which states were fixed
                         during initialization.
                - False - state variables are unfixed after initialization by
                          calling the relase_state method

            Returns:
                If hold_states is True, returns a dict containing flags for
                which states were fixed during initialization.
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")

        # Create solver
        opt = get_solver(solver, optarg)

        if state_args is None:
            state_args = {
                "flow_mol_comp": {
                    "N2": 1.0,
                    "CO2": 1.0,
                    "NO": 1.0,
                    "O2": 1.0,
                    "H2O": 1.0,
                    "SO2": 1.0,
                },
                "pressure": 1e5,
                "temperature": 495.0,
            }

        if state_vars_fixed is False:
            flags = fix_state_vars(self, state_args)
        # Check when the state vars are fixed already result in dof 0
        for b in self.values():
            if degrees_of_freedom(b) != 0:
                raise InitializationError(
                    f"{self.name} State vars fixed but degrees of freedom not 0"
                )
        # ---------------------------------------------------------------------
        # Solve 1st stage
        for k, b in self.items():
            deactivate_list = []
            if hasattr(b, "enthalpy_correlation"):
                deactivate_list.append(b.enthalpy_correlation)
            if hasattr(b, "volumetric_flow_calculation"):
                deactivate_list.append(b.volumetric_flow_calculation)
            if hasattr(b, "entropy_correlation"):
                deactivate_list.append(b.entropy_correlation)
            for c in deactivate_list:
                c.deactivate()

            if number_activated_constraints(b) > 0:
                with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                    res = opt.solve(b, tee=slc.tee)
            else:
                res = "skipped"
            init_log.info_high(
                "Initialization Step 1 {}.".format(idaeslog.condition(res))
            )

            for c in deactivate_list:
                c.activate()

            if number_activated_constraints(b) > 0:
                with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                    res = opt.solve(b, tee=slc.tee)
            else:
                res = "skipped"
            init_log.info_high(
                "Initialization Step 2 {}.".format(idaeslog.condition(res))
            )

        if not res == "skipped" and not check_optimal_termination(res):
            raise InitializationError(f"Solve failed {res}.")

        init_log.info(f"Initialisation Complete, {idaeslog.condition(res)}.")
        # ---------------------------------------------------------------------
        # If input block, return flags, else release state
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                self.release_state(flags)

    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        """
        Method to relase state variables fixed during initialisation.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        """
        # Unfix state variables
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        revert_state_vars(self, flags)
        init_log.info("{} State Released.".format(self.name))


@declare_process_block_class("FlueGasStateBlock", block_class=_FlueGasStateBlock)
class FlueGasStateBlockData(StateBlockData):
    """
    This is an example of a property package for calculating the thermophysical
    properties of flue gas using the ideal gas assumption.
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(FlueGasStateBlockData, self).build()
        comps = self.params.component_list
        # Add state variables
        self.flow_mol_comp = Var(
            comps,
            domain=Reals,
            initialize=1.0,
            bounds=(0, 1e6),
            doc="Component molar flowrate [mol/s]",
            units=pyunits.mol / pyunits.s,
        )
        self.pressure = Var(
            domain=Reals,
            initialize=1.01325e5,
            bounds=(1, 5e7),
            doc="State pressure [Pa]",
            units=pyunits.Pa,
        )
        self.temperature = Var(
            domain=Reals,
            initialize=500,
            bounds=(200, 1500),
            doc="State temperature [K]",
            units=pyunits.K,
        )

        # Add expressions for some basic oft-used quantiies
        self.flow_mol = Expression(expr=sum(self.flow_mol_comp[j] for j in comps))

        def rule_mole_frac(b, c):
            return b.flow_mol_comp[c] / b.flow_mol

        self.mole_frac_comp = Expression(
            comps, rule=rule_mole_frac, doc="mole fraction of component i"
        )

        self.flow_mass = Expression(
            expr=sum(self.flow_mol_comp[j] * self.params.mw_comp[j] for j in comps),
            doc="total mass flow",
        )

        def rule_mw_comp(b, j):
            return b.params.mw_comp[j]

        self.mw_comp = Expression(comps, rule=rule_mw_comp)

        def rule_mw(b):
            return sum(b.mw_comp[j] * b.mole_frac_comp[j] for j in comps)

        self.mw = Expression(rule=rule_mw)

        self.pressure_crit = Expression(
            expr=sum(
                self.params.pressure_crit[j] * self.mole_frac_comp[j] for j in comps
            )
        )
        self.temperature_crit = Expression(
            expr=sum(
                self.params.temperature_crit[j] * self.mole_frac_comp[j] for j in comps
            )
        )
        self.pressure_red = Expression(expr=self.pressure / self.pressure_crit)
        self.temperature_red = Expression(expr=self.temperature / self.temperature_crit)

        self.compress_fact = Expression(expr=1.0, doc="Vapor Compressibility Factor")

        def rule_dens_mol_phase(b, p):
            return (
                b.pressure
                / b.compress_fact
                / constants.Constants.gas_constant
                / b.temperature
            )

        self.dens_mol_phase = Expression(
            self.params.phase_list, rule=rule_dens_mol_phase, doc="Molar Density"
        )

        self.flow_vol = Expression(
            doc="Volumetric Flowrate", expr=self.flow_mol / self.dens_mol_phase["Vap"]
        )

    def _heat_cap_calc(self):
        # heat capacity J/mol-K
        self.cp_mol = Var(
            initialize=1000,
            doc="heat capacity [J/mol-K]",
            units=pyunits.J / pyunits.mol / pyunits.K,
        )

        def rule_cp_phase(b, p):
            # This property module only has one phase
            return self.cp_mol

        self.cp_mol_phase = Expression(self.params.phase_list, rule=rule_cp_phase)

        try:
            ft = sum(self.flow_mol_comp[j] for j in self.params.component_list)
            t = pyunits.convert(self.temperature, to_units=pyunits.kK)
            self.heat_cap_correlation = Constraint(
                expr=(
                    self.cp_mol * ft
                    == sum(
                        self.flow_mol_comp[j]
                        * (
                            self.params.cp_mol_ig_comp_coeff_A[j]
                            + self.params.cp_mol_ig_comp_coeff_B[j] * t
                            + self.params.cp_mol_ig_comp_coeff_C[j] * t**2
                            + self.params.cp_mol_ig_comp_coeff_D[j] * t**3
                            + self.params.cp_mol_ig_comp_coeff_E[j] / t**2
                        )
                        for j in self.params.component_list
                    )
                )
            )
        except AttributeError:
            self.del_component(self.cp_mol)
            self.del_component(self.heat_cap_correlation)

    def _enthalpy_calc(self):
        self.enth_mol = Var(
            doc="Specific Enthalpy [J/mol]", units=pyunits.J / pyunits.mol
        )

        def rule_enth_phase(b, p):
            # This property module only has one phase
            return self.enth_mol

        self.enth_mol_phase = Expression(self.params.phase_list, rule=rule_enth_phase)

        def enthalpy_correlation(b):
            ft = sum(self.flow_mol_comp[j] for j in self.params.component_list)
            t = pyunits.convert(self.temperature, to_units=pyunits.kK)
            return self.enth_mol * ft == sum(
                self.flow_mol_comp[j]
                * pyunits.convert(
                    self.params.cp_mol_ig_comp_coeff_A[j] * t
                    + self.params.cp_mol_ig_comp_coeff_B[j] * t**2 / 2
                    + self.params.cp_mol_ig_comp_coeff_C[j] * t**3 / 3
                    + self.params.cp_mol_ig_comp_coeff_D[j] * t**4 / 4
                    - self.params.cp_mol_ig_comp_coeff_E[j] / t
                    + self.params.cp_mol_ig_comp_coeff_F[j],
                    to_units=pyunits.J / pyunits.mol,
                )
                for j in self.params.component_list
            )
            # NOTE: the H term (from the Shomate Equation) is not
            # included here so that the reference state enthalpy is the
            # enthalpy of formation (not 0).

        try:
            self.enthalpy_correlation = Constraint(rule=enthalpy_correlation)
        except AttributeError:
            self.del_component(self.enth_mol_phase)
            self.del_component(self.enth_mol)
            self.del_component(self.enthalpy_correlation)

    def _entropy_calc(self):
        self.entr_mol = Var(
            doc="Specific Entropy [J/mol/K]", units=pyunits.J / pyunits.mol / pyunits.K
        )
        # Specific Entropy

        def rule_entr_phase(b, p):
            # This property module only has one phase
            return self.entr_mol

        self.entr_mol_phase = Expression(self.params.phase_list, rule=rule_entr_phase)

        def entropy_correlation(b):
            ft = sum(self.flow_mol_comp[j] for j in self.params.component_list)
            t = pyunits.convert(self.temperature, to_units=pyunits.kK)
            n = self.flow_mol_comp
            x = self.mole_frac_comp
            p = self.pressure
            r_gas = constants.Constants.gas_constant
            return (self.entr_mol + r_gas * log(p / 1e5)) * ft == sum(
                n[j]
                * (
                    self.params.cp_mol_ig_comp_coeff_A[j] * log(t)
                    + self.params.cp_mol_ig_comp_coeff_B[j] * t
                    + self.params.cp_mol_ig_comp_coeff_C[j] * t**2 / 2
                    + self.params.cp_mol_ig_comp_coeff_D[j] * t**3 / 3
                    - self.params.cp_mol_ig_comp_coeff_E[j] / t**2 / 2
                    + self.params.cp_mol_ig_comp_coeff_G[j]
                    + r_gas * log(x[j])
                )
                for j in self.params.component_list
            )

        try:
            self.entropy_correlation = Constraint(rule=entropy_correlation)
        except AttributeError:
            self.del_component(self.entr_mol_phase)
            self.del_component(self.entropy_correlation)

    def _therm_cond(self):
        comps = self.params.component_list
        self.therm_cond_comp = Var(
            comps,
            initialize=0.05,
            doc="thermal conductivity J/m-K-s",
            units=pyunits.J / pyunits.m / pyunits.K / pyunits.s,
        )
        self.therm_cond = Var(
            initialize=0.05,
            doc="thermal conductivity of gas mixture J/m-K-s",
            units=pyunits.J / pyunits.m / pyunits.K / pyunits.s,
        )
        self.visc_d_comp = Var(
            comps,
            initialize=2e-5,
            doc="dynamic viscocity of pure gas species",
            units=pyunits.kg / pyunits.m / pyunits.s,
        )
        self.visc_d = Var(
            initialize=2e-5,
            doc="viscosity of gas mixture kg/m-s",
            units=pyunits.kg / pyunits.m / pyunits.s,
        )

        try:

            def rule_therm_cond(b, c):
                t = pyunits.convert(b.temperature, to_units=pyunits.kK)
                return (
                    b.therm_cond_comp[c]
                    == (
                        (
                            (
                                b.params.cp_mol_ig_comp_coeff_A[c]
                                + b.params.cp_mol_ig_comp_coeff_B[c] * t
                                + b.params.cp_mol_ig_comp_coeff_C[c] * t**2
                                + b.params.cp_mol_ig_comp_coeff_D[c] * t**3
                                + b.params.cp_mol_ig_comp_coeff_E[c] / t**2
                            )
                            / b.params.mw_comp[c]
                        )
                        + 1.25
                        * (constants.Constants.gas_constant / b.params.mw_comp[c])
                    )
                    * b.visc_d_comp[c]
                )

            self.therm_cond_con = Constraint(comps, rule=rule_therm_cond)

            def rule_theta(b, c):
                return b.temperature / b.params.ep_Kappa[c]

            self.theta = Expression(comps, rule=rule_theta)

            def rule_omega(b, c):
                return (
                    1.5794145
                    + 0.00635771 * b.theta[c]
                    - 0.7314 * log(b.theta[c])
                    + 0.2417357 * log(b.theta[c]) ** 2
                    - 0.0347045 * log(b.theta[c]) ** 3
                )

            self.omega = Expression(comps, rule=rule_omega)

            # Pure gas viscocity - from Chapman-Enskog theory
            def rule_visc_d(b, c):
                return (
                    pyunits.convert(
                        b.visc_d_comp[c], to_units=pyunits.g / pyunits.cm / pyunits.s
                    )
                    == b.params.ce_param
                    * sqrt(
                        pyunits.convert(
                            b.params.mw_comp[c], to_units=pyunits.g / pyunits.mol
                        )
                        * b.temperature
                    )
                    / b.params.sigma[c] ** 2
                    / b.omega[c]
                )

            self.visc_d_con = Constraint(comps, rule=rule_visc_d)

            # section to calculate viscosity of gas mixture
            def rule_phi(b, i, j):
                return (
                    1
                    / 2.8284
                    * (1 + (b.params.mw_comp[i] / b.params.mw_comp[j])) ** (-0.5)
                    * (
                        1
                        + sqrt(b.visc_d_comp[i] / b.visc_d_comp[j])
                        * (b.params.mw_comp[j] / b.params.mw_comp[i]) ** 0.25
                    )
                    ** 2
                )

            self.phi_ij = Expression(comps, comps, rule=rule_phi)

            # viscosity of Gas mixture kg/m-s
            def rule_visc_d_mix(b):
                return b.visc_d == sum(
                    (b.mole_frac_comp[i] * b.visc_d_comp[i])
                    / sum(b.mole_frac_comp[j] * b.phi_ij[i, j] for j in comps)
                    for i in comps
                )

            self.vis_d_mix_con = Constraint(rule=rule_visc_d_mix)

            # thermal conductivity of gas mixture in kg/m-s
            def rule_therm_mix(b):
                return b.therm_cond == sum(
                    (b.mole_frac_comp[i] * b.therm_cond_comp[i])
                    / sum(b.mole_frac_comp[j] * b.phi_ij[i, j] for j in comps)
                    for i in comps
                )

            self.therm_mix_con = Constraint(rule=rule_therm_mix)

        except AttributeError:
            self.del_component(self.therm_cond_comp)
            self.del_component(self.therm_cond)
            self.del_component(self.visc_d_comp)
            self.del_component(self.visc_d)
            self.del_component(self.omega)
            self.del_component(self.theta)
            self.del_component(self.phi_ij)
            self.del_component(self.sigma)
            self.del_component(self.ep_Kappa)
            self.del_component(self.therm_cond_con)
            self.del_component(self.theta_con)
            self.del_component(self.omega_con)
            self.del_component(self.visc_d_con)
            self.del_component(self.phi_con)

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_terms(self, p, j):
        return self.flow_mol_comp[j]

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def get_enthalpy_flow_terms(self, p):
        if not self.is_property_constructed("enthalpy_flow_terms"):
            try:

                def rule_enthalpy_flow_terms(b, p):
                    return self.enth_mol_phase[p] * self.flow_mol

                self.enthalpy_flow_terms = Expression(
                    self.params.phase_list, rule=rule_enthalpy_flow_terms
                )
            except AttributeError:
                self.del_component(self.enthalpy_flow_terms)
        return self.enthalpy_flow_terms[p]

    def get_material_density_terms(self, p, j):
        return self.dens_mol_phase[p]

    def get_energy_density_terms(self, p):
        if not self.is_property_constructed("energy_density_terms"):
            try:

                def rule_energy_density_terms(b, p):
                    return (
                        self.enth_mol_phase[p] * self.dens_mol_phase[p] - self.pressure
                    )

                self.energy_density_terms = Expression(
                    self.params.phase_list, rule=rule_energy_density_terms
                )
            except AttributeError:
                self.del_component(self.energy_density_terms)
        return self.energy_density_terms[p]

    def define_state_vars(self):
        return {
            "flow_mol_comp": self.flow_mol_comp,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }

    def model_check(self):
        """
        Model checks for property block
        """
        # Check temperature bounds
        for v in self.compoent_object_data(Var, descend_into=True):
            if value(v) < v.lb:
                _log.error(f"{v} is below lower bound in {self.name}")
            if value(v) > v.ub:
                _log.error(f"{v} is above upper bound in {self.name}")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # Get some scale factors that are frequently used to calculate others
        sf_flow = iscale.get_scaling_factor(self.flow_mol)
        sf_mol_fraction = {}
        comps = self.params.component_list
        for i in comps:
            sf_mol_fraction[i] = iscale.get_scaling_factor(self.mole_frac_comp[i])
        # calculate flow_mol_comp scale factors
        for i, c in self.flow_mol_comp.items():
            iscale.set_scaling_factor(c, sf_flow * sf_mol_fraction[i])

        if self.is_property_constructed("energy_density_terms"):
            for i, c in self.energy_density_terms.items():
                sf1 = iscale.get_scaling_factor(self.enth_mol_phase[i])
                sf2 = iscale.get_scaling_factor(self.dens_mol_phase[i])
                iscale.set_scaling_factor(c, sf1 * sf2)

        if self.is_property_constructed("enthalpy_flow_terms"):
            for i, c in self.enthalpy_flow_terms.items():
                sf1 = iscale.get_scaling_factor(self.enth_mol_phase[i])
                sf2 = iscale.get_scaling_factor(self.flow_mol)
                iscale.set_scaling_factor(c, sf1 * sf2)

        if self.is_property_constructed("heat_cap_correlation"):
            iscale.constraint_scaling_transform(
                self.heat_cap_correlation,
                iscale.get_scaling_factor(self.cp_mol)
                * iscale.get_scaling_factor(self.flow_mol),
                overwrite=False,
            )
        if self.is_property_constructed("enthalpy_correlation"):
            for p, c in self.enthalpy_correlation.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.enth_mol)
                    * iscale.get_scaling_factor(self.flow_mol),
                    overwrite=False,
                )
        if self.is_property_constructed("entropy_correlation"):
            iscale.constraint_scaling_transform(
                self.entropy_correlation,
                iscale.get_scaling_factor(self.entr_mol)
                * iscale.get_scaling_factor(self.flow_mol),
                overwrite=False,
            )
        if self.is_property_constructed("vapor_pressure_correlation"):
            iscale.constraint_scaling_transform(
                self.vapor_pressure_correlation,
                log(iscale.get_scaling_factor(self.pressure_sat))
                * iscale.get_scaling_factor(self.flow_mol),
                overwrite=False,
            )
        if self.is_property_constructed("therm_cond_con"):
            for i, c in self.therm_cond_con.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.therm_cond_comp[i]),
                    overwrite=False,
                )
        if self.is_property_constructed("therm_mix_con"):
            iscale.constraint_scaling_transform(
                self.therm_mix_con,
                iscale.get_scaling_factor(self.therm_cond),
                overwrite=False,
            )
        if self.is_property_constructed("visc_d_con"):
            for i, c in self.visc_d_con.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.visc_d_comp[i]), overwrite=False
                )
            for i, c in self.vis_d_mix_con.items():
                iscale.constraint_scaling_transform(
                    self.vis_d_mix_con,
                    iscale.get_scaling_factor(self.visc_d),
                    overwrite=False,
                )
