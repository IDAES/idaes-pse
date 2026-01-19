##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Example ideal parameter block for the VLE calculations for a
Benzene-Toluene-o-Xylene system.
"""


# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    Expression,
    log,
    NonNegativeReals,
    Var,
    Set,
    Param,
    sqrt,
    log10,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.config import ConfigValue

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    Component,
    LiquidPhase,
    VaporPhase,
)
from idaes.core.util.constants import Constants as const
from idaes.core.util.initialization import fix_state_vars, solve_indexed_blocks
from idaes.core.initialization import InitializerBase
from idaes.core.util.misc import add_object_reference
from idaes.core.util.model_statistics import number_unfixed_variables
from idaes.core.util.misc import extract_data
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)


class HDAInitializer(InitializerBase):
    """
    Initializer for HDA Property package.

    """

    CONFIG = InitializerBase.CONFIG()
    CONFIG.declare(
        "solver",
        ConfigValue(default=None, domain=str, description="Initialization solver"),
    )
    CONFIG.declare(
        "solver_options",
        ConfigValue(default=None, description="Initialization solver options"),
    )

    def initialization_routine(self, blk):
        init_log = idaeslog.getInitLogger(
            blk.name, self.config.output_level, tag="properties"
        )
        solve_log = idaeslog.getSolveLogger(
            blk.name, self.config.output_level, tag="properties"
        )

        # Set solver
        solver = get_solver(self.config.solver, self.config.solver_options)

        # ---------------------------------------------------------------------
        # If present, initialize bubble and dew point calculations
        for k in blk.keys():
            if hasattr(blk[k], "eq_temperature_dew"):
                calculate_variable_from_constraint(
                    blk[k].temperature_dew, blk[k].eq_temperature_dew
                )

            if hasattr(blk[k], "eq_pressure_dew"):
                calculate_variable_from_constraint(
                    blk[k].pressure_dew, blk[k].eq_pressure_dew
                )

        init_log.info_high(
            "Initialization Step 1 - Dew and bubble points " "calculation completed."
        )

        # ---------------------------------------------------------------------
        # If flash, initialize T1 and Teq
        for k in blk.keys():
            if blk[k].config.has_phase_equilibrium and not blk[k].config.defined_state:
                blk[k]._t1.value = max(
                    blk[k].temperature.value, blk[k].temperature_bubble.value
                )
                blk[k]._teq.value = min(blk[k]._t1.value, blk[k].temperature_dew.value)

        init_log.info_high(
            "Initialization Step 2 - Equilibrium temperature " " calculation completed."
        )

        # ---------------------------------------------------------------------
        # Initialize flow rates and compositions
        free_vars = 0
        for k in blk.keys():
            free_vars += number_unfixed_variables(blk[k])
        if free_vars > 0:
            try:
                with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                    res = solve_indexed_blocks(solver, [blk], tee=slc.tee)
            except:
                res = None
        else:
            res = None

        init_log.info("Initialization Complete")

        return res


@declare_process_block_class("HDAParameterBlock")
class HDAParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    def build(self):
        """
        Callable method for Block construction.
        """
        super(HDAParameterData, self).build()

        self._state_block_class = IdealStateBlock

        self.benzene = Component()
        self.toluene = Component()
        self.methane = Component()
        self.hydrogen = Component()

        self.Liq = LiquidPhase()
        self.Vap = VaporPhase()

        # List of components in each phase (optional)
        self.phase_comp = {"Liq": self.component_list, "Vap": self.component_list}

        # List of phase equilibrium index
        self.phase_equilibrium_idx = Set(initialize=[1, 2, 3, 4])

        self.phase_equilibrium_list = {
            1: ["benzene", ("Vap", "Liq")],
            2: ["toluene", ("Vap", "Liq")],
            3: ["hydrogen", ("Vap", "Liq")],
            4: ["methane", ("Vap", "Liq")],
        }

        # Thermodynamic reference state
        self.pressure_ref = Param(
            mutable=True, default=101325, units=pyunits.Pa, doc="Reference pressure"
        )
        self.temperature_ref = Param(
            mutable=True, default=298.15, units=pyunits.K, doc="Reference temperature"
        )

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        pressure_crit_data = {
            "benzene": 48.9e5,
            "toluene": 41e5,
            "hydrogen": 12.9e5,
            "methane": 46e5,
        }

        self.pressure_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=True,
            units=pyunits.Pa,
            initialize=extract_data(pressure_crit_data),
            doc="Critical pressure",
        )

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        temperature_crit_data = {
            "benzene": 562.2,
            "toluene": 591.8,
            "hydrogen": 33.0,
            "methane": 190.4,
        }

        self.temperature_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=True,
            units=pyunits.K,
            initialize=extract_data(temperature_crit_data),
            doc="Critical temperature",
        )

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        mw_comp_data = {
            "benzene": 78.1136e-3,
            "toluene": 92.1405e-3,
            "hydrogen": 2.016e-3,
            "methane": 16.043e-3,
        }

        self.mw_comp = Param(
            self.component_list,
            mutable=True,
            units=pyunits.kg / pyunits.mol,
            initialize=extract_data(mw_comp_data),
            doc="molecular weight",
        )

        # Constants for liquid densities
        # Source: Perry's Chemical Engineers Handbook
        #         - Robert H. Perry (Cp_liq)
        dens_liq_data = {
            ("benzene", "1"): 1.0162,
            ("benzene", "2"): 0.2655,
            ("benzene", "3"): 562.16,
            ("benzene", "4"): 0.28212,
            ("toluene", "1"): 0.8488,
            ("toluene", "2"): 0.26655,
            ("toluene", "3"): 591.8,
            ("toluene", "4"): 0.2878,
            ("hydrogen", "1"): 5.414,
            ("hydrogen", "2"): 0.34893,
            ("hydrogen", "3"): 33.19,
            ("hydrogen", "4"): 0.2706,
            ("methane", "1"): 2.9214,
            ("methane", "2"): 0.28976,
            ("methane", "3"): 190.56,
            ("methane", "4"): 0.28881,
        }

        self.dens_liq_param_1 = Param(
            self.component_list,
            mutable=True,
            initialize={c: v for (c, j), v in dens_liq_data.items() if j == "1"},
            doc="Parameter 1 to compute liquid densities",
            units=pyunits.kmol * pyunits.m**-3,
        )

        self.dens_liq_param_2 = Param(
            self.component_list,
            mutable=True,
            initialize={c: v for (c, j), v in dens_liq_data.items() if j == "2"},
            doc="Parameter 2 to compute liquid densities",
            units=pyunits.dimensionless,
        )

        self.dens_liq_param_3 = Param(
            self.component_list,
            mutable=True,
            initialize={c: v for (c, j), v in dens_liq_data.items() if j == "3"},
            doc="Parameter 3 to compute liquid densities",
            units=pyunits.K,
        )

        self.dens_liq_param_4 = Param(
            self.component_list,
            mutable=True,
            initialize={c: v for (c, j), v in dens_liq_data.items() if j == "4"},
            doc="Parameter 4 to compute liquid densities",
            units=pyunits.dimensionless,
        )

        # Boiling point at standard pressure
        # Source: Perry's Chemical Engineers Handbook
        #         - Robert H. Perry (Cp_liq)
        bp_data = {
            ("benzene"): 353.25,
            ("toluene"): 383.95,
            ("hydrogen"): 20.45,
            ("methane"): 111.75,
        }

        self.temperature_boil = Param(
            self.component_list,
            mutable=True,
            units=pyunits.K,
            initialize=extract_data(bp_data),
            doc="Pure component boiling points at standard pressure",
        )

        # Constants for specific heat capacity, enthalpy
        # Sources: The Properties of Gases and Liquids (1987)
        #         4th edition, Chemical Engineering Series - Robert C. Reid
        #         Perry's Chemical Engineers Handbook
        #         - Robert H. Perry (Cp_liq)
        cp_ig_data = {
            ("Liq", "benzene", "1"): 1.29e5,
            ("Liq", "benzene", "2"): -1.7e2,
            ("Liq", "benzene", "3"): 6.48e-1,
            ("Liq", "benzene", "4"): 0,
            ("Liq", "benzene", "5"): 0,
            ("Vap", "benzene", "1"): -3.392e1,
            ("Vap", "benzene", "2"): 4.739e-1,
            ("Vap", "benzene", "3"): -3.017e-4,
            ("Vap", "benzene", "4"): 7.130e-8,
            ("Vap", "benzene", "5"): 0,
            ("Liq", "toluene", "1"): 1.40e5,
            ("Liq", "toluene", "2"): -1.52e2,
            ("Liq", "toluene", "3"): 6.95e-1,
            ("Liq", "toluene", "4"): 0,
            ("Liq", "toluene", "5"): 0,
            ("Vap", "toluene", "1"): -2.435e1,
            ("Vap", "toluene", "2"): 5.125e-1,
            ("Vap", "toluene", "3"): -2.765e-4,
            ("Vap", "toluene", "4"): 4.911e-8,
            ("Vap", "toluene", "5"): 0,
            ("Liq", "hydrogen", "1"): 0,  # 6.6653e1,
            ("Liq", "hydrogen", "2"): 0,  # 6.7659e3,
            ("Liq", "hydrogen", "3"): 0,  # -1.2363e2,
            ("Liq", "hydrogen", "4"): 0,  # 4.7827e2, # Eqn 2
            ("Liq", "hydrogen", "5"): 0,
            ("Vap", "hydrogen", "1"): 2.714e1,
            ("Vap", "hydrogen", "2"): 9.274e-3,
            ("Vap", "hydrogen", "3"): -1.381e-5,
            ("Vap", "hydrogen", "4"): 7.645e-9,
            ("Vap", "hydrogen", "5"): 0,
            ("Liq", "methane", "1"): 0,  # 6.5708e1,
            ("Liq", "methane", "2"): 0,  # 3.8883e4,
            ("Liq", "methane", "3"): 0,  # -2.5795e2,
            ("Liq", "methane", "4"): 0,  # 6.1407e2, # Eqn 2
            ("Liq", "methane", "5"): 0,
            ("Vap", "methane", "1"): 1.925e1,
            ("Vap", "methane", "2"): 5.213e-2,
            ("Vap", "methane", "3"): 1.197e-5,
            ("Vap", "methane", "4"): -1.132e-8,
            ("Vap", "methane", "5"): 0,
        }

        self.cp_ig_1 = Param(
            self.phase_list,
            self.component_list,
            mutable=True,
            initialize={(p, c): v for (p, c, j), v in cp_ig_data.items() if j == "1"},
            doc="Parameter 1 to compute Cp_comp",
            units=pyunits.J / pyunits.mol / pyunits.K,
        )

        self.cp_ig_2 = Param(
            self.phase_list,
            self.component_list,
            mutable=True,
            initialize={(p, c): v for (p, c, j), v in cp_ig_data.items() if j == "2"},
            doc="Parameter 2 to compute Cp_comp",
            units=pyunits.J / pyunits.mol / pyunits.K**2,
        )

        self.cp_ig_3 = Param(
            self.phase_list,
            self.component_list,
            mutable=True,
            initialize={(p, c): v for (p, c, j), v in cp_ig_data.items() if j == "3"},
            doc="Parameter 3 to compute Cp_comp",
            units=pyunits.J / pyunits.mol / pyunits.K**3,
        )

        self.cp_ig_4 = Param(
            self.phase_list,
            self.component_list,
            mutable=True,
            initialize={(p, c): v for (p, c, j), v in cp_ig_data.items() if j == "4"},
            doc="Parameter 4 to compute Cp_comp",
            units=pyunits.J / pyunits.mol / pyunits.K**4,
        )

        self.cp_ig_5 = Param(
            self.phase_list,
            self.component_list,
            mutable=True,
            initialize={(p, c): v for (p, c, j), v in cp_ig_data.items() if j == "5"},
            doc="Parameter 5 to compute Cp_comp",
            units=pyunits.J / pyunits.mol / pyunits.K**5,
        )

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        # fitted to Antoine form
        # H2, Methane from NIST webbook
        pressure_sat_coeff_data = {
            ("benzene", "A"): 4.202,
            ("benzene", "B"): 1322,
            ("benzene", "C"): -38.56,
            ("toluene", "A"): 4.216,
            ("toluene", "B"): 1435,
            ("toluene", "C"): -43.33,
            ("hydrogen", "A"): 3.543,
            ("hydrogen", "B"): 99.40,
            ("hydrogen", "C"): 7.726,
            ("methane", "A"): 3.990,
            ("methane", "B"): 443.0,
            ("methane", "C"): -0.49,
        }

        self.pressure_sat_coeff_A = Param(
            self.component_list,
            mutable=True,
            initialize={
                c: v for (c, j), v in pressure_sat_coeff_data.items() if j == "A"
            },
            doc="Parameter A to compute saturated pressure",
            units=pyunits.dimensionless,
        )

        self.pressure_sat_coeff_B = Param(
            self.component_list,
            mutable=True,
            initialize={
                c: v for (c, j), v in pressure_sat_coeff_data.items() if j == "B"
            },
            doc="Parameter B to compute saturated pressure",
            units=pyunits.K,
        )

        self.pressure_sat_coeff_C = Param(
            self.component_list,
            mutable=True,
            initialize={
                c: v for (c, j), v in pressure_sat_coeff_data.items() if j == "C"
            },
            doc="Parameter C to compute saturated pressure",
            units=pyunits.K,
        )

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        dh_vap = {"benzene": 3.387e4, "toluene": 3.8262e4, "hydrogen": 0, "methane": 0}

        self.dh_vap = Param(
            self.component_list,
            mutable=True,
            units=pyunits.J / pyunits.mol,
            initialize=extract_data(dh_vap),
            doc="heat of vaporization",
        )

        # Set default scaling factors
        self.set_default_scaling("flow_mol", 1e3)
        self.set_default_scaling("flow_mol_phase_comp", 1e3)
        self.set_default_scaling("flow_mol_phase", 1e3)
        self.set_default_scaling("material_flow_terms", 1e3)
        self.set_default_scaling("enthalpy_flow_terms", 1e-2)
        self.set_default_scaling("mole_frac_comp", 1e1)
        self.set_default_scaling("temperature", 1e-2)
        self.set_default_scaling("temperature_dew", 1e-2)
        self.set_default_scaling("temperature_bubble", 1e-2)
        self.set_default_scaling("pressure", 1e-5)
        self.set_default_scaling("pressure_sat", 1e-5)
        self.set_default_scaling("pressure_dew", 1e-5)
        self.set_default_scaling("pressure_bubble", 1e-5)
        self.set_default_scaling("mole_frac_phase_comp", 1e1)
        self.set_default_scaling("enth_mol_phase", 1e-3, index="Liq")
        self.set_default_scaling("enth_mol_phase", 1e-4, index="Vap")
        self.set_default_scaling("enth_mol", 1e-3)
        self.set_default_scaling("entr_mol_phase", 1e-2)
        self.set_default_scaling("entr_mol", 1e-2)

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {
                "flow_mol": {"method": None},
                "flow_mol_phase_comp": {"method": None},
                "mole_frac_comp": {"method": None},
                "temperature": {"method": None},
                "pressure": {"method": None},
                "flow_mol_phase": {"method": None},
                "dens_mol_phase": {"method": "_dens_mol_phase"},
                "pressure_sat": {"method": "_pressure_sat"},
                "mole_frac_phase_comp": {"method": "_mole_frac_phase"},
                "energy_internal_mol_phase_comp": {
                    "method": "_energy_internal_mol_phase_comp"
                },
                "energy_internal_mol_phase": {"method": "_energy_internal_mol_phase"},
                "enth_mol_phase_comp": {"method": "_enth_mol_phase_comp"},
                "enth_mol_phase": {"method": "_enth_mol_phase"},
                "entr_mol_phase_comp": {"method": "_entr_mol_phase_comp"},
                "entr_mol_phase": {"method": "_entr_mol_phase"},
                "temperature_bubble": {"method": "_temperature_bubble"},
                "temperature_dew": {"method": "_temperature_dew"},
                "pressure_bubble": {"method": "_pressure_bubble"},
                "pressure_dew": {"method": "_pressure_dew"},
                "fug_phase_comp": {"method": "_fug_phase_comp"},
            }
        )

        obj.define_custom_properties(
            {
                # Enthalpy of vaporization
                "dh_vap": {"method": "_dh_vap", "units": obj.derived_units.ENERGY_MOLE},
                # Entropy of vaporization
                "ds_vap": {
                    "method": "_ds_vap",
                    "units": obj.derived_units.ENTROPY_MOLE,
                },
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


class _IdealStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    default_initializer = HDAInitializer

    def fix_initialization_states(blk):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """

        # Fix state variables
        fix_state_vars(blk)

        # Also need to deactivate sum of mole fraction constraint
        for k in blk.values():
            if not k.config.defined_state and k.config.has_phase_equilibrium:
                k.equilibrium_constraint.deactivate()


@declare_process_block_class("IdealStateBlock", block_class=_IdealStateBlock)
class IdealStateBlockData(StateBlockData):
    """An example property package for ideal VLE."""

    def build(self):
        """Callable method for Block construction."""
        super().build()

        # Add state variables
        self.flow_mol_phase_comp = Var(
            self._params.phase_list,
            self._params.component_list,
            initialize=0.5,
            units=pyunits.mol / pyunits.s,
            bounds=(1e-12, 100),
            doc="Phase-component molar flow rates",
        )

        self.pressure = Var(
            initialize=101325,
            bounds=(100000, 1000000),
            units=pyunits.Pa,
            domain=NonNegativeReals,
            doc="State pressure",
        )
        self.temperature = Var(
            initialize=298.15,
            units=pyunits.K,
            bounds=(298, 1000),
            domain=NonNegativeReals,
            doc="State temperature",
        )

        # Add supporting variables
        def flow_mol_phase(b, p):
            return sum(b.flow_mol_phase_comp[p, j] for j in b._params.component_list)

        self.flow_mol_phase = Expression(
            self._params.phase_list, rule=flow_mol_phase, doc="Phase molar flow rates"
        )

        def flow_mol(b):
            return sum(
                b.flow_mol_phase_comp[p, j]
                for j in b._params.component_list
                for p in b._params.phase_list
            )

        self.flow_mol = Expression(rule=flow_mol, doc="Total molar flowrate")

        def mole_frac_phase_comp(b, p, j):
            return b.flow_mol_phase_comp[p, j] / b.flow_mol_phase[p]

        self.mole_frac_phase_comp = Expression(
            self._params.phase_list,
            self._params.component_list,
            rule=mole_frac_phase_comp,
            doc="Phase mole fractions",
        )

        def mole_frac_comp(b, j):
            return (
                sum(b.flow_mol_phase_comp[p, j] for p in b._params.phase_list)
                / b.flow_mol
            )

        self.mole_frac_comp = Expression(
            self._params.component_list,
            rule=mole_frac_comp,
            doc="Mixture mole fractions",
        )

        # Reaction Stoichiometry
        add_object_reference(
            self, "phase_equilibrium_list_ref", self._params.phase_equilibrium_list
        )

        if self.config.has_phase_equilibrium and self.config.defined_state is False:
            # Definition of equilibrium temperature for smooth VLE
            self._teq = Var(
                initialize=self.temperature.value,
                units=pyunits.K,
                doc="Temperature for calculating phase equilibrium",
            )
            self._t1 = Var(
                initialize=self.temperature.value,
                units=pyunits.K,
                doc="Intermediate temperature for calculating Teq",
            )

            self.eps_1 = Param(
                default=0.01,
                units=pyunits.K,
                mutable=True,
                doc="Smoothing parameter for Teq",
            )
            self.eps_2 = Param(
                default=0.0005,
                units=pyunits.K,
                mutable=True,
                doc="Smoothing parameter for Teq",
            )

            # PSE paper Eqn 13
            def rule_t1(b):
                return b._t1 == 0.5 * (
                    b.temperature
                    + b.temperature_bubble
                    + sqrt((b.temperature - b.temperature_bubble) ** 2 + b.eps_1**2)
                )

            self._t1_constraint = Constraint(rule=rule_t1)

            # PSE paper Eqn 14
            # TODO : Add option for supercritical extension
            def rule_teq(b):
                return b._teq == 0.5 * (
                    b._t1
                    + b.temperature_dew
                    - sqrt((b._t1 - b.temperature_dew) ** 2 + b.eps_2**2)
                )

            self._teq_constraint = Constraint(rule=rule_teq)

            def rule_tr_eq(b, i):
                return b._teq / b._params.temperature_crit[i]

            self._tr_eq = Expression(
                self._params.component_list,
                rule=rule_tr_eq,
                doc="Component reduced temperatures",
            )

            def rule_equilibrium(b, i):
                return b.fug_phase_comp["Liq", i] == b.fug_phase_comp["Vap", i]

            self.equilibrium_constraint = Constraint(
                self._params.component_list, rule=rule_equilibrium
            )

    # -----------------------------------------------------------------------------
    # Property Methods
    def _dens_mol_phase(self):
        self.dens_mol_phase = Var(
            self._params.phase_list,
            initialize=1.0,
            units=pyunits.mol * pyunits.m**-3,
            doc="Molar density",
        )

        def rule_dens_mol_phase(b, p):
            if p == "Vap":
                return b._dens_mol_vap()
            else:
                return b._dens_mol_liq()

        self.eq_dens_mol_phase = Constraint(
            self._params.phase_list, rule=rule_dens_mol_phase
        )

    def _energy_internal_mol_phase_comp(self):
        self.energy_internal_mol_phase_comp = Var(
            self._params.phase_list,
            self._params.component_list,
            units=pyunits.J / pyunits.mol,
            doc="Phase-component molar specific internal energies",
        )

        def rule_energy_internal_mol_phase_comp(b, p, j):
            if p == "Vap":
                return b.energy_internal_mol_phase_comp[p, j] == b.enth_mol_phase_comp[
                    p, j
                ] - const.gas_constant * (b.temperature - b._params.temperature_ref)
            else:
                return (
                    b.energy_internal_mol_phase_comp[p, j]
                    == b.enth_mol_phase_comp[p, j]
                )

        self.eq_energy_internal_mol_phase_comp = Constraint(
            self._params.phase_list,
            self._params.component_list,
            rule=rule_energy_internal_mol_phase_comp,
        )

    def _energy_internal_mol_phase(self):
        self.energy_internal_mol_phase = Var(
            self._params.phase_list,
            units=pyunits.J / pyunits.mol,
            doc="Phase molar specific internal energies",
        )

        def rule_energy_internal_mol_phase(b, p):
            return b.energy_internal_mol_phase[p] == sum(
                b.energy_internal_mol_phase_comp[p, i] * b.mole_frac_phase_comp[p, i]
                for i in b._params.component_list
            )

        self.eq_energy_internal_mol_phase = Constraint(
            self._params.phase_list, rule=rule_energy_internal_mol_phase
        )

    def _enth_mol_phase_comp(self):
        self.enth_mol_phase_comp = Var(
            self._params.phase_list,
            self._params.component_list,
            initialize=7e5,
            units=pyunits.J / pyunits.mol,
            doc="Phase-component molar specific enthalpies",
        )

        def rule_enth_mol_phase_comp(b, p, j):
            if p == "Vap":
                return b._enth_mol_comp_vap(j)
            else:
                return b._enth_mol_comp_liq(j)

        self.eq_enth_mol_phase_comp = Constraint(
            self._params.phase_list,
            self._params.component_list,
            rule=rule_enth_mol_phase_comp,
        )

    def _enth_mol_phase(self):
        self.enth_mol_phase = Var(
            self._params.phase_list,
            initialize=7e5,
            units=pyunits.J / pyunits.mol,
            doc="Phase molar specific enthalpies",
        )

        def rule_enth_mol_phase(b, p):
            return b.enth_mol_phase[p] == sum(
                b.enth_mol_phase_comp[p, i] * b.mole_frac_phase_comp[p, i]
                for i in b._params.component_list
            )

        self.eq_enth_mol_phase = Constraint(
            self._params.phase_list, rule=rule_enth_mol_phase
        )

    def _entr_mol_phase_comp(self):
        self.entr_mol_phase_comp = Var(
            self._params.phase_list,
            self._params.component_list,
            units=pyunits.J / pyunits.mol / pyunits.K,
            doc="Phase-component molar specific entropies",
        )

        def rule_entr_mol_phase_comp(b, p, j):
            if p == "Vap":
                return b._entr_mol_comp_vap(j)
            else:
                return b._entr_mol_comp_liq(j)

        self.eq_entr_mol_phase_comp = Constraint(
            self._params.phase_list,
            self._params.component_list,
            rule=rule_entr_mol_phase_comp,
        )

    def _entr_mol_phase(self):
        self.entr_mol_phase = Var(
            self._params.phase_list,
            units=pyunits.J / pyunits.mol / pyunits.K,
            doc="Phase molar specific enthropies",
        )

        def rule_entr_mol_phase(b, p):
            return b.entr_mol_phase[p] == sum(
                b.entr_mol_phase_comp[p, i] * b.mole_frac_phase_comp[p, i]
                for i in b._params.component_list
            )

        self.eq_entr_mol_phase = Constraint(
            self._params.phase_list, rule=rule_entr_mol_phase
        )

    # -----------------------------------------------------------------------------
    # General Methods
    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        if not self.is_property_constructed("material_flow_terms"):
            try:

                def rule_material_flow_terms(blk, p, j):
                    return blk.flow_mol_phase_comp[p, j]

                self.material_flow_terms = Expression(
                    self.params.phase_list,
                    self.params.component_list,
                    rule=rule_material_flow_terms,
                )
            except AttributeError:
                self.del_component(self.material_flow_terms)

        if j in self.params.component_list:
            return self.material_flow_terms[p, j]
        else:
            return 0

    def get_enthalpy_flow_terms(self, p):
        """Create enthalpy flow terms."""
        if not self.is_property_constructed("enthalpy_flow_terms"):
            try:

                def rule_enthalpy_flow_terms(blk, p):
                    return blk.flow_mol_phase[p] * blk.enth_mol_phase[p]

                self.enthalpy_flow_terms = Expression(
                    self.params.phase_list, rule=rule_enthalpy_flow_terms
                )
            except AttributeError:
                self.del_component(self.enthalpy_flow_terms)
        return self.enthalpy_flow_terms[p]

    def get_material_density_terms(self, p, j):
        """Create material density terms."""
        if not self.is_property_constructed("material_density_terms"):
            try:

                def rule_material_density_terms(b, p, j):
                    return self.dens_mol_phase[p] * self.mole_frac_phase_comp[p, j]

                self.material_density_terms = Expression(
                    self.params.phase_list,
                    self.params.component_list,
                    rule=rule_material_density_terms,
                )
            except AttributeError:
                self.del_component(self.material_density_terms)

        if j in self.params.component_list:
            return self.material_density_terms[p, j]
        else:
            return 0

    def get_enthalpy_density_terms(self, p):
        """Create energy density terms."""
        if not self.is_property_constructed("enthalpy_density_terms"):
            try:

                def rule_energy_density_terms(b, p):
                    return self.dens_mol_phase[p] * self.energy_internal_mol_phase[p]

                self.energy_density_terms = Expression(
                    self.params.phase_list, rule=rule_energy_density_terms
                )
            except AttributeError:
                self.del_component(self.energy_density_terms)
        return self.enthalpy_density_terms[p]

    def default_material_balance_type(self):
        return MaterialBalanceType.componentPhase

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_basis(b):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        """Define state vars."""
        return {
            "flow_mol_phase_comp": self.flow_mol_phase_comp,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }

    # Property package utility functions
    def calculate_bubble_point_temperature(self, clear_components=True):
        """ "To compute the bubble point temperature of the mixture."""

        if hasattr(self, "eq_temperature_bubble"):
            # Do not delete components if the block already has the components
            clear_components = False

        calculate_variable_from_constraint(
            self.temperature_bubble, self.eq_temperature_bubble
        )

        return self.temperature_bubble.value

        if clear_components is True:
            self.del_component(self.eq_temperature_bubble)
            self.del_component(self._p_sat_bubbleT)
            self.del_component(self.temperature_bubble)

    def calculate_dew_point_temperature(self, clear_components=True):
        """ "To compute the dew point temperature of the mixture."""

        if hasattr(self, "eq_temperature_dew"):
            # Do not delete components if the block already has the components
            clear_components = False

        calculate_variable_from_constraint(
            self.temperature_dew, self.eq_temperature_dew
        )

        return self.temperature_dew.value

        # Delete the var/constraint created in this method that are part of the
        # IdealStateBlock if the user desires
        if clear_components is True:
            self.del_component(self.eq_temperature_dew)
            self.del_component(self._p_sat_dewT)
            self.del_component(self.temperature_dew)

    def calculate_bubble_point_pressure(self, clear_components=True):
        """ "To compute the bubble point pressure of the mixture."""

        if hasattr(self, "eq_pressure_bubble"):
            # Do not delete components if the block already has the components
            clear_components = False

        calculate_variable_from_constraint(
            self.pressure_bubble, self.eq_pressure_bubble
        )

        return self.pressure_bubble.value

        # Delete the var/constraint created in this method that are part of the
        # IdealStateBlock if the user desires
        if clear_components is True:
            self.del_component(self.eq_pressure_bubble)
            self.del_component(self._p_sat_bubbleP)
            self.del_component(self.pressure_bubble)

    def calculate_dew_point_pressure(self, clear_components=True):
        """ "To compute the dew point pressure of the mixture."""

        if hasattr(self, "eq_pressure_dew"):
            # Do not delete components if the block already has the components
            clear_components = False

        calculate_variable_from_constraint(self.pressure_dew, self.eq_pressure_dew)

        return self.pressure_dew.value

        # Delete the var/constraint created in this method that are part of the
        # IdealStateBlock if the user desires
        if clear_components is True:
            self.del_component(self.eq_pressure_dew)
            self.del_component(self._p_sat_dewP)
            self.del_component(self.pressure_dew)

    # -----------------------------------------------------------------------------
    # Bubble and Dew Points
    # Ideal-Ideal properties allow for the simplifications below
    # Other methods require more complex equations with shadow compositions

    # For future work, propose the following:
    # Core class writes a set of constraints Phi_L_i == Phi_V_i
    # Phi_L_i and Phi_V_i make calls to submethods which add shadow compositions
    # as needed
    def _temperature_bubble(self):
        self.temperature_bubble = Param(
            initialize=33.0, units=pyunits.K, doc="Bubble point temperature"
        )

    def _temperature_dew(self):

        self.temperature_dew = Var(
            initialize=298.15, units=pyunits.K, doc="Dew point temperature"
        )

        def rule_psat_dew(b, j):
            return (
                1e5
                * pyunits.Pa
                * 10
                ** (
                    b._params.pressure_sat_coeff_A[j]
                    - b._params.pressure_sat_coeff_B[j]
                    / (b.temperature_dew + b._params.pressure_sat_coeff_C[j])
                )
            )

        try:
            # Try to build expression
            self._p_sat_dewT = Expression(
                self._params.component_list, rule=rule_psat_dew
            )

            def rule_temp_dew(b):
                return (
                    b.pressure
                    * sum(
                        b.mole_frac_comp[i] / b._p_sat_dewT[i]
                        for i in ["toluene", "benzene"]
                    )
                    - 1
                    == 0
                )

            self.eq_temperature_dew = Constraint(rule=rule_temp_dew)
        except AttributeError:
            # If expression fails, clean up so that DAE can try again later
            # Deleting only var/expression as expression construction will fail
            # first; if it passes then constraint construction will not fail.
            self.del_component(self.temperature_dew)
            self.del_component(self._p_sat_dewT)

    def _pressure_bubble(self):
        self.pressure_bubble = Param(
            initialize=1e8, units=pyunits.Pa, doc="Bubble point pressure"
        )

    def _pressure_dew(self):
        self.pressure_dew = Var(
            initialize=298.15, units=pyunits.Pa, doc="Dew point pressure"
        )

        def rule_psat_dew(b, j):
            return (
                1e5
                * pyunits.Pa
                * 10
                ** (
                    b._params.pressure_sat_coeff_A[j]
                    - b._params.pressure_sat_coeff_B[j]
                    / (b.temperature + b._params.pressure_sat_coeff_C[j])
                )
            )

        try:
            # Try to build expression
            self._p_sat_dewP = Expression(
                self._params.component_list, rule=rule_psat_dew
            )

            def rule_pressure_dew(b):
                return (
                    b.pressure_dew
                    * sum(
                        b.mole_frac_comp[i] / b._p_sat_dewP[i]
                        for i in ["toluene", "benzene"]
                    )
                    - 1
                    == 0
                )

            self.eq_pressure_dew = Constraint(rule=rule_pressure_dew)
        except AttributeError:
            # If expression fails, clean up so that DAE can try again later
            # Deleting only var/expression as expression construction will fail
            # first; if it passes then constraint construction will not fail.
            self.del_component(self.pressure_dew)
            self.del_component(self._p_sat_dewP)

    # -----------------------------------------------------------------------------
    # Liquid phase properties
    def _dens_mol_liq(b):
        return b.dens_mol_phase["Liq"] == 1e3 * sum(
            b.mole_frac_phase_comp["Liq", j]
            * b._params.dens_liq_param_1[j]
            / b._params.dens_liq_param_2[j]
            ** (
                1
                + (1 - b.temperature / b._params.dens_liq_param_3[j])
                ** b._params.dens_liq_param_4[j]
            )
            for j in ["benzene", "toluene"]
        )

    def _fug_phase_comp(self):
        def fug_phase_comp_rule(b, p, i):
            if p == "Liq":
                if i in ["hydrogen", "methane"]:
                    return b.mole_frac_phase_comp["Liq", i]
                else:
                    return b.pressure_sat[i] * b.mole_frac_phase_comp["Liq", i]
            else:
                if i in ["hydrogen", "methane"]:
                    return 1e-6
                else:
                    return b.mole_frac_phase_comp["Vap", i] * b.pressure

        self.fug_phase_comp = Expression(
            self._params.phase_list,
            self._params.component_list,
            rule=fug_phase_comp_rule,
        )

    def _pressure_sat(self):
        self.pressure_sat = Var(
            self._params.component_list,
            initialize=101325,
            units=pyunits.Pa,
            doc="Vapor pressure",
        )

        def rule_P_sat(b, j):
            return (
                (
                    log10(b.pressure_sat[j] / pyunits.Pa * 1e-5)
                    - b._params.pressure_sat_coeff_A[j]
                )
                * (b._teq + b._params.pressure_sat_coeff_C[j])
            ) == -b._params.pressure_sat_coeff_B[j]

        self.eq_pressure_sat = Constraint(self._params.component_list, rule=rule_P_sat)

    def _enth_mol_comp_liq(b, j):
        return b.enth_mol_phase_comp["Liq", j] * 1e3 == (
            (b._params.cp_ig_5["Liq", j] / 5)
            * (b.temperature**5 - b._params.temperature_ref**5)
            + (b._params.cp_ig_4["Liq", j] / 4)
            * (b.temperature**4 - b._params.temperature_ref**4)
            + (b._params.cp_ig_3["Liq", j] / 3)
            * (b.temperature**3 - b._params.temperature_ref**3)
            + (b._params.cp_ig_2["Liq", j] / 2)
            * (b.temperature**2 - b._params.temperature_ref**2)
            + b._params.cp_ig_1["Liq", j] * (b.temperature - b._params.temperature_ref)
        )

    def _entr_mol_comp_liq(b, j):
        return b.entr_mol_phase_comp["Liq", j] * 1e3 == (
            (
                (b._params.cp_ig_5["Liq", j] / 4)
                * (b.temperature**4 - b._params.temperature_ref**4)
                + (b._params.cp_ig_4["Liq", j] / 3)
                * (b.temperature**3 - b._params.temperature_ref**3)
                + (b._params.cp_ig_3["Liq", j] / 2)
                * (b.temperature**2 - b._params.temperature_ref**2)
                + b._params.cp_ig_2["Liq", j]
                * (b.temperature - b._params.temperature_ref)
                + b._params.cp_ig_1["Liq", j]
                * log(b.temperature / b._params.temperature_ref)
            )
            - const.gas_constant
            * log(
                b.mole_frac_phase_comp["Liq", j] * b.pressure / b._params.pressure_ref
            )
        )

    # -----------------------------------------------------------------------------
    # Vapour phase properties
    def _dens_mol_vap(b):
        return b.pressure == (
            b.dens_mol_phase["Vap"] * const.gas_constant * b.temperature
        )

    def _dh_vap(self):
        # heat of vaporization
        add_object_reference(self, "dh_vap", self._params.dh_vap)

    def _ds_vap(self):
        # entropy of vaporization = dh_Vap/T_boil
        # TODO : something more rigorous would be nice
        self.ds_vap = Var(
            self._params.component_list,
            initialize=86,
            units=pyunits.J / pyunits.mol / pyunits.K,
            doc="Entropy of vaporization",
        )

        def rule_ds_vap(b, j):
            return b.dh_vap[j] == (b.ds_vap[j] * b._params.temperature_boil[j])

        self.eq_ds_vap = Constraint(self._params.component_list, rule=rule_ds_vap)

    def _enth_mol_comp_vap(b, j):
        return b.enth_mol_phase_comp["Vap", j] == b.dh_vap[j] + (
            (b._params.cp_ig_5["Vap", j] / 5)
            * (b.temperature**5 - b._params.temperature_ref**5)
            + (b._params.cp_ig_4["Vap", j] / 4)
            * (b.temperature**4 - b._params.temperature_ref**4)
            + (b._params.cp_ig_3["Vap", j] / 3)
            * (b.temperature**3 - b._params.temperature_ref**3)
            + (b._params.cp_ig_2["Vap", j] / 2)
            * (b.temperature**2 - b._params.temperature_ref**2)
            + b._params.cp_ig_1["Vap", j] * (b.temperature - b._params.temperature_ref)
        )

    def _entr_mol_comp_vap(b, j):
        return b.entr_mol_phase_comp["Vap", j] == (
            b.ds_vap[j]
            + (
                (b._params.cp_ig_5["Vap", j] / 4)
                * (b.temperature**4 - b._params.temperature_ref**4)
                + (b._params.cp_ig_4["Vap", j] / 3)
                * (b.temperature**3 - b._params.temperature_ref**3)
                + (b._params.cp_ig_3["Vap", j] / 2)
                * (b.temperature**2 - b._params.temperature_ref**2)
                + b._params.cp_ig_2["Vap", j]
                * (b.temperature - b._params.temperature_ref)
                + b._params.cp_ig_1["Vap", j]
                * log(b.temperature / b._params.temperature_ref)
            )
            - const.gas_constant
            * log(
                b.mole_frac_phase_comp["Vap", j] * b.pressure / b._params.pressure_ref
            )
        )

    def calculate_scaling_factors(self):
        # Get default scale factors
        super().calculate_scaling_factors()

        is_two_phase = len(self._params.phase_list) == 2
        sf_flow = iscale.get_scaling_factor(self.flow_mol, default=1, warning=True)
        sf_T = iscale.get_scaling_factor(self.temperature, default=1, warning=True)
        sf_P = iscale.get_scaling_factor(self.pressure, default=1, warning=True)

        if self.is_property_constructed("_teq"):
            iscale.set_scaling_factor(self._teq, sf_T)
        if self.is_property_constructed("_teq_constraint"):
            iscale.constraint_scaling_transform(
                self._teq_constraint, sf_T, overwrite=False
            )

        if self.is_property_constructed("_t1"):
            iscale.set_scaling_factor(self._t1, sf_T)
        if self.is_property_constructed("_t1_constraint"):
            iscale.constraint_scaling_transform(
                self._t1_constraint, sf_T, overwrite=False
            )

        if self.is_property_constructed("_mole_frac_pdew"):
            iscale.set_scaling_factor(self._mole_frac_pdew, 1e3)
            iscale.constraint_scaling_transform(
                self._sum_mole_frac_pdew, 1e3, overwrite=False
            )

        if self.is_property_constructed("total_flow_balance"):
            iscale.constraint_scaling_transform(
                self.total_flow_balance, sf_flow, overwrite=False
            )

        if self.is_property_constructed("component_flow_balances"):
            for i, c in self.component_flow_balances.items():
                if is_two_phase:
                    s = iscale.get_scaling_factor(
                        self.mole_frac_comp[i], default=1, warning=True
                    )
                    s *= sf_flow
                    iscale.constraint_scaling_transform(c, s, overwrite=False)
                else:
                    s = iscale.get_scaling_factor(
                        self.mole_frac_comp[i], default=1, warning=True
                    )
                    iscale.constraint_scaling_transform(c, s, overwrite=False)

        if self.is_property_constructed("dens_mol_phase"):
            for c in self.eq_dens_mol_phase.values():
                iscale.constraint_scaling_transform(c, sf_P, overwrite=False)

        if self.is_property_constructed("dens_mass_phase"):
            for p, c in self.eq_dens_mass_phase.items():
                sf = iscale.get_scaling_factor(
                    self.dens_mass_phase[p], default=1, warning=True
                )
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if self.is_property_constructed("enth_mol_phase"):
            for p, c in self.eq_enth_mol_phase.items():
                sf = iscale.get_scaling_factor(
                    self.enth_mol_phase[p], default=1, warning=True
                )
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if self.is_property_constructed("enth_mol"):
            sf = iscale.get_scaling_factor(self.enth_mol, default=1, warning=True)
            sf *= sf_flow
            iscale.constraint_scaling_transform(self.eq_enth_mol, sf, overwrite=False)

        if self.is_property_constructed("entr_mol_phase"):
            for p, c in self.eq_entr_mol_phase.items():
                sf = iscale.get_scaling_factor(
                    self.entr_mol_phase[p], default=1, warning=True
                )
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if self.is_property_constructed("entr_mol"):
            sf = iscale.get_scaling_factor(self.entr_mol, default=1, warning=True)
            sf *= sf_flow
            iscale.constraint_scaling_transform(self.eq_entr_mol, sf, overwrite=False)

        if self.is_property_constructed("gibbs_mol_phase"):
            for p, c in self.eq_gibbs_mol_phase.items():
                sf = iscale.get_scaling_factor(
                    self.gibbs_mol_phase[p], default=1, warning=True
                )
                iscale.constraint_scaling_transform(c, sf, overwrite=False)
