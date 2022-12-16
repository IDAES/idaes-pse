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
General Cubic Equation of State property package with VLE calucations.
Cubic formulation and pure component property correlations from:
"The Properties of Gases and Liquids, 4th Edition", Reid, Prausnitz and Poling,
McGraw-Hill, 1987

Smooth Vapor-Liquid Equilibrium formulation from:
"A Smooth, Square Flash Formulation for Equation-Oriented Flowsheet
Optimization", Burgard et al., Proceedings of the 13 the International
Symposium on Process Systems Engineering – PSE 2018, July 1-5, 2018, San Diego

All results have been cross-referenced against other sources.
"""
# Import Python libraries
import math

# Import Pyomo libraries
from pyomo.environ import (
    check_optimal_termination,
    Constraint,
    exp,
    Expression,
    log,
    NonNegativeReals,
    sqrt,
    Param,
    PositiveReals,
    value,
    Var,
    units as pyunits,
)
from pyomo.common.config import ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    LiquidPhase,
    VaporPhase,
)
from idaes.core.util.initialization import (
    solve_indexed_blocks,
    fix_state_vars,
    revert_state_vars,
)
from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    InitializationError,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_activated_equalities,
)
from idaes.core.util.math import safe_log
from idaes.core.solvers import get_solver
from idaes.core.util.constants import Constants as const
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.models.properties.modular_properties.eos.ceos_common import (
    EoS_param,
    cubic_roots_available,
    CubicThermoExpressions,
    CubicType as CubicEoS,
)

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("CubicParameterBlock")
class CubicParameterData(PhysicalParameterBlock):
    """
    General Property Parameter Block Class
    """

    # Config block for the _IdealStateBlock
    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare(
        "valid_phase",
        ConfigValue(
            default=("Vap", "Liq"),
            domain=In(["Liq", "Vap", ("Vap", "Liq"), ("Liq", "Vap")]),
            description="Flag indicating the valid phase",
            doc="""Flag indicating the valid phase for a given set of
conditions, and thus corresponding constraints  should be included,
**default** - ('Vap', 'Liq').
**Valid values:** {
**'Liq'** - Liquid only,
**'Vap'** - Vapor only,
**('Vap', 'Liq')** - Vapor-liquid equilibrium,
**('Liq', 'Vap')** - Vapor-liquid equilibrium,}""",
        ),
    )

    def build(self):
        """
        Callable method for Block construction.
        """
        super(CubicParameterData, self).build()

        self._state_block_class = CubicStateBlock

        # Create Phase objects
        if (
            self.config.valid_phase == ("Liq", "Vap")
            or self.config.valid_phase == ("Vap", "Liq")
            or self.config.valid_phase == "Liq"
        ):
            self.Liq = LiquidPhase()

        if (
            self.config.valid_phase == ("Liq", "Vap")
            or self.config.valid_phase == ("Vap", "Liq")
            or self.config.valid_phase == "Vap"
        ):
            self.Vap = VaporPhase()

        self.set_default_scaling("flow_mol", 1e-3)
        self.set_default_scaling("mole_frac_comp", 10)
        self.set_default_scaling("temperature", 1e-1)
        self.set_default_scaling("temperature_dew", 1e-1)
        self.set_default_scaling("temperature_bubble", 1e-1)
        self.set_default_scaling("flow_mol_phase", 1e-2)
        self.set_default_scaling("dens_mol_phase", 1)
        self.set_default_scaling("dens_mass_phase", 1e-1)
        self.set_default_scaling("pressure", 1e-6)
        self.set_default_scaling("pressure_sat", 1e-6)
        self.set_default_scaling("pressure_dew", 1e-6)
        self.set_default_scaling("pressure_bubble", 1e-6)
        self.set_default_scaling("mole_frac_phase_comp", 10)
        self.set_default_scaling("enth_mol_phase", 1e-3, index="Liq")
        self.set_default_scaling("enth_mol_phase", 1e-4, index="Vap")
        self.set_default_scaling("enth_mol", 1e-3)
        self.set_default_scaling("entr_mol_phase", 1e-2)
        self.set_default_scaling("entr_mol", 1e-2)
        self.set_default_scaling("mw", 100)
        self.set_default_scaling("mw_phase", 100)
        self.set_default_scaling("fug_phase_comp", 1)
        self.set_default_scaling("fug_coeff_phase_comp", 1)
        self.set_default_scaling("gibbs_mol_phase", 1e-3)

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {
                "flow_mol": {"method": None, "units": "mol/s"},
                "mole_frac_comp": {"method": None, "units": "none"},
                "temperature": {"method": None, "units": "K"},
                "pressure": {"method": None, "units": "Pa"},
                "flow_mol_phase": {"method": None, "units": "mol/s"},
                "dens_mol_phase": {"method": "_dens_mol_phase", "units": "mol/m^3"},
                "dens_mass_phase": {"method": "_dens_mass_phase", "units": "kg/m^3"},
                "pressure_sat": {"method": "_pressure_sat", "units": "Pa"},
                "mole_frac_phase_comp": {
                    "method": "_mole_frac_phase",
                    "units": "no unit",
                },
                "energy_internal_mol_phase": {
                    "method": "_energy_internal_mol_phase",
                    "units": "J/mol",
                },
                "energy_internal_mol": {
                    "method": "_energy_internal_mol",
                    "units": "J/mol",
                },
                "enth_mol_phase": {"method": "_enth_mol_phase", "units": "J/mol"},
                "enth_mol": {"method": "_enth_mol", "units": "J/mol"},
                "entr_mol_phase": {"method": "_entr_mol_phase", "units": "J/mol"},
                "entr_mol": {"method": "_entr_mol", "units": "J/mol.K"},
                "fug_phase_comp": {"method": "_fug_phase", "units": "Pa"},
                "fug_coeff_phase_comp": {"method": "_fug_coeff_phase", "units": "-"},
                "gibbs_mol_phase": {"method": "_gibbs_mol_phase", "units": "J/mol"},
                "mw": {"method": "_mw", "units": "kg/mol"},
                "mw_phase": {"method": "_mw_phase", "units": "kg/mol"},
                "temperature_bubble": {"method": "_temperature_bubble", "units": "K"},
                "temperature_dew": {"method": "_temperature_dew", "units": "K"},
                "pressure_bubble": {"method": "_pressure_bubble", "units": "Pa"},
                "pressure_dew": {"method": "_pressure_dew", "units": "Pa"},
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


class _CubicStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(
        blk,
        state_args=None,
        state_vars_fixed=False,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Initialization routine for property package.
        Keyword Arguments:
            state_args : Dictionary with initial guesses for the state vars
                         chosen. Note that if this method is triggered
                         through the control volume, and if initial guesses
                         were not provided at the unit model level, the
                         control volume passes the inlet values as initial
                         guess. Expected keys in state_args dict are:
                         * flow_mol
                         * mole_frac_comp (dict with components as keys)
                         * pressure
                         * temperature
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            state_vars_fixed: Flag to denote if state vars have already been
                              fixed.
                              - True - states have already been fixed and
                                       initialization does not need to worry
                                       about fixing and unfixing variables.
                             - False - states have not been fixed. The state
                                       block will deal with fixing/unfixing.
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)
            hold_state : flag indicating whether the initialization routine
                         should unfix any state variables fixed during
                         initialization (default=False).
                         - True - states varaibles are not unfixed, and
                                 a dict of returned containing flags for
                                 which states were fixed during
                                 initialization.
                        - False - state variables are unfixed after
                                 initialization by calling the
                                 relase_state method
        Returns:
            If hold_states is True, returns a dict containing flags for
            which states were fixed during initialization.
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="properties")

        init_log.info("Starting initialization")

        # Deactivate the constraints specific for outlet block i.e.
        # when defined state is False
        for k in blk.keys():
            if blk[k].config.defined_state is False:
                blk[k].sum_mole_frac_out.deactivate()

        # Fix state variables if not already fixed
        if state_vars_fixed is False:
            flags = fix_state_vars(blk, state_args)

        else:
            # Check when the state vars are fixed already result in dof 0
            for k in blk.keys():
                if degrees_of_freedom(blk[k]) != 0:
                    raise Exception(
                        "State vars fixed but degrees of freedom "
                        "for state block is not zero during "
                        "initialization."
                    )

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # If present, initialize bubble and dew point calculations
        # Antoine equation
        def antoine_P(b, j, T):
            return pyunits.convert_value(
                value(
                    10
                    ** (
                        b.params.antoine_coeff_A[j]
                        - b.params.antoine_coeff_B[j]
                        / (T + b.params.antoine_coeff_C[j])
                    )
                ),
                from_units=pyunits.bar,
                to_units=pyunits.Pa,
            )

        # Bubble temperature initialization
        for k in blk.keys():
            if hasattr(blk[k], "_mole_frac_tbub"):
                Tbub0 = 0
                for j in blk[k].params.component_list:
                    Tbub0 += value(
                        blk[k].mole_frac_comp[j]
                        * (
                            blk[k].params.antoine_coeff_B[j]
                            / (
                                blk[k].params.antoine_coeff_A[j]
                                - math.log10(
                                    value(
                                        pyunits.convert(
                                            blk[k].pressure, to_units=pyunits.bar
                                        )
                                    )
                                )
                            )
                            - blk[k].params.antoine_coeff_C[j]
                        )
                    )

                err = 1
                counter = 0

                while err > 1e-2 and counter < 100:
                    f = value(
                        sum(
                            antoine_P(blk[k], j, Tbub0) * blk[k].mole_frac_comp[j]
                            for j in blk[k].params.component_list
                        )
                        - blk[k].pressure
                    )
                    df = value(
                        sum(
                            blk[k].mole_frac_comp[j]
                            * blk[k].params.antoine_coeff_B[j]
                            * math.log(10)
                            * antoine_P(blk[k], j, Tbub0)
                            / (Tbub0 + blk[k].params.antoine_coeff_C[j]) ** 2
                            for j in blk[k].params.component_list
                        )
                    )

                    if f / df > 20:
                        Tbub1 = Tbub0 - 20
                    elif f / df < -20:
                        Tbub1 = Tbub0 + 20
                    else:
                        Tbub1 = Tbub0 - f / df

                    err = abs(Tbub1 - Tbub0)
                    Tbub0 = Tbub1
                    counter += 1

                blk[k].temperature_bubble.value = Tbub0

                for j in blk[k].params.component_list:
                    blk[k]._mole_frac_tbub[j].value = value(
                        blk[k].mole_frac_comp[j]
                        * antoine_P(blk[k], j, Tbub0)
                        / blk[k].pressure
                    )

        # Dew temperature initialization
        for k in blk.keys():
            if hasattr(blk[k], "_mole_frac_tdew"):
                Tdew0 = 0
                for j in blk[k].params.component_list:
                    Tdew0 += value(
                        blk[k].mole_frac_comp[j]
                        * (
                            blk[k].params.antoine_coeff_B[j]
                            / (
                                blk[k].params.antoine_coeff_A[j]
                                - math.log10(
                                    value(
                                        pyunits.convert(
                                            blk[k].pressure, to_units=pyunits.bar
                                        )
                                    )
                                )
                            )
                            - blk[k].params.antoine_coeff_C[j]
                        )
                    )

                err = 1
                counter = 0

                while err > 1e-2 and counter < 100:
                    f = value(
                        blk[k].pressure
                        * sum(
                            blk[k].mole_frac_comp[j] / antoine_P(blk[k], j, Tdew0)
                            for j in blk[k].params.component_list
                        )
                        - 1
                    )
                    df = -value(
                        blk[k].pressure
                        * math.log(10)
                        * sum(
                            blk[k].mole_frac_comp[j]
                            * blk[k].params.antoine_coeff_B[j]
                            / (
                                (Tdew0 + blk[k].params.antoine_coeff_C[j]) ** 2
                                * antoine_P(blk[k], j, Tdew0)
                            )
                            for j in blk[k].params.component_list
                        )
                    )

                    if f / df > 20:
                        Tdew1 = Tdew0 - 20
                    elif f / df < -20:
                        Tdew1 = Tdew0 + 20
                    else:
                        Tdew1 = Tdew0 - f / df

                    err = abs(Tdew1 - Tdew0)
                    Tdew0 = Tdew1
                    counter += 1

                blk[k].temperature_dew.value = Tdew0

                for j in blk[k].params.component_list:
                    blk[k]._mole_frac_tdew[j].value = value(
                        blk[k].mole_frac_comp[j]
                        * blk[k].pressure
                        / antoine_P(blk[k], j, Tdew0)
                    )

        # Bubble pressure initialization
        for k in blk.keys():
            if hasattr(blk[k], "_mole_frac_pbub"):
                blk[k].pressure_bubble.value = value(
                    sum(
                        blk[k].mole_frac_comp[j]
                        * antoine_P(blk[k], j, blk[k].temperature)
                        for j in blk[k].params.component_list
                    )
                )

                for j in blk[k].params.component_list:
                    blk[k]._mole_frac_pbub[j].value = value(
                        blk[k].mole_frac_comp[j]
                        * antoine_P(blk[k], j, blk[k].temperature)
                        / blk[k].pressure_bubble
                    )

                blk[k].pressure_bubble.display()
                blk[k]._mole_frac_pbub.display()

        # Dew pressure initialization
        for k in blk.keys():
            if hasattr(blk[k], "_mole_frac_pdew"):
                blk[k].pressure_dew.value = value(
                    sum(
                        1
                        / (
                            blk[k].mole_frac_comp[j]
                            / antoine_P(blk[k], j, blk[k].temperature)
                        )
                        for j in blk[k].params.component_list
                    )
                )

                for j in blk[k].params.component_list:
                    blk[k]._mole_frac_pdew[j].value = value(
                        blk[k].mole_frac_comp[j]
                        * blk[k].pressure_bubble
                        / antoine_P(blk[k], j, blk[k].temperature)
                    )

        # Solve bubble and dew point constraints
        cons_count = 0
        for k in blk.keys():
            for c in blk[k].component_objects(Constraint):
                # Deactivate all property constraints
                if c.local_name not in (
                    "eq_pressure_dew",
                    "eq_pressure_bubble",
                    "eq_temperature_dew",
                    "eq_temperature_bubble",
                    "_sum_mole_frac_tbub",
                    "_sum_mole_frac_tdew",
                    "_sum_mole_frac_pbub",
                    "_sum_mole_frac_pdew",
                ):
                    c.deactivate()
            cons_count += number_activated_equalities(blk[k])

        if cons_count > 0:
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
        else:
            res = ""
        init_log.info(
            "Dew and bubble point init complete {}.".format(idaeslog.condition(res))
        )

        # ---------------------------------------------------------------------
        # If flash, initialize T1 and Teq
        for k in blk.keys():
            if (blk[k].params.config.valid_phase == ("Liq", "Vap")) or (
                blk[k].params.config.valid_phase == ("Vap", "Liq")
            ):
                blk[k]._t1.value = max(
                    blk[k].temperature.value, blk[k].temperature_bubble.value
                )
                blk[k]._teq.value = min(blk[k]._t1.value, blk[k].temperature_dew.value)
        init_log.info("Equilibrium temperature init complete.")

        # ---------------------------------------------------------------------
        # Initialize flow rates and compositions
        # TODO : This will need to be generalised more when we move to a
        # modular implementation
        for k in blk.keys():
            if blk[k].params.config.valid_phase == "Liq":
                blk[k].flow_mol_phase["Liq"].value = blk[k].flow_mol.value

                for j in blk[k].params.component_list:
                    blk[k].mole_frac_phase_comp["Liq", j].value = (
                        blk[k].mole_frac_comp[j].value
                    )

            elif blk[k].params.config.valid_phase == "Vap":
                blk[k].flow_mol_phase["Vap"].value = blk[k].flow_mol.value

                for j in blk[k].params.component_list:
                    blk[k].mole_frac_phase_comp["Vap", j].value = (
                        blk[k].mole_frac_comp[j].value
                    )

            else:
                if blk[k].temperature.value > blk[k].temperature_dew.value:
                    # Pure vapour
                    blk[k].flow_mol_phase["Vap"].value = blk[k].flow_mol.value
                    blk[k].flow_mol_phase["Liq"].value = 1e-5 * blk[k].flow_mol.value

                    for j in blk[k].params.component_list:
                        blk[k].mole_frac_phase_comp["Vap", j].value = (
                            blk[k].mole_frac_comp[j].value
                        )
                        blk[k].mole_frac_phase_comp["Liq", j].value = (
                            blk[k]._mole_frac_tdew[j].value
                        )
                elif blk[k].temperature.value < blk[k].temperature_bubble.value:
                    # Pure liquid
                    blk[k].flow_mol_phase["Vap"].value = 1e-5 * blk[k].flow_mol.value
                    blk[k].flow_mol_phase["Liq"].value = blk[k].flow_mol.value

                    for j in blk[k].params.component_list:
                        blk[k].mole_frac_phase_comp["Vap", j].value = (
                            blk[k]._mole_frac_tbub[j].value
                        )
                        blk[k].mole_frac_phase_comp["Liq", j].value = (
                            blk[k].mole_frac_comp[j].value
                        )
                else:
                    # Two-phase
                    # Estimate vapor fraction from distance from Tbub and Tdew
                    # Thanks to Rahul Gandhi for the method
                    vapRatio = value(
                        (blk[k].temperature - blk[k].temperature_bubble)
                        / (blk[k].temperature_dew - blk[k].temperature_bubble)
                    )

                    blk[k].flow_mol_phase["Vap"].value = value(
                        vapRatio * blk[k].flow_mol
                    )
                    blk[k].flow_mol_phase["Liq"].value = value(
                        (1 - vapRatio) * blk[k].flow_mol
                    )

                    # Initialize compositions using Rachford-Rice equation
                    for j in blk[k].params.component_list:
                        kfact = value(
                            antoine_P(blk[k], j, blk[k].temperature.value)
                            / blk[k].pressure
                        )
                        blk[k].mole_frac_phase_comp["Liq", j].value = value(
                            blk[k].mole_frac_comp[j] / (1 + vapRatio * (kfact - 1))
                        )
                        blk[k].mole_frac_phase_comp["Vap", j].value = value(
                            blk[k].mole_frac_phase_comp["Liq", j] * kfact
                        )

        # ---------------------------------------------------------------------
        # Solve phase equilibrium constraints
        for k in blk.keys():
            for c in blk[k].component_objects(Constraint):
                # Activate equilibrium constraints
                if c.local_name in (
                    "total_flow_balance",
                    "component_flow_balances",
                    "equilibrium_constraint",
                    "sum_mole_frac",
                    "_t1_constraint",
                    "_teq_constraint",
                ):
                    c.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = solve_indexed_blocks(opt, [blk], tee=slc.tee)
        init_log.info("Phase equilibrium init: {}.".format(idaeslog.condition(results)))

        # ---------------------------------------------------------------------
        # Initialize other properties
        for k in blk.keys():
            for c in blk[k].component_objects(Constraint):
                # Activate all constraints except sum_mole_frac_out
                if c.local_name not in ("sum_mole_frac_out"):
                    c.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = solve_indexed_blocks(opt, [blk], tee=slc.tee)
        init_log.info("Property init: {}.".format(idaeslog.condition(results)))

        # ---------------------------------------------------------------------
        if not check_optimal_termination(results):
            raise InitializationError(
                f"{blk.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags, outlvl=outlvl)

        init_log.info("Initialization complete.")

    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        """
        Method to relase state variables fixed during initialization.
        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        """
        for k in blk.keys():
            if not blk[k].config.defined_state:
                blk[k].sum_mole_frac_out.activate()

        if flags is None:
            return

        # Unfix state variables
        revert_state_vars(blk, flags)

        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        init_log.info_high("States released.")


@declare_process_block_class("CubicStateBlock", block_class=_CubicStateBlock)
class CubicStateBlockData(StateBlockData):
    """An general property package for cubic equations of state with VLE."""

    def build(self):
        """Callable method for Block construction."""
        super(CubicStateBlockData, self).build()

        self.expr_write = CubicThermoExpressions(self)
        # Add state variables
        self.flow_mol = Var(
            initialize=1.0,
            domain=NonNegativeReals,
            doc="Component molar flowrate [mol/s]",
            units=pyunits.mol / pyunits.s,
        )
        self.mole_frac_comp = Var(
            self.params.component_list,
            bounds=(0, None),
            initialize=1 / len(self.params.component_list),
            doc="Mixture mole fractions [-]",
        )
        self.pressure = Var(
            initialize=101325,
            domain=NonNegativeReals,
            doc="State pressure [Pa]",
            units=pyunits.Pa,
        )
        self.temperature = Var(
            initialize=298.15,
            domain=NonNegativeReals,
            doc="State temperature [K]",
            units=pyunits.K,
        )

        # Add supporting variables
        self.flow_mol_phase = Var(
            self.params.phase_list,
            initialize=0.5,
            domain=NonNegativeReals,
            doc="Phase molar flow rates [mol/s]",
            units=pyunits.mol / pyunits.s,
        )

        self.mole_frac_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=1 / len(self.params.component_list),
            bounds=(0, None),
            doc="Phase mole fractions [-]",
        )

        if self.params.config.valid_phase == "Liq":
            self._make_liq_phase_eq()
        elif self.params.config.valid_phase == "Vap":
            self._make_vap_phase_eq()
        elif (self.params.config.valid_phase == ("Liq", "Vap")) or (
            self.params.config.valid_phase == ("Vap", "Liq")
        ):
            self._make_flash_eq()
        else:
            raise BurntToast(
                "{} found unexpected value for valid_phases. "
                "Please contact the "
                "IDAES developers with this bug.".format(self.name)
            )

    def _make_liq_phase_eq(self):
        # Add equilibrium temperature - in this case the state temperature
        self._teq = Expression(expr=self.temperature)

        # Add supporting equations for Cubic EoS
        self.common_cubic()

        def rule_total_mass_balance(b):
            return b.flow_mol_phase["Liq"] == b.flow_mol

        self.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.mole_frac_comp[i] == b.mole_frac_phase_comp["Liq", i]

        self.component_flow_balances = Constraint(
            self.params.component_list, rule=rule_comp_mass_balance
        )

        if self.config.defined_state is False:
            # applied at outlet only
            self.sum_mole_frac_out = Constraint(
                expr=1
                == sum(self.mole_frac_comp[i] for i in self.params.component_list)
            )

    def _make_vap_phase_eq(self):
        # Add equilibrium temperature - in this case the state temperature
        self._teq = Expression(expr=self.temperature)

        # Add supporting equations for Cubic EoS
        self.common_cubic()

        def rule_total_mass_balance(b):
            return b.flow_mol_phase["Vap"] == b.flow_mol

        self.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.mole_frac_comp[i] == b.mole_frac_phase_comp["Vap", i]

        self.component_flow_balances = Constraint(
            self.params.component_list, rule=rule_comp_mass_balance
        )

        if self.config.defined_state is False:
            # applied at outlet only
            self.sum_mole_frac_out = Constraint(
                expr=1
                == sum(self.mole_frac_comp[i] for i in self.params.component_list)
            )

    def _make_flash_eq(self):
        """
        Implementation of smooth VLE formulation.
        See module header for reference.
        """

        def rule_total_mass_balance(b):
            return b.flow_mol_phase["Liq"] + b.flow_mol_phase["Vap"] == b.flow_mol

        self.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return (
                b.flow_mol * b.mole_frac_comp[i]
                == b.flow_mol_phase["Liq"] * b.mole_frac_phase_comp["Liq", i]
                + b.flow_mol_phase["Vap"] * b.mole_frac_phase_comp["Vap", i]
            )

        self.component_flow_balances = Constraint(
            self.params.component_list, rule=rule_comp_mass_balance
        )

        def rule_mole_frac(b):
            return (
                sum(b.mole_frac_phase_comp["Liq", i] for i in b.params.component_list)
                - sum(b.mole_frac_phase_comp["Vap", i] for i in b.params.component_list)
                == 0
            )

        self.sum_mole_frac = Constraint(rule=rule_mole_frac)

        if self.config.defined_state is False:
            # applied at outlet only
            self.sum_mole_frac_out = Constraint(
                expr=1
                == sum(self.mole_frac_comp[i] for i in self.params.component_list)
            )

        # Definition of equilibrium temperature for smooth VLE
        self._teq = Var(
            initialize=self.temperature.value,
            doc="Temperature for calculating phase equilibrium",
            units=pyunits.K,
        )
        self._t1 = Var(
            initialize=self.temperature.value,
            doc="Intermediate temperature for calculating Teq",
            units=pyunits.K,
        )

        self.eps_1 = Param(
            default=0.01,
            mutable=True,
            doc="Smoothing parameter for Teq",
            units=pyunits.K,
        )
        self.eps_2 = Param(
            default=0.0005,
            mutable=True,
            doc="Smoothing parameter for Teq",
            units=pyunits.K,
        )

        # Add supporting equations for Cubic EoS
        self.common_cubic()

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
            return b._teq / b.params.temperature_crit[i]

        self._tr_eq = Expression(
            self.params.component_list,
            rule=rule_tr_eq,
            doc="Component reduced temperatures [-]",
        )

        def rule_equilibrium(b, i):
            return b._log_equilibrium_cubic("Vap", i) == b._log_equilibrium_cubic(
                "Liq", i
            )

        self.equilibrium_constraint = Constraint(
            self.params.component_list, rule=rule_equilibrium
        )

    # -----------------------------------------------------------------------------
    # Property Methods
    def _dens_mol_phase(self):
        self.dens_mol_phase = Var(
            self.params.phase_list,
            doc="Molar density [mol/m^3]",
            units=pyunits.mol / pyunits.m**3,
        )

        def rule_dens_mol_phase(b, p):
            if p == "Vap":
                return b._dens_mol_vap()
            else:
                return b._dens_mol_liq()

        self.eq_dens_mol_phase = Constraint(
            self.params.phase_list, rule=rule_dens_mol_phase
        )

    def _dens_mass_phase(self):
        self.dens_mass_phase = Var(
            self.params.phase_list,
            doc="Mass density [kg/m^3]",
            units=pyunits.kg / pyunits.m**3,
        )

        def rule_dens_mass_phase(b, p):
            if p == "Vap":
                return b._dens_mass_vap()
            else:
                return b._dens_mass_liq()

        self.eq_dens_mass_phase = Constraint(
            self.params.phase_list, rule=rule_dens_mass_phase
        )

    def _energy_internal_mol_phase(self):
        self.energy_internal_mol_phase = Var(
            self.params.phase_list,
            doc="Phase molar specific internal energies [J/mol]",
            units=pyunits.J / pyunits.mol,
        )

        def rule_energy_internal_mol_phase(b, p):
            if p == "Vap":
                return b.energy_internal_mol_phase[p] == b._energy_internal_mol_vap()
            else:
                return b.energy_internal_mol_phase[p] == b._energy_internal_mol_liq()

        self.eq_energy_internal_mol_phase = Constraint(
            self.params.phase_list, rule=rule_energy_internal_mol_phase
        )

    def _energy_internal_mol(self):
        self.energy_internal_mol = Var(
            doc="Mixture molar specific internal energy [J/mol]",
            units=pyunits.J / pyunits.mol,
        )

        def rule_energy_internal_mol(b):
            return b.energy_internal_mol * b.flow_mol == sum(
                b.flow_mol_phase[p] * b.energy_internal_mol_phase[p]
                for p in b.params.phase_list
            )

        self.eq_energy_internal_mol = Constraint(rule=rule_energy_internal_mol)

    def _enth_mol_phase(self):
        self.enth_mol_phase = Var(
            self.params.phase_list,
            doc="Phase molar specific enthalpies [J/mol]",
            units=pyunits.J / pyunits.mol,
        )

        def rule_enth_mol_phase(b, p):
            if p == "Vap":
                return b.enth_mol_phase[p] == b._enth_mol_vap()
            else:
                return b.enth_mol_phase[p] == b._enth_mol_liq()

        self.eq_enth_mol_phase = Constraint(
            self.params.phase_list, rule=rule_enth_mol_phase
        )

    def _enth_mol(self):
        self.enth_mol = Var(
            doc="Mixture molar specific enthalpies [J/mol]",
            units=pyunits.J / pyunits.mol,
        )

        def rule_enth_mol(b):
            return b.enth_mol * b.flow_mol == sum(
                b.flow_mol_phase[p] * b.enth_mol_phase[p] for p in b.params.phase_list
            )

        self.eq_enth_mol = Constraint(rule=rule_enth_mol)

    def _entr_mol(self):
        self.entr_mol = Var(
            doc="Mixture molar specific entropies [J/mol.K]",
            units=pyunits.J / pyunits.mol / pyunits.K,
        )

        def rule_entr_mol(b):
            return b.entr_mol * b.flow_mol == sum(
                b.flow_mol_phase[p] * b.entr_mol_phase[p] for p in b.params.phase_list
            )

        self.eq_entr_mol = Constraint(rule=rule_entr_mol)

    def _entr_mol_phase(self):
        self.entr_mol_phase = Var(
            self.params.phase_list,
            doc="Phase molar specific entropies [J/mol.K]",
            units=pyunits.J / pyunits.mol / pyunits.K,
        )

        def rule_entr_mol_phase(b, p):
            if p == "Vap":
                return b.entr_mol_phase[p] == b._entr_mol_vap()
            else:
                return b.entr_mol_phase[p] == b._entr_mol_liq()

        self.eq_entr_mol_phase = Constraint(
            self.params.phase_list, rule=rule_entr_mol_phase
        )

    def _fug_phase(self):
        def rule_fug_phase(b, p, j):
            if p == "Vap":
                return b._fug_vap(j)
            else:
                return b._fug_liq(j)

        self.fug_phase_comp = Expression(
            self.params.phase_list, self.params.component_list, rule=rule_fug_phase
        )

    def _fug_coeff_phase(self):
        def rule_fug_coeff_phase(b, p, j):
            if p == "Vap":
                return b._fug_coeff_vap(j)
            else:
                return b._fug_coeff_liq(j)

        self.fug_coeff_phase_comp = Expression(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_fug_coeff_phase,
        )

    def _gibbs_mol_phase(self):
        self.gibbs_mol_phase = Var(
            self.params.phase_list,
            doc="Phase molar specific Gibbs energy [J/mol]",
            units=pyunits.J / pyunits.mol,
        )

        def rule_gibbs_mol_phase(b, p):
            return b.gibbs_mol_phase[p] == (
                b.enth_mol_phase[p] - b.temperature * b.entr_mol_phase[p]
            )

        self.eq_gibbs_mol_phase = Constraint(
            self.params.phase_list, rule=rule_gibbs_mol_phase
        )

    def _mw(self):
        def rule_mw(b):
            return sum(b.mw_phase[p] for p in b.params.phase_list)

        self.mw = Expression(rule=rule_mw)

    def _mw_phase(self):
        def rule_mw_phase(b, p):
            return sum(
                b.mole_frac_phase_comp[p, j] * b.params.mw_comp[j]
                for j in b.params.component_list
            )

        self.mw_phase = Expression(self.params.phase_list, rule=rule_mw_phase)

    # -----------------------------------------------------------------------------
    # General Methods
    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        if not self.is_property_constructed("material_flow_terms"):
            try:

                def rule_material_flow_terms(b, p, j):
                    return self.flow_mol_phase[p] * self.mole_frac_phase_comp[p, j]

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

                def rule_enthalpy_flow_terms(b, p):
                    return self.flow_mol_phase[p] * self.enth_mol_phase[p]

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

    def get_energy_density_terms(self, p):
        """Create energy density terms."""
        if not self.is_property_constructed("energy_density_terms"):
            try:

                def rule_energy_density_terms(b, p):
                    return self.dens_mol_phase[p] * self.energy_internal_mol_phase[p]

                self.energy_density_terms = Expression(
                    self.params.phase_list, rule=rule_energy_density_terms
                )
            except AttributeError:
                self.del_component(self.energy_density_terms)
        return self.energy_density_terms[p]

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_basis(b):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        """Define state vars."""
        return {
            "flow_mol": self.flow_mol,
            "mole_frac_comp": self.mole_frac_comp,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }

    def define_display_vars(b):
        return {
            "Molar Flowrate": b.flow_mol,
            "Mole Fractions": b.mole_frac_comp,
            "Temperature": b.temperature,
            "Pressure": b.pressure,
        }

    def model_check(blk):
        """Model checks for property block."""
        # Check temperature bounds
        if value(blk.temperature) < blk.temperature.lb:
            _log.error("{} Temperature set below lower bound.".format(blk.name))
        if value(blk.temperature) > blk.temperature.ub:
            _log.error("{} Temperature set above upper bound.".format(blk.name))

        # Check pressure bounds
        if value(blk.pressure) < blk.pressure.lb:
            _log.error("{} Pressure set below lower bound.".format(blk.name))
        if value(blk.pressure) > blk.pressure.ub:
            _log.error("{} Pressure set above upper bound.".format(blk.name))

    # -----------------------------------------------------------------------------
    # Bubble and Dew Points
    def _temperature_bubble(self):
        self.temperature_bubble = Var(
            doc="Bubble point temperature (K)", units=pyunits.K
        )

        self._mole_frac_tbub = Var(
            self.params.component_list,
            initialize=1 / len(self.params.component_list),
            bounds=(0, None),
            doc="Vapor mole fractions at bubble point",
        )

        self._sum_mole_frac_tbub = Constraint(
            expr=1e3
            == 1e3 * sum(self._mole_frac_tbub[j] for j in self.params.component_list)
        )

        def rule_bubble_temp(b, j):
            return log(b.mole_frac_comp[j]) + log(b.bubble_temp_liq(j)) == log(
                b._mole_frac_tbub[j]
            ) + log(b.bubble_temp_vap(j))

        self.eq_temperature_bubble = Constraint(
            self.params.component_list, rule=rule_bubble_temp
        )

    def _temperature_dew(self):
        self.temperature_dew = Var(doc="Dew point temperature (K)", units=pyunits.K)

        self._mole_frac_tdew = Var(
            self.params.component_list,
            initialize=1 / len(self.params.component_list),
            bounds=(0, None),
            doc="Liquid mole fractions at dew point",
        )

        self._sum_mole_frac_tdew = Constraint(
            expr=1e3
            == 1e3 * sum(self._mole_frac_tdew[j] for j in self.params.component_list)
        )

        def rule_dew_temp(b, j):
            return log(b._mole_frac_tdew[j]) + log(b.dew_temp_liq(j)) == log(
                b.mole_frac_comp[j]
            ) + log(b.dew_temp_vap(j))

        self.eq_temperature_dew = Constraint(
            self.params.component_list, rule=rule_dew_temp
        )

    def _pressure_bubble(self):
        self.pressure_bubble = Var(
            domain=PositiveReals, doc="Bubble point pressure (Pa)", units=pyunits.Pa
        )

        self._mole_frac_pbub = Var(
            self.params.component_list,
            initialize=1 / len(self.params.component_list),
            bounds=(0, None),
            doc="Vapor mole fractions at bubble point",
        )

        self._sum_mole_frac_pbub = Constraint(
            expr=1e3
            == 1e3 * sum(self._mole_frac_pbub[j] for j in self.params.component_list)
        )

        def rule_bubble_pres(b, j):
            return log(b.mole_frac_comp[j]) + log(b.bubble_pres_liq(j)) == log(
                b._mole_frac_pbub[j]
            ) + log(b.bubble_pres_vap(j))

        self.eq_pressure_bubble = Constraint(
            self.params.component_list, rule=rule_bubble_pres
        )

    def _pressure_dew(self):
        self.pressure_dew = Var(
            domain=PositiveReals, doc="Dew point pressure (Pa)", units=pyunits.Pa
        )

        self._mole_frac_pdew = Var(
            self.params.component_list,
            initialize=1 / len(self.params.component_list),
            bounds=(0, None),
            doc="Liquid mole fractions at dew point",
        )

        self._sum_mole_frac_pdew = Constraint(
            expr=1 == sum(self._mole_frac_pdew[j] for j in self.params.component_list)
        )

        def rule_dew_press(b, j):
            return log(b._mole_frac_pdew[j]) + log(b.dew_press_liq(j)) == log(
                b.mole_frac_comp[j]
            ) + log(b.dew_press_vap(j))

        self.eq_pressure_dew = Constraint(
            self.params.component_list, rule=rule_dew_press
        )

    # -----------------------------------------------------------------------------
    # Liquid phase properties
    def _vol_mol_liq(b):
        return b._vol_mol_cubic("Liq")

    def _dens_mol_liq(b):
        return b._dens_mol_cubic("Liq")

    def _dens_mass_liq(b):
        return b._dens_mass_cubic("Liq")

    def _fug_liq(b, j):
        return b._fug_cubic("Liq", j)

    def _fug_coeff_liq(b, j):
        return b._fug_coeff_cubic("Liq", j)

    def _energy_internal_mol_liq(b):
        return b._energy_internal_mol_cubic("Liq")

    def _enth_mol_liq(b):
        return b._enth_mol_cubic("Liq")

    def _enth_mol_liq_ig(b):
        return b._enth_mol_ig("Liq")

    def _entr_mol_liq(b):
        return b._entr_mol_cubic("Liq")

    def _entr_mol_liq_ig(b):
        return b._entr_mol_ig("Liq")

    def bubble_temp_liq(b, j):
        def a(k):
            return (
                b.omegaA
                * (
                    (const.gas_constant * b.params.temperature_crit[k]) ** 2
                    / b.params.pressure_crit[k]
                )
                * (
                    (
                        1
                        + b.fw[k]
                        * (
                            1
                            - sqrt(b.temperature_bubble / b.params.temperature_crit[k])
                        )
                    )
                    ** 2
                )
            )

        am = sum(
            sum(
                b.mole_frac_comp[i]
                * b.mole_frac_comp[j]
                * sqrt(a(i) * a(j))
                * (1 - b.params.kappa[i, j])
                for j in b.params.component_list
            )
            for i in b.params.component_list
        )
        bm = sum(b.mole_frac_comp[i] * b.b[i] for i in b.params.component_list)

        A = am * b.pressure / (const.gas_constant * b.temperature_bubble) ** 2
        B = bm * b.pressure / (const.gas_constant * b.temperature_bubble)

        delta = (
            2
            * sqrt(a(j))
            / am
            * sum(
                b.mole_frac_comp[i] * sqrt(a(i)) * (1 - b.params.kappa[j, i])
                for i in b.params.component_list
            )
        )

        Z = b.expr_write.z_liq(eos=b.params.cubic_type, A=A, B=B)

        return exp(
            (
                b.b[j] / bm * (Z - 1) * (B * b.EoS_p)
                - safe_log(Z - B, eps=1e-6) * (B * b.EoS_p)
                + A
                * (b.b[j] / bm - delta)
                * safe_log(
                    (2 * Z + B * (b.EoS_u + b.EoS_p))
                    / (2 * Z + B * (b.EoS_u - b.EoS_p)),
                    eps=1e-6,
                )
            )
            / (B * b.EoS_p)
        )

    def dew_temp_liq(b, j):
        def a(k):
            return (
                b.omegaA
                * (
                    (const.gas_constant * b.params.temperature_crit[k]) ** 2
                    / b.params.pressure_crit[k]
                )
                * (
                    (
                        1
                        + b.fw[k]
                        * (1 - sqrt(b.temperature_dew / b.params.temperature_crit[k]))
                    )
                    ** 2
                )
            )

        am = sum(
            sum(
                b._mole_frac_tdew[i]
                * b._mole_frac_tdew[j]
                * sqrt(a(i) * a(j))
                * (1 - b.params.kappa[i, j])
                for j in b.params.component_list
            )
            for i in b.params.component_list
        )
        bm = sum(b._mole_frac_tdew[i] * b.b[i] for i in b.params.component_list)

        A = am * b.pressure / (const.gas_constant * b.temperature_dew) ** 2
        B = bm * b.pressure / (const.gas_constant * b.temperature_dew)

        delta = (
            2
            * sqrt(a(j))
            / am
            * sum(
                b._mole_frac_tdew[i] * sqrt(a(i)) * (1 - b.params.kappa[j, i])
                for i in b.params.component_list
            )
        )

        Z = b.expr_write.z_liq(eos=b.params.cubic_type, A=A, B=B)

        return exp(
            (
                b.b[j] / bm * (Z - 1) * (B * b.EoS_p)
                - safe_log(Z - B, eps=1e-6) * (B * b.EoS_p)
                + A
                * (b.b[j] / bm - delta)
                * safe_log(
                    (2 * Z + B * (b.EoS_u + b.EoS_p))
                    / (2 * Z + B * (b.EoS_u - b.EoS_p)),
                    eps=1e-6,
                )
            )
            / (B * b.EoS_p)
        )

    def bubble_pres_liq(b, j):
        am = sum(
            sum(
                b.mole_frac_comp[i]
                * b.mole_frac_comp[j]
                * sqrt(b.a[i] * b.a[j])
                * (1 - b.params.kappa[i, j])
                for j in b.params.component_list
            )
            for i in b.params.component_list
        )
        bm = sum(b.mole_frac_comp[i] * b.b[i] for i in b.params.component_list)

        A = am * b.pressure_bubble / (const.gas_constant * b.temperature) ** 2
        B = bm * b.pressure_bubble / (const.gas_constant * b.temperature)

        delta = (
            2
            * sqrt(b.a[j])
            / am
            * sum(
                b.mole_frac_comp[i] * sqrt(b.a[i]) * (1 - b.params.kappa[j, i])
                for i in b.params.component_list
            )
        )

        Z = b.expr_write.z_liq(eos=b.params.cubic_type, A=A, B=B)

        return exp(
            (
                b.b[j] / bm * (Z - 1) * (B * b.EoS_p)
                - safe_log(Z - B, eps=1e-6) * (B * b.EoS_p)
                + A
                * (b.b[j] / bm - delta)
                * safe_log(
                    (2 * Z + B * (b.EoS_u + b.EoS_p))
                    / (2 * Z + B * (b.EoS_u - b.EoS_p)),
                    eps=1e-6,
                )
            )
            / (B * b.EoS_p)
        )

    def dew_press_liq(b, j):
        am = sum(
            sum(
                b._mole_frac_pdew[i]
                * b._mole_frac_pdew[j]
                * sqrt(b.a[i] * b.a[j])
                * (1 - b.params.kappa[i, j])
                for j in b.params.component_list
            )
            for i in b.params.component_list
        )
        bm = sum(b._mole_frac_pdew[i] * b.b[i] for i in b.params.component_list)

        A = am * b.pressure_dew / (const.gas_constant * b.temperature) ** 2
        B = bm * b.pressure_dew / (const.gas_constant * b.temperature)

        delta = (
            2
            * sqrt(b.a[j])
            / am
            * sum(
                b._mole_frac_pdew[i] * sqrt(b.a[i]) * (1 - b.params.kappa[j, i])
                for i in b.params.component_list
            )
        )

        Z = b.expr_write.z_liq(eos=b.params.cubic_type, A=A, B=B)

        return exp(
            (
                b.b[j] / bm * (Z - 1) * (B * b.EoS_p)
                - safe_log(Z - B, eps=1e-6) * (B * b.EoS_p)
                + A
                * (b.b[j] / bm - delta)
                * safe_log(
                    (2 * Z + B * (b.EoS_u + b.EoS_p))
                    / (2 * Z + B * (b.EoS_u - b.EoS_p)),
                    eps=1e-6,
                )
            )
            / (B * b.EoS_p)
        )

    # -----------------------------------------------------------------------------
    # Vapour phase properties
    def _vol_mol_vap(b):
        return b._vol_mol_cubic("Vap")

    def _dens_mol_vap(b):
        return b._dens_mol_cubic("Vap")

    def _dens_mass_vap(b):
        return b._dens_mass_cubic("Vap")

    def _fug_vap(b, j):
        return b._fug_cubic("Vap", j)

    def _fug_coeff_vap(b, j):
        return b._fug_coeff_cubic("Vap", j)

    def _energy_internal_mol_vap(b):
        return b._energy_internal_mol_cubic("Vap")

    def _enth_mol_vap(b):
        return b._enth_mol_cubic("Vap")

    def _enth_mol_vap_ig(b):
        return b._enth_mol_ig("Vap")

    def _entr_mol_vap(b):
        return b._entr_mol_cubic("Vap")

    def _entr_mol_vap_ig(b):
        return b._entr_mol_ig("Vap")

    def bubble_temp_vap(b, j):
        def a(k):
            return (
                b.omegaA
                * (
                    (const.gas_constant * b.params.temperature_crit[k]) ** 2
                    / b.params.pressure_crit[k]
                )
                * (
                    (
                        1
                        + b.fw[k]
                        * (
                            1
                            - sqrt(b.temperature_bubble / b.params.temperature_crit[k])
                        )
                    )
                    ** 2
                )
            )

        am = sum(
            sum(
                b._mole_frac_tbub[i]
                * b._mole_frac_tbub[j]
                * sqrt(a(i) * a(j))
                * (1 - b.params.kappa[i, j])
                for j in b.params.component_list
            )
            for i in b.params.component_list
        )
        bm = sum(b._mole_frac_tbub[i] * b.b[i] for i in b.params.component_list)

        A = am * b.pressure / (const.gas_constant * b.temperature_bubble) ** 2
        B = bm * b.pressure / (const.gas_constant * b.temperature_bubble)

        delta = (
            2
            * sqrt(a(j))
            / am
            * sum(
                b._mole_frac_tbub[i] * sqrt(a(i)) * (1 - b.params.kappa[j, i])
                for i in b.params.component_list
            )
        )

        Z = b.expr_write.z_vap(eos=b.params.cubic_type, A=A, B=B)

        return exp(
            (
                b.b[j] / bm * (Z - 1) * (B * b.EoS_p)
                - safe_log(Z - B, eps=1e-6) * (B * b.EoS_p)
                + A
                * (b.b[j] / bm - delta)
                * safe_log(
                    (2 * Z + B * (b.EoS_u + b.EoS_p))
                    / (2 * Z + B * (b.EoS_u - b.EoS_p)),
                    eps=1e-6,
                )
            )
            / (B * b.EoS_p)
        )

    def dew_temp_vap(b, j):
        def a(k):
            return (
                b.omegaA
                * (
                    (const.gas_constant * b.params.temperature_crit[k]) ** 2
                    / b.params.pressure_crit[k]
                )
                * (
                    (
                        1
                        + b.fw[k]
                        * (1 - sqrt(b.temperature_dew / b.params.temperature_crit[k]))
                    )
                    ** 2
                )
            )

        am = sum(
            sum(
                b.mole_frac_comp[i]
                * b.mole_frac_comp[j]
                * sqrt(a(i) * a(j))
                * (1 - b.params.kappa[i, j])
                for j in b.params.component_list
            )
            for i in b.params.component_list
        )
        bm = sum(b.mole_frac_comp[i] * b.b[i] for i in b.params.component_list)

        A = am * b.pressure / (const.gas_constant * b.temperature_dew) ** 2
        B = bm * b.pressure / (const.gas_constant * b.temperature_dew)

        delta = (
            2
            * sqrt(a(j))
            / am
            * sum(
                b.mole_frac_comp[i] * sqrt(a(i)) * (1 - b.params.kappa[j, i])
                for i in b.params.component_list
            )
        )

        Z = b.expr_write.z_vap(eos=b.params.cubic_type, A=A, B=B)

        return exp(
            (
                b.b[j] / bm * (Z - 1) * (B * b.EoS_p)
                - safe_log(Z - B, eps=1e-6) * (B * b.EoS_p)
                + A
                * (b.b[j] / bm - delta)
                * safe_log(
                    (2 * Z + B * (b.EoS_u + b.EoS_p))
                    / (2 * Z + B * (b.EoS_u - b.EoS_p)),
                    eps=1e-6,
                )
            )
            / (B * b.EoS_p)
        )

    def bubble_pres_vap(b, j):
        am = sum(
            sum(
                b._mole_frac_pbub[i]
                * b._mole_frac_pbub[j]
                * sqrt(b.a[i] * b.a[j])
                * (1 - b.params.kappa[i, j])
                for j in b.params.component_list
            )
            for i in b.params.component_list
        )
        bm = sum(b._mole_frac_pbub[i] * b.b[i] for i in b.params.component_list)

        A = am * b.pressure_bubble / (const.gas_constant * b.temperature) ** 2
        B = bm * b.pressure_bubble / (const.gas_constant * b.temperature)

        delta = (
            2
            * sqrt(b.a[j])
            / am
            * sum(
                b._mole_frac_pbub[i] * sqrt(b.a[i]) * (1 - b.params.kappa[j, i])
                for i in b.params.component_list
            )
        )

        Z = b.expr_write.z_vap(eos=b.params.cubic_type, A=A, B=B)

        return exp(
            (
                b.b[j] / bm * (Z - 1) * (B * b.EoS_p)
                - safe_log(Z - B, eps=1e-6) * (B * b.EoS_p)
                + A
                * (b.b[j] / bm - delta)
                * safe_log(
                    (2 * Z + B * (b.EoS_u + b.EoS_p))
                    / (2 * Z + B * (b.EoS_u - b.EoS_p)),
                    eps=1e-6,
                )
            )
            / (B * b.EoS_p)
        )

    def dew_press_vap(b, j):
        am = sum(
            sum(
                b.mole_frac_comp[i]
                * b.mole_frac_comp[j]
                * sqrt(b.a[i] * b.a[j])
                * (1 - b.params.kappa[i, j])
                for j in b.params.component_list
            )
            for i in b.params.component_list
        )
        bm = sum(b.mole_frac_comp[i] * b.b[i] for i in b.params.component_list)

        A = am * b.pressure_dew / (const.gas_constant * b.temperature) ** 2
        B = bm * b.pressure_dew / (const.gas_constant * b.temperature)

        delta = (
            2
            * sqrt(b.a[j])
            / am
            * sum(
                b.mole_frac_comp[i] * sqrt(b.a[i]) * (1 - b.params.kappa[j, i])
                for i in b.params.component_list
            )
        )

        Z = b.expr_write.z_vap(eos=b.params.cubic_type, A=A, B=B)

        return exp(
            (
                b.b[j] / bm * (Z - 1) * (B * b.EoS_p)
                - safe_log(Z - B, eps=1e-6) * (B * b.EoS_p)
                + A
                * (b.b[j] / bm - delta)
                * safe_log(
                    (2 * Z + B * (b.EoS_u + b.EoS_p))
                    / (2 * Z + B * (b.EoS_u - b.EoS_p)),
                    eps=1e-6,
                )
            )
            / (B * b.EoS_p)
        )

    # -----------------------------------------------------------------------------
    # Common Cubic Functions
    # All of these equations drawn from Properties of Gases and Liquids
    # Quantities appended with _eq represent calculations @ equilibrium temperature
    def common_cubic(blk):
        if hasattr(blk, "omegaA"):
            return

        blk.omegaA = EoS_param[blk.params.cubic_type]["omegaA"]

        blk.EoS_Bc = EoS_param[blk.params.cubic_type]["coeff_b"]
        blk.EoS_u = EoS_param[blk.params.cubic_type]["u"]
        blk.EoS_w = EoS_param[blk.params.cubic_type]["w"]
        blk.EoS_p = sqrt(blk.EoS_u**2 - 4 * blk.EoS_w)

        # Create expressions for coefficients
        def func_fw(b, j):
            if b.params.cubic_type == CubicEoS.PR:
                return (
                    0.37464
                    + 1.54226 * b.params.omega[j]
                    - 0.26992 * b.params.omega[j] ** 2
                )
            elif b.params.cubic_type == CubicEoS.SRK:
                return 0.48 + 1.574 * b.params.omega[j] - 0.176 * b.params.omega[j] ** 2
            else:
                raise BurntToast(
                    "{} received unrecognised cubic type. This should "
                    "never happen, so please contact the IDAES developers "
                    "with this bug.".format(b.name)
                )

        blk.fw = Param(
            blk.params.component_list, initialize=func_fw, doc="EoS S factor"
        )

        def func_b(b, j):
            return (
                b.EoS_Bc
                * const.gas_constant
                * b.params.temperature_crit[j]
                / b.params.pressure_crit[j]
            )

        blk.b = Param(
            blk.params.component_list,
            initialize=func_b,
            doc="Component b coefficient",
            units=pyunits.m**3 / pyunits.mol,
        )

        def func_a(b, j):
            return (
                b.omegaA
                * (
                    (const.gas_constant * b.params.temperature_crit[j]) ** 2
                    / b.params.pressure_crit[j]
                )
                * (
                    (
                        1
                        + b.fw[j]
                        * (1 - sqrt(b.temperature / b.params.temperature_crit[j]))
                    )
                    ** 2
                )
            )

        blk.a = Expression(
            blk.params.component_list, rule=func_a, doc="Component a coefficient"
        )

        def func_a_eq(b, j):
            return (
                b.omegaA
                * (
                    (const.gas_constant * b.params.temperature_crit[j]) ** 2
                    / b.params.pressure_crit[j]
                )
                * (
                    (1 + b.fw[j] * (1 - sqrt(b._teq / b.params.temperature_crit[j])))
                    ** 2
                )
            )

        blk._a_eq = Expression(
            blk.params.component_list,
            rule=func_a_eq,
            doc="Component a coefficient at Teq",
        )

        def rule_am(b, p):
            return sum(
                sum(
                    b.mole_frac_phase_comp[p, i]
                    * b.mole_frac_phase_comp[p, j]
                    * sqrt(b.a[i] * b.a[j])
                    * (1 - b.params.kappa[i, j])
                    for j in b.params.component_list
                )
                for i in b.params.component_list
            )

        blk.am = Expression(blk.params.phase_list, rule=rule_am)

        def rule_am_eq(b, p):
            return sum(
                sum(
                    b.mole_frac_phase_comp[p, i]
                    * b.mole_frac_phase_comp[p, j]
                    * sqrt(b._a_eq[i] * b._a_eq[j])
                    * (1 - b.params.kappa[i, j])
                    for j in b.params.component_list
                )
                for i in b.params.component_list
            )

        blk._am_eq = Expression(blk.params.phase_list, rule=rule_am_eq)

        def rule_bm(b, p):
            return sum(
                b.mole_frac_phase_comp[p, i] * b.b[i] for i in b.params.component_list
            )

        blk.bm = Expression(blk.params.phase_list, rule=rule_bm)

        def rule_A(b, p):
            return b.am[p] * b.pressure / (const.gas_constant * b.temperature) ** 2

        blk.A = Expression(blk.params.phase_list, rule=rule_A)

        def rule_B(b, p):
            return b.bm[p] * b.pressure / (const.gas_constant * b.temperature)

        blk.B = Expression(blk.params.phase_list, rule=rule_B)

        def rule_A_eq(b, p):
            return b._am_eq[p] * b.pressure / (const.gas_constant * b._teq) ** 2

        blk._A_eq = Expression(blk.params.phase_list, rule=rule_A_eq)

        def rule_B_eq(b, p):
            return b.bm[p] * b.pressure / (const.gas_constant * b._teq)

        blk._B_eq = Expression(blk.params.phase_list, rule=rule_B_eq)

        def rule_delta(b, p, i):
            # See pg. 145 in Properties of Gases and Liquids
            return (
                2
                * sqrt(blk.a[i])
                / b.am[p]
                * sum(
                    b.mole_frac_phase_comp[p, j]
                    * sqrt(blk.a[j])
                    * (1 - b.params.kappa[i, j])
                    for j in b.params.component_list
                )
            )

        blk.delta = Expression(
            blk.params.phase_list, blk.params.component_list, rule=rule_delta
        )

        def rule_delta_eq(b, p, i):
            # See pg. 145 in Properties of Gases and Liquids
            return (
                2
                * sqrt(blk._a_eq[i])
                / b._am_eq[p]
                * sum(
                    b.mole_frac_phase_comp[p, j]
                    * sqrt(blk._a_eq[j])
                    * (1 - b.params.kappa[i, j])
                    for j in b.params.component_list
                )
            )

        blk._delta_eq = Expression(
            blk.params.phase_list, blk.params.component_list, rule=rule_delta_eq
        )

        def rule_dadT(b, p):
            # See pg. 102 in Properties of Gases and Liquids
            return -(
                (const.gas_constant / 2)
                * sqrt(b.omegaA)
                * sum(
                    sum(
                        b.mole_frac_phase_comp[p, i]
                        * b.mole_frac_phase_comp[p, j]
                        * (1 - b.params.kappa[i, j])
                        * (
                            b.fw[j]
                            * sqrt(
                                b.a[i]
                                * b.params.temperature_crit[j]
                                / b.params.pressure_crit[j]
                            )
                            + b.fw[i]
                            * sqrt(
                                b.a[j]
                                * b.params.temperature_crit[i]
                                / b.params.pressure_crit[i]
                            )
                        )
                        for j in b.params.component_list
                    )
                    for i in b.params.component_list
                )
                / sqrt(b.temperature)
            )

        blk.dadT = Expression(blk.params.phase_list, rule=rule_dadT)

        def rule_compress_fact(b, p):
            if p == "Vap":
                return blk.expr_write.z_vap(eos=b.params.cubic_type, A=b.A[p], B=b.B[p])
            else:
                return blk.expr_write.z_liq(eos=b.params.cubic_type, A=b.A[p], B=b.B[p])

        blk.compress_fact_phase = Expression(
            blk.params.phase_list, rule=rule_compress_fact
        )

        def rule_compress_fact_eq(b, p):
            if p == "Vap":
                return blk.expr_write.z_vap(
                    eos=b.params.cubic_type, A=b._A_eq[p], B=b._B_eq[p]
                )
            else:
                return blk.expr_write.z_liq(
                    eos=b.params.cubic_type, A=b._A_eq[p], B=b._B_eq[p]
                )

        blk._compress_fact_eq = Expression(
            blk.params.phase_list, rule=rule_compress_fact_eq
        )

    def _vol_mol_cubic(b, p):
        return (
            b.pressure * b.vol_mol_phase[p]
            == b.compress_fact_phase[p] * const.gas_constant * b.temperature
        )

    def _dens_mol_cubic(b, p):
        return b.pressure == (
            b.dens_mol_phase[p]
            * b.compress_fact_phase[p]
            * const.gas_constant
            * b.temperature
        )

    def _dens_mass_cubic(b, p):
        return b.dens_mass_phase[p] == b.dens_mol_phase[p] * b.mw_phase[p]

    def _fug_cubic(b, p, j):
        return b.mole_frac_phase_comp[p, j] * b.pressure * b.fug_coeff_phase_comp[p, j]

    def _fug_coeff_cubic(b, p, j):
        # See pg. 145 in Properties of Gases and Liquids
        return exp(
            (
                b.b[j] / b.bm[p] * (b.compress_fact_phase[p] - 1) * (b.B[p] * b.EoS_p)
                - safe_log(b.compress_fact_phase[p] - b.B[p], eps=1e-6)
                * (b.B[p] * b.EoS_p)
                + b.A[p]
                * (b.b[j] / b.bm[p] - b.delta[p, j])
                * safe_log(
                    (2 * b.compress_fact_phase[p] + b.B[p] * (b.EoS_u + b.EoS_p))
                    / (2 * b.compress_fact_phase[p] + b.B[p] * (b.EoS_u - b.EoS_p)),
                    eps=1e-6,
                )
            )
            / (b.B[p] * b.EoS_p)
        )

    def _log_equilibrium_cubic(b, p, j):
        # See pg. 145 in Properties of Gases and Liquids
        return (
            b.b[j] / b.bm[p] * (b._compress_fact_eq[p] - 1) * (b._B_eq[p] * b.EoS_p)
            - safe_log(b._compress_fact_eq[p] - b._B_eq[p], eps=1e-6)
            * (b._B_eq[p] * b.EoS_p)
            + b._A_eq[p]
            * (b.b[j] / b.bm[p] - b._delta_eq[p, j])
            * safe_log(
                (2 * b._compress_fact_eq[p] + b._B_eq[p] * (b.EoS_u + b.EoS_p))
                / (2 * b._compress_fact_eq[p] + b._B_eq[p] * (b.EoS_u - b.EoS_p)),
                eps=1e-6,
            )
        ) / (b._B_eq[p] * b.EoS_p) + log(b.mole_frac_phase_comp[p, j])

    def _energy_internal_mol_cubic(b, p):
        # Derived from equation on pg. 120 in Properties of Gases and Liquids
        return (
            (b.temperature * b.dadT[p] - b.am[p])
            * safe_log(
                (2 * b.compress_fact_phase[p] + b.B[p] * (b.EoS_u + b.EoS_p))
                / (2 * b.compress_fact_phase[p] + b.B[p] * (b.EoS_u - b.EoS_p)),
                eps=1e-6,
            )
        ) / (b.bm[p] * b.EoS_p) + b._energy_internal_mol_ig(p)

    def _energy_internal_mol_ig(b, p):
        dU_form = {}
        for j in b.params.component_list:
            ele_comp = b.params.get_component(j).config.elemental_composition
            if ele_comp is None:
                raise ConfigurationError(
                    "{} calculation of internal energy requires elemental "
                    "composition of all species. Please set this using the "
                    "elemental_composition argument in the component "
                    "declaration ({}).".format(b.name, j)
                )

            delta_n = 0
            for e, s in ele_comp.items():
                # Check for any element which is vapor at standard state
                if e in ["He", "Ne", "Ar", "Kr", "Xe", "Ra"]:
                    delta_n += -s
                elif e in ["F", "Cl", "H", "N", "O"]:
                    delta_n += -s / 2  # These are diatomic at standard state

            if p == "Vap":
                delta_n += 1  # One mole of gaseous compound is formed
            dU_form[j] = delta_n * const.gas_constant * b.params.temperature_ref

        return sum(
            b.mole_frac_phase_comp[p, j]
            * (
                b._enth_mol_comp_ig(j)
                - const.gas_constant * (b.temperature - b.params.temperature_ref)
                + b.params.enth_mol_form_ref[j]
                + dU_form[j]
            )
            for j in b.params.component_list
        )

    def _enth_mol_cubic(b, p):
        # Derived from equation on pg. 120 in Properties of Gases and Liquids
        return (
            (b.temperature * b.dadT[p] - b.am[p])
            * safe_log(
                (2 * b.compress_fact_phase[p] + b.B[p] * (b.EoS_u + b.EoS_p))
                / (2 * b.compress_fact_phase[p] + b.B[p] * (b.EoS_u - b.EoS_p)),
                eps=1e-6,
            )
            + const.gas_constant
            * b.temperature
            * (b.compress_fact_phase[p] - 1)
            * b.bm[p]
            * b.EoS_p
        ) / (b.bm[p] * b.EoS_p) + b._enth_mol_ig(p)

    def _enth_mol_ig(b, p):
        return sum(
            b.mole_frac_phase_comp[p, j]
            * (b._enth_mol_comp_ig(j) + b.params.enth_mol_form_ref[j])
            for j in b.params.component_list
        )

    def _entr_mol_cubic(b, p):
        # See pg. 102 in Properties of Gases and Liquids
        return (
            const.gas_constant
            * safe_log(
                (b.compress_fact_phase[p] - b.B[p]) / b.compress_fact_phase[p], eps=1e-6
            )
            * b.bm[p]
            * b.EoS_p
            + const.gas_constant
            * safe_log(
                b.compress_fact_phase[p] * b.params.pressure_ref / b.pressure, eps=1e-6
            )
            * b.bm[p]
            * b.EoS_p
            + b.dadT[p]
            * safe_log(
                (2 * b.compress_fact_phase[p] + b.B[p] * (b.EoS_u + b.EoS_p))
                / (2 * b.compress_fact_phase[p] + b.B[p] * (b.EoS_u - b.EoS_p)),
                eps=1e-6,
            )
        ) / (b.bm[p] * b.EoS_p) + b._entr_mol_ig(p)

    def _entr_mol_ig(b, p):
        return sum(
            b.mole_frac_phase_comp[p, j]
            * (b._entr_mol_comp_ig(j) + b.params.entr_mol_form_ref[j])
            for j in b.params.component_list
        )

    # -----------------------------------------------------------------------------
    # Pure component properties
    def _enth_mol_comp_ig(b, j):
        return (
            (b.params.cp_mol_ig_comp_coeff_4[j] / 4)
            * (b.temperature**4 - b.params.temperature_ref**4)
            + (b.params.cp_mol_ig_comp_coeff_3[j] / 3)
            * (b.temperature**3 - b.params.temperature_ref**3)
            + (b.params.cp_mol_ig_comp_coeff_2[j] / 2)
            * (b.temperature**2 - b.params.temperature_ref**2)
            + b.params.cp_mol_ig_comp_coeff_1[j]
            * (b.temperature - b.params.temperature_ref)
        )

    def _entr_mol_comp_ig(b, j):
        return (
            (b.params.cp_mol_ig_comp_coeff_4[j] / 3)
            * (b.temperature**3 - b.params.temperature_ref**3)
            + (b.params.cp_mol_ig_comp_coeff_3[j] / 2)
            * (b.temperature**2 - b.params.temperature_ref**2)
            + b.params.cp_mol_ig_comp_coeff_2[j]
            * (b.temperature - b.params.temperature_ref)
            + b.params.cp_mol_ig_comp_coeff_1[j]
            * log(b.temperature / b.params.temperature_ref)
        )

    def calculate_scaling_factors(self):
        # Get default scale factors and do calculations from base classes
        super().calculate_scaling_factors()

        phases = self.params.config.valid_phase
        is_two_phase = len(phases) == 2  # possibly {Liq}, {Vap}, or {Liq, Vap}
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
            s = iscale.get_scaling_factor(self.flow_mol, default=1, warning=True)
            iscale.constraint_scaling_transform(
                self.total_flow_balance, s, overwrite=False
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
