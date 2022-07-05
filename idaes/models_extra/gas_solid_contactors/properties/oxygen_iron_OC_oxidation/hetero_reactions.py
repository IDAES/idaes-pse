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
Reaction property package for the oxidation of an iron-based OC.
Overall reactions:
(1) O2 + 4Fe3O4 => 6Fe2O3

Equations written in this model were primarily derived from:
A. Abad, J. Adánez, F. García-Labiano, L.F. de Diego, P. Gayán, J. Celaya,
Mapping of the range of operational conditions for cu-, Fe-, and Ni-based
oxygen carriers in chemical-looping combustion,
Chem. Eng. Sci. 62 (2007) 533–549.

"""

# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    exp,
    Param,
    Reals,
    Set,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.config import ConfigBlock, ConfigValue, Bool


# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    ReactionParameterBlock,
    ReactionBlockDataBase,
    ReactionBlockBase,
)
from idaes.core.util.misc import add_object_reference
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.util.model_statistics import (
    number_unfixed_variables_in_activated_equalities,
)
from idaes.core.util.config import (
    is_state_block,
    is_physical_parameter_block,
    is_reaction_parameter_block,
)
from idaes.core.util.constants import Constants
import idaes.logger as idaeslog
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

# Some more information about this module
__author__ = "Chinedu Okoli"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("HeteroReactionParameterBlock")
class ReactionParameterData(ReactionParameterBlock):
    """
    Property Parameter Block Class

    Contains parameters and indexing sets associated with properties for
    superheated steam.

    """

    # Create Class ConfigBlock
    CONFIG = ConfigBlock()
    CONFIG.declare(
        "gas_property_package",
        ConfigValue(
            description="Reference to associated PropertyPackageParameter "
            "object for the gas phase.",
            domain=is_physical_parameter_block,
        ),
    )
    CONFIG.declare(
        "solid_property_package",
        ConfigValue(
            description="Reference to associated PropertyPackageParameter "
            "object for the solid phase.",
            domain=is_physical_parameter_block,
        ),
    )
    CONFIG.declare(
        "default_arguments",
        ConfigBlock(
            description="Default arguments to use with Property Package", implicit=True
        ),
    )

    def build(self):
        """
        Callable method for Block construction.
        """
        super(ReactionParameterBlock, self).build()

        self._reaction_block_class = ReactionBlock

        self.default_scaling_factor = {}

        # Reaction Index
        self.rate_reaction_idx = Set(initialize=["R1"])

        # Reaction Stoichiometry
        self.rate_reaction_stoichiometry = {
            ("R1", "Vap", "O2"): -1,
            ("R1", "Vap", "N2"): 0,
            ("R1", "Vap", "CO2"): 0,
            ("R1", "Vap", "H2O"): 0,
            ("R1", "Sol", "Fe2O3"): 6,
            ("R1", "Sol", "Fe3O4"): -4,
            ("R1", "Sol", "Al2O3"): 0,
        }

        # Standard Heat of Reaction - J/mol_rxn
        dh_rxn_dict = {"R1": -469.4432e3}
        self.dh_rxn = Param(
            self.rate_reaction_idx,
            initialize=dh_rxn_dict,
            doc="Heat of reaction [J/mol]",
            units=pyunits.J / pyunits.mol,
        )

        # Smoothing factor
        self.eps = Param(
            mutable=True,
            default=1e-8,
            doc="Smoothing Factor",
            units=pyunits.mol / pyunits.m**3,
        )
        # Reaction rate scale factor
        self._scale_factor_rxn = Param(
            mutable=True,
            default=1,
            doc="Scale Factor for reaction eqn." "Used to help initialization routine",
        )

        # -------------------------------------------------------------------------
        """ Reaction properties that can be estimated"""

        # Particle grain radius within OC particle
        self.grain_radius = Var(
            domain=Reals,
            initialize=2.6e-7,
            doc="Representative particle grain" "radius within OC particle [m]",
            units=pyunits.m,
        )
        self.grain_radius.fix()

        # Molar density OC particle
        self.dens_mol_sol = Var(
            domain=Reals,
            initialize=22472,
            doc="Molar density of OC particle [mol/m^3]",
            units=pyunits.mol / pyunits.m**3,
        )
        self.dens_mol_sol.fix()

        # Available volume for reaction - from EPAT report (1-ep)'
        self.a_vol = Var(
            domain=Reals,
            initialize=0.28,
            doc="Available reaction vol. per vol. of OC",
            units=pyunits.dimensionless,
        )
        self.a_vol.fix()

        # Activation Energy
        self.energy_activation = Var(
            self.rate_reaction_idx,
            domain=Reals,
            initialize=1.4e4,
            doc="Activation energy [J/mol]",
            units=pyunits.J / pyunits.mol,
        )
        self.energy_activation.fix()

        # Reaction order
        self.rxn_order = Var(
            self.rate_reaction_idx,
            domain=Reals,
            initialize=1.0,
            doc="Reaction order in gas species [-]",
            units=pyunits.dimensionless,
        )
        self.rxn_order.fix()

        # Pre-exponential factor
        self.k0_rxn = Var(
            self.rate_reaction_idx,
            domain=Reals,
            initialize=3.1e-4,
            doc="Pre-exponential factor" "[mol^(1-N_reaction)m^(3*N_reaction -2)/s]",
            units=pyunits.m / pyunits.s,
        )
        self.k0_rxn.fix()

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "k_rxn": {"method": "_k_rxn"},
                "OC_conv": {"method": "_OC_conv"},
                "OC_conv_temp": {"method": "_OC_conv_temp"},
                "reaction_rate": {"method": "_reaction_rate"},
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


class _ReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Reaction Blocks as a
    whole, rather than individual elements of indexed Reaction Blocks.
    """

    def initialize(blk, outlvl=idaeslog.NOTSET, optarg={"tol": 1e-8}, solver="ipopt"):
        """
        Initialisation routine for reaction package.

        Keyword Arguments:
            outlvl : sets output level of initialization routine
                 * 0 = Use default idaes.init logger setting
                 * 1 = Maximum output
                 * 2 = Include solver output
                 * 3 = Return solver state for each step in subroutines
                 * 4 = Return solver state for each step in routine
                 * 5 = Final initialization status and exceptions
                 * 6 = No output
            optarg : solver options dictionary object (default=None)
            solver : str indicating whcih solver to use during
                     initialization (default = "ipopt")
        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="reactions")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="reactions")

        init_log.info_high("Starting initialization")

        # TODO - Update in the future as needed
        # Get a single representative block for getting config arguments
        for k in blk.keys():
            break

        # Fix state variables if not already fixed
        # Fix state variables of the primary (solid) state block
        state_var_flags = fix_state_vars(blk[k].config.solid_state_block)

        # Fix values of secondary (gas) state block variables if not fixed,
        # as well as the solid density variable.
        # This is done to keep the initialization problem square
        Cflag = {}  # Gas concentration flag
        Dflag = {}  # Solid density flag

        for k in blk.keys():
            for j in blk[k].gas_state_ref._params.component_list:
                if blk[k].gas_state_ref.dens_mol_comp[j].fixed is True:
                    Cflag[k, j] = True
                else:
                    Cflag[k, j] = False
                    blk[k].gas_state_ref.dens_mol_comp[j].fix(
                        blk[k].gas_state_ref.dens_mol_comp[j].value
                    )
            if blk[k].solid_state_ref.dens_mass_skeletal.fixed is True:
                Dflag[k] = True
            else:
                Dflag[k] = False
                blk[k].solid_state_ref.dens_mass_skeletal.fix(
                    blk[k].solid_state_ref.dens_mass_skeletal.value
                )

        # Initialise values
        for k in blk.keys():
            if hasattr(blk[k], "OC_conv_eqn"):
                calculate_variable_from_constraint(blk[k].OC_conv, blk[k].OC_conv_eqn)

            if hasattr(blk[k], "OC_conv_temp_eqn"):
                calculate_variable_from_constraint(
                    blk[k].OC_conv_temp, blk[k].OC_conv_temp_eqn
                )

            for j in blk[k]._params.rate_reaction_idx:
                if hasattr(blk[k], "rate_constant_eqn"):
                    calculate_variable_from_constraint(
                        blk[k].k_rxn[j], blk[k].rate_constant_eqn[j]
                    )

                if hasattr(blk[k], "gen_rate_expression"):
                    calculate_variable_from_constraint(
                        blk[k].reaction_rate[j], blk[k].gen_rate_expression[j]
                    )

        # Solve property block if non-empty
        free_vars = 0
        for k in blk.keys():
            free_vars += number_unfixed_variables_in_activated_equalities(blk[k])

        if free_vars > 0:
            # Create solver
            opt = get_solver(solver, optarg)
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
        else:
            res = ""
        init_log.info_high(
            "reactions initialization complete {}.".format(idaeslog.condition(res))
        )

        # ---------------------------------------------------------------------
        # Revert state vars and other variables to pre-initialization states
        # Revert state variables of the primary (solid) state block
        revert_state_vars(blk[k].config.solid_state_block, state_var_flags)

        for k in blk.keys():
            for j in blk[k].gas_state_ref._params.component_list:
                if Cflag[k, j] is False:
                    blk[k].gas_state_ref.dens_mol_comp[j].unfix()
            if Dflag[k] is False:
                blk[k].solid_state_ref.dens_mass_skeletal.unfix()

        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="reactions")
        init_log.info_high("States released.")


@declare_process_block_class("ReactionBlock", block_class=_ReactionBlock)
class ReactionBlockData(ReactionBlockDataBase):
    """
    Heterogeneous reaction package for methane reacting with Fe2O3 based OC
    """

    # Create Class ConfigBlock
    CONFIG = ConfigBlock()
    CONFIG.declare(
        "parameters",
        ConfigValue(
            domain=is_reaction_parameter_block,
            description="""
            A reference to an instance of the Reaction Parameter
            Block associated with this property package.
            """,
        ),
    )
    CONFIG.declare(
        "solid_state_block",
        ConfigValue(
            domain=is_state_block,
            description="""
            A reference to an instance of a StateBlock for the
            solid phase with which this reaction block should be associated.
            """,
        ),
    )
    CONFIG.declare(
        "gas_state_block",
        ConfigValue(
            domain=is_state_block,
            description="""
            A reference to an instance of a StateBlock for the
            gas phase with which this reaction block should be associated.
            """,
        ),
    )
    CONFIG.declare(
        "has_equilibrium",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Equilibrium reaction construction flag",
            doc="""
        Indicates whether terms for equilibrium controlled reactions
        should be constructed,
        **default** - True.
        **Valid values:** {
        **True** - include equilibrium reaction terms,
        **False** - exclude equilibrium reaction terms.}
        """,
        ),
    )

    def build(self):
        """
        Callable method for Block construction
        """
        super(ReactionBlockDataBase, self).build()

        # Object references to the corresponding state blocks and parameters
        add_object_reference(self, "_params", self.config.parameters)

        add_object_reference(
            self, "solid_state_ref", self.config.solid_state_block[self.index()]
        )
        add_object_reference(
            self, "gas_state_ref", self.config.gas_state_block[self.index()]
        )

        # Object reference for parameters if needed by CV1D
        # Reaction stoichiometry
        add_object_reference(
            self,
            "rate_reaction_stoichiometry",
            self.config.parameters.rate_reaction_stoichiometry,
        )

        # Heat of reaction
        add_object_reference(self, "dh_rxn", self.config.parameters.dh_rxn)

    # Rate constant method
    def _k_rxn(self):
        self.k_rxn = Var(
            self._params.rate_reaction_idx,
            domain=Reals,
            initialize=1,
            doc="Rate constant [mol^(1-rxn_order) * m^(3*rxn_order -2)/s]",
            units=pyunits.m / pyunits.s,
        )

        def rate_constant_eqn(b, j):
            if j == "R1":
                return self.k_rxn[j] == (
                    self._params.k0_rxn[j]
                    * exp(
                        -self._params.energy_activation[j]
                        / (
                            pyunits.convert(
                                Constants.gas_constant,  # J/mol/K
                                to_units=pyunits.J / pyunits.mol / pyunits.K,
                            )
                            * self.solid_state_ref.temperature
                        )
                    )
                )
            else:
                return Constraint.Skip

        try:
            # Try to build constraint
            self.rate_constant_eqn = Constraint(
                self._params.rate_reaction_idx, rule=rate_constant_eqn
            )
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.k_rxn)
            self.del_component(self.rate_constant_eqn)
            raise

    # Conversion of oxygen carrier
    def _OC_conv(self):
        self.OC_conv = Var(
            domain=Reals,
            initialize=0.0,
            doc="Fraction of metal oxide converted",
            units=pyunits.dimensionless,
        )

        def OC_conv_eqn(b):
            return (
                b.OC_conv
                * (
                    b.solid_state_ref.mass_frac_comp["Fe2O3"]
                    + (
                        b.solid_state_ref._params.mw_comp["Fe2O3"]
                        / b.solid_state_ref._params.mw_comp["Fe3O4"]
                    )
                    * (
                        b._params.rate_reaction_stoichiometry["R1", "Sol", "Fe2O3"]
                        / -b._params.rate_reaction_stoichiometry["R1", "Sol", "Fe3O4"]
                    )
                    * b.solid_state_ref.mass_frac_comp["Fe3O4"]
                )
                == b.solid_state_ref.mass_frac_comp["Fe2O3"]
            )

        try:
            # Try to build constraint
            self.OC_conv_eqn = Constraint(rule=OC_conv_eqn)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.OC_conv)
            self.del_component(self.OC_conv_eqn)

    # Conversion of oxygen carrier reformulated
    def _OC_conv_temp(self):
        self.OC_conv_temp = Var(
            domain=Reals,
            initialize=1.0,
            doc="Reformulation term for" "X to help eqn scaling",
            units=pyunits.dimensionless,
        )

        def OC_conv_temp_eqn(b):
            return b.OC_conv_temp**3 == (1 - b.OC_conv) ** 2

        try:
            # Try to build constraint
            self.OC_conv_temp_eqn = Constraint(rule=OC_conv_temp_eqn)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.OC_conv_temp)
            self.del_component(self.OC_conv_temp_eqn)

    # General rate of reaction method
    def _reaction_rate(self):
        self.reaction_rate = Var(
            self._params.rate_reaction_idx,
            domain=Reals,
            initialize=0,
            doc="Gen. rate of reaction [mol_rxn/m3.s]",
            units=pyunits.mol / pyunits.m**3 / pyunits.s,
        )

        def rate_rule(b, r):
            return b.reaction_rate[r] == b._params._scale_factor_rxn * (
                b.solid_state_ref.mass_frac_comp["Fe3O4"]
                * (1 - b.solid_state_ref.particle_porosity)
                * b.solid_state_ref.dens_mass_skeletal
                * (b._params.a_vol / (b.solid_state_ref._params.mw_comp["Fe3O4"]))
                * 3
                * -b._params.rate_reaction_stoichiometry["R1", "Sol", "Fe3O4"]
                * b.k_rxn[r]
                * (
                    (
                        (b.gas_state_ref.dens_mol_comp["O2"] ** 2 + b._params.eps**2)
                        ** 0.5
                    )
                    ** b._params.rxn_order[r]
                )
                * b.OC_conv_temp
                / (b._params.dens_mol_sol * b._params.grain_radius)
                / (-b._params.rate_reaction_stoichiometry["R1", "Sol", "Fe3O4"])
            )

        try:
            # Try to build constraint
            self.gen_rate_expression = Constraint(
                self._params.rate_reaction_idx, rule=rate_rule
            )
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.reaction_rate)
            self.del_component(self.gen_rate_expression)
            raise

    def get_reaction_rate_basis(b):
        return MaterialFlowBasis.molar

    def model_check(blk):
        """
        Model checks for property block
        """
        # Check temperature bounds
        if value(blk.temperature) < blk.temperature.lb:
            _log.error("{} Temperature set below lower bound.".format(blk.name))
        if value(blk.temperature) > blk.temperature.ub:
            _log.error("{} Temperature set above upper bound.".format(blk.name))

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # Set default scaling
        def _set_default_factor(v, s):
            for i in v:
                if iscale.get_scaling_factor(v[i]) is None:
                    iscale.set_scaling_factor(v[i], s)

        _set_default_factor(self.k_rxn, 1e6)
        _set_default_factor(self.OC_conv, 1e6)
        _set_default_factor(self.OC_conv_temp, 1e3)
        _set_default_factor(self.reaction_rate, 1e4)

        if self.is_property_constructed("OC_conv_eqn"):
            iscale.constraint_scaling_transform(
                self.OC_conv_eqn,
                iscale.get_scaling_factor(self.OC_conv),
                overwrite=False,
            )

        if self.is_property_constructed("OC_conv_temp_eqn"):
            iscale.constraint_scaling_transform(
                self.OC_conv_temp_eqn,
                iscale.get_scaling_factor(self.OC_conv_temp),
                overwrite=False,
            )

        if self.is_property_constructed("rate_constant_eqn"):
            for i, c in self.rate_constant_eqn.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.k_rxn[i]), overwrite=False
                )

        if self.is_property_constructed("gen_rate_expression"):
            for i, c in self.gen_rate_expression.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.reaction_rate[i]), overwrite=False
                )
