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
Property package for the heterogeneous reaction of CH4 with an iron-based OC.
Overall reducer reactions for Methane combustion:
    (1) CH4 + 12Fe2O3 => 8Fe3O4 + CO2 + 2H2O
"""

# Changes the divide behavior to not do integer division
from __future__ import division

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           exp,
                           Param,
                           PositiveReals,
                           Reals,
                           Set,
                           value,
                           Var)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.opt import SolverFactory
from pyomo.common.config import ConfigValue

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        ReactionParameterBlock,
                        ReactionBlockDataBase,
                        ReactionBlockBase)
from idaes.core.util.misc import add_object_reference
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            solve_indexed_blocks)
from idaes.core.util.model_statistics import number_unfixed_variables
from idaes.core.util.config import is_state_block
import idaes.logger as idaeslog

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
    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(ReactionParameterData, self).build()

        self.reaction_block_class = ReactionBlock

        # List of valid phases in property package
        self.phase_list = Set(initialize=['Vap, Sol'])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['CO2', 'H2O', 'CH4', 'Fe2O3',
                                              'Fe3O4', 'Al2O3'])

        self.gas_component_list = Set(initialize=['CO2', 'H2O', 'CH4'])

        self.sol_component_list = Set(initialize=['Fe2O3', 'Fe3O4', 'Al2O3'])

        # Reaction Index
        self.rate_reaction_idx = Set(initialize=["R1"])

        # Reaction Stoichiometry
        self.rate_reaction_stoichiometry = {("R1", "Vap", "CH4"): -1,
                                            ("R1", "Vap", "CO2"): 1,
                                            ("R1", "Vap", "H2O"): 2,
                                            ("R1", "Sol", "Fe2O3"): -12,
                                            ("R1", "Sol", "Fe3O4"): 8,
                                            ("R1", "Sol", "Al2O3"): 0}

        # Gas Constant
        self.gas_const = Param(within=PositiveReals,
                               mutable=False,
                               default=8.314459848e-3,
                               doc='Gas Constant [kJ/mol.K]')

        # Particle grain radius within OC particle
        self.grain_radius = Param(within=PositiveReals,
                                  mutable=True,
                                  default=2.6e-7,
                                  doc='Representative particle grain'
                                  'radius within OC particle [m]')

        # Molar density OC particle
        self.dens_mol_sol = Param(within=PositiveReals,
                                  mutable=True,
                                  default=32811,
                                  doc='Molar density of OC particle [mol/m^3]')

        # Available volume for reaction - from EPAT report (1-ep)'
        self.a_vol = Param(default=0.28,
                           mutable=True,
                           doc='Available reaction vol. per vol. of OC')

        # Activation Energy
        self.energy_activation = Param(self.rate_reaction_idx,
                                       default=4.9e1,
                                       mutable=True,
                                       doc='Activation energy [kJ/mol]')

        # Reaction order
        self.rxn_order = Param(self.rate_reaction_idx,
                               default=1.3,
                               mutable=True,
                               doc='Reaction order in gas species [-]')

        # Reaction stoichiometric coefficient
        self.rxn_stoich_coeff = Param(self.rate_reaction_idx,
                                      default=12,
                                      mutable=True,
                                      doc='Reaction stoichiometric'
                                      'coefficient [-]')

        # Pre-exponential factor
        self.k0_rxn = Param(self.rate_reaction_idx,
                            default=8e-4,
                            mutable=True,
                            doc='Pre-exponential factor'
                                '[mol^(1-N_reaction)m^(3*N_reaction -2)/s]')

        # TODO - Generalize this for r equations to compute automatically
        # Standard Heat of Reaction - kJ/mol_rxn
        dh_rxn_dict = {"R1": 136.5843}
        self.dh_rxn = Param(self.rate_reaction_idx,
                            initialize=dh_rxn_dict,
                            doc="Heat of reaction [kJ/mol]")

        self._eps = Param(default=1e-8, doc='Smoothing Factor')
        self._scale_factor_rxn = Param(mutable=True, default=1,
                                       doc='Scale Factor for reaction eqn.'
                                       'Used to help initialization routine')

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
                'k_rxn': {'method': '_k_rxn',
                          'units': 'mol^(1-N_reaction)m^(3*N_reaction -2)/s]'},
                'OC_conv': {'method': "_OC_conv", 'units': None},
                'OC_conv_temp': {'method': "_OC_conv_temp", 'units': None},
                'reaction_rate': {'method': "_reaction_rate",
                                  'units': 'mol_rxn/m3.s'}
                })
        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'kg',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'kJ'})


class _ReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Reaction Blocks as a
    whole, rather than individual elements of indexed Reaction Blocks.
    """
    def initialize(blk, outlvl=idaeslog.NOTSET,
                   optarg={'tol': 1e-8}, solver='ipopt'):
        '''
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
        '''
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="reactions")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="reactions")

        init_log.info('Starting initialization')

        # TODO - Update in the future as needed
        # Get a single representative block for getting config arguments
        for k in blk.keys():
            break

        # Fix state variables if not already fixed
        state_var_flags = fix_state_vars(blk[k].config.state_block)

        # Fix values of secondary (gas) state block variables if not fixed
        # This is done to keep the initialization problem square
        Cflag = {}  # Gas concentration flag

        for k in blk.keys():
            for j in blk[k]._params.gas_component_list:
                if blk[k].gas_state_ref.dens_mole_comp_vap[j].fixed is True:
                    Cflag[k, j] = True
                else:
                    Cflag[k, j] = False
                    blk[k].gas_state_ref.dens_mole_comp_vap[j].fix(
                            blk[k].gas_state_ref.dens_mole_comp_vap[j].value)

        # Set solver options
        opt = SolverFactory(solver)
        opt.options = optarg

        # Initialise values
        for k in blk.keys():
            if hasattr(blk[k], "OC_conv_eqn"):
                calculate_variable_from_constraint(
                            blk[k].OC_conv,
                            blk[k].OC_conv_eqn)

            if hasattr(blk[k], "OC_conv_temp_eqn"):
                calculate_variable_from_constraint(
                            blk[k].OC_conv_temp,
                            blk[k].OC_conv_temp_eqn)

            for j in blk[k]._params.rate_reaction_idx:
                if hasattr(blk[k], "rate_constant_eqn"):
                    calculate_variable_from_constraint(
                                blk[k].k_rxn[j],
                                blk[k].rate_constant_eqn[j])

                if hasattr(blk[k], "gen_rate_expression"):
                    calculate_variable_from_constraint(
                            blk[k].reaction_rate[j],
                            blk[k].gen_rate_expression[j])

        # Solve property block if non-empty
        free_vars = 0
        for k in blk.keys():
            free_vars += number_unfixed_variables(blk[k])

        if free_vars > 0:
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
        else:
            res = ""
        init_log.info("reactions initialization complete {}.".format(
            idaeslog.condition(res))
                        )

        # ---------------------------------------------------------------------
        # Revert state vars and other variables to pre-initialization states
        revert_state_vars(blk[k].config.state_block, state_var_flags)

        for k in blk.keys():
            for j in blk[k]._params.gas_component_list:
                if Cflag[k, j] is False:
                    blk[k].gas_state_ref.dens_mole_comp_vap[j].unfix()

        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="reactions")
        init_log.info_high('States released.')

@declare_process_block_class("ReactionBlock",
                             block_class=_ReactionBlock)
class ReactionBlockData(ReactionBlockDataBase):
    """
    Reaction package for methane reacting with Fe2O3 based OC
    """
    CONFIG = ReactionBlockDataBase.CONFIG()
    CONFIG.declare("gas_state_block", ConfigValue(
            default=None,
            domain=is_state_block,
            description="Additional state block for heterogeneous reactions",
            doc="""Additional State block object used to define gas property
calculations
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a StateBlock object"""))

    def build(self):
        """
        Callable method for Block construction
        """
        super(ReactionBlockData, self).build()

        # Object reference to the corresponding gas state block
        add_object_reference(self,
                             "gas_state_ref",
                             self.config.gas_state_block[self.index()])

        # Object reference for parameters if needed by CV1D
        # Reaction stoichiometry
        add_object_reference(
                self,
                "rate_reaction_stoichiometry",
                self.config.parameters.rate_reaction_stoichiometry)

        # Heat of reaction
        add_object_reference(
                self,
                "dh_rxn",
                self.config.parameters.dh_rxn)

    # Rate constant method
    def _k_rxn(self):
        self.k_rxn = Var(self._params.rate_reaction_idx,
                         domain=Reals,
                         initialize=1,
                         doc='Rate constant '
                         '[mol^(1-N_reaction)m^(3*N_reaction -2)/s]')

        def rate_constant_eqn(b, j):
            if j == 'R1':
                return 1e6 * self.k_rxn[j] == \
                        1e6 * (self._params.k0_rxn[j] *
                               exp(-self._params.energy_activation[j] /
                                   (self._params.gas_const *
                                    self.state_ref.temperature)))
            else:
                return Constraint.Skip
        try:
            # Try to build constraint
            self.rate_constant_eqn = Constraint(self._params.rate_reaction_idx,
                                                rule=rate_constant_eqn)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.k_rxn)
            self.del_component(self.rate_constant_eqn)
            raise

    # Conversion of oxygen carrier
    def _OC_conv(self):
        self.OC_conv = Var(domain=Reals, initialize=0.0,
                           doc='Fraction of metal oxide converted')

        def OC_conv_eqn(b):
            return 1e6 * b.OC_conv * \
                   (b.state_ref.mass_frac['Fe3O4'] +
                    (b.state_ref._params.mw['Fe3O4'] /
                        b.state_ref._params.mw['Fe2O3']) *
                    (b._params.rate_reaction_stoichiometry
                       ['R1', 'Sol', 'Fe3O4']
                       / -b._params.rate_reaction_stoichiometry
                       ['R1', 'Sol', 'Fe2O3']) *
                    b.state_ref.mass_frac['Fe2O3']) == \
                   1e6 * b.state_ref.mass_frac['Fe3O4']
        try:
            # Try to build constraint
            self.OC_conv_eqn = Constraint(rule=OC_conv_eqn)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.OC_conv)
            self.del_component(self.OC_conv_eqn)

    # Conversion of oxygen carrier reformulated
    def _OC_conv_temp(self):
        self.OC_conv_temp = Var(domain=Reals, initialize=1.0,
                                doc='Reformulation term for'
                                    'X to help eqn scaling')

        def OC_conv_temp_eqn(b):
            return 1e3*b.OC_conv_temp**3 == 1e3*(1-b.OC_conv)**2
        try:
            # Try to build constraint
            self.OC_conv_temp_eqn = Constraint(rule=OC_conv_temp_eqn)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.OC_conv_temp)
            self.del_component(self.OC_conv_temp_eqn)

    # General rate of reaction method
    def _reaction_rate(self):
        self.reaction_rate = Var(self._params.rate_reaction_idx,
                                 domain=Reals,
                                 initialize=0,
                                 doc="Gen. rate of reaction [mol_rxn/m3.s]")

        def rate_rule(b, r):
            return b.reaction_rate[r]*1e4 == b._params._scale_factor_rxn*1e4*(
                b.state_ref.mass_frac['Fe2O3'] *
                b.state_ref._params.dens_mass_sol_fresh *
                (b._params.a_vol/(b.state_ref._params.mw['Fe2O3'])) *
                3*b._params.rxn_stoich_coeff[r]*b.k_rxn[r] *
                (((b.gas_state_ref.dens_mole_comp_vap['CH4']**2 +
                  b._params._eps**2)**0.5) **
                 b._params.rxn_order[r]) *
                b.OC_conv_temp/(b._params.dens_mol_sol *
                                b._params.grain_radius) /
                (-b._params.rate_reaction_stoichiometry['R1', 'Sol', 'Fe2O3']))
        try:
            # Try to build constraint
            self.gen_rate_expression = Constraint(
                                            self._params.rate_reaction_idx,
                                            rule=rate_rule)
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
            _log.error('{} Temperature set below lower bound.'
                       .format(blk.name))
        if value(blk.temperature) > blk.temperature.ub:
            _log.error('{} Temperature set above upper bound.'
                       .format(blk.name))
