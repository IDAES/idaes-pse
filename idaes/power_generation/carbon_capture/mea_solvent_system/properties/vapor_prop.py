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
First Generation (GEN 1) vapor property -package

This property package provides the necessary constraints or expressions for the
vapor phase properties of  amine-based (MEA) scrubbing of CO2 acid gas.
The GEN 1 MEA model uses the Enhancement factor calculation.

For absorption process
   Vapor components: Carbondioxide (CO2), Water (H2O), Oxygen(O2), Nitrogen(N2)
For stripping process
   Vapor components: Carbondioxide (CO2), Water (H2O)
"""

# Import Pyomo libraries
from pyomo.environ import (Constraint, Expression, Param, SolverFactory,
                           Reference, PositiveReals, Reals, NonNegativeReals,
                           value, Var, sqrt, units as pyunits)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.config import ConfigValue, In

# Import IDAES libraries
from idaes.core.util.constants import Constants as CONST
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock,
                        Component,
                        VaporPhase)
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            solve_indexed_blocks)
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_unfixed_variables)
from idaes.power_generation.carbon_capture.mea_solvent_system.unit_models.column \
    import ProcessType
from idaes.core.util import get_solver
import idaes.logger as idaeslog

__author__ = "Paul Akula, John Eslick"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("VaporParameterBlock")
class PhysicalParameterData(PhysicalParameterBlock):
    """
    Vapor Phase Property Parameter Block Class

    Contains parameters and indexing sets associated with
    vapor phase properties for amine-based scrubbing process.

    """
    # Remove zero flow components in the vapor phase
    # using config to set the vapor components according to process type

    CONFIG = PhysicalParameterBlock.CONFIG()
    CONFIG.declare("process_type", ConfigValue(
        default=ProcessType.absorber,
        domain=In(ProcessType),
        description="Flag indicating the type of  process",
        doc="""Flag indicating either absorption or stripping process.
            **default** - ProcessType.absorber.
            **Valid values:** {
            **ProcessType.absorber** - absorption process,
            **ProcessType.stripper** - stripping process.}"""))

    def build(self):
        '''
        Callable method for Block construction.
        '''
        '''
            Create Component objects
            components created using the component class are added
            to the component_list object which is used by the framework
            '''

        super(PhysicalParameterData, self).build()

        if (self.config.process_type == ProcessType.stripper):
            self.CO2 = Component()
            self.H2O = Component()
        elif (self.config.process_type == ProcessType.absorber):
            self.CO2 = Component()
            self.H2O = Component()
            self.O2 = Component()
            self.N2 = Component()


        self._state_block_class = VaporStateBlock

        # Create Phase object
        self.Vap = VaporPhase()

        # Thermodynamic reference state
        self.pressure_ref = Param(within=PositiveReals,
                                  mutable=True,
                                  default=101325,
                                  units=pyunits.Pa,
                                  doc='Reference pressure ')
        self.temperature_ref = Param(default=298.15,
                                     units=pyunits.K,
                                     doc='Thermodynamic Reference Temperature [K]')

        # Mol. weights of vapor component - units = kg/mol.
        mw_comp_dict_all = {'CO2': 0.04401,
                            'H2O': 0.01802, 'O2': 0.032, 'N2': 0.02801}
        mw_comp_dict = {}
        for i in self.component_list:
            mw_comp_dict[i] = mw_comp_dict_all[i]

        self.mw_comp = Param(
            self.component_list,
            mutable=False,
            initialize=mw_comp_dict,
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weights of vapor components")

        # from Van Ness, J. Smith (appendix C):  CO2,N2,O2,WATER
        # cpig/R = A + BT + CT^-2
        # unit depends on unit of gas constant(J/mol.K)
        cp_param_dict_all = {
            ('O2', 1): 3.639 * CONST.gas_constant,
            ('O2', 2): 0.506E-3 * CONST.gas_constant,
            ('O2', 3): -0.227E5 * CONST.gas_constant,
            ('CO2', 1): 5.457 * CONST.gas_constant,
            ('CO2', 2): 1.045E-3 * CONST.gas_constant,
            ('CO2', 3): -1.157E5 * CONST.gas_constant,
            ('H2O', 1): 3.47 * CONST.gas_constant,
            ('H2O', 2): 1.45E-3 * CONST.gas_constant,
            ('H2O', 3): 0.121E5 * CONST.gas_constant,
            ('N2', 1): 3.28 * CONST.gas_constant,
            ('N2', 2): 0.593E-3 * CONST.gas_constant,
            ('N2', 3): 0.04E5 * CONST.gas_constant
        }
        cp_param_dict = {}
        for i in self.component_list:
            for j in [1, 2, 3]:
                cp_param_dict[i, j] = cp_param_dict_all[i, j]

        self.cp_param = Param(self.component_list,
                              range(1, 4),
                              mutable=False,
                              initialize=cp_param_dict,
                              units=pyunits.J / (pyunits.mol * pyunits.K),
                              doc="Ideal gas heat capacity parameters")

        # Viscosity constants
        #CO2 & H2O
        # calculated from  :C1*T^(C2)/(1+C3/T)
        # Reference: Perry and Green Handbook; McGraw Hill,8th edition 2008
        #O2 & N2
        # calculated from Sutherland Formular
        # C1*(C2 + C3)/(T+C3)*(T/C2)^1.5:
        # constants C1,C2,C3 in sutherlands' formular are:
        #C1 = vis_d_ref
        #C2 = temperature_ref
        # C3 = sutherland constant
        visc_d_param_dict_all = {
            ('N2', 1): 0.01781e-3,
            ('N2', 2): 300.55,
            ('N2', 3): 111,
            ('O2', 1): 0.02018e-3,
            ('O2', 2): 292.25,
            ('O2', 3): 127,
            ('CO2', 1): 2.148e-6,
            ('CO2', 2): 0.46,
            ('CO2', 3): 290,
            ('H2O', 1): 1.7096e-8,
            ('H2O', 2): 1.1146,
            ('H2O', 3): 0.0
        }
        visc_d_param_dict = {}
        for i in self.component_list:
            for j in [1, 2, 3]:
                visc_d_param_dict[i, j] = visc_d_param_dict_all[i, j]
        self.visc_d_param = Param(self.component_list,
                                  range(1, 4),
                                  mutable=True,
                                  initialize=visc_d_param_dict,
                                  units=pyunits.Pa * pyunits.s,
                                  doc="Dynamic viscosity constants")

        # Thermal conductivity constants -
        # Reference: Perry and Green Handbook; McGraw Hill, 8th edition 2008
        therm_cond_param_dict_all = {('N2', 1): 0.000331,
                                     ('N2', 2): 0.7722,
                                     ('N2', 3): 16.323,
                                     ('N2', 4): 373.72,
                                     ('CO2', 1): 3.69,
                                     ('CO2', 2): -0.3838,
                                     ('CO2', 3): 964,
                                     ('CO2', 4): 1.86e6,
                                     ('H2O', 1): 6.204e-6,
                                     ('H2O', 2): 1.3973,
                                     ('H2O', 3): 0,
                                     ('H2O', 4): 0,
                                     ('O2', 1): 0.00045,
                                     ('O2', 2): 0.7456,
                                     ('O2', 3): 56.699,
                                     ('O2', 4): 0.0}
        therm_cond_param_dict = {}
        for i in self.component_list:
            for j in [1, 2, 3, 4]:
                therm_cond_param_dict[i, j] = therm_cond_param_dict_all[i, j]
        self.therm_cond_param = Param(self.component_list,
                                      range(1, 5),
                                      mutable=True,
                                      initialize=therm_cond_param_dict,
                                      units=pyunits.W / pyunits.m / pyunits.K,
                                      doc="Thermal conductivity constants")

        # Diffusion Coefficient(binary) constants -
        # Diffusion volumes in Fuller-Schettler-Giddings correlation
        # for estimating binary diffusivities of components in vapor phase
        # Reference: Table 3.1 pp 71 Seader Henley (2006)
        diffus_binary_param_dict_all = {
            'N2': 18.5,
            'CO2': 26.7,
            'H2O': 13.1,
            'O2': 16.3}
        diffus_binary_param_dict = {}
        for i in self.component_list:
            diffus_binary_param_dict[i] = diffus_binary_param_dict_all[i]

        self.diffus_binary_param = \
            Param(self.component_list,
                  initialize=diffus_binary_param_dict,
                  units=pyunits.m**2 / pyunits.s,
                  doc="Diffusion volume parameter for binary diffusivity")

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
            'flow_mol': {'method': None, 'units': 'mol/s'},
            'pressure': {'method': None, 'units': 'Pa'},
            'temperature': {'method': None, 'units': 'K'},
            'mole_frac_comp': {'method': None, 'units': None},
            'flow_mol_comp': {'method': '_flow_mol_comp', 'units': 'mol/s'},
            'mw': {'method': '_mw', 'units': 'kg/mol'},
            'conc_mol': {'method': '_conc_mol', 'units': 'mol/m^3'},
            'conc_mol_comp': {'method': '_conc_mol_comp', 'units': 'mol/m^3'},
            'dens_mass': {'method': '_dens_mass', 'units': 'kg/m^3'},
            'cp_mol': {'method': '_cp_mol', 'units': 'J/mol.K'},
            'cp_mol_mean': {'method': '_cp_mol_mean', 'units': 'J/mol.K'},
            'cp_mol_comp': {'method': '_cp_mol_comp', 'units': 'J/mol.K'},
            'cp_mol_comp_mean': {'method': '_cp_mol_comp_mean', 'units': 'J/mol.K'},
            'enth_mean': {'method': '_enth_mean', 'units': 'J/s'},
            'enth_vap_density': {'method': '_enth_vap_density', 'units': 'J/m^3'},
            'diffus': {'method': '_diffus', 'units': 'm^2/s'},
            'visc_d': {'method': '_visc_d', 'units': 'kg/m.s'},
            'visc_d_comp': {'method': '_visc_d_comp', 'units': 'kg/m.s'},
            'therm_cond': {'method': '_therm_cond', 'units': 'J/m.K.s'},
            'therm_cond_comp': {'method': '_therm_cond_comp', 'units': 'J/m.K.s'}})

        obj.add_default_units({'time': pyunits.s,
                               'length': pyunits.m,
                               'mass': pyunits.kg,
                               'amount': pyunits.mol,
                               'temperature': pyunits.K})


class VaporStateBlockMethods(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """
    def initialize(blk, state_args=None,
                   state_vars_fixed=False,
                   hold_state=False, outlvl=idaeslog.NOTSET,
                   solver=None, optarg=None):
        """
        Initialization routine for property package.

        Keyword Arguments:
          state_args : Dictionary with initial guesses for the state vars
                       chosen. Note that if this method is triggered through
                       the control volume, and if initial guesses were not
                       provided at the unit model level, the control volume
                       passes the inlet values as initial guess.Keys for the
                       state_args dictionary are: flow_mol, temperature,
                       pressure and mole_frac_comp.
          outlvl : sets output level of initialization routine
          optarg : solver options dictionary object (default=None, use
                     default solver options)
          solver : str indicating which solver to use during
                   initialization (default = None)
          hold_state :
                  flag indicating whether the initialization routine
                  should unfix any state variables fixed during initialization
                  (default=False).

                  valid options:
                    True :
                      states varaibles are not unfixed, and a dict of returned
                      containing flags for which states were fixed during
                      initialization.
                    False :
                      state variables are unfixed after initialization by
                      calling the relase_state method
        Returns:
          If hold_states is True, returns a dict containing flags for which
          states were fixed during initialization.

        """

        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="properties")

        init_log.info('Starting Vapor phase properties initialization')

        # Deactivate the constraints specific for non-inlet blocks i.e.
        # when defined state is False
        for k in blk.keys():
            if blk[k].config.defined_state is False:
                blk[k].sum_component_eqn.deactivate()

        # Fix state variables if not already fixed
        if state_vars_fixed is False:
            flags = fix_state_vars(blk, state_args)
        for k in blk.keys():
            dof = degrees_of_freedom(blk[k])
            if dof != 0:
                raise RuntimeError(
                    "{} - degrees of freedom for state block is not zero "
                    "during initialization. DoF = {}".format(blk.name, dof))

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialise values
        for k in blk.keys():
            for j in blk[k].component_list:
                if hasattr(blk[k], "cp_mol_comp_eqn"):
                    calculate_variable_from_constraint(blk[k].cp_mol_comp[j],
                                                       blk[k].cp_mol_comp_eqn[j])

                if hasattr(blk[k], "flow_mol_comp_eqn"):
                    calculate_variable_from_constraint(blk[k].flow_mol_comp[j],
                                                       blk[k].flow_mol_comp_eqn[j])

                if hasattr(blk[k], "cp_mol_comp_mean_eqn"):
                    calculate_variable_from_constraint(blk[k].cp_mol_comp_mean[j],
                                                       blk[k].cp_mol_comp_mean_eqn[j])

            if hasattr(blk[k], "cp_mol_eqn"):
                calculate_variable_from_constraint(blk[k].cp_mol,
                                                   blk[k].cp_mol_eqn)

            if hasattr(blk[k], "cp_mol_mean_eqn"):
                calculate_variable_from_constraint(blk[k].cp_mol_mean,
                                                   blk[k].cp_mol_mean_eqn)

            if hasattr(blk[k], "enth_mean_eqn"):
                calculate_variable_from_constraint(blk[k].enth_mean,
                                                   blk[k].enth_mean_eqn)

        # Solve property block if non-empty
        free_vars = 0
        for k in blk.keys():
            free_vars += number_unfixed_variables(blk[k])
        if free_vars > 0:
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
        else:
            res = ""

        init_log.info("Vapor properties initialization complete {}.".format(
            idaeslog.condition(res)))

        # ----------------------------------------------------------------------
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)

    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        """
        Method to release state variables fixed during initialisation.

        Keyword Arguments:
            flags : dict containing information of which state variables were fixed
                    during initialization, and should now be unfixed. This dict is
                    returned by initialize if hold_state=True.
            outlvl : sets output level of of logging
        """
        if flags is None:
            return

        # Unfix state variables
        revert_state_vars(blk, flags)

        # Activate state variable related constraints
        for k in blk.keys():
            if blk[k].config.defined_state is False:
                blk[k].sum_component_eqn.activate()

        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        init_log.info_high('States released.')


@declare_process_block_class("VaporStateBlock",
                             block_class=VaporStateBlockMethods)
class VaporStateBlockData(StateBlockData):
    """
    Vapor phase property package of amine-based scrubbing of CO2
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(VaporStateBlockData, self).build()

        # Object reference for molecular weight if needed by CV1D
        # Molecular weights
        self.mw_comp = Reference(self.params.mw_comp)

        self.flow_mol = Var(initialize=1.0,
                            domain=NonNegativeReals,
                            units=pyunits.mol / pyunits.s,
                            doc='Total molar flowrate')

        self.mole_frac_comp = Var(self.component_list,
                                  domain=NonNegativeReals,
                                  bounds=(0, 1),
                                  units=None,
                                  initialize=1 / len(self.component_list),
                                  doc='Component mole fractions [-]')

        self.pressure = Var(initialize=101325,
                            domain=NonNegativeReals,
                            units=pyunits.Pa,
                            doc='Pressure [Pa]')

        self.temperature = Var(initialize=298.15,
                               domain=NonNegativeReals,
                               units=pyunits.K,
                               doc='Temperature [K]')

        # Sum mole fractions if not inlet block
        if self.config.defined_state is False:
            def sum_component_eqn(b):
                return  b.flow_mol == sum(b.flow_mol_comp[j]
                                        for j in b._params.component_list)
            self.sum_component_eqn = Constraint(rule=sum_component_eqn)

    def _flow_mol_comp(self):

        def rule_flow_mol_comp(b, j):
            return b.flow_mol_comp[j] == b.mole_frac_comp[j] * b.flow_mol

        try:
            self.flow_mol_comp = Var(self.component_list,
                                     initialize=1.0,
                                     domain=NonNegativeReals,
                                     units=pyunits.mol / pyunits.s,
                                     doc='Component molar flowrate')
            self.flow_mol_comp_eq = Constraint(self.component_list,
                                               rule=rule_flow_mol_comp,
                                               doc="Component molar flow in vapor phase"
                                                   " [mol/s]")
        except AttributeError:
            self.del_component(self.flow_mol_comp)
            self.del_component(self.flow_mol_comp_eq)
            raise

    def _mw(self):
        def rule_mw(b):
            return sum(b.mole_frac_comp[j] * b.mw_comp[j]
                       for j in b.component_list)

        try:
            self.mw = Expression(rule=rule_mw,
                                 doc="Average molecular weight of  vapor phase"
                                 "[kg/mol]")
        except AttributeError:
            self.del_component(self.mw)
            raise

    def _conc_mol(self):
        def rule_conc_mol(b):
            return b.pressure / (CONST.gas_constant * b.temperature)

        try:
            self.conc_mol = Expression(rule=rule_conc_mol, doc="concentration [mol/m3]")
        except AttributeError:
            self.del_component(self.conc_mol)
            raise

    def _conc_mol_comp(self):
        # Vapor phase component concentration
        def rule_conc_mol_comp(b, i):
            return b.conc_mol * b.mole_frac_comp[i]

        try:
            self.conc_mol_comp = Expression(self.component_list,
                                            rule=rule_conc_mol_comp,
                                            doc="concentration of "
                                                "vapor components [mol/m3]")
        except AttributeError:
            self.del_component(self.conc_mol_comp)
            raise

    def _dens_mass(self):
        # dens_massity
        def rule_dens_mass(b):
            return b.mw * b.conc_mol

        try:
            self.dens_mass = Expression(rule=rule_dens_mass, doc="density [kg/m3]")
        except AttributeError:
            self.del_component(self.dens_mass)
            raise

    def _cp_mol_comp(self):
        # Pure component vapour heat capacities
        self.cp_mol_comp = Var(self.component_list,
                               domain=Reals,
                               initialize=1.0,
                               units=pyunits.J / pyunits.K / pyunits.mol,
                               doc="Pure component vapour heat capacities ")

        def rule_cp_mol_comp(b, j):
            return b.cp_mol_comp[j] == (
                   b._params.cp_param[j, 1] +
                   b._params.cp_param[j, 2] * (b.temperature) +
                   b._params.cp_param[j, 3] * (b.temperature)**-2)

        try:
            # Try to build constraint
            self.cp_mol_comp_eqn = Constraint(self.component_list,
                                              rule=rule_cp_mol_comp)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol_comp)
            self.del_component(self.cp_mol_comp_eqn)
            raise

    def _cp_mol(self):
        # vapour heat capacities
        self.cp_mol = Var(domain=Reals,
                          initialize=1.0,
                          units=pyunits.J / pyunits.K / pyunits.mol,
                          doc="vapour heat capacities ")

        def rule_cp_mol(b):
            return b.cp_mol == sum(b.cp_mol_comp[j] * b.mole_frac_comp[j]
                                   for j in b.component_list)

        try:
            # Try to build constraint
            self.cp_mol_eqn = Constraint(rule=rule_cp_mol)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol)
            self.del_component(self.cp_mol_eqn)
            raise

    def _cp_mol_comp_mean(self):
        # average Pure component vapour heat capacities btw T and T_ref
        self.cp_mol_comp_mean = Var(self.component_list,
                                    domain=Reals,
                                    initialize=1.0,
                                    units=pyunits.J / pyunits.K / pyunits.mol,
                                    doc="avearge pure component vapour heat capacities "
                                    "between T and T_ref ")

        def rule_cp_mol_comp_mean(b, j):
            tau = b.temperature / b._params.temperature_ref
            return b.cp_mol_comp_mean[j] == (
                   b._params.cp_param[j, 1] +
                   b._params.cp_param[j, 2] * b._params.temperature_ref * 0.5 *
                   (tau + 1) +
                   b._params.cp_param[j, 3] / (tau * b._params.temperature_ref**2))

        try:
            # Try to build constraint
            self.cp_mol_comp_mean__eqn = Constraint(self.component_list,
                                                    rule=rule_cp_mol_comp_mean)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol_comp_mean)
            self.del_component(self.cp_mol_comp_mean_eqn)
            raise

    def _cp_mol_mean(self):
        # Average vapour heat capacities btw T and T_ref
        self.cp_mol_mean = Var(domain=Reals,
                               initialize=1.0,
                               doc="Mean vapour heat capacities "
                               "[J/mol.K]")

        def rule_cp_mol_mean(b):
            return b.cp_mol_mean == sum(b.cp_mol_comp_mean[j] * b.mole_frac_comp[j]
                                        for j in b.component_list)

        try:
            # Try to build constraint
            self.cp_mol_mean_eqn = Constraint(rule=rule_cp_mol_mean)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol_mean)
            self.del_component(self.cp_mol_mean_eqn)
            raise

    def _enth_mean(self):
        # Average vapour Enthalpy btw T and T_ref
        self.enth_mean = Var(domain=Reals,
                                 initialize=1.0,
                                 units=pyunits.J / pyunits.s,
                                 doc="Mean Vapor Enthalpy flow btw T and Tref "
                                 "[J/s]")

        def rule_enth_mean(b):
            return b.enth_mean == b.cp_mol_mean * b.flow_mol * b.temperature

        try:
            # Try to build constraint
            self.enth_mean_eqn = Constraint(rule=rule_enth_mean)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.enth_mean)
            self.del_component(self.enth_mean_eqn)
            raise

    def _visc_d_comp(self):
        '''
        Dynamic viscosity of vapor components
        Sutherland formular for N2,O2
        DIPPR method for H2O,CO2
        '''
        DIPPR_list = ['H2O', 'CO2']
        Sutherland_list = ['N2', 'O2']

        def rule_visc_d_comp(b, j):
            if j in DIPPR_list:
                return b._params.visc_d_param[j, 1] *\
                       b.temperature**b._params.visc_d_param[j, 2] /\
                       (1 + b._params.visc_d_param[j, 3] / b.temperature)
            elif j in Sutherland_list:
                return b._params.visc_d_param[j, 1] *\
                       (b._params.visc_d_param[j, 2] + b._params.visc_d_param[j, 3]) /\
                       (b.temperature + b._params.visc_d_param[j, 3]) *\
                       (b.temperature / b._params.visc_d_param[j, 2])**1.5

        try:
            self.visc_d_comp = Expression(self.component_list,
                                          rule=rule_visc_d_comp,
                                          doc="dynamic viscosity of "
                                              "vapor components [Pa.s or kg/m.s]")
        except AttributeError:
            self.del_component(self.visc_d_comp)
            raise

    def _visc_d(self):
        # Vapor dynamic viscosity (Wilke,1950)
        bin_set = []
        for i in self.component_list:
            for j in self.component_list:
                bin_set.append((i, j))

        self.thetha_ij = Expression(
            bin_set, doc='diffusivity interaction parameter')
        comp = [i for i in self.component_list]
        o = dict()
        for (i, j) in enumerate(comp, 1):
            o[i] = j

        for i in range(1, len(comp)):
            for j in range(i + 1, len(comp) + 1):
                self.thetha_ij[o[i], o[j]] =\
                    (1 + 2 * sqrt(self.visc_d_comp[o[i]] / self.visc_d_comp[o[j]]) *
                     (self.mw_comp[o[j]] / self.mw_comp[o[i]])**0.25 +
                     self.visc_d_comp[o[i]] / self.visc_d_comp[o[j]] *
                     (self.mw_comp[o[j]] / self.mw_comp[o[i]])**0.5) /\
                    (8 + 8 * self.mw_comp[o[i]] / self.mw_comp[o[j]])**0.5

                self.thetha_ij[o[j], o[i]] =\
                    self.visc_d_comp[o[j]] / self.visc_d_comp[o[i]] *\
                    self.mw_comp[o[i]] / self.mw_comp[o[j]] *\
                    self.thetha_ij[o[i], o[j]]

        for i in self.component_list:
            for j in self.component_list:
                if i == j:
                    self.thetha_ij[i, j] = 1

        mu_vap = sum(self.mole_frac_comp[i] * self.visc_d_comp[i] /
                     sum(self.mole_frac_comp[j] * self.thetha_ij[i, j]
                         for j in self.component_list)
                     for i in self.component_list)

        try:
            self.visc_d = Expression(expr=mu_vap,
                                     doc="Vapor dynamic viscosity [Pa.s or kg/m.s]")
        except AttributeError:
            self.del_component(self.visc_d)
            raise

    def _therm_cond_comp(self):
        # Thermal conductivity of vapor components

        def rule_therm_cond_comp(b, i):
            return b._params.therm_cond_param[i, 1] *\
                   (b.temperature**b._params.therm_cond_param[i, 2])/ \
                   ((1 + (b._params.therm_cond_param[i, 3] / b.temperature)) +
                    (b._params.therm_cond_param[i, 4] / (b.temperature**2)))

        try:
            self.therm_cond_comp = Expression(self.component_list,
                                              rule=rule_therm_cond_comp,
                                              doc='Vapor  component thermal'
                                                  'conductivity [J/(m.K.s)]')
        except AttributeError:
            self.del_component(self.therm_cond_comp)
            raise

    def _therm_cond(self):
        """
        Thermal conductivity of vapor phase
        Wassiljewa-Mason-Saxena mixing rule(low pressure)
        """
        xy = dict()  # used to label the components e.g 1->CO2,2->N2
        for (i, j) in enumerate(self.component_list, 1):
            xy[i] = j

        k_vap = 0
        for i in range(1, len(self.component_list) + 1):
            sumij = 0
            for j in range(1, len(self.component_list) + 1):
                Aij = (1 + (self.visc_d_comp[xy[i]] / self.visc_d_comp[xy[j]])**0.5 *
                      (self.mw_comp[xy[j]] / self.mw_comp[xy[i]])**0.25)**2 *\
                      (8 * (1 + self.mw_comp[xy[i]] / self.mw_comp[xy[j]]))**-0.5
                sumij += self.mole_frac_comp[xy[j]] * Aij
            k_vap += self.mole_frac_comp[xy[i]] * self.therm_cond_comp[xy[i]] / sumij

        try:
            self.therm_cond = Expression(expr=k_vap,
                                         doc='Vapor thermal'
                                             'conductivity [J/(m.K.s)]')
        except AttributeError:
            self.del_component(self.therm_cond)
            raise

    def _diffus(self):
        """
        Diffusivity of vapor phase using Fuller method
        """
        binary_set = []
        for i in self.component_list:
            for j in self.component_list:
                if i != j and (j, i) not in binary_set:
                    binary_set.append((i, j))

        # Binary diffusivities
        def rule_diffus_binary(b, i, j):
            return 1.013e-2 * b.temperature**1.75 / b.pressure *\
                sqrt(1e-3 * (1 / b.mw_comp[i] + 1 / b.mw_comp[j])) /\
                (b._params.diffus_binary_param[i]**(0.333333) +
                 b._params.diffus_binary_param[j]**(0.333333))**2

        try:
            self.diffus_binary = Expression(binary_set,
                                            rule=rule_diffus_binary,
                                            doc='binary diffusion Coefficient[m^2/s]')
        except AttributeError:
            self.del_component(self.diffus_binary)
            raise

        def rule_diffus(blk, i):
            return (1 - self.mole_frac_comp[i]) / (
                sum(self.mole_frac_comp[j] / self.diffus_binary[i, j]
                    for j in self.component_list if (i, j) in binary_set) +
                sum(self.mole_frac_comp[j] / self.diffus_binary[j, i]
                    for j in self.component_list if (j, i) in binary_set))

        try:
            self.diffus = Expression(self.component_list,
                                     rule=rule_diffus,
                                     doc='diffusivity of component i in vapor [m^2/s]')
        except AttributeError:
            self.del_component(self.diffus)
            raise

    def _enth_vap_density(self):
        #  molar enthalpy holdup per unit volume
        self.enth_vap_density = Var(
            domain=Reals,
            initialize=1.0,
            doc=' enthalpy holdup [J/m3]')
        try:
            # Try to build constraint
            self.eq_enth_vap_density =\
                Constraint(expr=self.enth_vap_density ==
                           self.enth_mean * self.conc_mol)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.enth_vap_density)
            self.del_component(self.eq_enth_vap_density)
            raise

    # ==========================================================================

    def get_material_flow_terms(self, p, j):
        return self.flow_mol_comp[j]

    def get_enthalpy_flow_terms(self, p):
        return self.enth_mean

    def get_material_density_terms(self, p, j):
        return self.conc_mol_comp[j]

    def get_energy_density_terms(self, p):
        return self.enth_vap_density

    def define_state_vars(self):
        return {"flow_mol": self.flow_mol,
                "temperature": self.temperature,
                "pressure": self.pressure,
                "mole_frac_comp": self.mole_frac_comp}

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def model_check(blk):
        """
        Model checks for property block
        """
        # Check temperature bounds
        if value(blk.temperature) < blk.temperature.lb:
            _log.error('{} Temperature set below lower bound.'.format(blk.name))
        if value(blk.temperature) > blk.temperature.ub:
            _log.error('{} Temperature set above upper bound.'.format(blk.name))

        # Check pressure bounds
        if value(blk.pressure) < blk.pressure.lb:
            _log.error('{} Pressure set below lower bound.'.format(blk.name))
        if value(blk.pressure) > blk.pressure.ub:
            _log.error('{} Pressure set above upper bound.'.format(blk.name))
