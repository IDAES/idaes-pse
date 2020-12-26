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
GEN1 vapor property -package
This property package provides the necessary constraints or expressions for the
vapor phase properties of  amine-based(MEA) scrubbing of CO2 acid gas.
The first generation (GEN1) MEA model uses the Enhancement factor calculation.

For absorption process
   vapor components: Carbondioxide (CO2), Water (H2O),Oxygen(O2),Nitrogen(N2)
For stripping process
   vapor components: Carbondioxide (CO2), Water (H2O)
"""

# Changes the divide behavior to avoid integer division
from __future__ import division

# Import Pyomo libraries
from pyomo.environ import (Constraint, Expression, Param,
                          PositiveReals, Reals, NonNegativeReals,
                          value, Var, sqrt, units as pyunits)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.config import  ConfigValue, In
from pyomo.opt import SolverFactory

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
from idaes.core.util.misc import add_object_reference
from idaes.core.util.model_statistics  import (degrees_of_freedom,
                                              number_unfixed_variables)
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
        default='Absorber',
        domain=In(['Absorber', 'Stripper','Regenerator']),
        description="Flag indicating the type of  process",
        doc="""Flag indicating either absoprtion or stripping process.
             Hence, Absorber process has O2 and N2 in the vapor phase,while
             the stripping process does not """))

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

        if (self.config.process_type == 'Regenerator' or
           self.config.process_type == 'Stripper'):
                   self.CO2 = Component()
                   self.H2O = Component()
        elif self.config.process_type == 'Absorber':
                   self.CO2 = Component()
                   self.H2O = Component()
                   self.O2  = Component()
                   self.N2  = Component()

        self._state_block_class = VaporStateBlock

        # Create Phase object
        self.Vap = VaporPhase()

        # Thermodynamic reference state
        self.pressure_ref = Param(within=PositiveReals,
                                  mutable=True,
                                  default=101325,
                                  doc='Reference pressure [Pa]')
        self.temperature_ref = Param(default=298.15,
                    doc='Thermodynamic Reference Temperature [K]')

        # Mol. weights of vapor component - units = kg/mol.
        mw_comp_dict_all = {'CO2':0.04401,'H2O':0.01802,'O2':0.032,'N2':0.02801}
        mw_comp_dict ={}
        for i in self.component_list:
          mw_comp_dict[i] = mw_comp_dict_all[i]

        self.mw = Param(
                self.component_list,
                mutable=False,
                initialize= mw_comp_dict,
                doc="Molecular weights of vapor components [kg/mol]")

        #from Van Ness, J. Smith (appendix C):  CO2,N2,O2,WATER
        # cpig/R = A + BT + CT^-2
        #unit depends on unit of gas constant(J/mol.K)
        cp_param_dict_all = {
                        ('O2', 1): 3.639*CONST.gas_constant,
                        ('O2', 2): 0.506E-3*CONST.gas_constant,
                        ('O2', 3): -0.227E5*CONST.gas_constant,
                        ('CO2', 1): 5.457*CONST.gas_constant,
                        ('CO2', 2): 1.045E-3*CONST.gas_constant,
                        ('CO2', 3): -1.157E5*CONST.gas_constant,
                        ('H2O', 1): 3.47*CONST.gas_constant,
                        ('H2O', 2): 1.45E-3*CONST.gas_constant,
                        ('H2O', 3): 0.121E5*CONST.gas_constant,
                        ('N2', 1): 3.28*CONST.gas_constant,
                        ('N2', 2): 0.593E-3*CONST.gas_constant,
                        ('N2', 3): 0.04E5*CONST.gas_constant
                         }
        cp_param_dict ={}
        for i in self.component_list:
          for j in[1,2,3]:
            cp_param_dict[i,j] = cp_param_dict_all[i,j]

        self.cp_param = Param(self.component_list,
                               range(1, 4),
                               mutable=False,
                               initialize=cp_param_dict,
                               doc="ideal gas heat capacity parameters")


        #Viscosity constants
        #CO2 & H2O
        # calculated from  :C1*T^(C2)/(1+C3/T)
        # Reference: Perry and Green Handbook; McGraw Hill,8th edition 2008
        #O2 & N2
        #calculated from Sutherland Formular
        # C1*(C2 + C3)/(T+C3)*(T/C2)^1.5:
        # constants C1,C2,C3 in sutherlands' formular are:
        #C1 = vis_d_ref
        #C2 = temperature_ref
        #C3 = sutherland constant
        visc_d_param_dict_all = {
                         ('N2',1):0.01781e-3,
                         ('N2',2):300.55,
                         ('N2',3):111,
                         ('O2',1):0.02018e-3,
                         ('O2',2):292.25,
                         ('O2',3):127,
                         ('CO2',1):2.148e-6,
                         ('CO2',2):0.46,
                         ('CO2',3):290,
                         ('H2O',1):1.7096e-8,
                         ('H2O',2):1.1146,
                         ('H2O',3):0.0
                         }
        visc_d_param_dict ={}
        for i in self.component_list:
          for j in[1,2,3]:
            visc_d_param_dict[i,j] = visc_d_param_dict_all[i,j]
        self.visc_d_param = Param(self.component_list,
                               range(1, 4),
                               mutable=True,
                               initialize=visc_d_param_dict,
                               doc="Dynamic viscosity constants")

        #Thermal conductivity constants -
        #Reference: Perry and Green Handbook; McGraw Hill, 8th edition 2008
        therm_cond_param_dict_all = {('N2',1):0.000331,
                                 ('N2',2):0.7722,
                                 ('N2',3):16.323,
                                 ('N2',4):373.72,
                                 ('CO2',1):3.69,
                                 ('CO2',2):-0.3838,
                                 ('CO2',3):964,
                                 ('CO2',4):1.86e6,
                                 ('H2O',1):6.204e-6,
                                 ('H2O',2):1.3973,
                                 ('H2O',3):0,
                                 ('H2O',4):0,
                                 ('O2',1):0.00045,
                                 ('O2',2):0.7456,
                                 ('O2',3):56.699,
                                 ('O2',4):0.0}
        therm_cond_param_dict ={}
        for i in self.component_list:
          for j in[1,2,3,4]:
            therm_cond_param_dict[i,j] = therm_cond_param_dict_all[i,j]
        self.therm_cond_param = Param(self.component_list,
                               range(1, 5),
                               mutable=True,
                               initialize=therm_cond_param_dict,
                               doc="Thermal conductivity constants")

        #Diffusion Coefficient(binary) constants -
        #Diffusion volumes in Fuller-Schettler-Giddings correlation
        #for estimating binary diffusivities of components in vapor phase
        #Reference: Table 3.1 pp 71 Seader Henley (2006)
        diffus_binary_param_dict_all = {
                                 'N2' :18.5,
                                 'CO2':26.7,
                                 'H2O':13.1,
                                 'O2' :16.3}
        diffus_binary_param_dict ={}
        for i in self.component_list:
          diffus_binary_param_dict[i] = diffus_binary_param_dict_all[i]

        self.diffus_binary_param = \
            Param(self.component_list,
                  initialize=diffus_binary_param_dict,
                  doc="Diffusion volume parameter for binary diffusivity")

        # Transport  parameters
        #reference:  Billet and Schultes, 1999
        packing_dict = {
        'mellapak_250Y':{
                'a':250,      # specific surface Area of packing (m2/m3)
                'S':0.017,    # Channel Angle (m)
                'eps':0.97,   # Porosity (m3/m3)
                'LpA':237,   # ratio of packing length to cross sectional area
                'Cv':0.357,  # vapor packing specific Constant
                'Cl':0.5,    # liquid packing specific Constant
                'C1':5,      #Stichmair-Bravo-Fair pressure drop parameter
                'C2':3,      #Stichmair-Bravo-Fair pressure drop parameter
                'C3':0.45}}   #Stichmair-Bravo-Fair pressure drop paramete

        self.eps_p = Param(initialize=packing_dict['mellapak_250Y']['eps'],
                           doc="Packing void space m3/m3")
        self.LpA   = Param(initialize=packing_dict['mellapak_250Y']['LpA'],
                           doc = "Packing specific wetted perimeter m/m2")
        self.S     = Param(initialize=packing_dict['mellapak_250Y']['S'],
                           doc="Parameter related to packing geometry m")
        self.a     = Param(initialize=packing_dict['mellapak_250Y']['a'],
                           doc="Packing specific surface area m2/m3")
        self.dia_hydraulic = Param(initialize=
                    4*packing_dict['mellapak_250Y']['eps']/
                      packing_dict['mellapak_250Y']['a'],
                           doc="Hydraulic diameter [m]")

        #interfacial area model parameters
        #reference: Tsai 2010,. Chinen et al, 2018
        self.ae_para = Var(initialize=1.42)
        self.ae_parb = Var(initialize=0.12)
        self.ae_para.fix()
        self.ae_parb.fix()

       #holdup model parameters: epsilon
       #reference :Chinen et al, 2018
        self.eps_para = Var(initialize=11.45)
        self.eps_parb = Var(initialize=0.6471)
        self.eps_para.fix()
        self.eps_parb.fix()

        #specific constants for volumetric mass transfer coefficients
        #reference:  Billet and Schultes, 1999
        self.Cv = Var(initialize=packing_dict['mellapak_250Y']['Cv'])
        self.Cl = Var(initialize=packing_dict['mellapak_250Y']['Cl'])
        self.Cv.fix()
        self.Cl.fix()

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
          'flow_mol'        :
                            {'method': None, 'units': 'mol/s'},
          'pressure'        :
                            {'method': None, 'units': 'Pa'},
          'temperature'     :
                            {'method': None, 'units': 'K'},
          'mole_frac'       :
                            {'method': None, 'units': None},
          'mass_frac'       :
                            {'method': None, 'units': None},
          'flow_mol_comp'   :
                            {'method': '_flow_mol_comp', 'units': 'mol/s'},
          'mw_ave'          :
                            {'method': '_mw_ave', 'units': 'kg/mol'},
          'conc'            :
                            {'method': '_conc','units': 'mol/m^3'},
          'conc_comp'       :
                            {'method': '_conc_comp','units': 'mol/m^3'},
          'dens_mass'            :
                            {'method': '_dens_mass', 'units': 'kg/m^3'},
          'cp_mol'          :
                            {'method': '_cp_mol', 'units': 'J/mol.K'},
          'cp_mol_mean'     :
                            {'method': '_cp_mol_mean', 'units': 'J/mol.K'},
          'cp_mol_comp'     :
                            {'method': '_cp_mol_comp','units': 'J/mol.K'},
          'cp_mol_comp_mean':
                            {'method': '_cp_mol_comp_mean','units': 'J/mol.K'},
          'enth_mol_mean'   :
                            {'method': '_enth_mol_mean','units': 'J/s'},
          'enth_mol_vap_density':
                            {'method': '_enth_mol_vap_density','units': 'J/m^3'},
          'diffus'          :
                            {'method': '_diffus', 'units': 'm^2/s'},
          'visc_d'          :
                            {'method': '_visc_d', 'units': 'kg/m.s'},
          'visc_d_comp'     :
                            {'method': '_visc_d_comp', 'units': 'kg/m.s'},
          'therm_cond'      :
                            {'method': '_therm_cond', 'units': 'J/m.K.s'},
          'therm_cond_comp' :
                            {'method': '_therm_cond_comp', 'units': 'J/m.K.s'}
                            })

        obj.add_default_units({'time':   pyunits.s,
                               'length': pyunits.m,
                               'mass':   pyunits.kg,
                               'amount': pyunits.mol,
                               'temperature': pyunits.K,
                               'pressure':pyunits.Pa,
                               'energy': pyunits.J,
                               'holdup': pyunits.mol})

class VaporStateBlockMethods(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """
    def initialize(blk, state_args={},
                   state_vars_fixed=False,
                   hold_state=False, outlvl=idaeslog.NOTSET,
                   solver='ipopt', optarg={'tol': 1e-8}):
        '''
        Initialization routine for property package.
        Keyword Arguments:
            state_args : Dictionary with initial guesses for the state vars
                         chosen. Note that if this method is triggered
                         through the control volume, and if initial guesses
                         were not provided at the unit model level, the
                         control volume passes the inlet values as initial
                         guess.
                         Keys for the state_args dictionary are:
                         flow_mol, temperature, pressure and mole_frac_comp
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating whcih solver to use during
                     initialization (default = "ipopt")
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
        '''

        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="properties")

        init_log.info('\nStarting Vapor phase properties initialization')

        #Fix state variables if not already fixed
        if state_vars_fixed is False:
            flags = fix_state_vars(blk, state_args)
        else:
            # Check  dof  for a square problem
            for k in blk.keys():
                if degrees_of_freedom(blk[k]) != 0:
                    raise Exception("State vars fixed but degrees of freedom "
                                    "for state block is not zero during "
                                    "initialization.")

        # Set solver options
        opt = SolverFactory(solver)
        opt.options = optarg

        # ---------------------------------------------------------------------
        # Initialise values
        for k in blk.keys():
            for j in blk[k]._params.component_list:
                if hasattr(blk[k], "cp_mol_comp_eqn"):
                    calculate_variable_from_constraint(blk[k].cp_mol_comp[j],
                                                      blk[k].cp_mol_comp_eqn[j])

                if hasattr(blk[k], "flow_mol_comp_eqn"):
                    calculate_variable_from_constraint(blk[k].flow_mol_comp[j],
                                                      blk[k].flow_mol_comp_eqn[j])

                if hasattr(blk[k], "cp_mol_comp_mean_eqn"):
                    calculate_variable_from_constraint(
                            blk[k].cp_mol_comp_mean[j],
                            blk[k].cp_mol_comp_mean_eqn[j])

            if hasattr(blk[k], "cp_mol_eqn"):
                calculate_variable_from_constraint(blk[k].cp_mol,
                                                  blk[k].cp_mol_eqn)

            if hasattr(blk[k], "cp_mol_mean_eqn"):
                calculate_variable_from_constraint(
                        blk[k].cp_mol_mean,
                        blk[k].cp_mol_mean_eqn)

            if hasattr(blk[k], "enth_mol_mean_eqn"):
                calculate_variable_from_constraint(
                        blk[k].enth_mol_mean,
                        blk[k].enth_mol_mean_eqn)

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

        #----------------------------------------------------------------------
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)

    def release_state(blk, flags, outlvl=0):
        '''
        Method to release state variables fixed during initialisation.
        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        '''
        if flags is None:
            return

        # Unfix state variables
        revert_state_vars(blk, flags)

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
        add_object_reference(self, "mw",
                             self.config.parameters.mw)

        self._make_state_vars()

    def _make_state_vars(self):
        """List the necessary state variable objects."""
        self.flow_mol = Var(initialize=1.0,
                            domain=NonNegativeReals,
                            doc='Total molar flowrate [mol/s]')

        self.mole_frac = Var(self._params.component_list,
                             domain=NonNegativeReals,
                             bounds=(0, 1),
                             initialize=1 / len(self._params.component_list),
                             doc='Component mole fractions [-]')
        self.pressure = Var(initialize=101325,
                            domain=NonNegativeReals,
                            doc='Pressure [Pa]')
        self.temperature = Var(initialize=298.15,
                               domain=NonNegativeReals,
                               doc='Temperature [K]')

    def _flow_mol_comp(self):
        self.flow_mol_comp = Var(self._params.component_list,
                    initialize=1.0,
                    domain=NonNegativeReals,
                    doc='Component molar flowrate [mol/s]')

        def rule_flow_mol_comp(b,j):
            return  b.flow_mol_comp[j] == b.mole_frac[j]*b.flow_mol

        self.flow_mol_comp_eq = Constraint(self._params.component_list,
                                rule = rule_flow_mol_comp,
                                doc ="Component molar flow in vapor phase"
                                     " [mol/s]" )

    def _mw_ave(self):
        def rule_mw_ave(b):
            return  sum(b.mole_frac[j]*b._params.mw[j]
                       for j in b._params.component_list)

        self.mw_ave = Expression(rule = rule_mw_ave,
                                doc ="Average molecular weight of  vapor phase"
                                  "components [kg/mol]" )

    def _conc(self):
        def rule_conc(b):
            return b.pressure/(CONST.gas_constant*b.temperature)

        self.conc = Expression(rule = rule_conc, doc="concentration [mol/m3]" )

    def _conc_comp(self):
        # Vapor phase component concentration
        def rule_conc_comp(b,i):
            return b.conc*b.mole_frac[i]

        self.conc_comp = Expression(self._params.component_list,
                        rule = rule_conc_comp, doc="concentration of "
                        "vapor components [mol/m3]" )

    def _dens_mass(self):
        # dens_massity
        def rule_dens_mass(b):
            return  b.mw_ave*b.conc

        self.dens_mass = Expression(rule=rule_dens_mass,doc="density [kg/m3]")

    def _cp_mol_comp(self):
        # Pure component vapour heat capacities
        self.cp_mol_comp = Var(self._params.component_list,
                               domain=Reals,
                               initialize=1.0,
                               doc="Pure component vapour heat capacities "
                               "[J/mol.K]")

        def rule_cp_mol_comp(b, j):
            return b.cp_mol_comp[j] == (
                        b._params.cp_param[j, 1] +
                        b._params.cp_param[j, 2]*(b.temperature) +
                        b._params.cp_param[j, 3]*(b.temperature)**-2)

        try:
            # Try to build constraint
            self.cp_mol_comp_eqn = Constraint(self._params.component_list,
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
                               doc="vapour heat capacities "
                               "[J/mol.K]")

        def rule_cp_mol(b):
            return b.cp_mol == sum(b.cp_mol_comp[j]*b.mole_frac[j]
                                   for j in b._params.component_list)

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
        self.cp_mol_comp_mean = Var(self._params.component_list,
                               domain=Reals,
                               initialize=1.0,
                               doc="avearge pure component vapour heat capacities "
                               "between T and T_ref [J/mol.K]")

        def rule_cp_mol_comp_mean(b, j):
            tau =  b.temperature/b._params.temperature_ref
            return b.cp_mol_comp_mean[j] == (
                   b._params.cp_param[j, 1] +
                   b._params.cp_param[j, 2]*b._params.temperature_ref*0.5*
                   (tau + 1) +
                   b._params.cp_param[j, 3]/(tau*b._params.temperature_ref**2))

        try:
            # Try to build constraint
            self.cp_mol_comp_mean__eqn = Constraint(self._params.component_list,
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
            return b.cp_mol_mean == sum(b.cp_mol_comp_mean[j]*b.mole_frac[j]
                                   for j in b._params.component_list)

        try:
            # Try to build constraint
            self.cp_mol_mean_eqn = Constraint(rule=rule_cp_mol_mean)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol_mean)
            self.del_component(self.cp_mol_mean_eqn)
            raise

    def _enth_mol_mean(self):
        # Average vapour Enthalpy btw T and T_ref
        self.enth_mol_mean = Var(domain=Reals,
                               initialize=1.0,
                               doc="Mean Vapor Enthalpy btw T and Tref "
                               "[J/s]")
        def rule_enth_mol_mean(b):
            return b.enth_mol_mean == b.cp_mol_mean*b.flow_mol*b.temperature

        try:
            # Try to build constraint
            self.enth_mol_mean_eqn = Constraint(rule=rule_enth_mol_mean)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.enth_mol_mean)
            self.del_component(self.enth_mol_mean_eqn)
            raise

    def _visc_d_comp(self):
        '''
        Dynamic viscosity of vapor components
        Sutherland formular for N2,O2
        DIPPR method for H2O,CO2
        '''
        DIPPR_list = ['H2O','CO2']
        Sutherland_list = ['N2','O2']

        def rule_visc_d_comp(b,j):
            if j in DIPPR_list:
                return b._params.visc_d_param[j,1]*\
                     b.temperature**b._params.visc_d_param[j,2]/\
                     (1 + b._params.visc_d_param[j,3]/b.temperature)
            elif j in Sutherland_list:
                return b._params.visc_d_param[j,1]*\
                    (b._params.visc_d_param[j,2] + b._params.visc_d_param[j,3])/\
                     (b.temperature + b._params.visc_d_param[j,3])*\
                     (b.temperature/b._params.visc_d_param[j,2])**1.5

        self.visc_d_comp = Expression(self._params.component_list,
                        rule = rule_visc_d_comp, doc="dynamic viscosity of "
                        "vapor components [Pa.s or kg/m.s]" )

    def _visc_d(self):
        #Vapor dynamic viscosity (Wilke,1950)
        bin_set =[]
        for i in self._params.component_list:
            for j in self._params.component_list:
                    bin_set.append((i,j))

        self.thetha_ij = Expression(bin_set,doc='diffusivity interaction parameter')
        comp = [i for i in self._params.component_list]
        o = dict()
        for (i,j) in enumerate(comp,1):
            o[i] =j

        for i in range(1,len(comp)):
            for j in range(i+1,len(comp)+1):
                self.thetha_ij[o[i],o[j]] =\
                   (1+ 2*sqrt(self.visc_d_comp[o[i]]/self.visc_d_comp[o[j]])*\
                   (self._params.mw[o[j]]/self._params.mw[o[i]])**0.25+
                   self.visc_d_comp[o[i]]/self.visc_d_comp[o[j]]*\
                   (self._params.mw[o[j]]/self._params.mw[o[i]])**0.5)/\
                   (8+8*self._params.mw[o[i]]/self._params.mw[o[j]])**0.5

                self.thetha_ij[o[j],o[i]] =\
                    self.visc_d_comp[o[j]]/self.visc_d_comp[o[i]]*\
                    self._params.mw[o[i]]/self._params.mw[o[j]]*\
                    self.thetha_ij[o[i],o[j]]

        for i in self._params.component_list:
            for j in self._params.component_list:
                  if i == j:
                        self.thetha_ij[i,j] = 1

        mu_vap = sum(self.mole_frac[i]*self.visc_d_comp[i]/\
                 sum(self.mole_frac[j]*self.thetha_ij[i,j]\
                 for j in self._params.component_list)\
                 for i in self._params.component_list)

        self.visc_d = Expression(expr = mu_vap,
                            doc="Vapor dynamic viscosity [Pa.s or kg/m.s]")

    def _therm_cond_comp(self):
        # Thermal conductivity of vapor components

        def rule_therm_cond_comp(b,i):
            return b._params.therm_cond_param[i,1] \
                        *(b.temperature**b._params.therm_cond_param[i,2]) \
                        / ((1 + (b._params.therm_cond_param[i,3]/b.temperature)) \
                        + (b._params.therm_cond_param[i,4]/(b.temperature**2)))
        self.therm_cond_comp = Expression(self._params.component_list,
                                      rule=rule_therm_cond_comp,
                                      doc = 'Vapor  component thermal'
                                      'conductivity [J/(m.K.s)]')

    def _therm_cond(self):
        """
        Thermal conductivity of vapor phase
        Wassiljewa-Mason-Saxena mixing rule(low pressure)
        """
        xy=dict() #used to label the components e.g 1->CO2,2->N2
        for (i,j) in enumerate(self._params.component_list,1):
            xy[i]=j

        k_vap = 0
        for i in range(1,len(self._params.component_list)+1):
          sumij =0
          for j in range(1,len(self._params.component_list)+1):
            Aij = (1 + (self.visc_d_comp[xy[i]]/self.visc_d_comp[xy[j]])**0.5\
                *(self._params.mw[xy[j]]/self._params.mw[xy[i]])**0.25)**2\
                *(8*(1+self._params.mw[xy[i]]/self._params.mw[xy[j]]))**-0.5
            sumij += self.mole_frac[xy[j]]*Aij
          k_vap += self.mole_frac[xy[i]]*self.therm_cond_comp[xy[i]]/sumij

        self.therm_cond = Expression(expr = k_vap,
                                      doc = 'Vapor thermal'
                                      'conductivity [J/(m.K.s)]')

    def _diffus(self):
        """
        Diffusivity of vapor phase using Fuller method
        """
        binary_set = []
        for i in self._params.component_list:
            for j in self._params.component_list:
                if i!=j and (j,i) not in binary_set:
                    binary_set.append((i,j))

        # Binary diffusivities
        def rule_diffus_binary(b,i,j):
            return 1.013e-2*b.temperature**1.75/b.pressure*\
                sqrt(1e-3*(1/b._params.mw[i]+ 1/b._params.mw[j]))/\
                (b._params.diffus_binary_param[i]**(0.333333) +
                 b._params.diffus_binary_param[j]**(0.333333))**2

        self.diffus_binary = Expression(binary_set,rule=rule_diffus_binary,
                                    doc='binary diffusion Coefficient[m^2/s]')

        def rule_diffus(blk, i):
            return (1 - self.mole_frac[i])/(
                  sum(self.mole_frac[j]/self.diffus_binary[i, j]
                  for j in self._params.component_list if (i, j) in binary_set)
                  + sum(self.mole_frac[j]/self.diffus_binary[j, i]
                 for j in self._params.component_list if (j, i) in binary_set))

        self.diffus = Expression(self._params.component_list,
                                 rule=rule_diffus,
                             doc='diffusivity of component i in vapor [m^2/s]')

    def _enth_mol_vap_density(self):
        #  molar enthalpy holdup per unit volume
        self.enth_mol_vap_density = Var(
                            domain=Reals,
                            initialize=1.0,
                            doc=' enthalpy holdup [J/m3]')
        try:
            # Try to build constraint
            self.eq_enth_mol_vap_density = Constraint(expr=
                self.enth_mol_vap_density == self.enth_mol_mean*self.conc)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.enth_mol_vap_density)
            self.del_component(self.eq_enth_mol_vap_density)
            raise

    #==========================================================================

    def get_material_flow_terms(b, p, j):
        return b.flow_mol_comp[j]

    def get_enthalpy_flow_terms(b,p):
        return b.enth_mol_mean

    def get_material_density_terms(b, p, j):
        return b.conc_comp[j]

    def get_energy_density_terms(b, p):
        return b.enth_mol_vap_density


    def define_state_vars(b):
        return {"flow_mol": b.flow_mol,
                "temperature": b.temperature,
                "pressure": b.pressure,
                "mole_frac": b.mole_frac}

    def get_material_flow_basis(b):
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

        # Check pressure bounds
        if value(blk.pressure) < blk.pressure.lb:
            _log.error('{} Pressure set below lower bound.'.format(blk.name))
        if value(blk.pressure) > blk.pressure.ub:
            _log.error('{} Pressure set above upper bound.'.format(blk.name))


def evaluate_property_package():
        from pyomo.environ import ConcreteModel, SolverFactory
        from idaes.core import FlowsheetBlock
        from idaes.core.util.model_statistics  import degrees_of_freedom

        # See if ipopt is available and set up solver
        solver = SolverFactory('ipopt')
        solver.options = {'tol': 1e-6,
                          'mu_init': 1e-8,
                          'bound_push': 1e-8}

        m = ConcreteModel()
        # Create a flowsheet object to test inlet state blocks
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Vapor properties and state block For stripping process
        m.fs.properties = VaporParameterBlock(
                          default={'process_type':'Stripper'})
        m.fs.state_block = m.fs.properties.state_block_class(
                           default={"parameters": m.fs.properties})


        #fix state variables
        m.fs.state_block.flow_mol.fix(22)
        m.fs.state_block.temperature.fix(313)
        m.fs.state_block.pressure.fix(101325)
        m.fs.state_block.mole_frac["CO2"].fix(0.8)
        m.fs.state_block.mole_frac["H2O"].fix(0.2)

        #add properties
        hasattr(m.fs.state_block, "flow_mol_comp")
        hasattr(m.fs.state_block, "mw_ave")
        hasattr(m.fs.state_block, "conc")
        hasattr(m.fs.state_block, "conc_comp")
        hasattr(m.fs.state_block, "dens_mass")
        hasattr(m.fs.state_block, "cp_mol")
        hasattr(m.fs.state_block, "cp_mol_mean")
        hasattr(m.fs.state_block, "cp_mol_comp")
        hasattr(m.fs.state_block, "cp_mol_comp_mean")
        hasattr(m.fs.state_block, "enth_mol_mean")
        hasattr(m.fs.state_block, "diffus")
        hasattr(m.fs.state_block, "visc_d_comp")
        hasattr(m.fs.state_block, "visc_d")
        hasattr(m.fs.state_block, "therm_cond")
        hasattr(m.fs.state_block, "therm_cond_comp")
        hasattr(m.fs.state_block, "enth_mol_vap_density")

        #initialize variables
        m.fs.state_block.initialize()
        #determine property valus
        opt = SolverFactory('ipopt')
        init_log = idaeslog.getInitLogger(m.fs.state_block.name, 0, tag='unit')
        solve_log = idaeslog.getSolveLogger(m.fs.state_block.name,0,tag="unit")

        init_log.info("\n===============Property Calculations================")
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                opt.solve(m.fs.state_block, tee=slc.tee)

        print('degrees of freedom = {}\n'.format(degrees_of_freedom(m)))

        print('======================Variables===============================')
        for k in m.fs.state_block.component_objects(Var,descend_into =True):
            print('\n{}'.format(k.parent_component().doc))
            for i,j in k.iteritems():
                print('{} = {:0.4g}'.format(
                               str(j).replace('fs.state_block.',''),value(j)))

        print('\n=====================Expressions============================')
        for k in m.fs.state_block.component_objects(
            Expression,descend_into =True):
            print('\n{}'.format(k.parent_component().doc))
            for i,j in k.iteritems():
                print('{} = {:0.4g}'.format(
                             str(j).replace('fs.state_block.',''),value(j)))

if __name__ == '__main__':
  evaluate_property_package()

