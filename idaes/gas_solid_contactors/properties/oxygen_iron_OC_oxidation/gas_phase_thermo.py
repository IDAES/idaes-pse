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
This package provides the necessary constraints for gas phase properties for
the oxidation of a chemical looping oxygen carrier.
Components - Oxygen (O2), Nitrogen (N2), Carbon Dioxide (CO2), Water (H2O)

Equations written in this model were derived from:
(1) B.E. Poling, J.M. Prausnitz, J.P. O'connell, The Properties of Gases and
Liquids, Mcgraw-Hill, New York, 2001.
(2) National Institute of Standards and Technology, NIST Chemistry WebBook,
https://webbook.nist.gov/chemistry/ (accessed March 10, 2018).

"""

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           Param,
                           PositiveReals,
                           Reals,
                           value,
                           Var,
                           units as pyunits)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
# from pyomo.opt import SolverFactory

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        Component,
                        VaporPhase)
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            solve_indexed_blocks)
from idaes.core.util.misc import add_object_reference
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_unfixed_variables_in_activated_equalities)
import idaes.logger as idaeslog

# Some more information about this module
__author__ = "Chinedu Okoli"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("GasPhaseParameterBlock")
class PhysicalParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class
    Contains parameters and indexing sets associated with properties for
    oxidation of oxygen carrier with oxygen.
    """

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(PhysicalParameterData, self).build()

        self._state_block_class = GasPhaseStateBlock

        # Create Phase object
        self.Vap = VaporPhase()

        # Create Component objects
        self.N2 = Component()
        self.O2 = Component()
        self.CO2 = Component()
        self.H2O = Component()

        # Gas Constant
        self.gas_const = Param(within=PositiveReals,
                               default=8.314459848e-3,
                               doc='Gas Constant [kJ/mol.K]',
                               units=pyunits.kJ/pyunits.mol/pyunits.K)

        # Mol. weights of gas - units = kg/mol. ref: NIST webbook
        mw_comp_dict = {'O2': 0.032, 'N2': 0.028, 'CO2': 0.044, 'H2O': 0.018}
        self.mw_comp = Param(
                self.component_list,
                mutable=False,
                initialize=mw_comp_dict,
                doc="Molecular weights of gas components [kg/mol]",
                units=pyunits.kg/pyunits.mol)

        # Std. heat of formation of comp. - units = kJ/(mol comp) - ref: NIST
        enth_mol_form_comp_dict = {'O2': 0, 'N2': 0, 'CO2': -393.5224,
                                   'H2O': -241.8264}
        self.enth_mol_form_comp = Param(
                self.component_list,
                mutable=False,
                initialize=enth_mol_form_comp_dict,
                doc="Component molar heats of formation [kJ/mol]",
                units=pyunits.kJ/pyunits.mol)

        # Ideal gas spec. heat capacity parameters(Shomate) of
        # components - ref: NIST webbook. Shomate equations from NIST.
        # Parameters A-E are used for cp calcs while A-H are used for enthalpy
        # calc.
        # 1e3*cp_comp = A + B*T + C*T^2 + D*T^3 + E/(T^2)
        # where T = Temperature (K)/1000, and cp_comp = (kJ/mol.K)
        # H_comp = H - H(298.15) = A*T + B*T^2/2 + C*T^3/3 +
        # D*T^4/4 - E/T + F - H where T = Temp (K)/1000 and H_comp = (kJ/mol)
        cp_param_dict = {
                        ('O2', 1): 30.03235,
                        ('O2', 2): 8.772972,
                        ('O2', 3): -3.988133,
                        ('O2', 4): 0.788313,
                        ('O2', 5): -0.741599,
                        ('O2', 6): -11.32468,
                        ('O2', 7): 236.1663,
                        ('O2', 8): 0.0000,
                        ('N2', 1): 19.50583,
                        ('N2', 2): 19.88705,
                        ('N2', 3): -8.598535,
                        ('N2', 4): 1.369784,
                        ('N2', 5): 0.527601,
                        ('N2', 6): -4.935202,
                        ('N2', 7): 212.3900,
                        ('N2', 8): 0.0000,
                        ('CO2', 1): 24.9973500,
                        ('CO2', 2): 55.1869600,
                        ('CO2', 3): -33.6913700,
                        ('CO2', 4): 7.9483870,
                        ('CO2', 5): -0.1366380,
                        ('CO2', 6): -403.6075000,
                        ('CO2', 7): 228.2431000,
                        ('CO2', 8): -393.5224000,
                        ('H2O', 1): 30.0920000,
                        ('H2O', 2): 6.8325140,
                        ('H2O', 3): 6.7934350,
                        ('H2O', 4): -2.5344800,
                        ('H2O', 5): 0.0821390,
                        ('H2O', 6): -250.8810000,
                        ('H2O', 7): 223.3967000,
                        ('H2O', 8): -241.8264000
                        }
        self.cp_param = Param(self.component_list,
                              range(1, 10),
                              mutable=False,
                              initialize=cp_param_dict,
                              doc="Shomate equation heat capacity parameters")

        # Viscosity constants:
        # Reference: Perry and Green Handbook; McGraw Hill, 2008
        visc_d_param_dict = {('O2', 1): 1.101e-6, ('O2', 2): 0.5634,
                             ('O2', 3): 96.3, ('O2', 4): 0,
                             ('N2', 1): 6.5592e-7, ('N2', 2): 0.6081,
                             ('N2', 3): 54.714, ('N2', 4): 0,
                             ('CO2', 1): 2.148e-6, ('CO2', 2): 0.46,
                             ('CO2', 3): 290, ('CO2', 4): 0,
                             ('H2O', 1): 1.7096e-8, ('H2O', 2): 1.1146,
                             ('H2O', 3): 0, ('H2O', 4): 0}
        self.visc_d_param = Param(self.component_list,
                                  range(1, 10),
                                  mutable=True,
                                  initialize=visc_d_param_dict,
                                  doc="Dynamic viscosity constants")

        # Thermal conductivity constants:
        # Reference: Perry and Green Handbook; McGraw Hill, 2008
        therm_cond_param_dict = {('N2', 1): 3.3143e-4, ('N2', 2): 0.7722,
                                 ('N2', 3): 16.323, ('N2', 4): 0,
                                 ('O2', 1): 4.4994e-4, ('O2', 2): 0.7456,
                                 ('O2', 3): 56.699, ('O2', 4): 0,
                                 ('CO2', 1): 3.69, ('CO2', 2): -0.3838,
                                 ('CO2', 3): 964, ('CO2', 4): 1.86e6,
                                 ('H2O', 1): 6.204e-6, ('H2O', 2): 1.3973,
                                 ('H2O', 3): 0, ('H2O', 4): 0}
        self.therm_cond_param = Param(self.component_list,
                                      range(1, 10),
                                      mutable=True,
                                      initialize=therm_cond_param_dict,
                                      doc="Thermal conductivity constants")

        # Component diffusion volumes:
        # Ref: (1) Prop gas & liquids (2) Fuller et al. IECR, 58(5), 19, 1966
        diff_vol_param_dict = {'O2': 16.6, 'N2': 17.9,
                               'CO2': 26.9, 'H2O': 13.1}
        self.diff_vol_param = Param(self.component_list,
                                    mutable=True,
                                    initialize=diff_vol_param_dict,
                                    doc="Component diffusion volumes")

        self._eps = Param(initialize=1e-8,
                          doc="Smooth abs reformulation parameter")

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
                'flow_mol': {'method': None, 'units': 'mol/s'},
                'pressure': {'method': None, 'units': 'bar'},
                'temperature': {'method': None, 'units': 'K'},
                'mole_frac_comp': {'method': None, 'units': None},
                'mw': {'method': '_mw', 'units': 'kg/mol'},
                'cp_mol': {'method': '_cp_mol', 'units': 'kJ/mol.K'},
                'cp_mol_comp': {'method': '_cp_mol_comp',
                                'units': 'kJ/mol.K'},
                'cp_mass': {'method': '_cp_mass', 'units': 'kJ/kg.K'},
                'dens_mol': {'method': '_dens_mol',
                             'units': 'mol/m^3'},
                'dens_mol_comp': {'method': '_dens_mol_comp',
                                  'units': 'mol/m^3'},
                'dens_mass': {'method': '_dens_mass',
                              'units': 'kg/m^3'},
                'enth_mol': {'method': '_enth_mol', 'units': 'kJ/mol'},
                'enth_mol_comp': {'method': '_enth_mol_comp',
                                  'units': 'kJ/mol'},
                'visc_d': {'method': '_visc_d', 'units': 'kg/m.s'},
                'therm_cond': {'method': '_therm_cond', 'units': 'kJ/m.K.s'},
                'diffusion_comp': {'method': '_diffusion_comp',
                                   'units': 'cm2/s'}})

        obj.add_default_units({'time': pyunits.s,
                               'length': pyunits.m,
                               'mass': pyunits.kg,
                               'amount': pyunits.mol,
                               'temperature': pyunits.K})
        # def add_default_units(self, u): u (dict): Key=property, Value=units
        # def add_properties(self, p): p (dict): Key=property, Value=PropertyMetadata or equiv. dict
        # def get_derived_units(self, units):
        # obj.get_derived_units("power") = pyunits.kJ * pyunits.s ** -1


class _GasPhaseStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to State Blocks as a
    whole, rather than individual elements of indexed State Blocks.
    """
    def initialize(blk, state_args=None, hold_state=False,
                   state_vars_fixed=False, outlvl=idaeslog.NOTSET,
                   solver="ipopt", optarg={"tol": 1e-8}):
        """
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
            solver : str indicating which solver to use during
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
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="properties")

        init_log.info_high('Starting initialization')

        # Deactivate the constraints specific for non-inlet blocks i.e.
        # when defined state is False
        for k in blk.keys():
            if blk[k].config.defined_state is False:
                blk[k].sum_component_eqn.deactivate()

        # Fix state variables if not already fixed
        if state_vars_fixed is False:
            flags = fix_state_vars(blk, state_args)
        else:
            # Check when the state vars are fixed already result in dof 0
            for k in blk.keys():
                if degrees_of_freedom(blk[k]) != 0:
                    raise Exception("State vars fixed but degrees of freedom "
                                    "for state block is not zero during "
                                    "initialization.")

        # ---------------------------------------------------------------------
        # Initialise values
        for k in blk.keys():

            if hasattr(blk[k], "mw_eqn"):
                calculate_variable_from_constraint(
                            blk[k].mw,
                            blk[k].mw_eqn)

            if hasattr(blk[k], "ideal_gas"):
                calculate_variable_from_constraint(
                            blk[k].dens_mol,
                            blk[k].ideal_gas)

            if hasattr(blk[k], "dens_mass_basis"):
                calculate_variable_from_constraint(
                            blk[k].dens_mass,
                            blk[k].dens_mass_basis)

            if hasattr(blk[k], "mixture_heat_capacity_eqn"):
                calculate_variable_from_constraint(
                            blk[k].cp_mol,
                            blk[k].mixture_heat_capacity_eqn)

            if hasattr(blk[k], "cp_mass_basis"):
                calculate_variable_from_constraint(
                            blk[k].cp_mass,
                            blk[k].cp_mass_basis)

            if hasattr(blk[k], "visc_d_constraint"):
                calculate_variable_from_constraint(
                            blk[k].visc_d,
                            blk[k].visc_d_constraint)

            if hasattr(blk[k], "therm_cond_constraint"):
                calculate_variable_from_constraint(
                            blk[k].therm_cond,
                            blk[k].therm_cond_constraint)

            if hasattr(blk[k], "mixture_enthalpy_eqn"):
                calculate_variable_from_constraint(
                            blk[k].enth_mol,
                            blk[k].mixture_enthalpy_eqn)

            for j in blk[k]._params.component_list:

                if hasattr(blk[k], "comp_conc_eqn"):
                    calculate_variable_from_constraint(
                                blk[k].dens_mol_comp[j],
                                blk[k].comp_conc_eqn[j])

                if hasattr(blk[k], "diffusion_comp_constraint"):
                    calculate_variable_from_constraint(
                                blk[k].diffusion_comp[j],
                                blk[k].diffusion_comp_constraint[j])

                if hasattr(blk[k], "cp_shomate_eqn"):
                    calculate_variable_from_constraint(
                        blk[k].cp_mol_comp[j],
                        blk[k].cp_shomate_eqn[j])

                if hasattr(blk[k], "enthalpy_shomate_eqn"):
                    calculate_variable_from_constraint(
                            blk[k].enth_mol_comp[j],
                            blk[k].enthalpy_shomate_eqn[j])

        # Solve property block if non-empty
        free_vars = 0
        for k in blk.keys():
            free_vars += number_unfixed_variables_in_activated_equalities(
                blk[k])

        if free_vars > 0:
            # Create solver
            opt = get_solver(solver, optarg)
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
        else:
            res = ""
        init_log.info_high("Initialization complete {}.".format(
            idaeslog.condition(res))
        )

        # ---------------------------------------------------------------------
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)

    def release_state(blk, flags, outlvl=0):
        """
        Method to relase state variables fixed during initialization.
        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of logging
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


@declare_process_block_class("GasPhaseStateBlock",
                             block_class=_GasPhaseStateBlock)
class GasPhaseStateBlockData(StateBlockData):
    """
    Property package for gas phase properties of methane combustion in CLC FR
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(GasPhaseStateBlockData, self).build()

        # Object reference for molecular weight if needed by CV1D
        # Molecular weights
        add_object_reference(self, "mw_comp",
                             self.config.parameters.mw_comp)

        """List the necessary state variable objects."""
        self.flow_mol = Var(initialize=1.0,
                            domain=Reals,
                            doc='Component molar flowrate [mol/s]',
                            units=pyunits.mol/pyunits.s)
        self.mole_frac_comp = Var(
                self._params.component_list,
                domain=Reals,
                initialize=1 / len(self._params.component_list),
                doc='State component mole fractions [-]',
                units=pyunits.mol/pyunits.mol)
        self.pressure = Var(initialize=1.01325,
                            domain=Reals,
                            doc='State pressure [bar]',
                            units=pyunits.bar)
        self.temperature = Var(initialize=298.15,
                               domain=Reals,
                               doc='State temperature [K]',
                               units=pyunits.K)

        # Create standard constraints
        # Sum mole fractions if not inlet block
        if self.config.defined_state is False:
            def sum_component_eqn(b):
                return 1e2 == 1e2 * sum(b.mole_frac_comp[j]
                                        for j in b._params.component_list)
            self.sum_component_eqn = Constraint(rule=sum_component_eqn)

    def _mw(self):
        # Molecular weight of gas mixture
        self.mw = Var(domain=Reals,
                      initialize=1.0,
                      doc="Molecular weight of gas mixture [kg/mol]",
                      units=pyunits.kg/pyunits.mol)

        def mw_eqn(b):
            return (b.mw ==
                    sum(b.mole_frac_comp[j]*b._params.mw_comp[j]
                        for j in b._params.component_list))
        try:
            # Try to build constraint
            self.mw_eqn = Constraint(rule=mw_eqn)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.mw)
            self.del_component(self.mw_eqn)
            raise

    def _dens_mol(self):
        # Molar density
        self.dens_mol = Var(domain=Reals,
                            initialize=1.0,
                            doc="Molar density/concentration [mol/m3]",
                            units=pyunits.mol/pyunits.m**3)

        def ideal_gas(b):
            return (b.dens_mol*b._params.gas_const*b.temperature*1e-2 ==
                    b.pressure)
        try:
            # Try to build constraint
            self.ideal_gas = Constraint(rule=ideal_gas)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.dens_mol)
            self.del_component(self.ideal_gas)
            raise

    def _dens_mol_comp(self):
        # Mixture heat capacities
        self.dens_mol_comp = Var(self._params.component_list,
                                 domain=Reals,
                                 initialize=1.0,
                                 doc='Component molar concentration'
                                 '[mol/m3]',
                                 units=pyunits.mol/pyunits.m**3)

        def comp_conc_eqn(b, j):
            return (b.dens_mol_comp[j] ==
                    b.dens_mol*b.mole_frac_comp[j])
        try:
            # Try to build constraint
            self.comp_conc_eqn = Constraint(self._params.component_list,
                                            rule=comp_conc_eqn)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.dens_mol_comp)
            self.del_component(self.comp_conc_eqn)
            raise

    def _dens_mass(self):
        # Mass density
        self.dens_mass = Var(domain=Reals,
                             initialize=1.0,
                             doc="Mass density [kg/m3]",
                             units=pyunits.kg/pyunits.m**3)

        def dens_mass_basis(b):
            return b.dens_mass == b.mw*b.dens_mol
        try:
            # Try to build constraint
            self.dens_mass_basis = Constraint(rule=dens_mass_basis)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.dens_mass)
            self.del_component(self.dens_mass_basis)
            raise

    def _visc_d(self):
        # Mixture dynamic viscosity
        self.visc_d = Var(domain=Reals,
                          initialize=1e-5,
                          doc="Mixture dynamic viscosity [kg/m.s]",
                          units=pyunits.kg/pyunits.m/pyunits.s)

        def visc_d_comp(i):
            return self._params.visc_d_param[i, 1] * \
                    (self.temperature**self._params.visc_d_param[i, 2]) \
                    / ((1 + (self._params.visc_d_param[i, 3]/self.temperature))
                        + (self._params.visc_d_param[i, 4] /
                           (self.temperature**2)))

        def visc_d_constraint(b):
            return 1e6*b.visc_d == 1e6*sum(b.mole_frac_comp[i]*visc_d_comp(i)
                                           / (sum(b.mole_frac_comp[j]
                                                  * (b._params.mw_comp[j] /
                                                     b._params.mw_comp[i])**0.5
                                                  for j in
                                                  b._params.component_list))
                                           for i in
                                           b._params.component_list)
        try:
            # Try to build constraint
            self.visc_d_constraint = Constraint(rule=visc_d_constraint)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.visc_d)
            self.del_component(self.visc_d_constraint)
            raise

    def _diffusion_comp(self):
        # Component diffusion in a gas mixture - units of cm2/s to help scaling
        self.diffusion_comp = Var(self._params.component_list,
                                  domain=Reals,
                                  initialize=1e-5,
                                  doc='Component diffusion in a gas mixture'
                                  '[cm2/s]',
                                  units=pyunits.cm**2/pyunits.s)

        def D_bin(i, j):
            # 1e3 used to multiply MW to convert from kg/mol to kg/kmol
            return ((1.43e-3*(self.temperature**1.75) *
                     ((1e3 * self._params.mw_comp[i] +
                       1e3 * self._params.mw_comp[j])
                     / (2 * (1e3 * self._params.mw_comp[i]) *
                        (1e3*self._params.mw_comp[j])))**0.5)
                    / ((self.pressure)
                        * ((self._params.diff_vol_param[i]**(1/3))
                            + (self._params.diff_vol_param[j]**(1/3)))**2))

        def diffusion_comp_constraint(b, i):
            return (b.diffusion_comp[i]
                    * sum(b.mole_frac_comp[j]/D_bin(i, j)
                    for j in b._params.component_list if i != j)
                    == (1-b.mole_frac_comp[i]))
        try:
            # Try to build constraint
            self.diffusion_comp_constraint \
                = Constraint(self._params.component_list,
                             rule=diffusion_comp_constraint)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.diffusion_comp)
            self.del_component(self.diffusion_comp_constraint)
            raise

    def _therm_cond(self):
        # Thermal conductivity of gas
        self.therm_cond = Var(domain=Reals,
                              initialize=1e-5,
                              doc="Thermal conductivity of gas [kJ/m.K.s]",
                              units=pyunits.kJ/pyunits.m/pyunits.K/pyunits.s)

        def therm_cond_comp(i):
            return self._params.therm_cond_param[i, 1] \
                    * (self.temperature**self._params.therm_cond_param[i, 2]) \
                    / ((1 + (self._params.therm_cond_param[i, 3] /
                             self.temperature))
                        + (self._params.therm_cond_param[i, 4] /
                       (self.temperature**2)))

        def A_bin(i, j):
            return (1 + ((therm_cond_comp(j)/therm_cond_comp(i))**0.5)
                    * ((self._params.mw_comp[j] /
                        self._params.mw_comp[i])**0.25))**2 \
                    / (8*(1+(self._params.mw_comp[j] /
                             self._params.mw_comp[i])))**0.5

        def therm_cond_constraint(b):
            # The 1e-3 term is used as a conversion factor to a kJ basis
            return 1e6*b.therm_cond == 1e6*(1e-3) * \
                                       sum(b.mole_frac_comp[i]
                                           * therm_cond_comp(i)
                                           / (sum(b.mole_frac_comp[j] *
                                                  A_bin(i, j)**0.5
                                                  for j in
                                                  b._params.component_list))
                                           for i in
                                           b._params.component_list)
        try:
            # Try to build constraint
            self.therm_cond_constraint = Constraint(rule=therm_cond_constraint)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.therm_cond)
            self.del_component(self.therm_cond_constraint)
            raise

    def _cp_mol_comp(self):
        # Pure component vapour heat capacities
        self.cp_mol_comp = Var(self._params.component_list,
                               domain=Reals,
                               initialize=1.0,
                               doc="Pure component vapour heat capacities "
                               "[kJ/mol.K]",
                               units=pyunits.kJ/pyunits.mol/pyunits.K)

        def pure_component_cp_mol(b, j):
            return b.cp_mol_comp[j] == 1e-3*(
                        b._params.cp_param[j, 1] +
                        b._params.cp_param[j, 2]*(b.temperature*1e-3) +
                        b._params.cp_param[j, 3]*(b.temperature*1e-3)**2 +
                        b._params.cp_param[j, 4]*(b.temperature*1e-3)**3 +
                        b._params.cp_param[j, 5]/((b.temperature*1e-3)**2))
        try:
            # Try to build constraint
            self.cp_shomate_eqn = Constraint(self._params.component_list,
                                             rule=pure_component_cp_mol)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol_comp)
            self.del_component(self.cp_shomate_eqn)
            raise

    def _cp_mol(self):
        # Mixture heat capacities
        self.cp_mol = Var(domain=Reals,
                          initialize=1.0,
                          doc="Mixture heat capacity [kJ/mol.K]",
                          units=pyunits.kJ/pyunits.mol/pyunits.K)

        def cp_mol(b):
            return b.cp_mol == sum(b.cp_mol_comp[j]*b.mole_frac_comp[j]
                                   for j in b._params.component_list)
        try:
            # Try to build constraint
            self.mixture_heat_capacity_eqn = Constraint(
                                            rule=cp_mol)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol)
            self.del_component(self.mixture_heat_capacity_eqn)
            raise

    def _cp_mass(self):
        # Mixture heat capacities
        self.cp_mass = Var(domain=Reals,
                           initialize=1.0,
                           doc="Mixture heat capacity, mass-basis [kJ/kg.K]",
                           units=pyunits.kJ/pyunits.kg/pyunits.K)

        def cp_mass(b):
            return b.cp_mass*b.mw == b.cp_mol
        try:
            # Try to build constraint
            self.cp_mass_basis = Constraint(rule=cp_mass)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mass)
            self.del_component(self.cp_mass_basis)
            raise

    def _enth_mol_comp(self):
        # Pure component vapour enthalpies
        self.enth_mol_comp = Var(
                self._params.component_list,
                domain=Reals,
                initialize=1.0,
                doc="Pure component enthalpies [kJ/mol]",
                units=pyunits.kJ/pyunits.mol)

        def pure_comp_enthalpy(b, j):
            return b.enth_mol_comp[j] == (
                    b._params.cp_param[j, 1]*(b.temperature*1e-3) +
                    b._params.cp_param[j, 2]*((b.temperature*1e-3)**2)/2 +
                    b._params.cp_param[j, 3]*((b.temperature*1e-3)**3)/3 +
                    b._params.cp_param[j, 4]*((b.temperature*1e-3)**4)/4 -
                    b._params.cp_param[j, 5]/(b.temperature*1e-3) +
                    b._params.cp_param[j, 6] -
                    b._params.cp_param[j, 8])
        try:
            # Try to build constraint
            self.enthalpy_shomate_eqn = Constraint(self._params.component_list,
                                                   rule=pure_comp_enthalpy)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.enth_mol_comp)
            self.del_component(self.enthalpy_shomate_eqn)
            raise

    def _enth_mol(self):
        # Mixture molar enthalpy
        self.enth_mol = Var(
                            domain=Reals,
                            initialize=1.0,
                            doc='Mixture specific enthalpy [kJ/mol]',
                            units=pyunits.kJ/pyunits.mol)
        try:
            # Try to build constraint
            self.mixture_enthalpy_eqn = Constraint(expr=(
                        self.enth_mol == sum(self.mole_frac_comp[j] *
                                             self.enth_mol_comp[j]
                                             for j in
                                             self._params.component_list)))
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.enth_mol)
            self.del_component(self.mixture_enthalpy_eqn)
            raise

    def get_material_flow_terms(b, p, j):
        return b.flow_mol*b.mole_frac_comp[j]

    def get_enthalpy_flow_terms(b, p):
        return b.flow_mol*b.enth_mol

    def get_material_density_terms(b, p, j):
        return b.dens_mol_comp[j]

    def get_energy_density_terms(b, p):
        return b.dens_mol*b.enth_mol

    def define_state_vars(b):
        return {"flow_mol": b.flow_mol,
                "temperature": b.temperature,
                "pressure": b.pressure,
                "mole_frac_comp": b.mole_frac_comp}

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

    def default_material_balance_type(blk):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(blk):
        return EnergyBalanceType.enthalpyTotal
