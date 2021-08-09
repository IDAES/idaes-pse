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
This package provides the necessary constraints for gas phase properties for
the CLC of methane
Components - Methane (CH4), Carbon Dioxide (CO2), Water (H2O)
"""

# Changes the divide behavior to not do integer division
from __future__ import division

# Import Python libraries
# import logging

# Import Pyomo libraries
from pyomo.environ import Constraint, Param, \
                          PositiveReals, Reals, \
                          Set, value, Var
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.opt import SolverFactory

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock)
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            solve_indexed_blocks)
from idaes.core.util.misc import add_object_reference
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_unfixed_variables)
import idaes.logger as idaeslog

# Some more information about this module
__author__ = "Chinedu Okoli"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Gas_Phase_Thermo_ParameterBlock")
class PhysicalParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class

    Contains parameters and indexing sets associated with properties for
    methane CLC.

    """
    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(PhysicalParameterData, self).build()

        self.state_block_class = Gas_Phase_Thermo_StateBlock

        # List of valid phases in property package
        self.phase_list = Set(initialize=['Vap'])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['CH4', 'CO2', 'H2O'])

        # List of components in each phase (optional)
        self.phase_comp = {"Vap": self.component_list}

        # Thermodynamic reference state
        self.pressure_ref = Param(within=PositiveReals,
                                  mutable=True,
                                  default=1.01325,
                                  doc='Reference pressure [bar]')
        self.temperature_ref = Param(within=PositiveReals,
                                     mutable=True,
                                     default=298.15,
                                     doc='Thermodynamic Reference'
                                     'Temperature [K]')

        # Gas Constant
        self.gas_const = Param(within=PositiveReals,
                               default=8.314459848e-3,
                               doc='Gas Constant [kJ/mol.K]')

        # Mol. weights of gas - units = kg/mol. ref: NIST webbook
        mw_comp_dict = {'CH4': 0.016, 'CO2': 0.044, 'H2O': 0.018}
        self.mw = Param(
                self.component_list,
                mutable=False,
                initialize=mw_comp_dict,
                doc="Molecular weights of gas components [kg/mol]")

        # Std. heat of formation of comp. - units = kJ/(mol comp) - ref: NIST
        enth_mol_form_comp_dict = {'CH4': -74.8731, 'CO2': -393.5224,
                                   'H2O': -241.8264}
        self.enth_mol_form_comp = Param(
                self.component_list,
                mutable=False,
                initialize=enth_mol_form_comp_dict,
                doc="Component molar heats of formation [kJ/mol]")

        # Ideal gas spec. heat capacity parameters(Shomate) of
        # components - ref: NIST webbook. Shomate equations from NIST.
        # Parameters A-E are used for cp calcs while A-H are used for enthalpy
        # calc.
        # 1e3*cp_comp = A + B*T + C*T^2 + D*T^3 + E/(T^2)
        # where T = Temperature (K)/1000, and cp_comp = (kJ/mol.K)
        # H_comp = H - H(298.15) = A*T + B*T^2/2 + C*T^3/3 +
        # D*T^4/4 - E/T + F - H where T = Temp (K)/1000 and H_comp = (kJ/mol)
        cp_param_dict = {
                        ('CH4', 1): -0.7030290,
                        ('CH4', 2): 108.4773000,
                        ('CH4', 3): -42.5215700,
                        ('CH4', 4): 5.8627880,
                        ('CH4', 5): 0.6785650,
                        ('CH4', 6): -76.8437600,
                        ('CH4', 7): 158.7163000,
                        ('CH4', 8): -74.8731000,
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
        visc_d_param_dict = {('CH4', 1): 5.2546e-7, ('CH4', 2): 0.59006,
                             ('CH4', 3): 105.67, ('CH4', 4): 0,
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
        therm_cond_param_dict = {('CH4', 1): 8.3983e-6, ('CH4', 2): 1.4268,
                                 ('CH4', 3): -49.654, ('CH4', 4): 0,
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
        diff_vol_param_dict = {'CH4': 24.42, 'CO2': 26.9, 'H2O': 13.1}
        self.diff_vol_param = Param(self.component_list,
                                    mutable=True,
                                    initialize=diff_vol_param_dict,
                                    doc="Component diffusion volumes")

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
                'flow_mol': {'method': None, 'units': 'mol/s'},
                'pressure': {'method': None, 'units': 'bar'},
                'temperature': {'method': None, 'units': 'K'},
                'mole_frac': {'method': None, 'units': None},
                'mw_gas': {'method': '_mw_gas', 'units': 'kg/mol'},
                'cp_mol': {'method': '_cp_mol', 'units': 'kJ/mol.K'},
                'cp_mol_comp': {'method': '_cp_mol_comp',
                                'units': 'kJ/mol.K'},
                'cp_mass': {'method': '_cp_mass', 'units': 'kJ/kg.K'},
                'dens_mole_vap': {'method': '_dens_mole_vap',
                                  'units': 'mol/m^3'},
                'dens_mole_comp_vap': {'method': '_dens_mole_comp_vap',
                                       'units': 'mol/m^3'},
                'dens_mass_vap': {'method': '_dens_mass_vap',
                                  'units': 'kg/m^3'},
                'enth_mol': {'method': '_enth_mol', 'units': 'kJ/mol'},
                'enth_mol_vap_comp': {'method': '_enth_mol_vap_comp',
                                      'units': 'kJ/mol'},
                'visc_d': {'method': '_visc_d', 'units': 'kg/m.s'},
                'therm_cond': {'method': '_therm_cond', 'units': 'kJ/m.K.s'},
                'diffusion_comp': {'method': '_diffusion_comp',
                                   'units': 'cm2/s'}})

        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'kg',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'kJ',
                               'holdup': 'mol'})


class _Gas_Phase_Thermo_StateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
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
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="properties")

        init_log.info('Starting initialization')

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

        # Set solver options
        opt = SolverFactory(solver)
        opt.options = optarg

        # ---------------------------------------------------------------------
        # Initialise values
        for k in blk.keys():
            for j in blk[k]._params.component_list:
                if hasattr(blk[k], "diffusion_comp_constraint"):
                    calculate_variable_from_constraint(
                                blk[k].diffusion_comp[j],
                                blk[k].diffusion_comp_constraint[j])

                if hasattr(blk[k], "cp_shomate_eqn"):
                    calculate_variable_from_constraint(blk[k].cp_mol_comp[j],
                                                       blk[k].cp_shomate_eqn[j]
                                                       )

                if hasattr(blk[k], "enthalpy_shomate_eqn"):
                    calculate_variable_from_constraint(
                            blk[k].enth_mol_vap_comp[j],
                            blk[k].enthalpy_shomate_eqn[j])

            if hasattr(blk[k], "mw_gas_eqn"):
                calculate_variable_from_constraint(
                            blk[k].mw_gas,
                            blk[k].mw_gas_eqn)

            if hasattr(blk[k], "ideal_gas"):
                calculate_variable_from_constraint(
                            blk[k].dens_mole_vap,
                            blk[k].ideal_gas)

            if hasattr(blk[k], "comp_conc_eqn"):
                calculate_variable_from_constraint(
                            blk[k].dens_mole_comp_vap[j],
                            blk[k].comp_conc_eqn[j])

            if hasattr(blk[k], "dens_mass_basis"):
                calculate_variable_from_constraint(
                            blk[k].dens_mass_vap,
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

            if hasattr(blk[k], "comp_conc_eqn"):
                calculate_variable_from_constraint(
                            blk[k].dens_mole_comp_vap[j],
                            blk[k].comp_conc_eqn[j])

        # Solve property block if non-empty
        free_vars = 0
        for k in blk.keys():
            free_vars += number_unfixed_variables(blk[k])

        if free_vars > 0:
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
        else:
            res = ""
        init_log.info("Gas properties initialization complete {}.".format(
            idaeslog.condition(res))
        )

        # ---------------------------------------------------------------------
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)

        init_log.info("Initialization complete.")

    def release_state(blk, flags, outlvl=0):
        """
        Method to relase state variables fixed during initialization.
        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
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

@declare_process_block_class("Gas_Phase_Thermo_StateBlock",
                             block_class=_Gas_Phase_Thermo_StateBlock)
class Gas_Phase_Thermo_StateBlockData(StateBlockData):
    """
    Property package for gas phase properties of methane combustion in CLC FR
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(Gas_Phase_Thermo_StateBlockData, self).build()

        # Object reference for molecular weight if needed by CV1D
        # Molecular weights
        add_object_reference(self, "mw",
                             self.config.parameters.mw)

        """List the necessary state variable objects."""
        self.flow_mol = Var(initialize=1.0,
                            domain=Reals,
                            doc='Component molar flowrate [mol/s]')
        self.mole_frac = Var(self._params.component_list,
                             domain=Reals,
                             initialize=1 / len(self._params.component_list),
                             doc='State component mole fractions [-]')
        self.pressure = Var(initialize=1.01325,
                            domain=Reals,
                            doc='State pressure [bar]')
        self.temperature = Var(initialize=298.15,
                               domain=Reals,
                               doc='State temperature [K]')

        # Create standard constraints
        # Sum mole fractions if not inlet block
        if self.config.defined_state is False:
            def sum_component_eqn(b):
                return 1e2 == 1e2 * sum(b.mole_frac[j]
                                        for j in b._params.component_list)
            self.sum_component_eqn = Constraint(rule=sum_component_eqn)

    def _mw_gas(self):
        # Molecular weight of gas mixture
        self.mw_gas = Var(domain=Reals,
                          initialize=1.0,
                          doc="Molecular weight of gas mixture [kg/mol]")

        def mw_gas_eqn(b):
            return (b.mw_gas ==
                    sum(b.mole_frac[j]*b._params.mw[j]
                        for j in b._params.component_list))
        try:
            # Try to build constraint
            self.mw_gas_eqn = Constraint(rule=mw_gas_eqn)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.mw_gas)
            self.del_component(self.mw_gas_eqn)
            raise

    def _dens_mole_vap(self):
        # Molar density
        self.dens_mole_vap = Var(domain=Reals,
                                 initialize=1.0,
                                 doc="Molar density/concentration [mol/m3]")

        def ideal_gas(b):
            return (b.dens_mole_vap*b._params.gas_const*b.temperature*1e-2 ==
                    b.pressure)
        try:
            # Try to build constraint
            self.ideal_gas = Constraint(rule=ideal_gas)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.dens_mole_vap)
            self.del_component(self.ideal_gas)
            raise

    def _dens_mole_comp_vap(self):
        # Mixture heat capacities
        self.dens_mole_comp_vap = Var(self._params.component_list,
                                      domain=Reals,
                                      initialize=1.0,
                                      doc='Component molar concentration'
                                      '[mol/m3]')

        def comp_conc_eqn(b, j):
            return b.dens_mole_comp_vap[j] == b.dens_mole_vap*b.mole_frac[j]
        try:
            # Try to build constraint
            self.comp_conc_eqn = Constraint(self._params.component_list,
                                            rule=comp_conc_eqn)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.dens_mole_comp_vap)
            self.del_component(self.comp_conc_eqn)
            raise

    def _dens_mass_vap(self):
        # Mass density
        self.dens_mass_vap = Var(domain=Reals,
                                 initialize=1.0,
                                 doc="Mass density [kg/m3]")

        def dens_mass_basis(b):
            return b.dens_mass_vap == b.mw_gas*b.dens_mole_vap
        try:
            # Try to build constraint
            self.dens_mass_basis = Constraint(rule=dens_mass_basis)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.dens_mass_vap)
            self.del_component(self.dens_mass_basis)
            raise

    def _visc_d(self):
        # Mixture dynamic viscosity
        self.visc_d = Var(domain=Reals,
                          initialize=1e-5,
                          doc="Mixture dynamic viscosity [kg/m.s]")

        def visc_d_comp(i):
            return self._params.visc_d_param[i, 1] * \
                    (self.temperature**self._params.visc_d_param[i, 2]) \
                    / ((1 + (self._params.visc_d_param[i, 3]/self.temperature))
                        + (self._params.visc_d_param[i, 4] /
                           (self.temperature**2)))

        def visc_d_constraint(b):
            return 1e6*b.visc_d == 1e6*sum(b.mole_frac[i]*visc_d_comp(i)
                                           / (sum(b.mole_frac[j]
                                                  * (b._params.mw[j] /
                                                     b._params.mw[i])**0.5
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
                                  '[cm2/s]')

        def D_bin(i, j):
            # 1e3 used to multiply MW to convert from kg/mol to kg/kmol
            return ((1.43e-3*(self.temperature**1.75) *
                     ((1e3 * self._params.mw[i] + 1e3 * self._params.mw[j])
                     / (2 * (1e3 * self._params.mw[i]) *
                        (1e3*self._params.mw[j])))**0.5)
                    / ((self.pressure)
                        * ((self._params.diff_vol_param[i]**(1/3))
                            + (self._params.diff_vol_param[j]**(1/3)))**2))

        def diffusion_comp_constraint(b, i):
            return (b.diffusion_comp[i]
                    * sum(b.mole_frac[j]/D_bin(i, j)
                    for j in b._params.component_list if i != j)
                    == (1-b.mole_frac[i]))
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
                              doc="Thermal conductivity of gas [kJ/m.K.s]")

        def therm_cond_comp(i):
            return self._params.therm_cond_param[i, 1] \
                    * (self.temperature**self._params.therm_cond_param[i, 2]) \
                    / ((1 + (self._params.therm_cond_param[i, 3] /
                             self.temperature))
                        + (self._params.therm_cond_param[i, 4] /
                       (self.temperature**2)))

        def A_bin(i, j):
            return (1 + ((therm_cond_comp(j)/therm_cond_comp(i))**0.5)
                    * ((self._params.mw[j]/self._params.mw[i])**0.25))**2 \
                    / (8*(1+(self._params.mw[j]/self._params.mw[i])))**0.5

        def therm_cond_constraint(b):
            # The 1e-3 term is used as a conversion factor to a kJ basis
            return 1e6*b.therm_cond == 1e6*(1e-3) * \
                                       sum(b.mole_frac[i]
                                           * therm_cond_comp(i)
                                           / (sum(b.mole_frac[j] *
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
                               "[kJ/mol.K]")

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
                          doc="Mixture heat capacity [kJ/mol.K]")

        def cp_mol(b):
            return b.cp_mol == sum(b.cp_mol_comp[j]*b.mole_frac[j]
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
                           doc="Mixture heat capacity, mass-basis [kJ/kg.K]")

        def cp_mass(b):
            return b.cp_mass*b.mw_gas == b.cp_mol
        try:
            # Try to build constraint
            self.cp_mass_basis = Constraint(rule=cp_mass)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mass)
            self.del_component(self.cp_mass_basis)
            raise

    def _enth_mol_vap_comp(self):
        # Pure component vapour enthalpies
        self.enth_mol_vap_comp = Var(
                self._params.component_list,
                domain=Reals,
                initialize=1.0,
                doc="Pure component enthalpies [kJ/mol]")

        def pure_comp_enthalpy(b, j):
            return b.enth_mol_vap_comp[j] == (
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
                            doc='Mixture specific enthalpy [kJ/mol]')
        try:
            # Try to build constraint
            self.mixture_enthalpy_eqn = Constraint(expr=(
                        self.enth_mol == sum(self.mole_frac[j] *
                                             self.enth_mol_vap_comp[j]
                                             for j in
                                             self._params.component_list)))
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.enth_mol)
            self.del_component(self.mixture_enthalpy_eqn)
            raise

    def get_material_flow_terms(b, p, j):
        return b.flow_mol*b.mole_frac[j]

    def get_enthalpy_flow_terms(b, p):
        return b.flow_mol*b.enth_mol

    def get_material_density_terms(b, p, j):
        return b.dens_mole_comp_vap[j]

    def get_energy_density_terms(b, p):
        return b.dens_mole_vap*b.enth_mol

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
