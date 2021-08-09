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
Example property package for testing TGA model gas phase.

Assumes ideal gas behavior.
"""

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           NonNegativeReals,
                           Param,
                           PositiveReals,
                           Reals,
                           Set,
                           Var)

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import fix_state_vars, revert_state_vars

# Some more inforation about this module
__author__ = "Andrew Lee"


# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("AirParameterBlock")
class AirParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class
    """
    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(AirParameterData, self).build()

        self.state_block_class = AirStateBlock

        # List of valid phases in property package
        self.phase_list = Set(initialize=['Vap'])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['N2', 'O2'])

        # Heat capacity coefficients
        # "The Properties of Gases and Liquids, 4th Edition",
        # Reid, Prausnitz and Poling, McGraw-Hill, 1987
        self.cp_mol_ig_comp_coeff = Param(
                self.component_list,
                ["A", "B", "C", "D"],
                mutable=False,
                initialize={("N2", "A"): 3.115e1,
                            ("N2", "B"): -1.357e-2,
                            ("N2", "C"): 2.680e-5,
                            ("N2", "D"): -1.168e-8,
                            ("O2", "A"): 2.811e1,
                            ("O2", "B"): -3.680e-6,
                            ("O2", "C"): 1.746e-5,
                            ("O2", "D"): -1.065e-8},
                doc="Molar heat capacity coefficients [J/mol.K]")

        # Thermodynamic reference state
        self.pressure_ref = Param(within=PositiveReals,
                                  mutable=True,
                                  default=101325.0,
                                  doc='Reference pressure [Pa]')
        self.temperature_ref = Param(within=PositiveReals,
                                     mutable=True,
                                     default=298.15,
                                     doc='Reference temperature [K]')

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
                'flow_mol': {'method': None, 'units': 'mol/s'},
                'pressure': {'method': None, 'units': 'Pa'},
                'temperature': {'method': None, 'units': 'K'},
                'mole_frac_comp': {'method': None, 'units': '-'},
                'enth_mol': {'method': None, 'units': 'J/mol'},
                'dens_mol': {'method': None, 'units': 'mol/m^3'}})
        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'kg',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'mol'})


class _StateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """
    def initialize(blk, state_args={}, state_vars_fixed=False,
                   hold_state=False, outlvl=0,
                   solver='ipopt', optarg={'tol': 1e-8}):
        '''
        Initialization routine for property package.

        Keyword Arguments:
        state_args : Dictionary with initial guesses for the state vars
                     chosen. Note that if this method is triggered
                     through the control volume, and if initial guesses
                     were not provied at the unit model level, the
                     control volume passes the inlet values as initial
                     guess.The keys for the state_args dictionary are:

                     flow_mol_comp : value at which to initialize component
                                     flows (default=None)
                     pressure : value at which to initialize pressure
                                (default=None)
                     temperature : value at which to initialize temperature
                                  (default=None)
            outlvl : sets output level of initialization routine

                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = include solver output infomation (tee=True)
            state_vars_fixed: Flag to denote if state vars have already been
                              fixed.
                              - True - states have already been fixed by the
                                       control volume 1D. Control volume 0D
                                       does not fix the state vars, so will
                                       be False if this state block is used
                                       with 0D blocks.
                             - False - states have not been fixed. The state
                                       block will deal with fixing/unfixing.
            optarg : solver options dictionary object (default=None)
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')
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
        # Deactivate the constraints specific for outlet block
        for k in blk.keys():
            if blk[k].config.defined_state is False:
                blk[k].sum_mole_frac_comp_eqn.deactivate()

        if state_vars_fixed is False:
            # Fix state variables if not already fixed
            flags = fix_state_vars(blk, state_args)

        else:
            # Check when the state vars are fixed already result in dof 0
            for k in blk.keys():
                if degrees_of_freedom(blk[k]) != 0:
                    raise Exception("State vars fixed but degrees of freedom "
                                    "for state block is not zero during "
                                    "initialization.")

        # Reactivate constraints
        for k in blk.keys():
            if not blk[k].config.defined_state:
                blk[k].sum_mole_frac_comp_eqn.activate()

        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)

        if outlvl > 0:
            _log.info('{} Initialization Complete.'.format(blk.name))

    def release_state(blk, flags, outlvl=0):
        '''
        Method to relase state variables fixed during initialization.

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

        if outlvl > 0:
            if outlvl > 0:
                _log.info('{} State Released.'.format(blk.name))


@declare_process_block_class("AirStateBlock",
                             block_class=_StateBlock)
class AirStateBlockData(StateBlockData):
    """
    An example property package for ideal gas properties of air
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(AirStateBlockData, self).build()

        # Create state variables
        self.flow_mol = Var(initialize=1.0,
                            domain=NonNegativeReals,
                            doc='Total molar flowrate [mol/s]')
        self.pressure = Var(domain=Reals,
                            initialize=101325.0,
                            bounds=(1e3, 1e6),
                            doc='State pressure [Pa]')
        self.temperature = Var(domain=Reals,
                               initialize=298.15,
                               bounds=(298.15, 323.15),
                               doc='State temperature [K]')
        self.mole_frac_comp = Var(self._params.component_list,
                                  domain=NonNegativeReals,
                                  initialize=0.5,
                                  doc='Component mole fractions [-]')

        if self.config.defined_state is False:
            self.sum_mole_frac_comp_eqn = Constraint(
                    expr=1 == sum(self.mole_frac_comp[j]
                                  for j in self._params.component_list))

        self.enth_mol = Var(domain=Reals,
                            initialize=1000,
                            doc='Molar enthalpy of mixture [J/mol]')

        self.enth_mol_eqn = Constraint(
                expr=self.enth_mol == sum(
                        (self._params.cp_mol_ig_comp_coeff[j, "D"]/4) *
                        (self.temperature**4-self._params.temperature_ref**4) +
                        (self._params.cp_mol_ig_comp_coeff[j, "C"]/3) *
                        (self.temperature**3-self._params.temperature_ref**3) +
                        (self._params.cp_mol_ig_comp_coeff[j, "B"]/2) *
                        (self.temperature**2-self._params.temperature_ref**2) +
                        self._params.cp_mol_ig_comp_coeff[j, "A"] *
                        (self.temperature-self._params.temperature_ref)
                        for j in self._params.component_list))

        self.dens_mol = Var(domain=Reals,
                            initialize=1,
                            doc='Molar density of mixture [mol/m^3]')

        self.dens_mol_eqn = Constraint(
                expr=self.dens_mol*self.pressure == 8.314*self.temperature)

    def get_material_flow_terms(b, p, j):
        return b.flow_mol*b.mole_frac_comp[j]

    def get_enthalpy_flow_terms(b, p):
        return (b.flow_mol*b.enth_mol)

    def get_material_density_terms(b, p, j):
        return b.dens_mol*b.mole_frac_comp[j]

    def get_energy_density_terms(b, p):
        # Subtracting R*(T-Tref) to convert from enthalpy to internal energy
        # Cv = Cp - R
        return b.dens_mol*(b.enth_mol - 8.314*(b.temperature -
                                               b._params.temperature_ref))

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def define_state_vars(b):
        return {"flow_mol": b.flow_mol,
                "mole_frac_comp": b.mole_frac_comp,
                "temperature": b.temperature,
                "pressure": b.pressure}

    def get_material_flow_basis(b):
        return MaterialFlowBasis.molar
