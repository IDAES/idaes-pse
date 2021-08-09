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
Example property package for testing TGA model solid phase.

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


@declare_process_block_class("ZeoliteParameterBlock")
class ZeoliteParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class
    """
    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(ZeoliteParameterData, self).build()

        self.state_block_class = ZeoliteStateBlock

        # List of valid phases in property package
        self.phase_list = Set(initialize=['Sol'])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['Zeo', 'N2'])

        # Heat capacity
        self.cp_mass_sol = Param(initialize=840,
                                 doc="Specific heat capacity [J/kg.K]")

        self.dens_mass_sol = Param(initialize=2100,
                                   doc="Particle density [kg/m^3]")

        self.diameter_particle = Param(initialize=9e-6,
                                       doc="Particle diameter [m]")

        self.area_mass = Param(initialize=8e5,
                               doc="Specific surface area [m^2/kg]")

        self.mw_comp = Param(self.component_list,
                             initialize={"Zeo": 0.162143, "N2": 2.8e-2},
                             doc="Component molecular weights [kg/mol]")

        # 720 kg/m^3 bulk density
        # pore volume 0.3 cc/g
        # dHads H2O = 4186800 J/kg

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
                'flow_mass': {'method': None, 'units': 'kg/s'},
                'pressure': {'method': None, 'units': 'Pa'},
                'temperature': {'method': None, 'units': 'K'},
                'loading': {'method': None, 'units': 'mol/kg'},
                'enth_mass': {'method': None, 'units': 'J/kg'},
                'dens_mass': {'method': None, 'units': 'kg/m^3'}})
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
                blk[k].zeolite_loading_eqn.deactivate()

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
                blk[k].zeolite_loading_eqn.activate()

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


@declare_process_block_class("ZeoliteStateBlock",
                             block_class=_StateBlock)
class ZeoliteStateBlockData(StateBlockData):
    """
    An example property package for a simple zeolite sorbent
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(ZeoliteStateBlockData, self).build()

        # Create state variables
        self.flow_mass = Var(initialize=1.0,
                             domain=NonNegativeReals,
                             doc='Solids mass flow rate [kg/s]')
        self.pressure = Var(domain=Reals,
                            initialize=101325.0,
                            bounds=(1e3, 1e6),
                            doc='State pressure [Pa]')
        self.temperature = Var(domain=Reals,
                               initialize=298.15,
                               bounds=(298.15, 323.15),
                               doc='State temperature [K]')
        self.loading = Var(self._params.component_list,
                           domain=NonNegativeReals,
                           initialize=0,
                           doc='Loading of adsorbed species [mol/kg]')

        if self.config.defined_state is False:
            self.zeolite_loading_eqn = Constraint(
                    expr=1 == self.loading["Zeo"])

        # TODO: Should account for enthalpy of adsorbed speices, but leave this until later
        self.enth_mass = Var(domain=Reals,
                             initialize=1000,
                             doc='Specific enthalpy of solid [J/kg]')

        self.enth_mass_eqn = Constraint(
                expr=self.enth_mass == self._params.cp_mass_sol*(
                        self.temperature-self._params.temperature_ref))

        self.dens_mass = Var(domain=Reals,
                             initialize=1,
                             doc='Molar density of mixture [kg/m^3]')

        self.dens_mass_eqn = Constraint(
                expr=self.dens_mass == self._params.dens_mass_sol +
                sum(self.loading[j]*self._params.mw_comp[j]
                    for j in self._params.component_list))

    def get_material_flow_terms(b, p, j):
        if j == "Zeo":
            return b.flow_mass
        else:
            return b.flow_mass*b.loading[j]*b._params.mw_comp[j]

    def get_enthalpy_flow_terms(b, p):
        return b.flow_mass*b.enth_mass

    def get_material_density_terms(b, p, j):
        if j == "Zeo":
            return b._params.dens_mass_sol
        else:
            return b._params.dens_mass_sol*b.loading[j]*b._params.mw_comp[j]

    # TODO: As enth_mass does not include adsorbate, only consider base solid here
    def get_energy_density_terms(b, p):
        return b._params.dens_mass_sol*b.enth_mass

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def define_state_vars(b):
        return {"flow_mass": b.flow_mass,
                "loading": b.loading,
                "temperature": b.temperature,
                "pressure": b.pressure}

    def get_material_flow_basis(b):
        return MaterialFlowBasis.mass
