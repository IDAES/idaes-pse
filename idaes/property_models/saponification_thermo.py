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
Example property package for the saponification of Ethyl Acetate with NaOH
Assumes dilute solutions with properties of H2O.
"""

# Chages the divide behavior to not do integer division
from __future__ import division

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           NonNegativeReals,
                           Param,
                           PositiveReals,
                           Reals,
                           Set,
                           value,
                           Var)
from pyomo.opt import SolverFactory

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        PhysicalParameterBase,
                        StateBlockDataBase,
                        StateBlockBase)
from idaes.core.util.misc import add_object_reference

# Some more inforation about this module
__author__ = "Andrew Lee"


# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("PhysicalParameterBlock")
class PhysicalParameterData(PhysicalParameterBase):
    """
    Property Parameter Block Class

    Contains parameters and indexing sets associated with properties for
    superheated steam.

    """
    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(PhysicalParameterData, self).build()

        self.state_block_class = StateBlock

        # List of valid phases in property package
        self.phase_list = Set(initialize=['Liq'])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['H2O', 'NaOH',
                                              'EthylAcetate',
                                              'SodiumAcetate',
                                              'Ethanol'])

        # Heat capacity of water
        self.cp_mol = Param(mutable=False,
                            initialize=75.327,
                            doc="Molar heat capacity of water [J/mol.K]")

        # Density of water
        self.dens_mol = Param(mutable=False,
                              initialize=55388.0,
                              doc="Molar density of water [mol/m^3]")

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
                'flow_vol': {'method': None, 'units': 'm^3/s'},
                'pressure': {'method': None, 'units': 'Pa'},
                'temperature': {'method': None, 'units': 'K'},
                'conc_mol_comp': {'method': None, 'units': 'mol/m^3'},
                'cp_mol': {'method': None, 'units': 'J/mol.K'},
                'dens_mol': {'method': None, 'units': 'mol/m^3'}})
        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'g',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'mol'})


class _StateBlock(StateBlockBase):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """
    def initialize(blk, flow_vol=None, temperature=None, pressure=None,
                   conc_mol_comp=None,
                   hold_state=False, outlvl=0,
                   solver='ipopt', optarg={'tol': 1e-8}):
        '''
        Initialisation routine for property package.

        Keyword Arguments:
            flow_mol_comp : value at which to initialize component flows
                             (default=None)
            pressure : value at which to initialize pressure (default=None)
            temperature : value at which to initialize temperature
                          (default=None)
            outlvl : sets output level of initialisation routine

                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = include solver output infomation (tee=True)

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
        # Fix state variables if not already fixed
        Fflag = {}
        Pflag = {}
        Tflag = {}
        Cflag = {}

        for k in blk.keys():
            if blk[k].flow_vol.fixed is True:
                Fflag[k] = True
            else:
                Fflag[k] = False
                if flow_vol is None:
                    blk[k].flow_vol.fix(1.0)
                else:
                    blk[k].flow_vol.fix(flow_vol)

            for j in blk[k].component_list_ref:
                if blk[k].conc_mol_comp[j].fixed is True:
                    Cflag[k, j] = True
                else:
                    Cflag[k, j] = False
                    if conc_mol_comp is None:
                        if j == "H2O":
                            blk[k].conc_mol_comp[j].fix(55388.0)
                        else:
                            blk[k].conc_mol_comp[j].fix(100.0)
                    else:
                        blk[k].conc_mol_comp[j].fix(conc_mol_comp[j])

            if blk[k].pressure.fixed is True:
                Pflag[k] = True
            else:
                Pflag[k] = False
                if pressure is None:
                    blk[k].pressure.fix(101325.0)
                else:
                    blk[k].pressure.fix(pressure)

            if blk[k].temperature.fixed is True:
                Tflag[k] = True
            else:
                Tflag[k] = False
                if temperature is None:
                    blk[k].temperature.fix(298.15)
                else:
                    blk[k].temperature.fix(temperature)

        opt = SolverFactory(solver)
        opt.options = optarg

        # ---------------------------------------------------------------------
        # If input block, return flags, else release state
        flags = {"Fflag": Fflag, "Pflag": Pflag,
                 "Tflag": Tflag, "Cflag": Cflag}

        if outlvl > 0:
            if outlvl > 0:
                _log.info('{} Initialisation Complete.'.format(blk.name))

        if hold_state is True:
            return flags
        else:
            blk.release_state(flags)

    def release_state(blk, flags, outlvl=0):
        '''
        Method to relase state variables fixed during initialisation.

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
        for k in blk.keys():
            if flags['Fflag'][k] is False:
                blk[k].flow_vol.unfix()
            for j in blk[k].component_list_ref:
                if flags['Cflag'][k, j] is False:
                    blk[k].conc_mol_comp[j].unfix()
            if flags['Pflag'][k] is False:
                blk[k].pressure.unfix()
            if flags['Tflag'][k] is False:
                blk[k].temperature.unfix()

        if outlvl > 0:
            if outlvl > 0:
                _log.info('{} State Released.'.format(blk.name))


@declare_process_block_class("StateBlock",
                             block_class=_StateBlock)
class StateBlockData(StateBlockDataBase):
    """
    An example property package for properties for saponification of ethyl
    acetate
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(StateBlockData, self).build()

        # Create references to package parameters
        # List of valid phases in property package
        add_object_reference(self,
                             "phase_list_ref",
                             self.config.parameters.phase_list)

        # Component list - a list of component identifiers
        add_object_reference(self,
                             "component_list_ref",
                             self.config.parameters.component_list)

        # Heat capacity - no _ref ending as this is the actual property
        add_object_reference(self,
                             "cp_mol",
                             self.config.parameters.cp_mol)

        # Density - no _ref ending as this is the actual property
        add_object_reference(self,
                             "dens_mol",
                             self.config.parameters.dens_mol)

        # Thermodynamic reference state
        add_object_reference(self, "pressure_ref_ref",
                             self.config.parameters.pressure_ref)
        add_object_reference(self, "temperature_ref_ref",
                             self.config.parameters.temperature_ref)

        # Create state variables
        self.flow_vol = Var(initialize=1.0,
                            domain=NonNegativeReals,
                            doc='Total volumentric flowrate [m^3/s]')
        self.pressure = Var(domain=Reals,
                            initialize=101325.0,
                            bounds=(1e3, 1e6),
                            doc='State pressure [Pa]')
        self.temperature = Var(domain=Reals,
                               initialize=298.15,
                               bounds=(298.15, 323.15),
                               doc='State temperature [K]')
        self.conc_mol_comp = Var(self.component_list_ref,
                                 domain=NonNegativeReals,
                                 initialize=100.0,
                                 doc='Component molar concentrations '
                                     '[mol/m^3]')

        if self.config.defined_state is False:
            self.conc_water_eqn = Constraint(expr=self.conc_mol_comp["H2O"] ==
                                             self.dens_mol)

    def get_material_flow_terms(b, p, j):
        return b.flow_vol*b.conc_mol_comp[j]

    def get_enthalpy_flow_terms(b, p):
        return (b.flow_vol*b.dens_mol*b.cp_mol *
                (b.temperature - b.temperature_ref_ref))

    def define_state_vars(b):
        return {"flow_vol": b.flow_vol,
                "conc_mol_comp": b.conc_mol_comp,
                "temperature": b.temperature,
                "pressure": b.pressure}

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
