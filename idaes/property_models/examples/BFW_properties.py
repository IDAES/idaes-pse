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
"""Basic property package for water liquid.

Please verify the valid ranges for the temperature
before using this property package.

All unit is SI and mole basis
"""

# Chages the divide behavior to not do integer division
from __future__ import division

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import (Constraint, log, Param, value,
                           PositiveReals, RangeSet, Reals, Set, Var)
from pyomo.opt import SolverFactory, TerminationCondition

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock)
from idaes.core.util.misc import add_object_reference

# Some more inforation about this module
__author__ = "Jaffer Ghouse"


# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("BFWParameterBlock")
class PhysicalParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class.

    Contains parameters and indexing sets associated with properties for
    superheated steam.

    """
    def build(self):
        """Callable method for Block construction."""
        super(PhysicalParameterData, self).build()
        self.state_block_class = BFWStateBlock
        self._make_params()

    def _make_params(self):

        # List of valid phases in property package
        self.phase_list = Set(initialize=['Liq'])

        # Component list - a Set of component identifiers
        self.component_list = Set(initialize=['H2O'])

        # List of components in each phase (optional)
        self.phase_component_list = {"Liq": self.component_list}

        # Thermodynamic reference state
        self.temperature_ref = Param(within=PositiveReals,
                                     mutable=True,
                                     default=298.15,
                                     doc='Reference temperature [K]')

        # Specific heat capacity coefficients Temp range - 273.16 to 533.15 K
        # Cp in J/kmol/K converted to J/mol/K
        # Specific heat capacity coefficients
        cp_param = {1: 2.7637E5,
                    2: -2.0901E3,
                    3: 8.125,
                    4: -1.4116E-2,
                    5: 9.3701E-6}
        self.cp_param = Param(RangeSet(5), initialize=cp_param,
                              doc="specific heat parameters")

        # Viscosity coefficients Temp range - 273.16 K to 646.15 K
        # Viscosity is in centipoise converted to Pa.s
        # Source: DIPPR (From Aspen Properties)
        visc_d_phase_param = {1: -45.9352,
                              2: 3703.6,
                              3: 5.866,
                              4: -5.879E-29,
                              5: 10}
        self.visc_d_phase_param = Param(RangeSet(5),
                                        initialize=visc_d_phase_param,
                                        doc="dynamic viscosity parameters")

        # Thermal conductivoty coefficients Temp range - 273.16 K to 633.15 K
        # Thermal conductivity is in W/m-K
        # Source: DIPPR (From Aspen Properties)
        therm_cond_phase_param = {1: -0.432,
                                  2: 5.7255E-3,
                                  3: -8.078E-6,
                                  4: 1.861E-9}
        self.therm_cond_phase_param = Param(RangeSet(4),
                                            initialize=therm_cond_phase_param,
                                            doc="thermal conductivity W/m/K")

        self.mw = Param(initialize=18.015E-3, doc="molecular weight Kg/mol")

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({'mw': {'method': None, 'units': 'kg/mol'},
                            'flow_mol': {'method': None, 'units': 'mol/s'},
                            'pressure': {'method': None, 'units': 'Pa'},
                            'temperature': {'method': None, 'units': 'K'},
                            'vapor_frac': {'method': None, 'units': None},
                            'dens_mol_phase': {'method': None,
                                               'units': 'mol/m^3'},
                            'enth_mol_phase': {'method': None,
                                               'units': 'J/mol'},
                            'cp_mol_phase': {'method': None,
                                             'units': 'J/mol.K'},
                            'visc_d_phase': {'method': None, 'units': 'Pa.s'},
                            'therm_cond_phase': {'method': None,
                                                 'units': 'W/m.K'}})
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

    def initialize(blk, flow_mol=None, temperature=None, pressure=None,
                   vapor_frac=None, outlvl=0, hold_state=False,
                   state_vars_fixed=False, solver='ipopt',
                   optarg={'tol': 1e-8}):
        """
        Declare initialisation routine.

        Keyword Arguments:
            state_args = to be used if state block initialized independent of
                         control volume initialize
            outlvl : sets output level of initialisation routine

                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = include solver output infomation (tee=True)

            optarg : solver options dictionary object (default=None)
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')
            state_vars_fixed: Flag to denote if state vars have already been
                              fixed for the state block.
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
        # Fix state variables if not already fixed by the control volume block
        if state_vars_fixed is False:
            Fflag = {}
            Pflag = {}
            Tflag = {}
            vfflag = {}

            for k in blk.keys():
                if blk[k].flow_mol.fixed is True:
                    Fflag[k] = True
                else:
                    Fflag[k] = False
                    if flow_mol is None:
                        blk[k].flow_mol.fix(1.0)
                    else:
                        blk[k].flow_mol.fix(flow_mol)
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
                        blk[k].temperature.fix(300.0)
                    else:
                        blk[k].temperature.fix(temperature)

                if blk[k].vapor_frac.fixed is True:
                    vfflag[k] = True
                else:
                    vfflag[k] = False
                    if vapor_frac is None:
                        blk[k].vapor_frac.fix(300.0)
                    else:
                        blk[k].vapor_frac.fix(temperature)

            flags = {"Fflag": Fflag, "Pflag": Pflag,
                     "Tflag": Tflag, "vfflag": vfflag}

        # Set solver options
        if outlvl > 1:
            stee = True
        else:
            stee = False

        opt = SolverFactory(solver)
        opt.options = optarg

        # ---------------------------------------------------------------------
        # Solve property correlation
        for k in blk.keys():
            results = opt.solve(blk[k], tee=stee)

        if outlvl > 0:
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                _log.info('{} Initialisation Step 1 Complete.'
                          .format(blk.name))
            else:
                _log.warning('{} Initialisation Step 1 Failed.'
                             .format(blk.name))

        # ---------------------------------------------------------------------
        if state_vars_fixed is False:
            # release state vars fixed during initialization if control
            # volume didn't fix the state vars
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)

        if outlvl > 0:
            if outlvl > 0:
                _log.info('{} Initialisation Complete.'.format(blk.name))

    def release_state(blk, flags, outlvl=0):
        # Method to release states only if explicitly called

        if flags is None:
            return

        # Unfix state variables
        for k in blk.keys():
            if flags['Fflag'][k] is False:
                blk[k].flow_mol.unfix()
            if flags['Pflag'][k] is False:
                blk[k].pressure.unfix()
            if flags['Tflag'][k] is False:
                blk[k].temperature.unfix()
            if flags['vfflag'][k] is False:
                blk[k].vapor_frac.unfix()

        if outlvl > 0:
            if outlvl > 0:
                _log.info('{} State Released.'.format(blk.name))

@declare_process_block_class("BFWStateBlock",
                             block_class=_StateBlock)
class StateTestBlockData(StateBlockData):
    """An example property package for boiler feed water properties."""

    def build(self):
        """Callable method for Block construction."""
        super(StateTestBlockData, self).build()
        self._make_params()
        self._make_state_vars()
        self._make_prop_vars()
        self._make_constraints()

    def _make_params(self):
        """Make references to the necessary parameters contained."""
        # List of valid phases in property package
        add_object_reference(self, "phase_list",
                             self.config.parameters.phase_list)

        # Component list - a list of component identifiers
        add_object_reference(self, "component_list",
                             self.config.parameters.component_list)

        # Thermodynamic reference state
        add_object_reference(self, "temperature_ref",
                             self.config.parameters.temperature_ref)

        # Specific Enthalpy Coefficients
        add_object_reference(self, "cp_param",
                             self.config.parameters.cp_param)

        add_object_reference(self, "visc_d_phase_param",
                             self.config.parameters.visc_d_phase_param)

        add_object_reference(self, "therm_cond_phase_param",
                             self.config.parameters.therm_cond_phase_param)

        add_object_reference(self, "mw", self.config.parameters.mw)

    def _make_state_vars(self):
        """Declare the necessary state variable objects."""
        self.flow_mol = Var(domain=Reals,
                            initialize=0.5,
                            doc='Component mole flowrate [mol/s]')
        self.pressure = Var(domain=Reals,
                            initialize=1.01325E5,
                            doc='State pressure [Pa]')
        self.temperature = Var(domain=Reals,
                               initialize=298.15,
                               doc='State temperature [K]')

        self.vapor_frac = Var(initialize=0,
                              doc="Fraction of water in the vapor phase \
                                    [unitless]", bounds=(0, 1))

    def _make_prop_vars(self):
            """Make additional variables for calcuations."""
            self.dens_mol_phase = Param(self.phase_list,
                                        initialize=43E3,
                                        doc="molar density mol/m3")

            self.enth_mol_phase = Var(self.phase_list,
                                      doc='Specific Enthalpy [J/mol]')

            self.cp_mol_phase = Var(self.phase_list,
                                    doc="Specific heat capacity [J/mol/K]")

            self.visc_d_phase = Var(self.phase_list,
                                    initialize=0.001,
                                    doc="Dynamic viscosity [Pa.s]")

            self.therm_cond_phase = Var(self.phase_list,
                                        initialize=1,
                                        doc="thermal conductivity [W/m/K]")

    def _make_constraints(self):
            """Create property constraints."""
            # Specific heat capacity
            def cp_mol_phase_liq_correlation(self, p):
                return self.cp_mol_phase[p] == \
                    1E-3 * (self.cp_param[1] +
                            (self.cp_param[2] * self.temperature) +
                            (self.cp_param[3] * self.temperature**2) +
                            (self.cp_param[4] * self.temperature**3) +
                            (self.cp_param[5] * self.temperature**4))
            self.cp_mol_phase_liq_correlation = Constraint(
                self.phase_list, rule=cp_mol_phase_liq_correlation)

            # Specific Enthalpy
            def enthalpy_correlation(self, p):
                return self.enth_mol_phase[p] * 1E3 == \
                    ((self.cp_param[1] * self.temperature) +
                     (self.cp_param[2] * 0.5 * self.temperature**2) +
                     (self.cp_param[3] * 0.33 * self.temperature**3) +
                     (self.cp_param[4] * 0.25 * self.temperature**4) +
                     (self.cp_param[5] * 0.2 * self.temperature**5)) - \
                    ((self.cp_param[1] * self.temperature_ref) +
                     (self.cp_param[2] * 0.5 * self.temperature_ref**2) +
                     (self.cp_param[3] * 0.33 * self.temperature_ref**3) +
                     (self.cp_param[4] * 0.25 * self.temperature_ref**4) +
                     (self.cp_param[5] * 0.2 * self.temperature_ref**5))
            self.enthalpy_correlation = Constraint(self.phase_list,
                                                   rule=enthalpy_correlation)

            # Dynamic viscosity (1E3 factor to convert from CP to Pa.s)
            def visc_correlation(self, p):
                return self.temperature * log(self.visc_d_phase[p] * 1E3) == \
                    self.visc_d_phase_param[1] * self.temperature + \
                    self.visc_d_phase_param[2] + self.visc_d_phase_param[3] * \
                    self.temperature * log(self.temperature) + \
                    self.visc_d_phase_param[4] *\
                    self.temperature**(self.visc_d_phase_param[5] + 1)
            self.visc_correlation = Constraint(self.phase_list,
                                               rule=visc_correlation)

            # Thermal conductivity (1E3 factor is to convert from W/m-K)
            def therm_cond_phase_correlation(self, p):
                return self.therm_cond_phase[p] ==\
                    self.therm_cond_phase_param[1] + \
                    self.therm_cond_phase_param[2] * self.temperature + \
                    self.therm_cond_phase_param[3] * self.temperature**2 +\
                    self.therm_cond_phase_param[4] * self.temperature**3
            self.therm_cond_phase_correlation = Constraint(
                self.phase_list, rule=therm_cond_phase_correlation)

    def get_material_flow_terms(b, p, j):
        """Define material flow terms for control volume."""
        return b.flow_mol

    def get_enthalpy_flow_terms(b, p):
        """Define enthalpy flow terms for control volume."""
        return b.flow_mol * b.enth_mol_phase[p]

    def get_material_density_terms(b, p, j):
        """Define material density terms for control volume."""
        return b.dens_mol_phase[p]

    def get_enthalpy_density_terms(b, p):
        """Define enthalpy density terms for control volume."""
        return b.enth_mol_phase[p] * b.dens_mol_phase[p]

    def define_state_vars(b):
        """Define state variables for ports."""
        return {"flow_mol": b.flow_mol,
                "temperature": b.temperature,
                "pressure": b.pressure,
                "vapor_frac": b.vapor_frac}

    def model_check(blk):
        """Model checks for property block."""
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
