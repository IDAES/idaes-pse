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
Example property package for the saponification of Ethyl Acetate with NaOH
Assumes dilute solutions with properties of H2O.
"""

# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    NonNegativeReals,
    Param,
    PositiveReals,
    Reals,
    units,
    value,
    Var,
)

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    LiquidPhase,
    Component,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
import idaes.logger as idaeslog

# Some more inforation about this module
__author__ = "Andrew Lee"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("SaponificationParameterBlock")
class PhysicalParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class

    Contains parameters and indexing sets associated with properties for
    superheated steam.

    """

    def build(self):
        """
        Callable method for Block construction.
        """
        super(PhysicalParameterData, self).build()

        self._state_block_class = SaponificationStateBlock

        # Add Phase objects
        self.Liq = LiquidPhase()

        # Add Component objects
        self.H2O = Component()
        self.NaOH = Component()
        self.EthylAcetate = Component()
        self.SodiumAcetate = Component()
        self.Ethanol = Component()

        # Heat capacity of water
        self.cp_mol = Param(
            mutable=False,
            initialize=75.327,
            doc="Molar heat capacity of water [J/mol.K]",
            units=units.J / units.mol / units.K,
        )

        # Density of water
        self.dens_mol = Param(
            mutable=False,
            initialize=55388.0,
            doc="Molar density of water [mol/m^3]",
            units=units.mol / units.m**3,
        )

        # Thermodynamic reference state
        self.pressure_ref = Param(
            within=PositiveReals,
            mutable=True,
            default=101325.0,
            doc="Reference pressure [Pa]",
            units=units.Pa,
        )
        self.temperature_ref = Param(
            within=PositiveReals,
            mutable=True,
            default=298.15,
            doc="Reference temperature [K]",
            units=units.K,
        )

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_vol": {"method": None, "units": "m^3/s"},
                "pressure": {"method": None, "units": "Pa"},
                "temperature": {"method": None, "units": "K"},
                "conc_mol_comp": {"method": None, "units": "mol/m^3"},
                "dens_mol": {"method": None, "units": "mol/m^3"},
            }
        )
        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class _StateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(
        blk,
        state_args=None,
        state_vars_fixed=False,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
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
            state_vars_fixed: Flag to denote if state vars have already been
                              fixed.
                              - True - states have already been fixed and
                                       initialization does not need to worry
                                       about fixing and unfixing variables.
                             - False - states have not been fixed. The state
                                       block will deal with fixing/unfixing.
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)
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
        # Deactivate the constraints specific for outlet block i.e.
        # when defined state is False
        # This is needed as fixing state vars fixes conc_mol_comp["H2O"],
        # which is also specified by the conc_water_eqn constraint
        for k in blk.keys():
            if blk[k].config.defined_state is False:
                blk[k].conc_water_eqn.deactivate()

        if state_vars_fixed is False:
            # Fix state variables if not already fixed
            flags = fix_state_vars(blk, state_args)

        else:
            # Check when the state vars are fixed already result in dof 0
            for k in blk.keys():
                if degrees_of_freedom(blk[k]) != 0:
                    raise Exception(
                        "State vars fixed but degrees of freedom "
                        "for state block is not zero during "
                        "initialization."
                    )

        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)

        init_log.info("Initialization Complete.")

    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        """
        Method to relase state variables fixed during initialization.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")

        # Reactivate conc_water_eqn
        for k in blk.keys():
            if not blk[k].config.defined_state:
                blk[k].conc_water_eqn.activate()

        if flags is None:
            return
        # Unfix state variables
        revert_state_vars(blk, flags)
        init_log.info("State Released.")


@declare_process_block_class("SaponificationStateBlock", block_class=_StateBlock)
class SaponificationStateBlockData(StateBlockData):
    """
    An example property package for properties for saponification of ethyl
    acetate
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(SaponificationStateBlockData, self).build()

        # Create state variables
        self.flow_vol = Var(
            initialize=1.0,
            domain=NonNegativeReals,
            doc="Total volumentric flowrate [m^3/s]",
            units=units.m**3 / units.s,
        )
        self.pressure = Var(
            domain=Reals,
            initialize=101325.0,
            bounds=(1e3, 1e6),
            doc="State pressure [Pa]",
            units=units.Pa,
        )
        self.temperature = Var(
            domain=Reals,
            initialize=298.15,
            bounds=(298.15, 323.15),
            doc="State temperature [K]",
            units=units.K,
        )
        self.conc_mol_comp = Var(
            self.params.component_list,
            domain=NonNegativeReals,
            initialize=100.0,
            doc="Component molar concentrations " "[mol/m^3]",
            units=units.mol / units.m**3,
        )

        if self.config.defined_state is False:
            self.conc_water_eqn = Constraint(
                expr=self.conc_mol_comp["H2O"] == self.params.dens_mol
            )

    def get_material_flow_terms(b, p, j):
        return b.flow_vol * b.conc_mol_comp[j]

    def get_enthalpy_flow_terms(b, p):
        return (
            b.flow_vol
            * b.params.dens_mol
            * b.params.cp_mol
            * (b.temperature - b.params.temperature_ref)
        )

    def get_material_density_terms(b, p, j):
        return b.conc_mol_comp[j]

    def get_energy_density_terms(b, p):
        return (
            b.params.dens_mol
            * b.params.cp_mol
            * (b.temperature - b.params.temperature_ref)
        )

    def default_material_balance_type(self):
        return MaterialBalanceType.componentPhase

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def define_state_vars(b):
        return {
            "flow_vol": b.flow_vol,
            "conc_mol_comp": b.conc_mol_comp,
            "temperature": b.temperature,
            "pressure": b.pressure,
        }

    def define_display_vars(b):
        return {
            "Volumetric Flowrate": b.flow_vol,
            "Molar Concentration": b.conc_mol_comp,
            "Temperature": b.temperature,
            "Pressure": b.pressure,
        }

    def get_material_flow_basis(b):
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

        # Check pressure bounds
        if value(blk.pressure) < blk.pressure.lb:
            _log.error("{} Pressure set below lower bound.".format(blk.name))
        if value(blk.pressure) > blk.pressure.ub:
            _log.error("{} Pressure set above upper bound.".format(blk.name))
