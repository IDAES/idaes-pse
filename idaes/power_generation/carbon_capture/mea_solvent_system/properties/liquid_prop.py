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
First Generation (GEN 1) Liquid Property Package.

This property package provides the necessary constraints or expressions for the
liquid phase properties of amine-based (MEA) scrubbing of CO2 acid gas
The GEN 1 MEA model uses the Enhancement factor calculation. MEA is taken to
be non-volatile.

Liquid phase components: Carbondioxide (CO2), Monoethanolamine (MEA), Water (H2O)

"""

# Import Pyomo libraries
from pyomo.environ import (Constraint, Expression, Param, SolverFactory,
                           Reference, PositiveReals, Reals, NonNegativeReals,
                           Set, value, Var, exp, log, units as pyunits)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.config import ConfigValue, In


# Import IDAES libraries
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock,
                        Component,
                        LiquidPhase)
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            solve_indexed_blocks)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.power_generation.carbon_capture.mea_solvent_system.unit_models.column \
    import ProcessType
from idaes.core.util import get_solver
import idaes.logger as idaeslog

__author__ = "Paul Akula, John Eslick"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("LiquidParameterBlock")
class PhysicalParameterData(PhysicalParameterBlock):
    """
    Liquid Phase Property Parameter Block Class

    Contains parameters and indexing sets associated with
    liquid phase properties for amine-based scrubbing process.

    """
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
        super(PhysicalParameterData, self).build()

        self._state_block_class = LiquidStateBlock

        # Create Phase object
        self.Liq = LiquidPhase()

        # Create Component objects
        self.MEA = Component()
        self.CO2 = Component()
        self.H2O = Component()

        # component list for true species
        self.component_list_true = Set(
            initialize=['CO2', 'H2O', 'MEA', 'MEA+', 'MEACOO-', 'HCO3-'])
        # component list for solvent species
        self.component_list_solvent = Set(initialize=['H2O', 'MEA'])
        # list of all diffusing components aside the excess solvent
        self.component_list_diffus = Set(
            initialize=['CO2', 'MEA', 'MEA+', 'MEACOO-'])
        # list of diffusing component into/from the vapor phase
        self.component_list_d = Set(initialize=['CO2', 'H2O'])

        # Reaction index
        self.reaction_idx = Set(initialize=['R1', 'R2'])

        # Thermodynamic reference state
        self.pressure_ref = Param(within=PositiveReals,
                                  mutable=True,
                                  default=101325,
                                  units=pyunits.Pa,
                                  doc='Reference pressure')
        self.temperature_ref = Param(default=298.15,
                                     units=pyunits.K,
                                     doc='Thermodynamic Reference Temperature')

        # Critical Temperature dictionary, unit K (Aprion 2005)
        self.temperature_crit = Param(
            self.component_list,
            mutable=False,
            initialize={'CO2': 304.18, 'H2O': 647.13, 'MEA': 614.45},
            units=pyunits.K,
            doc="Critical Temperature")

        # Mol. weights of vapor component - units = kg/mol.
        self.mw_comp = Param(
            self.component_list,
            mutable=False,
            initialize={'CO2': 0.04401, 'H2O': 0.01802, 'MEA': 0.06108},
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weights of liquid components")

        '''
        Specific Heat capacity of Liquid,kJ/kgK
        Data from Hilliard thesis(1998)

        Cp (kJ/kgC) =  C1 + C2*t + C3*t^2 + C4*t^3 +C5*t^4
        t(C) = T(K) - a, a = 273.15
        C1 = C1
        C2*t = C2*(T-a)
        C3*t^2 = C3*(T^2 -2aT +a^2)
        C4*t^3 = C4*(T^3 -3aT^2 + 3a^2T -a^3)
        C5*t^4 = C5*(T^4 -4aT^3 + 6a^2T^2 -4a^3T + a^4)

        CpkJ/kgK = K1 + K2*T + K3*T^2 + K4*T^3 +K5*T^4
        required for mean cp
        K1 = C1 -aC2 + a^2C3 - a^3C4 + a^4C5
        K2 = C2 -2aC3 + 3a^2C4 -4a^3C5
        K3 = C3 -3aC4 + 6a^2C5
        K4 = C4 -4aC5
        K5 = C5
        '''
        _T = 273.15
        mea_c1 = 2.6161
        mea_c2 = 3.706e-3
        mea_c3 = 3.787e-6
        mea_c4 = 0.0
        mea_c5 = 0.0
        h2o_c1 = 4.2107
        h2o_c2 = -1.696e-3
        h2o_c3 = 2.568e-5
        h2o_c4 = -1.095e-7
        h2o_c5 = 3.038e-10

        mea_K1 = mea_c1 - _T * mea_c2 + _T**2 * \
            mea_c3 - _T**3 * mea_c4 + _T**4 * mea_c5
        mea_K2 = mea_c2 - 2 * _T * mea_c3 + 3 * _T**2 * mea_c4 - 4 * _T**3 * mea_c5
        mea_K3 = mea_c3 - 3 * _T * mea_c4 + 6 * _T**2 * mea_c5
        mea_K4 = mea_c4 - 4 * _T * mea_c5
        mea_K5 = mea_c5
        h2o_K1 = h2o_c1 - _T * h2o_c2 + _T**2 * \
            h2o_c3 - _T**3 * h2o_c4 + _T**4 * h2o_c5
        h2o_K2 = h2o_c2 - 2 * _T * h2o_c3 + 3 * _T**2 * h2o_c4 - 4 * _T**3 * h2o_c5
        h2o_K3 = h2o_c3 - 3 * _T * h2o_c4 + 6 * _T**2 * h2o_c5
        h2o_K4 = h2o_c4 - 4 * _T * h2o_c5
        h2o_K5 = h2o_c5

        cp_param_dict = {
            ('MEA', 1): mea_K1,
            ('MEA', 2): mea_K2,
            ('MEA', 3): mea_K3,
            ('MEA', 4): mea_K4,
            ('MEA', 5): mea_K5,
            ('H2O', 1): h2o_K1,
            ('H2O', 2): h2o_K2,
            ('H2O', 3): h2o_K3,
            ('H2O', 4): h2o_K4,
            ('H2O', 5): h2o_K5
        }
        self.cp_param = Param(self.component_list_solvent,
                              range(1, 6),
                              mutable=False,
                              initialize=cp_param_dict,
                              units=pyunits.kJ / (pyunits.kg * pyunits.K),
                              doc="Heat capacity parameters for solvent")

        # Viscosity constants for Aqueous MEA
        # Reference: Morgan et.al (2015)
        visc_d_param_dict = {
            1: -0.0838,
            2: 2.8817,
            3: 33.651,
            4: 1817.0,
            5: 0.00847,
            6: 0.0103,
            7: -2.3890
        }
        self.visc_d_param = Param(range(1, 8),
                                  mutable=False,
                                  initialize=visc_d_param_dict,
                                  units=pyunits.Pa * pyunits.s,
                                  doc="Dynamic viscosity constants")

        # Pure solvent molar volume parameters (MEA and H2O),
        # parameter unit in ml/mol or cm^3/mol
        # Reference: Morgan et.al (2015)
        vol_mol_comp_param_dict = {
            ('MEA', 1): -5.35162e-7,
            ('MEA', 2): -4.51417e-4,
            ('MEA', 3): 1.19451,
            ('H2O', 1): -3.2484e-6,
            ('H2O', 2): 0.00165,
            ('H2O', 3): 0.793,
        }
        self.vol_mol_comp_param = Param(self.component_list_solvent,
                                        range(1, 4),
                                        mutable=False,
                                        initialize=vol_mol_comp_param_dict,
                                        units=(pyunits.cm**3) / pyunits.mol,
                                        doc="Pure solvent molar volume constants")

        # Liquid molar volume parameters  Aqueous MEA,
        # parameter unit in ml/mol
        # Reference: Morgan et.al (2015)
        vol_mol_param_dict = {
            1: 10.2074,
            2: -2.2642,
            3: 3.0059,
            4: 207,
            5: -563.3701
        }
        self.vol_mol_param = Param(range(1, 6),
                                   mutable=False,
                                   initialize=vol_mol_param_dict,
                                   units=(pyunits.cm**3) / pyunits.mol,
                                   doc="Liquid molar volume constants")

        # Thermal conductivity parameters for water
        # k_H2O = C1 + C2*T/T_ref + C3(T/T_ref)^2, T_ref =298.15
        # DIPPR
        therm_cond_param_dict = {1: -0.9003,
                                 2: 2.5006,
                                 3: -0.9938
                                 }
        self.therm_cond_param = Param(range(1, 4),
                                      mutable=False,
                                      initialize=therm_cond_param_dict,
                                      units=pyunits.W / pyunits.m / pyunits.K,
                                      doc="Thermal conductivity constants")

        # Diffusion Coefficient(binary) constants, m^2/s
        # CO2 ---> Ying and Eimer (2012)
        # MEA ---> Snijder et al. (1993)
        # MEA+(MEAH+) ----> Hoff et al.(2004)
        # MEACOO-(MEACOO-)----> Hoff et al.(2004)
        diffus_param_dict = {
            ('CO2', 1): 2.35e-6,
            ('CO2', 2): 2.9837e-8,
            ('CO2', 3): -9.7078e-9,
            ('CO2', 4): -2119,
            ('CO2', 5): -20.132,
            ('MEA', 1): -13.275,
            ('MEA', 2): -2198.3,
            ('MEA', 3): -7.8142e-5,
            ('MEA', 4): 0.0,
            ('MEA', 5): 0.0,
            ('MEA+', 1): -22.64,
            ('MEA+', 2): -1000.0,
            ('MEA+', 3): -0.7,
            ('MEA+', 4): 0.0,
            ('MEA+', 5): 0.0,
            ('MEACOO-', 1): -22.64,
            ('MEACOO-', 2): -1000.0,
            ('MEACOO-', 3): -0.7,
            ('MEACOO-', 4): 0.0,
            ('MEACOO-', 5): 0.0
        }
        self.diffus_param = \
            Param(self.component_list_diffus, range(1, 6),
                  initialize=diffus_param_dict,
                  units=pyunits.m**2 / pyunits.s,
                  doc=" Diffusivity parameters")

        # Concentration-based Equilibruim constants -
        # Reference: Morgan et al 2017
        k_eq_conc_param_dict = {('R1', 1): 233.4,
                                ('R1', 2): -3410,
                                ('R1', 3): -36.8,
                                ('R2', 1): 176.72,
                                ('R2', 2): -2909,
                                ('R2', 3): -28.46
                                }
        self. k_eq_conc_param = Param(self.reaction_idx,
                                      range(1, 4),
                                      mutable=False,
                                      initialize=k_eq_conc_param_dict,
                                      units=pyunits.m**3 / pyunits.kmol,
                                      doc="Equilibrium constant parameters")

        # Surface tension constants, N/m
        # CO2 ---> Morgan et al. (2015)
        # MEA & H2O ---> Aprion(2005)
        # surface tension F parameter --->Morgan et al. (2015)
        surf_tens_param_dict = {
            ('CO2', 1): -5.987,
            ('CO2', 2): 3.7699,
            ('CO2', 3): -0.43164,
            ('CO2', 4): 0.018155,
            ('CO2', 5): -0.01207,
            ('CO2', 6): 0.002119,
            ('MEA', 1): 0.09945,
            ('MEA', 2): 1.067,
            ('MEA', 3): 0.0,
            ('MEA', 4): 0.0,
            ('MEA', 5): 0.0,
            ('MEA', 6): 0.0,
            ('H2O', 1): 0.18548,
            ('H2O', 2): 2.717,
            ('H2O', 3): -3.554,
            ('H2O', 4): 2.047,
            ('H2O', 5): 0.0,
            ('H2O', 6): 0.0
        }
        self.surf_tens_param = \
            Param(self.component_list, range(1, 7),
                  initialize=surf_tens_param_dict,
                  units=pyunits.N / pyunits.m,
                  doc=" Surface tension parameters")

        surf_tens_F_param_dict = {
            1: 2.4558,
            2: -1.5311,
            3: 3.4994,
            4: -5.6398,
            5: 10.2109,
            6: 2.3122,
            7: 4.5608,
            8: -2.3924,
            9: 5.3324,
            10: -12.0494
        }
        self.surf_tens_F_param = \
            Param(range(1, 11),
                  initialize=surf_tens_F_param_dict,
                  doc="Surface tension F parameters")

        # saturated vapor pressure parameters
        pressure_sat_param_dict = {
            ('MEA', 1): 172.78,
            ('MEA', 2): -13492,
            ('MEA', 3): -21.914,
            ('MEA', 4): 1.38e-5,
            ('H2O', 1): 72.55,
            ('H2O', 2): -7206.70,
            ('H2O', 3): -7.1385,
            ('H2O', 4): 4.05e-6
        }
        self.pressure_sat_param = \
            Param(self.component_list_solvent, range(1, 5),
                  initialize=pressure_sat_param_dict,
                  units=pyunits.Pa,
                  doc="Vapor pressure  parameters")

        # Heat of vaporization parameters of water
        heat_vaporization_param_dict = {
            1: 56.6,
            2: 0.61204,
            3: -0.6267,
            4: 0.3988}

        self.heat_vaporization_param = \
            Param(range(1, 5),
                  initialize=heat_vaporization_param_dict,
                  units=pyunits.kJ / pyunits.mol,
                  doc="Heat of vaporization parameters of water")

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
            'flow_mol':
            {'method': None, 'units': 'mol/s'},
            'pressure':
            {'method': None, 'units': 'Pa'},
            'temperature':
            {'method': None, 'units': 'K'},
            'mole_frac_comp':
            {'method': None, 'units': None},
            'flow_mol_comp':
            {'method': '_flow_mol_comp', 'units': 'mol/s'},
            'mw':
            {'method': '_mw', 'units': 'kg/mol'},
            'mass_frac_comp':
            {'method': '_mass_frac_comp', 'units': None},
            'mass_frac_co2_free':
            {'method': '_mass_frac_co2_free', 'units': None},
            'vol_mol':
            {'method': '_vol_mol', 'units': 'm^3/mol'},
            'vol_mol_comp':
            {'method': '_vol_mol_comp', 'units': 'cm^3/mol'},
            'conc_mol':
            {'method': '_conc_mol', 'units': 'mol/m^3'},
            'conc_mol_comp':
            {'method': '_conc_mol_comp', 'units': 'mol/m^3'},
            'lnk_eq_conc':
            {'method': '_lnk_eq_conc', 'units': 'm^3/mol'},
            'dens_mass':
            {'method': '_dens_mass', 'units': 'kg/m^3'},
            'cp_mass_comp':
            {'method': '_cp_mass_comp', 'units': 'J/kg.K'},
            'cp_mol':
            {'method': '_cp_mol', 'units': 'J/mol.K'},
            'cp_mol_mean':
            {'method': '_cp_mol_mean', 'units': 'J/mol.K'},
            'cp_mol_comp':
            {'method': '_cp_mol_comp', 'units': 'J/mol.K'},
            'cp_mol_comp_mean':
            {'method': '_cp_mol_comp_mean', 'units': 'J/mol.K'},
            'cp_mass_comp_mean':
            {'method': '_cp_mass_comp_mean', 'units': 'J/kg.K'},
            'enth_mean':
            {'method': '_enth_mean', 'units': 'J/s'},
            'enth_liq_density':
            {'method': '_enth_liq_density', 'units': 'J/m^3'},
            'diffus':
            {'method': '_diffus', 'units': 'm^2/s'},
            'visc_d':
            {'method': '_visc_d', 'units': 'kg/m.s'},
            'thermal_cond':
            {'method': '_thermal_cond', 'units': 'W/m.K'},
            'surf_tens':
            {'method': '_surf_tens', 'units': 'N/m'},
            'pressure_sat':
            {'method': '_pressure_sat', 'units': 'Pa'},
            'k2_rxn':
            {'method': '_k2_rxn', 'units': 'm3/mol.s'},
            'habs':
            {'method': '_habs', 'units': 'J/mol'},
            'hvap':
            {'method': '_hvap', 'units': 'J/mol'},
            'henry_N2O_analogy':
            {'method': '_henry_N2O_analogy', 'units': 'Pa.m3/mol'}
        })

        obj.add_default_units({'time': pyunits.s,
                               'length': pyunits.m,
                               'mass': pyunits.kg,
                               'amount': pyunits.mol,
                               'temperature': pyunits.K})


class LiquidStateBlockMethods(StateBlock):
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

        init_log.info('Starting Liquid phase properties initialization')

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
            for j in blk[k].params.component_list_solvent:
                if hasattr(blk[k], "cp_mass_comp_eqn"):
                    calculate_variable_from_constraint(blk[k].cp_mass_comp[j],
                                                       blk[k].cp_mass_comp_eqn[j])
                if hasattr(blk[k], "cp_mol_comp_eqn"):
                    calculate_variable_from_constraint(blk[k].cp_mol_comp[j],
                                                       blk[k].cp_mol_comp_eqn[j])

                if hasattr(blk[k], "flow_mol_comp_eqn"):
                    calculate_variable_from_constraint(blk[k].flow_mol_comp[j],
                                                       blk[k].flow_mol_comp_eqn[j])

                if hasattr(blk[k], "cp_mass_comp_mean_eqn"):
                    calculate_variable_from_constraint(
                        blk[k].cp_mass_comp_mean[j],
                        blk[k].cp_mass_comp_mean_eqn[j])

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

            if hasattr(blk[k], "enth_mean_eqn"):
                calculate_variable_from_constraint(
                    blk[k].enth_mean,
                    blk[k].enth_mean_eqn)

        # Initialise values of true species
        for k in blk.keys():
            if hasattr(blk[k], "speciation_model"):
                blk[k].conc_mol_comp_true['CO2'].value = 0.01539
                blk[k].conc_mol_comp_true['MEA'].value = 2032
                blk[k].conc_mol_comp_true['H2O'].value = 3.209e+04
                blk[k].conc_mol_comp_true['MEA+'].value = 2007
                blk[k].conc_mol_comp_true['HCO3-'].value = 24.43
                blk[k].conc_mol_comp_true['MEACOO-'].value = 1983

        # Solve property block if non-empty
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solve_indexed_blocks(opt, [blk], tee=slc.tee)

        init_log.info("Liquid properties initialization complete {}.".format(
            idaeslog.condition(res)))

        # ---------------------------------------------------------------------
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


@declare_process_block_class("LiquidStateBlock",
                             block_class=LiquidStateBlockMethods)
class LiquidStateBlockData(StateBlockData):
    """
    Liquid phase property package of amine-based scrubbing of CO2
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(LiquidStateBlockData, self).build()

        # Object reference for molecular weight if needed by CV1D
        # Molecular weights
        self.mw_comp = Reference(self.params.mw_comp)

        self._make_state_vars()
        self._make_speciation_model()

    def _make_speciation_model(self):
        self.conc_mol_comp_true = Var(self.params.component_list_true,
                                      domain=NonNegativeReals,
                                      initialize=10,
                                      units=pyunits.mol / pyunits.m**3,
                                      doc="Concentration of true species ")

        def rule_logC(blk, i):
            return log(blk.conc_mol_comp_true[i])

        self.logC = Expression(self.params.component_list_true,
                               rule=rule_logC,
                               doc='Logarithm of true species concentration')

        def rule_speciation_model(blk, r):
            xa = blk.conc_mol_comp
            xt = blk.conc_mol_comp_true
            CT = blk.logC
            lnk = blk.lnk_eq_conc

            # chemical equilibruim reaction 1
            if r == 1:
                return lnk['R1'] == \
                    CT['MEA+'] + CT['MEACOO-'] - 2 * CT['MEA'] - CT['CO2']
            # chemical equilibruim reaction 2
            elif r == 2:
                return lnk['R2'] == \
                    CT['MEA+'] + CT['HCO3-'] - \
                    CT['MEA'] - CT['CO2'] - CT['H2O']
            # MEA balance
            elif r == 3:
                return xa['MEA'] == xt['MEA'] + xt['MEA+'] + xt['MEACOO-']
            # CO2 balance
            elif r == 4:
                return xa['CO2'] == xt['CO2'] + xt['HCO3-'] + xt['MEACOO-']
            # H2O balance
            elif r == 5:
                return xa['H2O'] == xt['H2O'] + xt['HCO3-']
            # charge balance
            elif r == 6:
                return xt['MEA+'] == xt['HCO3-'] + xt['MEACOO-']

        self.speciation_model = Constraint([1, 2, 3, 4, 5, 6],
                                           rule=rule_speciation_model)

    def _make_state_vars(self):

        self.flow_mol = Var(initialize=1.0,
                            domain=NonNegativeReals,
                            units=pyunits.mol / pyunits.s,
                            doc='Total molar flowrate')
        self.mole_frac_comp = Var(self.component_list,
                                  domain=NonNegativeReals,
                                  bounds=(0, 1),
                                  initialize=1 /
                                  len(self.component_list),
                                  units=None,
                                  doc='Component mole fractions')
        self.pressure = Var(initialize=101325,
                            domain=NonNegativeReals,
                            units=pyunits.Pa,
                            doc='Pressure [Pa]')
        self.temperature = Var(initialize=298.15,
                               domain=NonNegativeReals,
                               units=pyunits.K,
                               doc='Temperature')

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
                                     doc='Component molar flowrate [mol/s]')
            self.flow_mol_comp_eq = Constraint(self.component_list,
                                               rule=rule_flow_mol_comp,
                                               doc="Component molar flow in"
                                                   "liquid phase [mol/s]")
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.flow_mol_comp)
            self.del_component(self.flow_mol_comp_eq)
            raise

    def _mw(self):
        def rule_mw(blk):
            return sum(blk.mole_frac_comp[j] * blk.mw_comp[j]
                       for j in blk.component_list)

        try:
            self.mw = Expression(rule=rule_mw,
                                 doc="Average molecular weight of  liquid phase"
                                 "components")
        except AttributeError:
            self.del_component(self.mw)
            raise

    def _mass_frac_comp(self):
        def rule_mass_frac_comp(blk, j):
            return blk.mole_frac_comp[j] * blk.mw_comp[j] / blk.mw

        try:
            self.mass_frac_comp = Expression(self.component_list,
                                             rule=rule_mass_frac_comp,
                                             doc="Mass fraction"
                                             "[-]")
        except AttributeError:
            self.del_component(self.mass_frac_comp)
            raise

    def _mass_frac_co2_free(self):
        '''
        mass fraction on CO2 free-basis
        '''
        def rule_mass_frac_co2_free(blk, j):
            return blk.mole_frac_comp[j] * blk.mw_comp[j] /\
                sum(blk.mole_frac_comp[i] * blk.mw_comp[i]
                    for i in blk.params.component_list_solvent)

        try:
            self.mass_frac_co2_free = Expression(self.params.component_list_solvent,
                                                 rule=rule_mass_frac_co2_free,
                                                 doc="Mass fraction on CO2 free basis"
                                                 "[-]")
        except AttributeError:
            self.del_component(self.mass_frac_co2_free)
            raise

    def _vol_mol_comp(self):
        # UNIT: cm**3/mol,
        def rule_vol_mol_comp(blk, j):
            T = blk.temperature
            mw = blk.mw_comp
            c = blk._params.vol_mol_comp_param
            return mw[j] * 1e3 / (c[j, 1] * T**2 + c[j, 2] * T + c[j, 3])

        try:
            self.vol_mol_comp = Expression(self.params.component_list_solvent,
                                           rule=rule_vol_mol_comp,
                                           doc="Pure solvent molar volume")
        except AttributeError:
            self.del_component(self.vol_mol_comp)
            raise

    def _vol_mol(self):
        # Liquid phase molar volume m^3/mol
        def rule_vol_mol(blk):
            v_mea = self.vol_mol_comp['MEA']
            v_h20 = self.vol_mol_comp['H2O']
            a = blk._params.vol_mol_param[1]
            b = blk._params.vol_mol_param[2]
            c = blk._params.vol_mol_param[3]
            d = blk._params.vol_mol_param[4]
            e = blk._params.vol_mol_param[5]
            v_soln = (blk.mole_frac_comp['MEA'] * v_mea + blk.mole_frac_comp['H2O']
                      * v_h20 + a * blk.mole_frac_comp['CO2'] +
                      (b + c * blk.mole_frac_comp['MEA'])
                      * blk.mole_frac_comp['MEA'] * blk.mole_frac_comp['H2O'] +
                      (d + e * blk.mole_frac_comp['MEA']) *
                      blk.mole_frac_comp['MEA'] * blk.mole_frac_comp['CO2']) * 1e-6
            return v_soln

        try:
            self.vol_mol = Expression(rule=rule_vol_mol,
                                      doc="Liquid molar volume")
        except AttributeError:
            self.del_component(self.vol_mol)
            raise

    def _conc_mol(self):
        # Liquid phase concentration
        def rule_conc_mol(b):
            return 1 / self.vol_mol

        try:
            self.conc_mol = Expression(rule=rule_conc_mol,
                                       doc="Concentration [mol/m3]")
        except AttributeError:
            self.del_component(self.conc_mol)
            raise

    def _conc_mol_comp(self):
        def rule_conc_mol_comp(blk, i):
            return blk.conc_mol * blk.mole_frac_comp[i]

        try:
            self.conc_mol_comp = Expression(self.component_list,
                                            rule=rule_conc_mol_comp,
                                            doc="Concentration of liquid components")
        except AttributeError:
            self.del_component(self.conc_mol_comp)
            raise

    def _lnk_eq_conc(self):
        '''
        Natural logarithm of concentration based equilibruim constant [m**3/mol]
        for chemical equilibruim in the liquid phase.
        '''
        def rule_lnk_eq_conc(blk, rxn):
            return (self._params.k_eq_conc_param[rxn, 1] +
                    self._params.k_eq_conc_param[rxn, 2] / blk.temperature +
                    self._params.k_eq_conc_param[rxn, 3] * log(blk.temperature) +
                    log(1e-3))

        try:
            self.lnk_eq_conc = Expression(self._params.reaction_idx,
                                          rule=rule_lnk_eq_conc,
                                          doc="Natural log of concentration"
                                          "based equilibruim constant [m3/mol]")
        except AttributeError:
            self.del_component(self.lnk_eq_conc)
            raise

    def _dens_mass(self):
        def rule_dens_mass(blk):
            return blk.mw * blk.conc_mol

        try:
            self.dens_mass = Expression(rule=rule_dens_mass,
                                        doc="mass density [kg/m^3]")
        except AttributeError:
            self.del_component(self.dens_mass)
            raise

    def _cp_mass_comp(self):
        self.cp_mass_comp = Var(self.params.component_list_solvent,
                                domain=Reals,
                                initialize=1.0,
                                units=pyunits.J / pyunits.K / pyunits.kg,
                                doc="Pure component specific heat capacities ")

        def rule_cp_mass_comp(blk, j):
            return blk.cp_mass_comp[j] == 1e3 * (
                blk._params.cp_param[j, 1] +
                blk._params.cp_param[j, 2] * (blk.temperature) +
                blk._params.cp_param[j, 3] * (blk.temperature)**2 +
                blk._params.cp_param[j, 4] * (blk.temperature)**3 +
                blk._params.cp_param[j, 5] * (blk.temperature)**4)

        try:
            # Try to build constraint
            self.cp_mass_comp_eqn = Constraint(self.params.component_list_solvent,
                                               rule=rule_cp_mass_comp)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mass_comp)
            self.del_component(self.cp_mass_comp_eqn)
            raise

    def _cp_mol_comp(self):
        # Pure component liquid heat capacities J/mol.K
        self.cp_mol_comp = Var(self.params.component_list_solvent,
                               domain=Reals,
                               initialize=1.0,
                               units=pyunits.J / pyunits.K / pyunits.mol,
                               doc="Pure component liquid heat capacities ")

        def rule_cp_mol_comp(blk, j):
            return blk.cp_mol_comp[j] == blk.cp_mass_comp[j] * blk.mw_comp[j]

        try:
            # Try to build constraint
            self.cp_mol_comp_eqn = Constraint(self.params.component_list_solvent,
                                              rule=rule_cp_mol_comp)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol_comp)
            self.del_component(self.cp_mol_comp_eqn)
            raise

    def _cp_mol(self):
        self.cp_mol = Var(domain=Reals,
                          initialize=1.0,
                          units=pyunits.J / pyunits.K / pyunits.mol,
                          doc="Liquid molar heat capacity ")

        def rule_cp_mol(blk):
            return blk.cp_mol == sum(blk.cp_mol_comp[j] * blk.mole_frac_comp[j]
                                     for j in blk.params.component_list_solvent)

        try:
            # Try to build constraint
            self.cp_mol_eqn = Constraint(rule=rule_cp_mol)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol)
            self.del_component(self.cp_mol_eqn)
            raise

    def _cp_mass_comp_mean(self):
        # mean Pure component liquid heat capacities
        self.cp_mass_comp_mean = Var(self.params.component_list_solvent,
                                     domain=Reals,
                                     initialize=1.0,
                                     units=pyunits.J / pyunits.K / pyunits.kg,
                                     doc="Mean Pure component specific heat capacities "
                                         "between T and T_ref ")

        def rule_cp_mass_comp_mean(blk, j):
            t = blk.temperature / blk._params.temperature_ref
            t0 = blk._params.temperature_ref
            return blk.cp_mass_comp_mean[j] == 1e3 * (
                blk._params.cp_param[j, 1] +
                blk._params.cp_param[j, 2] / 2 * t0 * (t + 1) +
                blk._params.cp_param[j, 3] / 3 * (t0**2) * (t**2 + t + 1) +
                blk._params.cp_param[j, 4] / 4 * (t0**3) * (t + 1) * (t**2 + 1) +
                blk._params.cp_param[j, 5] / 5 * (t0**4) * (t**4 + t**3 + t**2 + t + 1))

        try:
            # Try to build constraint
            self.cp_mass_comp_mean_eqn =\
                Constraint(self.params.component_list_solvent,
                           rule=rule_cp_mass_comp_mean)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mass_comp_mean)
            self.del_component(self.cp_mass_comp_mean_eqn)
            raise

    def _cp_mol_comp_mean(self):
        # average Pure component liquid heat capacities btw T and T_ref
        self.cp_mol_comp_mean = Var(self.params.component_list_solvent,
                                    domain=Reals,
                                    initialize=1.0,
                                    units=pyunits.J / pyunits.K / pyunits.mol,
                                    doc="Mean pure component liq. heat capacities "
                                        "between T and T_ref ")

        def rule_cp_mol_comp_mean(blk, j):
            return blk.cp_mol_comp_mean[j] == (blk.cp_mass_comp_mean[j] *
                                               blk.mw_comp[j])

        try:
            # Try to build constraint
            self.cp_mol_comp_mean__eqn = Constraint(self.params.component_list_solvent,
                                                    rule=rule_cp_mol_comp_mean)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol_comp_mean)
            self.del_component(self.cp_mol_comp_mean_eqn)
            raise

    def _cp_mol_mean(self):
        # Average liquid heat capacities btw T and T_ref
        self.cp_mol_mean = Var(domain=Reals,
                               initialize=1.0,
                               units=pyunits.J / pyunits.K / pyunits.mol,
                               doc="Mean molar liquid heat capacities ")

        def rule_cp_mol_mean(blk):
            return blk.cp_mol_mean == sum(blk.cp_mol_comp_mean[j] *
                                          blk.mole_frac_comp[j]
                                          for j in blk.params.component_list_solvent)

        try:
            # Try to build constraint
            self.cp_mol_mean_eqn = Constraint(rule=rule_cp_mol_mean)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol_mean)
            self.del_component(self.cp_mol_mean_eqn)
            raise

    def _enth_mean(self):
        # Average Liquid Phase Enthalpy  btw T and T_ref
        self.enth_mean = Var(domain=Reals,
                                 initialize=1.0,
                                 units=pyunits.J / pyunits.s,
                                 doc="Mean Liquid Enthalpy flow  btw T and Tref ")

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

    def _enth_liq_density(self):
        #  molar enthalpy holdup per unit volume
        self.enth_liq_density = Var(
            domain=Reals,
            initialize=1.0,
            units=pyunits.J / pyunits.m**3,
            doc=' Enthalpy holdup')
        try:
            # Try to build constraint
            self.eq_enth_liq_density = Constraint(
                expr=self.enth_liq_density == self.enth_mean * self.conc_mol)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.enth_liq_density)
            self.del_component(self.eq_enth_liq_density)
            raise

    def _visc_d(self):
        '''
        Dynamic viscosity of liquid
        '''
        def rule_visc_d(blk):
            r = blk.mass_frac_comp['MEA'] * 100
            T = blk.temperature
            alpha = blk.mole_frac_comp['CO2'] / blk.mole_frac_comp['MEA']
            mu_H2O = 1.002e-3 *\
                10**((1.3272 *
                      (293.15 - T - 0.001053 * (T - 293.15)**2)) / (T - 168.15))
            a = blk._params.visc_d_param[1]
            b = blk._params.visc_d_param[2]
            c = blk._params.visc_d_param[3]
            d = blk._params.visc_d_param[4]
            e = blk._params.visc_d_param[5]
            f = blk._params.visc_d_param[6]
            g = blk._params.visc_d_param[7]

            return mu_H2O * exp(r * (T * (a * r + b) + c * r + d) *
                                (alpha * (e * r + f * T + g) + 1) / T**2)

        try:
            self.visc_d = Expression(rule=rule_visc_d,
                                     doc="Dynamic viscosity [Pa.s]")

        except AttributeError:
            self.del_component(self.visc_d)
            raise

    def _thermal_cond(self):
        '''
        thermal conductivity of fluid[W/m/K]
        source:
        pure components :Sato-Riedel model
        mixing rule: Vredeveld mixing rule
        arguments
        T = temperature of the liquid phase
        Tb_MEA= boiling point of MEA [K]
        Tc_MEA = critical Temperature of MEA [K]
        Tr = reduced Temperature of the liquid
        Tbr_MEA = reduced boiling point temperature of MEA
        K_MEA = thermal conductivity of pure MEA
        K_H2O = thermal conductivity of pure water
        '''
        def rule_thermal_cond(blk):
            Tb_MEA = 443
            Tc_MEA = 614.2
            T = blk.temperature
            Tr = T / Tc_MEA
            Tbr_MEA = Tb_MEA / Tc_MEA
            K_MEA = 1.1053152 / (61.08**0.5) * (3 + 20 * (1 - Tr)
                                                ** (2 / 3)) / (3 + 20 *
                                                               (1 - Tbr_MEA)**(2 / 3))
            K_H2O = 0.6065 * (-1.48445 + 4.12292 * T / 298.15 -
                              1.63866 * (T / 298.15)**2)

            return ((blk.mole_frac_comp['H2O'] * K_H2O**-2 +
                     blk.mole_frac_comp['MEA'] * K_MEA**-2)**-1)**0.5

        try:
            self.thermal_cond = Expression(rule=rule_thermal_cond,
                                           doc='Thermal conductivity of fluid [W/m.K]')
        except AttributeError:
            self.del_component(self.thermal_cond)
            raise

    def _surf_tens(self):
        '''
        Surface tension,N/m
        '''
        def rule_surf_tens_comp(blk, j):
            r = blk.mass_frac_co2_free['MEA']
            T = blk.temperature
            Tc = blk._params.temperature_crit
            if j in blk.params.component_list_solvent:
                return blk._params.surf_tens_param[j, 1] * (1 - T / Tc[j]) **\
                    (blk._params.surf_tens_param[j, 2] +
                     blk._params.surf_tens_param[j, 3] * T / Tc[j] +
                     blk._params.surf_tens_param[j, 4] * (T / Tc[j])**2)
            elif j == 'CO2':
                return blk._params.surf_tens_param[j, 1] * r**2 + \
                    blk._params.surf_tens_param[j, 2] * r + \
                    blk._params.surf_tens_param[j, 3] +\
                    T * (blk._params.surf_tens_param[j, 4] * r**2 +
                         blk._params.surf_tens_param[j, 5] * r +
                         blk._params.surf_tens_param[j, 6])

        try:
            self.surf_tens_comp = Expression(self.component_list,
                                             rule=rule_surf_tens_comp,
                                             doc="Surface tension of components")
        except AttributeError:
            self.del_component(self.surf_tens_comp)
            raise

        def rule_surf_tens(blk):
            r = blk.mass_frac_co2_free['MEA']
            alpha = self.mole_frac_comp['CO2'] / blk.mole_frac_comp['MEA']
            Fa = blk._params.surf_tens_F_param[1]
            Fb = blk._params.surf_tens_F_param[2]
            Fc = blk._params.surf_tens_F_param[3]
            Fd = blk._params.surf_tens_F_param[4]
            Fe = blk._params.surf_tens_F_param[5]
            Ff = blk._params.surf_tens_F_param[6]
            Fg = blk._params.surf_tens_F_param[7]
            Fh = blk._params.surf_tens_F_param[8]
            Fi = blk._params.surf_tens_F_param[9]
            Fj = blk._params.surf_tens_F_param[10]
            return blk.surf_tens_comp['H2O'] + \
                (blk.surf_tens_comp['CO2'] - blk.surf_tens_comp['H2O']) *\
                blk.mole_frac_comp['CO2'] * (Fa + Fb * alpha + Fc * alpha**2 +
                                             Fd * r + Fe * r**2) +\
                (blk.surf_tens_comp['MEA'] - blk.surf_tens_comp['H2O']) *\
                (Ff + Fg * alpha + Fh * alpha**2 + Fi * r + Fj * r**2) *\
                blk.mole_frac_comp['MEA']

        try:
            self.surf_tens = Expression(rule=rule_surf_tens,
                                        doc="Surface tension [N/m] ")
        except AttributeError:
            self.del_component(self.surf_tens)
            raise

    def _pressure_sat(self):
        '''
        saturated vapor pressure of solvents
        '''
        def rule_pressure_sat_comp(blk, j):
            T = blk.temperature
            return (exp(blk._params.pressure_sat_param[j, 1] +
                        blk._params.pressure_sat_param[j, 2] / T +
                        blk._params.pressure_sat_param[j, 3] * log(T) +
                        blk._params.pressure_sat_param[j, 4] * T * T))

        try:
            self.pressure_sat = Expression(self.params.component_list_solvent,
                                           rule=rule_pressure_sat_comp,
                                           doc="Vapor pressure [Pa]")
        except AttributeError:
            self.del_component(self.pressure_sat)
            raise

    def _diffus(self):
        """
        Diffusivity of diffusing components, m**2/s
        """
        def rule_diffus(blk, i):
            if i == 'CO2':
                C_MEA = blk.conc_mol_comp['MEA'] * 1e-3
                return (blk._params.diffus_param[i, 1] +
                        blk._params.diffus_param[i, 2] * C_MEA +
                        blk._params.diffus_param[i, 3] * C_MEA**2) * exp(
                    (blk._params.diffus_param[i, 4] +
                     blk._params.diffus_param[i, 5]) / blk.temperature)
            elif i == 'MEA':
                C_MEA = blk.conc_mol_comp['MEA']
                return exp(blk._params.diffus_param[i, 1] +
                           blk._params.diffus_param[i, 2] / blk.temperature +
                           blk._params.diffus_param[i, 3] * C_MEA)
            elif i == "MEACOO-" or i == 'MEA+':
                return exp(blk._params.diffus_param[i, 1] +
                           blk._params.diffus_param[i, 2] / blk.temperature +
                           blk._params.diffus_param[i, 3] * log(blk.visc_d))

        try:
            self.diffus = Expression(self.params.component_list_diffus,
                                     rule=rule_diffus,
                                     doc='Diffusivity of component i in liquid')
        except AttributeError:
            self.del_component(self.diffus)
            raise

    def _k2_rxn(self):
        '''
        concentration based rate constant for overall reaction rate
        Second order rate contant(Luo et al., 2015)
        '''
        def rule_k2_rxn(blk):
            T = blk.temperature
            C_MEA = blk.conc_mol_comp['MEA'] * 1e-3
            C_H2O = blk.conc_mol_comp['H2O'] * 1e-3
            return (2.003e10 * exp(-4742 / T) * C_MEA +
                    4.147e6 * exp(-3110 / T) * C_H2O) * 1e-3

        try:
            self.k2_rxn = Expression(rule=rule_k2_rxn,
                                     doc="Second order rate contant [m3/(mol.s)]")
        except AttributeError:
            self.del_component(self.k2_rxn)
            raise

    def _habs(self):
        '''
        Heat of absorption of CO2 J/mol
        '''
        def rule_habs(blk):
            if blk.config.parameters.config.process_type == ProcessType.stripper:
                return -97000
            elif blk.config.parameters.config.process_type == ProcessType.absorber:
                return -84000
        try:
            self.habs = Expression(rule=rule_habs,
                                   doc="CO2 heat of absorption")
        except AttributeError:
            self.del_component(self.habs)
            raise

    def _hvap(self):
        '''
        Heat of vaporization of H2O
        '''
        def rule_hvap(blk):
            c = blk._params.heat_vaporization_param
            T = blk.temperature
            Tc = blk._params.temperature_crit['H2O']
            return 1000 * c[1] * (1 - T / Tc)**(c[2] + c[3] *
                                                T / Tc + c[4] * (T / Tc)**2)
        try:
            self.hvap = Expression(rule=rule_hvap,
                                   doc="Heat of vaporization of H2O [J/mol]")
        except AttributeError:
            self.del_component(self.hvap)
            raise

    def _henry_N2O_analogy(self):
        '''
        Henry's constant N2O Analogy Jiru et.al (2012)
        '''
        def rule_henry_N2O_analogy(blk):
            T = blk.temperature
            t = blk.temperature - 273.15
            wt_MEA = blk.mass_frac_co2_free['MEA']
            wt_H2O = blk.mass_frac_co2_free['H2O']
            H_N2O_MEA = 2.448e5 * exp(-1348 / T)
            H_CO2_H2O = 3.52e6 * exp(-2113 / T)
            H_N2O_H2O = 8.449e6 * exp(-2283 / T)
            H_CO2_MEA = H_N2O_MEA * (H_CO2_H2O / H_N2O_H2O)
            lwm = 1.70981 + 0.03972 * t - 4.3e-4 * t**2 - 2.20377 * wt_H2O

            return (exp(wt_MEA * log(H_CO2_MEA) + wt_H2O * log(H_CO2_H2O) +
                        wt_MEA * wt_H2O * lwm))

        try:
            self.henry_N2O_analogy = Expression(
                rule=rule_henry_N2O_analogy,
                doc="Henry's constant [Pa.m3/mol]")
        except AttributeError:
            self.del_component(self.henry_N2O_analogy)
            raise
    # --------------------------------------------------------------------------

    def get_material_flow_terms(self, p, j):
        return self.flow_mol_comp[j]

    def get_enthalpy_flow_terms(self, p):
        return self.enth_mean

    def get_material_density_terms(self, p, j):
        return self.conc_mol_comp[j]

    def get_energy_density_terms(self, p):
        return self.enth_liq_density

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
