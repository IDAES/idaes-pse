"""
Flowsheet combining absorber and stripper sections for Monoethanolamine solvent carbon capture system
"""

# Python imports
import copy
import math
import logging
import pkgutil

# Pyomo imports
import pyomo.environ as pyo
from pyomo.network import Arc, Port
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.common.config import ConfigValue, ConfigDict, In, Bool, ListOf, PositiveInt

# IDAES imports
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlockData, declare_process_block_class

from idaes.models.unit_models import Mixer, MomentumMixingType, Translator, StreamScaler

from idaes.models_extra.column_models.solvent_reboiler import SolventReboiler
from idaes.models_extra.column_models.solvent_condenser import SolventCondenser
from idaes.models_extra.column_models.MEAsolvent_column import _fix_vars, _restore_fixedness

from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale

from mea_stripper_reformulated import MEAStripperFlowsheet
from mea_absorber_reformulated import MEAAbsorberFlowsheet
from mea_properties import state_bounds_default
from MEA_NGCC_integrated_system_cost_module import get_ngcc_solvent_cost

def initialize_uninitialized_vars(blk):
    for var in blk.component_data_objects(ctype=pyo.Var, descend_into=True):
        if var.value is None:
            var.set_value(0)

@declare_process_block_class("MEACombinedFlowsheet")
class MEACombinedFlowsheetData(FlowsheetBlockData):
    CONFIG = FlowsheetBlockData.CONFIG()
    CONFIG.declare(
        "stripper_finite_element_set",
        ConfigValue(
            domain=ListOf(float),
            description="List containing coordinates finite element faces "
                        "of stripper column. Coordinates must start with zero, be "
                        "strictly increasing, and end with one",
        ),
    )
    CONFIG.declare(
        "absorber_finite_element_set",
        ConfigValue(
            domain=ListOf(float),
            description="List containing coordinates finite element faces "
                        "of absorber column. Coordinates must start with zero, be "
                        "strictly increasing, and end with one",
        ),
    )
    CONFIG.declare(
        "absorber_surrogate_enhancement_factor_model",
        ConfigValue(
            default=None,
            description="Placeholder",
            doc="""Placeholder""",
        ),
    )
    CONFIG.declare(
        "stripper_surrogate_enhancement_factor_model",
        ConfigValue(
            default=None,
            description="Placeholder",
            doc="""Placeholder""",
        ),
    )
    CONFIG.declare(
        "property_package_bounds",
        state_bounds_default()
    )
    def build(self):
        super().build()
        self._add_units()
        self._add_arcs()
        self._add_performance_math()
        self._set_design_inputs()

    def _add_units(self):
        self.stripper_section = MEAStripperFlowsheet(
            column_finite_element_set=self.config.stripper_finite_element_set,
            property_package_bounds=self.config.property_package_bounds,
            surrogate_enhancement_factor_model=self.config.stripper_surrogate_enhancement_factor_model,
        )
        self.absorber_section = MEAAbsorberFlowsheet(
            column_finite_element_set=self.config.absorber_finite_element_set,
            property_package_bounds=self.config.property_package_bounds,
            surrogate_enhancement_factor_model=self.config.absorber_surrogate_enhancement_factor_model,
        )
        self.flue_gas_feed_scaler = StreamScaler(
            property_package=self.absorber_section.vapor_properties
        )

        self.clean_gas_scaler = StreamScaler(
            property_package=self.absorber_section.vapor_properties
        )
        self.makeup_scaler = StreamScaler(
            property_package=self.absorber_section.liquid_properties_no_ions
        )
        self.distillate_scaler = StreamScaler(
            property_package=self.stripper_section.vapor_properties
        )

        # Use Translator block in order to manage specification of MEA flowrate in total recycle loop
        # Use the liquid property package without ions to reduce number of equations generated
        self.lean_solvent_translator = Translator(
            inlet_property_package=self.absorber_section.liquid_properties_no_ions,
            outlet_property_package=self.absorber_section.liquid_properties_no_ions,
            outlet_state_defined=False
        )

    def _add_arcs(self):
        # Define Arcs (streams)
        self.lean_solvent_stream01 = Arc(
            source=self.stripper_section.lean_solvent,
            destination=self.lean_solvent_translator.inlet
        )
        self.lean_solvent_stream02 = Arc(
            source=self.lean_solvent_translator.outlet,
            destination=self.absorber_section.lean_solvent
        )
        self.rich_solvent_stream = Arc(
            source=self.absorber_section.rich_solvent,
            destination=self.stripper_section.rich_solvent
        )
        self.flue_gas_feed_scaler_stream = Arc(
            source=self.flue_gas_feed_scaler.outlet,
            destination=self.absorber_section.flue_gas_feed
        )
        self.clean_gas_scaler_stream = Arc(
            source=self.absorber_section.clean_gas,
            destination=self.clean_gas_scaler.inlet
        )
        self.makeup_scaler_stream = Arc(
            source=self.makeup_scaler.outlet,
            destination=self.absorber_section.makeup,
        )
        self.distillate_scaler_stream = Arc(
            source=self.stripper_section.distillate,
            destination=self.distillate_scaler.inlet
        )

        # Transform Arcs
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)

        # Add sub-flowsheet level ports for easy connectivity
        self.flue_gas_feed = Port(extends=self.flue_gas_feed_scaler.inlet)
        self.clean_gas = Port(extends=self.clean_gas_scaler.outlet)
        self.makeup = Port(extends=self.makeup_scaler.inlet)
        self.distillate = Port(extends=self.distillate_scaler.outlet)

    def _add_performance_math(self):
        self.h2o_mea_ratio = pyo.Var(self.time, initialize=7.91)

        @self.Constraint(self.time)
        def h2o_mea_ratio_eqn(b, t):
            return (
                    b.h2o_mea_ratio[t]
                    * b.absorber_section.makeup_mixer.outlet.mole_frac_comp[t, "MEA"]
                    == b.absorber_section.makeup_mixer.outlet.mole_frac_comp[t, "H2O"]
            )

        self.lean_loading = pyo.Var(self.time)

        @self.Constraint(self.time)
        def lean_loading_eqn(fs, t):
            return (
                    fs.lean_loading[t]
                    * fs.absorber_section.absorber.liquid_inlet.mole_frac_comp[t, "MEA"]
                    == fs.absorber_section.absorber.liquid_inlet.mole_frac_comp[t, "CO2"]
            )

        self.rich_loading = pyo.Var(self.time)

        @self.Constraint(self.time)
        def rich_loading_eqn(fs, t):
            return (
                    fs.rich_loading[t]
                    * fs.absorber_section.absorber.liquid_outlet.mole_frac_comp[t, "MEA"]
                    == fs.absorber_section.absorber.liquid_outlet.mole_frac_comp[t, "CO2"]
            )

        @self.Constraint(self.time)
        def stripper_reflux_pressure(b, t):
            return (
                    b.stripper_section.reflux_mixer.rich_solvent.pressure[t]
                    == b.stripper_section.reflux_mixer.reflux.pressure[t]
            )
        
        self.number_trains = pyo.Var(initialize=4)
        self.number_trains.fix()

        @self.Constraint()
        def flue_gas_feed_multiplier_eqn(b):
            return b.flue_gas_feed_scaler.multiplier == 1 / b.number_trains
        
        @self.Constraint()
        def clean_gas_multiplier_eqn(b):
            return b.clean_gas_scaler.multiplier == b.number_trains
        
        @self.Constraint()
        def makeup_multiplier_eqn(b):
            return b.makeup_scaler.multiplier == 1 / b.number_trains
        
        @self.Constraint()
        def distillate_multiplier_eqn(b):
            return b.distillate_scaler.multiplier == b.number_trains
        
        # Specify translator block
        @self.lean_solvent_translator.Constraint(self.time)
        def temperature_eqn(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature
        
        @self.lean_solvent_translator.Constraint(self.time)
        def pressure_eqn(b, t):
            return b.properties_in[t].pressure == b.properties_out[t].pressure
        
        @self.lean_solvent_translator.Constraint(self.time, self.lean_solvent_translator.properties_in.component_list)
        def flow_mol_comp_eqn(b, t, j):
            return b.properties_in[t].flow_mol_comp[j] == b.properties_out[t].flow_mol_comp[j]
        # Need to deactivate this constraint in order to specify amount of MEA in recycle loop
        self.lean_solvent_translator.flow_mol_comp_eqn[:,"MEA"].deactivate()

        self.mea_recirculation_rate = pyo.Var(self.time, initialize=746, units=pyo.units.mol/pyo.units.s)
        @self.lean_solvent_translator.Constraint(self.time)
        def mea_recirculation_eqn(b, t):
            return b.properties_out[t].flow_mol_comp["MEA"] == self.mea_recirculation_rate[t]
        self.mea_recirculation_rate.fix()

    def _set_design_inputs(self):
        self.flue_gas_feed.flow_mol.fix(37116.4)  # mol/sec
        self.flue_gas_feed.temperature.fix(313.15)  # K
        self.flue_gas_feed.pressure.fix(105000)  # Pa
        self.flue_gas_feed.mole_frac_comp[0, "CO2"].fix(0.04226)
        self.flue_gas_feed.mole_frac_comp[0, "H2O"].fix(0.05480)
        self.flue_gas_feed.mole_frac_comp[0, "N2"].fix(0.76942 + 0.00920) # Added term is to make up for Argon, which is neglected
        self.flue_gas_feed.mole_frac_comp[0, "O2"].fix(0.12430)

        self.makeup.flow_mol.fix(7500)  # mol/sec
        self.makeup.temperature.fix(313.15)  # K
        self.makeup.pressure.fix(183700)  # Pa
        # Pressure determined by pressure equality with lean solvent stream
        self.makeup.mole_frac_comp[0, "CO2"].fix(1e-12)
        self.makeup.mole_frac_comp[0, "H2O"].fix(1.0)
        self.makeup.mole_frac_comp[0, "MEA"].fix(1e-12)
    
    def calculate_scaling_factors(self):
        for t in self.time:
            iscale.constraint_scaling_transform(self.stripper_reflux_pressure[t], 1e-5)

            sT = iscale.get_scaling_factor(self.lean_solvent_translator.properties_in[t].temperature)
            iscale.constraint_scaling_transform(self.lean_solvent_translator.temperature_eqn[t], sT)
            sP = iscale.get_scaling_factor(self.lean_solvent_translator.properties_in[t].pressure)
            iscale.constraint_scaling_transform(self.lean_solvent_translator.pressure_eqn[t], sP)
            for j in self.lean_solvent_translator.properties_in.component_list:
                sF = iscale.get_scaling_factor(self.lean_solvent_translator.properties_in[t].flow_mol_comp[j])
                iscale.constraint_scaling_transform(self.lean_solvent_translator.flow_mol_comp_eqn[t, j], sF)
            sF = iscale.set_and_get_scaling_factor(self.mea_recirculation_rate[t], 1/750)
            iscale.constraint_scaling_transform(self.lean_solvent_translator.mea_recirculation_eqn[t], sF)

        



    # def print_column_design_parameters(self):
    #     print("\n ******* Printing some results *******")
    #     print("\nColumn diameter: ", pyo.value(self.stripper.diameter_column), "m")
    #     print("Column height: ", pyo.value(self.stripper.length_column), "m")
    #     print("Column volume: ", pyo.value(self.volume_column), "m3")
    #     print("Column volume with heads: ", pyo.value(self.volume_column_withheads), "m3")
    #     print("Column height to diameter ratio: ", pyo.value(self.HDratio))
    #     print("\nSolvent inlet molar flowrate: ", pyo.value(self.stripper.liquid_inlet.flow_mol[0]), "mol/s")
    #     print("L/G ratio: ", pyo.value(self.LGratio))
    #
    # def strip_statevar_bounds_for_stripper_initialization(self):
    #     # Strip variable bounds
    #     for x in self.stripper.liquid_phase.length_domain:
    #         for j in self.stripper.config.liquid_phase.property_package.true_phase_component_set:
    #             self.stripper.liquid_phase.properties[0, x].flow_mol_phase_comp_true[j].setlb(None)
    #             self.stripper.liquid_phase.properties[0, x].flow_mol_phase_comp_true[j].domain = pyo.Reals
    #
    #             self.stripper.liquid_phase.properties[0, x].mole_frac_phase_comp_true[j].setlb(None)
    #
    #         for j in self.stripper.config.liquid_phase.property_package.apparent_species_set:
    #             self.stripper.liquid_phase.properties[0, x].mole_frac_comp[j].setlb(None)
    #
    #         for j in self.stripper.config.liquid_phase.property_package.apparent_phase_component_set:
    #             self.stripper.liquid_phase.properties[0, x].mole_frac_phase_comp.setlb(None)
    #
    #         self.stripper.conc_interface_MEA[0, x].setub(None)
    #         self.stripper.log_conc_interface_MEA[0, x].setub(None)
    #
    #     for x in self.stripper.vapor_phase.length_domain:
    #         for j in self.stripper.config.vapor_phase.property_package.component_list:
    #             self.stripper.vapor_phase.properties[0, x].mole_frac_comp[j].setlb(None)
    #
    #         for j in self.stripper.config.vapor_phase.property_package._phase_component_set:
    #             self.stripper.vapor_phase.properties[0, x].mole_frac_phase_comp[j].setlb(None)

    def initialize_build(
            self,
            outlvl=idaeslog.NOTSET,
            solver="ipopt",
            optarg=None,
            rich_solvent_guess=None,
            lean_solvent_guess=None,
            stripper_boilup_guess=None,
    ):
        solver_obj = get_solver(solver, optarg)
        init_log = idaeslog.getInitLogger(self.name, outlvl)
        solve_log = idaeslog.getSolveLogger(self.name, outlvl)

        def safe_solve(blk):
            assert degrees_of_freedom(blk) == 0
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = solver_obj.solve(blk, tee=slc.tee)
            pyo.assert_optimal_termination(results)

        # Set absorber section design and operating conditions
        # Flue gas inlet
        flue_gas_feed_flags = _fix_vars(self.flue_gas_feed.vars.values())
        # Makeup H2O
        # Flowrate will be calculated by H2O:MEA ratio
        makeup_flags = _fix_vars(self.makeup.vars.values())

        self.flue_gas_feed_scaler.multiplier.value = 1 / self.number_trains.value
        self.clean_gas_scaler.multiplier.value = self.number_trains.value
        self.makeup_scaler.multiplier.value = 1 / self.number_trains.value
        self.distillate_scaler.multiplier.value = self.number_trains.value

        self.flue_gas_feed_scaler.initialize_build(
            outlvl=outlvl, 
            solver=solver,
            optarg=optarg,
        )
        propagate_state(self.flue_gas_feed_scaler_stream, overwrite_fixed=True)
        self.makeup_scaler.initialize_build(
            outlvl=outlvl, 
            solver=solver,
            optarg=optarg,
        )
        propagate_state(self.makeup_scaler_stream, overwrite_fixed=True)

        if lean_solvent_guess is None:
            lean_pump_inlet = self.absorber_section.lean_solvent_pump.inlet
            lean_solvent_guess = ComponentMap()
            lean_solvent_guess[lean_pump_inlet.flow_mol] = 23000/self.number_trains.value  # mol/sec
            lean_solvent_guess[lean_pump_inlet.temperature] = 396  # K
            lean_solvent_guess[lean_pump_inlet.pressure] = 183700  # Pa
            lean_solvent_guess[lean_pump_inlet.mole_frac_comp[:, "CO2"]] = 0.015178
            lean_solvent_guess[lean_pump_inlet.mole_frac_comp[:, "H2O"]] = 0.840049
            lean_solvent_guess[lean_pump_inlet.mole_frac_comp[:, "MEA"]] = 0.144773

        for comp, guess in lean_solvent_guess.items():
            comp.fix(guess)

        # Initialize sub-flowsheets
        self.absorber_section.initialize_build(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg
        )
        propagate_state(self.clean_gas_scaler_stream)
        self.clean_gas_scaler.initialize_build(
            outlvl=outlvl, 
            solver=solver,
            optarg=optarg,
        )

        if rich_solvent_guess is None:
            rich_solvent = self.stripper_section.reflux_mixer.rich_solvent
            rich_solvent_guess = ComponentMap()
            rich_solvent_guess[rich_solvent.flow_mol] = 28000/self.number_trains.value  # mol/sec
            rich_solvent_guess[rich_solvent.temperature] = 366.7  # K
            rich_solvent_guess[rich_solvent.pressure] = 183700  # Pa
            rich_solvent_guess[rich_solvent.mole_frac_comp[:, "CO2"]] = 0.058298
            rich_solvent_guess[rich_solvent.mole_frac_comp[:, "H2O"]] = 0.823349
            rich_solvent_guess[rich_solvent.mole_frac_comp[:, "MEA"]] =  0.118353

        for comp, guess in rich_solvent_guess.items():
            comp.fix(guess)

        self.stripper_reflux_pressure.deactivate()
        self.lean_solvent_stream02_expanded.deactivate()
        self.rich_solvent_stream_expanded.deactivate()

        self.stripper_section.initialize_build(
            outlvl=outlvl,
            solver=solver, 
            optarg=optarg,
            boilup_guess=stripper_boilup_guess
        )
        propagate_state(self.distillate_scaler_stream)
        self.distillate_scaler.initialize_build(
            outlvl=outlvl, 
            solver=solver,
            optarg=optarg,
        )
        propagate_state(self.lean_solvent_stream01)
        self.lean_solvent_translator.initialize_build(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
        )

        init_log.info("Switching specification of reboiler and makeup")
        self.absorber_section.flue_gas_feed.unfix()
        self.absorber_section.makeup.unfix()
        self.stripper_section.reboiler.heat_duty.unfix()
        self.stripper_section.reboiler.bottoms.temperature.fix()

        self.makeup.flow_mol.unfix()
        self.h2o_mea_ratio.fix(7.91)

        safe_solve(self)

        init_log.info("Coupling rich solvent stream")
        self.rich_solvent_stream_expanded.activate()

        self.stripper_section.rich_solvent.unfix()

        safe_solve(self)

        init_log.info("Coupling lean solvent stream")
        self.lean_solvent_stream02_expanded.activate()
        propagate_state(self.lean_solvent_stream02)

        self.absorber_section.lean_solvent.unfix()
        safe_solve(self)

        init_log.info("Final Steps")

        self.absorber_section.rich_solvent_pump.deltaP.unfix()
        self.stripper_reflux_pressure.activate()

        # TODO: Due to the low assumed pressure drop in the HX unit and the lack of pumping cost to elevate
        # solvent to the top of the columns, there is actually no need for the lean solvent pump
        # This should be addressed in future, but for now set deltaP = 0
        self.absorber_section.lean_solvent_pump.deltaP.fix(0)
        safe_solve(self)

    def add_costing(self, export_economic_results=True, overwrite_economic_results=True):
        """
        Method to add costing equations using IDAES Power Plant Costing 
        Framework. Relevant process variables are collected from various 
        process results, and passed to an external economic module where 
        costing blocks are generated as childs of `self`. Note that `self` 
        should be a flowsheet object, which is defined in the script 
        `run_combined_flowsheet`.
        
        Results are exported to CSV unless set to `False`. New results will
        automatically overwrite (delete) old files by default, if `False` the 
        new files will not be generated if old files are found.
        """

        # temporarily fix unfixed process model variables
        vars_temporarily_fixed = []
        for var in self.parent_block().component_data_objects(pyo.Var, descend_into=True):
            if not var.is_fixed():
                vars_temporarily_fixed.append(var)
                var.fix()

        # define hard-coded and scaled variables that don't exist in the model
        # some of these initialize from existing model variables
        # some of these scale by the number of trains when initializing
        # some of these scale with hard-coded reference values (expressions)
        self.costing_setup = pyo.Block()

        # reference total flue gas flow for all trains
        self.costing_setup.ref_fg_feed_flow = pyo.Param(
            initialize=38446.3716099917, mutable=False,
            units=pyo.units.mol/pyo.units.s
            )
        self.costing_setup.fg_feed_ratio = pyo.Var(
            self.time, initialize=1, units=pyo.units.dimensionless
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_fg_feed_ratio(b, t):
            return (b.fg_feed_ratio[t] ==
                    (
                        b.parent_block().flue_gas_feed.flow_mol[t] /
                        b.ref_fg_feed_flow
                    )
                    )

        # reference total stripper condenser vapor outlet for one train
        self.costing_setup.ref_cond_vap_flow = pyo.Param(
            initialize=740.772547565646, mutable=False,
            units=pyo.units.mol/pyo.units.s
            )
        self.costing_setup.cond_vap_ratio = pyo.Var(
            self.time, initialize=1, units=pyo.units.dimensionless
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_cond_vap_ratio(b, t):
            return (b.cond_vap_ratio[t] ==
                    (
                        b.parent_block().stripper_section.condenser
                        .vapor_outlet.flow_mol[t] /
                        b.ref_cond_vap_flow
                        )
                    )

        # NGCC components
        self.costing_setup.feedwater_flowrate = pyo.Var(
            self.time, initialize=478535.935,
            units=pyo.units.kg/pyo.units.hr
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_feedwater_flowrate(b, t):
            return (b.feedwater_flowrate[t] ==
                    (
                        478535.935 * pyo.units.kg/pyo.units.hr *
                        b.fg_feed_ratio[t]
                        )
                    )

        self.costing_setup.fuelgas_flowrate = pyo.Var(
            self.time, initialize=93272.0,
            units=pyo.units.kg/pyo.units.hr
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_fuelgas_flowrate(b, t):
            return (b.fuelgas_flowrate[t] ==
                    (
                        93272.0 * pyo.units.kg/pyo.units.hr *
                        b.fg_feed_ratio[t]
                        )
                    )

        self.costing_setup.fluegas_flowrate = pyo.Var(
            self.time, initialize=1282.72285,
            units=pyo.units.m**3/pyo.units.s
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_fluegas_flowrate(b, t):
            return (b.fluegas_flowrate[t] ==
                    (
                        1282.72285 * pyo.units.m**3/pyo.units.s *
                        b.fg_feed_ratio[t]
                        )
                    )

        # Solvent System components
        # cooling componnents
        # self.costing_setup.absorber_IC1_temperature_in = pyo.Var(
        #     self.time, initialize=49.4774364 + 273.15, units=pyo.units.K
        #     )
        # self.costing_setup.absorber_IC1_temperature_out = pyo.Var(
        #     self.time, initialize=40.13 + 273.15, units=pyo.units.K
        #     )
        # self.costing_setup.absorber_IC1_duty = pyo.Var(
        #     self.time, initialize=-13656639.9, units=pyo.units.W
        #     )
        # self.costing_setup.absorber_IC2_temperature_in = pyo.Var(
        #     self.time, initialize=44.4484669 + 273.15, units=pyo.units.K
        #     )
        # self.costing_setup.absorber_IC2_temperature_out = pyo.Var(
        #     self.time, initialize=43.32 + 273.15, units=pyo.units.K
        #     )
        # self.costing_setup.absorber_IC2_duty = pyo.Var(
        #     self.time, initialize=-1202141.89, units=pyo.units.W
        #     )
        # self.costing_setup.lean_solvent_cooler_temperature_in = pyo.Var(
        #     self.time, initialize=53.4528554 + 273.15, units=pyo.units.K
        #     )
        # self.costing_setup.lean_solvent_cooler_temperature_out = pyo.Var(
        #     self.time, initialize=39.9666668 + 273.15, units=pyo.units.K
        #     )
        # self.costing_setup.lean_solvent_cooler_duty = pyo.Var(
        #     self.time, initialize=-25453567.2, units=pyo.units.W
        #     )

        # turbine components
        self.costing_setup.flue_gas_blower_load = pyo.Var(
            self.time, initialize=8037475.28, units=pyo.units.W
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_flue_gas_blower_load(b, t):
            return (b.flue_gas_blower_load[t] ==
                    (
                        8037475.28 * pyo.units.W *
                        b.fg_feed_ratio[t]
                        )
                    )

        # constant from upstream NGCC/HRSG section, don't need to scale this
        self.costing_setup.IP_LP_crossover_steam_fraction = pyo.Var(
            self.time, initialize=0.662806512, units=pyo.units.dimensionless
            )

        self.costing_setup.flue_gas_flowrate_blower_discharge = pyo.Var(
            self.time, initialize=607.06744, units=pyo.units.m**3/pyo.units.s
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_flue_gas_flowrate_blower_discharge(b, t):
            return (b.flue_gas_flowrate_blower_discharge[t] ==
                    (
                        607.06744 * pyo.units.m**3/pyo.units.s *
                        b.fg_feed_ratio[t]
                        )
                    )

        # wash section components
        # don't need to scale heights
        self.costing_setup.wash_section_packing_height = pyo.Var(
            initialize=6.0, units=pyo.units.m
            )

        # DCC components
        self.costing_setup.dcc_packing_height = pyo.Var(
            initialize=18.0, units=pyo.units.m
            )

        self.costing_setup.dcc_diameter = pyo.Var(
            initialize=16.0, units=pyo.units.m
            )

        # area scales linearly with flow ratio, so diameter scales as sqrt
        @self.costing_setup.Constraint()
        def calculate_dcc_diameter(b):
            return (b.dcc_diameter ==
                    (
                        16.0 * pyo.units.m *
                        pyo.sqrt(b.fg_feed_ratio[0])
                        )
                    )

        self.costing_setup.dcc_duty = pyo.Var(
            self.time, initialize=-81409348.0, units=pyo.units.W
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_dcc_duty(b, t):
            return (b.dcc_duty[t] ==
                    (
                        -81409348.0 * pyo.units.W *
                        b.fg_feed_ratio[t]
                        )
                    )

        # don't need to scale temperatures
        self.costing_setup.dcc_temperature_in = pyo.Var(
            self.time, initialize=54.5472838 + 273.15, units=pyo.units.K
            )
        self.costing_setup.dcc_temperature_out = pyo.Var(
            self.time, initialize=21.0 + 273.15, units=pyo.units.K
            )
        self.costing_setup.dcc_pump_load = pyo.Var(
            self.time, initialize=13234.4015, units=pyo.units.W
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_dcc_pump_load(b, t):
            return (b.dcc_pump_load[t] ==
                    (
                        13234.4015 * pyo.units.W *
                        b.fg_feed_ratio[t]
                        )
                    )

        self.costing_setup.dcc_pump_flowrate = pyo.Var(
            self.time, initialize=0.534461996, units=pyo.units.m**3/pyo.units.s
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_dcc_pump_flowrate(b, t):
            return (b.dcc_pump_flowrate[t] ==
                    (
                        0.534461996 * pyo.units.m**3/pyo.units.s *
                        b.fg_feed_ratio[t]
                        )
                    )

        # pump-around components - flow (kg/hr) / rho (kg/m3) = flow (m3/s)
        # self.costing_setup.pumparound1_flowrate = pyo.Var(
        #     self.time, initialize=(1975177.54/1058.44532) / 3600,
        #     units=pyo.units.m**3/pyo.units.s
        #     )
        # self.costing_setup.pumparound2_flowrate = pyo.Var(
        #     self.time, initialize=(1975177.54/1085.34636) / 3600,
        #     units=pyo.units.m**3/pyo.units.s)

        # Compressor components
        self.costing_setup.CO2_compressor_auxiliary_load = pyo.Var(
            self.time, initialize=self.number_trains.value * 11020626.6,
            units=pyo.units.W
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_CO2_compressor_auxiliary_load(b, t):
            return (b.CO2_compressor_auxiliary_load[t] ==
                    (
                        b.parent_block().number_trains.value *
                        11020626.6 * pyo.units.W *
                        b.cond_vap_ratio[t]
                        )
                    )

        self.costing_setup.CO2_compressor_intercooling_duty = pyo.Var(
            self.time, initialize=self.number_trains.value * -16058955.0,
            units=pyo.units.W
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_CO2_compressor_intercooling_duty(b, t):
            return (b.CO2_compressor_intercooling_duty[t] ==
                    (
                        b.parent_block().number_trains.value *
                        -16058955.0 * pyo.units.W *
                        b.cond_vap_ratio[t]
                        )
                    )

        # define variables/expressions that calculate based on other variables
        # NGCC components
        self.costing_setup.stackgas_flowrate = pyo.Var(
            self.time, initialize=1, units=pyo.units.m**3/pyo.units.s
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_stackgas_flowrate(b, t):
            return (
                b.stackgas_flowrate[t] == pyo.units.convert(
                    b.parent_block().number_trains.value *
                    b.parent_block().absorber_section.stack_gas_flow_vol[t],
                    to_units=pyo.units.m**3/pyo.units.s
                    )
                )

        # emissions components
        self.costing_setup.CO2_capture_rate = pyo.Var(
            self.time, initialize=1, units=pyo.units.lb/pyo.units.hr
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_CO2_capture_rate(b, t):
            return (
                b.CO2_capture_rate[t] == pyo.units.convert(
                    b.parent_block().number_trains.value *
                    (
                        b.parent_block().absorber_section.flue_gas_flow_mass_CO2[t] -
                        b.parent_block().absorber_section.stack_gas_flow_mass_CO2[t]
                        ),
                    to_units=pyo.units.lb/pyo.units.hr
                    )
                )

        self.costing_setup.CO2_emissions_no_cap = pyo.Var(
            self.time, initialize=1, units=pyo.units.lb/pyo.units.hr
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_CO2_emissions_no_cap(b, t):
            return (
                b.CO2_emissions_no_cap[t] == pyo.units.convert(
                    b.parent_block().number_trains.value *
                    b.parent_block().absorber_section.flue_gas_flow_mass_CO2[t],
                    to_units=pyo.units.lb/pyo.units.hr
                    )
                )

        self.costing_setup.CO2_emissions_with_cap = pyo.Var(
            self.time, initialize=1, units=pyo.units.lb/pyo.units.hr
            )

        @self.costing_setup.Constraint(self.time)
        def calculate_CO2_emissions_with_cap(b, t):
            return (
                b.CO2_emissions_with_cap[t] == pyo.units.convert(
                    b.parent_block().number_trains.value *
                    b.parent_block().absorber_section.stack_gas_flow_mass_CO2[t],
                    to_units=pyo.units.lb/pyo.units.hr
                    )
                )

        # solvent makeup components
        self.costing_setup.MW_MEA = pyo.Param(initialize=61.08, mutable=False,
                                              units=pyo.units.g/pyo.units.mol)
        self.costing_setup.MW_H2O = pyo.Param(initialize=18.01528, mutable=False,
                                              units=pyo.units.g/pyo.units.mol)
        self.costing_setup.MW_CO2 = pyo.Param(initialize=44, mutable=False,
                                              units=pyo.units.g/pyo.units.mol)
        self.costing_setup.wMEA_noCO2 = pyo.Param(initialize=0.3, mutable=False,
                                                  units=pyo.units.dimensionless)

        self.costing_setup.lean_solvent_xMEA = pyo.Var(
            self.time, initialize=0.1, units=pyo.units.dimensionless
            )

        # mole fraction MEA in lean solvent
        @self.costing_setup.Constraint(self.time)
        def calculate_lean_solvent_xMEA(b, t):
            return (
                b.lean_solvent_xMEA[t] == (1 + b.parent_block().lean_loading[t] +
                                        (b.MW_MEA/b.MW_H2O) * (1/b.wMEA_noCO2 - 1))
                ** (-1)
                )

        self.costing_setup.lean_solvent_xCO2 = pyo.Var(
            self.time, initialize=0.1, units=pyo.units.dimensionless)

        # mole fraction CO2 in lean solvent
        @self.costing_setup.Constraint(self.time)
        def calculate_lean_solvent_xCO2(b, t):
            return (
                b.lean_solvent_xCO2[t] == b.parent_block().lean_loading[t]
                * b.lean_solvent_xMEA[t]
                )

        self.costing_setup.lean_solvent_xH2O = pyo.Var(
            self.time, initialize=0.1, units=pyo.units.dimensionless)

        # mole fraction H2O in lean solvent
        @self.costing_setup.Constraint(self.time)
        def calculate_lean_solvent_xH2O(b, t):
            return (
                b.lean_solvent_xH2O[t] == 1 - b.lean_solvent_xMEA[t] -
                b.lean_solvent_xCO2[t]
                )

        self.costing_setup.lean_solvent_wMEA = pyo.Var(
            self.time, initialize=0.1, units=pyo.units.dimensionless)

        # mass fraction MEA in lean solvent
        @self.costing_setup.Constraint(self.time)
        def calculate_lean_solvent_wMEA(b, t):
            return (
                b.lean_solvent_wMEA[t] == b.lean_solvent_xMEA[t] * b.MW_MEA /
                (
                    b.lean_solvent_xMEA[t] * b.MW_MEA +
                    b.lean_solvent_xCO2[t] * b.MW_CO2 +
                    b.lean_solvent_xH2O[t] * b.MW_H2O
                )
                )

        self.costing_setup.solvent_makeup_flow = pyo.Var(
            self.time, initialize=1, units=pyo.units.kg/pyo.units.s)

        # MEA makeup = expected flow from loading - return flow from stripper
        @self.costing_setup.Constraint(self.time)
        def calculate_solvent_makeup_flow(b, t):
            return (
                b.solvent_makeup_flow[t] ==
                pyo.units.convert(
                    (1/b.parent_block().lean_loading[t]) *  # we want mol MEA / mol CO2
                    b.MW_MEA/b.MW_CO2 *  # this gives mass MEA / mass CO2
                    b.parent_block().stripper_section.distillate_flow_mass_CO2[t],
                    to_units=pyo.units.kg/pyo.units.s) -
                pyo.units.convert(
                    b.parent_block().stripper_section.lean_solvent_flow_mass[t] *
                    b.lean_solvent_wMEA[t],
                    to_units=pyo.units.kg/pyo.units.s
                    )
                    )

        self.costing_setup.solvent_makeup_loading = pyo.Var(
            self.time, initialize=0.1, units=pyo.units.kg/pyo.units.tonne)

        @self.costing_setup.Constraint(self.time)
        def calculate_solvent_makeup_loading(b, t):
            return (
                b.solvent_makeup_loading[t] ==
                pyo.units.convert(
                    b.solvent_makeup_flow[t],
                    to_units=pyo.units.kg/pyo.units.s
                    ) /
                pyo.units.convert(
                    b.CO2_capture_rate[t],
                    to_units=pyo.units.tonne/pyo.units.s
                    )
                )

        self.costing_setup.solvent_fill_init = pyo.Var(
            self.time, initialize=1, units=pyo.units.kg)

        @self.costing_setup.Constraint(self.time)
        def calculate_solvent_fill_init(b, t):
            return (
                b.solvent_fill_init[t] ==
                pyo.units.convert(
                    b.parent_block().stripper_section.lean_solvent_flow_mass[t] *
                    b.lean_solvent_wMEA[t],
                    to_units=pyo.units.kg/pyo.units.hr
                    ) *
                b.parent_block().number_trains.value * 5 * pyo.units.hr
                )

        # pre-solve setup variables and constraints
        print("\nCosting setup variables and constraints built, solving...\n")
        results = get_solver().solve(self.costing_setup, tee=True)

        # unfix frozen process model variables
        for var in vars_temporarily_fixed:
            var.unfix()

        # define process variable lists for economic inputs
        print("\nCosting setup solved, building main costing blocks...\n")
        args_ngcc = [
            self.costing_setup.feedwater_flowrate[0],
            self.costing_setup.fuelgas_flowrate[0],
            self.costing_setup.fluegas_flowrate[0],
            self.costing_setup.stackgas_flowrate[0]
            ]

        # 0 values in lieu of components that don't exist at the moment
        # may update later as components are added
        args_solvent = [
            self.costing_setup.CO2_capture_rate[0],
            self.number_trains,
            self.absorber_section.absorber.length_column,
            self.absorber_section.absorber.diameter_column,
            self.stripper_section.stripper.length_column,
            self.stripper_section.stripper.diameter_column,
            self.costing_setup.solvent_makeup_loading[0],
            self.costing_setup.CO2_emissions_no_cap[0],
            self.costing_setup.CO2_emissions_with_cap[0],
            self.costing_setup.solvent_fill_init[0],
            self.stripper_section.reboiler.heat_duty[0],
            self.absorber_section.lean_rich_heat_exchanger.area,
            0,  # self.costing_setup.absorber_IC1_temperature_in[0],
            0,  # self.costing_setup.absorber_IC1_temperature_out[0],
            0,  # self.costing_setup.absorber_IC1_duty[0],
            0,  # self.costing_setup.absorber_IC2_temperature_in[0],
            0,  # self.costing_setup.absorber_IC2_temperature_out[0],
            0,  # self.costing_setup.absorber_IC2_duty[0],
            0,  # self.costing_setup.lean_solvent_cooler_temperature_in[0],
            0,  # self.costing_setup.lean_solvent_cooler_temperature_out[0],
            0,  # self.costing_setup.lean_solvent_cooler_duty[0],
            self.stripper_section.condenser.vapor_phase.properties_in[0].temperature,
            self.stripper_section.condenser.vapor_phase.properties_out[0].temperature,
            self.stripper_section.condenser.heat_duty[0],
            self.stripper_section.reboiler.liquid_phase.properties_out[0].temperature,
            self.costing_setup.flue_gas_blower_load[0],
            self.absorber_section.rich_solvent_pump.work_mechanical[0],
            self.absorber_section.lean_solvent_pump.work_mechanical[0],
            self.costing_setup.IP_LP_crossover_steam_fraction[0],
            self.costing_setup.wash_section_packing_height,
            self.costing_setup.dcc_packing_height,
            self.costing_setup.dcc_diameter,
            self.costing_setup.dcc_duty[0],
            self.costing_setup.dcc_temperature_in[0],
            self.costing_setup.dcc_temperature_out[0],
            self.costing_setup.dcc_pump_load[0],
            self.costing_setup.dcc_pump_flowrate[0],
            self.absorber_section.rich_solvent_flow_vol[0],
            self.stripper_section.lean_solvent_flow_vol[0],
            self.absorber_section.stack_gas_flow_mass[0],
            self.costing_setup.flue_gas_flowrate_blower_discharge[0],
            0,  # self.costing_setup.pumparound1_flowrate[0],
            0,  # self.costing_setup.pumparound2_flowrate[0]
            ]

        args_compr = [
            self.costing_setup.CO2_compressor_auxiliary_load[0],
            self.costing_setup.CO2_compressor_intercooling_duty[0]
            ]

        print()
        get_ngcc_solvent_cost(self, args_ngcc, args_solvent, args_compr,
                              export_economic_results=True,
                              overwrite_economic_results=True)
