"""
Absorber sub-flowsheet for Monoethanolamine solvent carbon capture system
"""
# Python imports
import math
import logging

# Pyomo imports
import pyomo.environ as pyo
from pyomo.network import Arc, Port
from pyomo.common.config import ConfigValue, ListOf

# IDAES imports
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlockData, declare_process_block_class

from idaes.models.unit_models import Mixer, MomentumMixingType, Pump, HeatExchanger
from idaes.models_extra.column_models.MEAsolvent_column import MEAColumn
# from intercooler_absorber import MEAAbsorberIntercooler

from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale

# External IDAES imports (workspace)

from mea_properties import (
    MEALiquidParameterBlock,
    FlueGasParameterBlock,
    scale_mea_liquid_params,
    scale_mea_vapor_params,
    # switch_liquid_to_parmest_params,
    state_bounds_default
)

logging.getLogger('pyomo.repn.plugins.nl_writer').setLevel(logging.ERROR)

@declare_process_block_class("MEAAbsorberFlowsheet")
class MEAAbsorberFlowsheetData(FlowsheetBlockData):
    CONFIG = FlowsheetBlockData.CONFIG()
    CONFIG.declare(
        "column_finite_element_set",
        ConfigValue(
            domain=ListOf(float),
            description="List containing coordinates finite element faces "
            "of absorber column. Coordinates must start with zero, be "
            "strictly increasing, and end with one",
        ),
    )
    CONFIG.declare(
        "property_package_bounds",
        state_bounds_default()
    )
    CONFIG.declare(
        "surrogate_enhancement_factor_model",
        ConfigValue(
            default=None,
            description="Placeholder",
            doc="""Placeholder""",
        ),
    )
    CONFIG.declare(
        "has_absorber_intercooling",
        ConfigValue(
            default=False,
            description="Implements intercooler models within the absorber colulm",
        ),
    )
    def build(self):
        super().build()
        self._add_properties()
        self._add_units()
        self._add_arcs()
        self._add_performance_math()
        self._set_absorber_design_inputs()
        self._scaling()

    def _add_properties(self):
        if len(self.config.property_package_bounds) == 0:
            bounds = None
        else:
            bounds = self.config.property_package_bounds
        self.vapor_properties = FlueGasParameterBlock(state_bounds=bounds)
        self.liquid_properties = MEALiquidParameterBlock(ions=True, state_bounds=bounds)
        self.liquid_properties_no_ions = MEALiquidParameterBlock(ions=False, state_bounds=bounds)

    def _add_units(self):
        fe_set = self.config.column_finite_element_set
        if self.config.has_absorber_intercooling:
            self.absorber = MEAAbsorberIntercooler(
                column_finite_element_set=fe_set,
                liquid_properties=self.liquid_properties,
                vapor_properties=self.vapor_properties,
                # surrogate_enhancement_factor_model=self.config.surrogate_enhancement_factor_model,
            )
        else:
            self.absorber = MEAColumn(
                finite_elements=len(fe_set) - 1,
                length_domain_set=fe_set,
                vapor_phase={"property_package": self.vapor_properties},
                liquid_phase={"property_package": self.liquid_properties},
                # corrected_hx_coeff_eqn=False,
                # surrogate_enhancement_factor_model=self.config.surrogate_enhancement_factor_model,
            )
        
        self.rich_solvent_pump = Pump(property_package=self.liquid_properties)
        
        self.lean_solvent_pump = Pump(property_package=self.liquid_properties)

        self.lean_rich_heat_exchanger = HeatExchanger(
            hot_side={
                "property_package": self.liquid_properties,
                "has_pressure_change": True,
            },
            cold_side={
                "property_package": self.liquid_properties,
                "has_pressure_change": True,
            },
        )

        self.makeup_mixer = Mixer(
            property_package=self.liquid_properties_no_ions,
            inlet_list=["lean_solvent", "h2o_makeup"],
            momentum_mixing_type=MomentumMixingType.minimize,
        )
        
    def _add_arcs(self):
        # Define Arcs (streams)
        self.hx_hot_feed = Arc(
            source=self.lean_solvent_pump.outlet,
            destination=self.lean_rich_heat_exchanger.hot_side_inlet,
        )

        self.absorber_liquid_product = Arc(
            source=self.absorber.liquid_outlet,
            destination=self.rich_solvent_pump.inlet,
        )
        self.hx_cold_feed = Arc(
            source=self.rich_solvent_pump.outlet,
            destination=self.lean_rich_heat_exchanger.cold_side_inlet,
        )
        self.solvent_to_makeup = Arc(
            source=self.lean_rich_heat_exchanger.hot_side_outlet,
            destination=self.makeup_mixer.lean_solvent,
        )
        self.absorber_liquid_feed = Arc(
            source=self.makeup_mixer.outlet,
            destination=self.absorber.liquid_inlet,
        )

        # Transform Arcs
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)

        # Add sub-flowsheet level ports for easy connectivity
        self.flue_gas_feed = Port(extends=self.absorber.vapor_inlet)
        self.clean_gas = Port(extends=self.absorber.vapor_outlet)
        self.rich_solvent = Port(extends=self.lean_rich_heat_exchanger.cold_side_outlet)
        self.lean_solvent = Port(extends=self.lean_solvent_pump.inlet)
        self.makeup = Port(extends=self.makeup_mixer.h2o_makeup)
        
    def _add_performance_math(self):
        self.absorber.co2_capture = pyo.Var(
            self.time,
            initialize=98,
            doc='''CO2 Capture Rate [%]'''
        )

        @self.absorber.Constraint(self.time, doc='''Correlation for CO2 Capture Rate''')
        def co2_capture_eqn(b, t):
            return (b.co2_capture[t] * (b.vapor_inlet.mole_frac_comp[0, "CO2"]
                                        * b.vapor_inlet.flow_mol[0]) == (
                            (b.vapor_inlet.mole_frac_comp[0, "CO2"] * b.vapor_inlet.flow_mol[0] -
                             b.vapor_outlet.mole_frac_comp[0, "CO2"] * b.vapor_outlet.flow_mol[0]) * 100))

        self.HDratio = pyo.Expression(
            expr=self.absorber.length_column / self.absorber.diameter_column,
            doc="Column height to column diameter ratio (-)"
        )
        self.LGratio = pyo.Expression(
            expr=self.absorber.liquid_inlet.flow_mol[0] / self.absorber.vapor_inlet.flow_mol[0],
            doc="Inlet liquid (solvent) flowrate to inlet gas flowrate ratio (-)"
        )
        self.volume_column = pyo.Expression(
            expr=math.pi * self.absorber.length_column * (self.absorber.diameter_column / 2) ** 2,
            doc="Volume of column (cubic m)"
        )
        self.volume_column_withheads = pyo.Expression(
            expr=(
                math.pi * (
                    self.absorber.diameter_column ** 2 * self.absorber.length_column
                ) / 4
                + math.pi / 3 * self.absorber.diameter_column ** 3
            ),
            doc="Volume of column with heads (cubic m)",
        )
        
        if not self.config.has_absorber_intercooling:
            @self.Expression(self.time, doc="Flue gas CO2 mass flow rate")
            def flue_gas_flow_mass_CO2(b, t):
                return(
                    b.absorber.vapor_phase.properties[t, 0].flow_mass_comp["CO2"]
                )
            @self.Expression(self.time, doc="Volumetric flow rate of CO2-poor gas to stack")
            def stack_gas_flow_vol(b, t):
                return (
                    b.absorber.vapor_phase.properties[t, 1.0].flow_vol
                )
            @self.Expression(self.time, doc="Stack mass flow rate")
            def stack_gas_flow_mass(b, t):
                return(
                    b.absorber.vapor_phase.properties[t, 1].flow_mass
                )
            @self.Expression(self.time, doc="Stack CO2 mass flow rate")
            def stack_gas_flow_mass_CO2(b, t):
                return(
                    b.absorber.vapor_phase.properties[t, 1].flow_mass_comp["CO2"]
                )
            @self.Expression(self.time, doc="Rich solvent volumentric flow rate")
            def rich_solvent_flow_vol(b, t):
                return(
                    b.absorber.liquid_phase.properties[t, 0].flow_vol
                )
        else:
            @self.Expression(self.time, doc="Flue gas CO2 mass flow rate")
            def flue_gas_flow_mass_CO2(b, t):
                return(
                    b.absorber.absorber_bed_3.vapor_phase.properties[t, 0].flow_mass_comp["CO2"]
                )
            @self.Expression(self.time, doc="Volumetric flow rate of CO2-poor gas to stack")
            def stack_gas_flow_vol(b, t):
                return (
                    b.absorber.absorber_bed_1.vapor_phase.properties[t, 1.0].flow_vol
                )
            @self.Expression(self.time, doc="Stack mass flow rate")
            def stack_gas_flow_mass(b, t):
                return(
                    b.absorber.absorber_bed_1.vapor_phase.properties[t, 1].flow_mass
                )
            @self.Expression(self.time, doc="Stack CO2 mass flow rate")
            def stack_gas_flow_mass_CO2(b, t):
                return(
                    b.absorber.absorber_bed_1.vapor_phase.properties[t, 1].flow_mass_comp["CO2"]
                )
            @self.Expression(self.time, doc="Rich solvent volumentric flow rate")
            def rich_solvent_flow_vol(b, t):
                return(
                    b.absorber.absorber_bed_3.liquid_phase.properties[t, 0].flow_vol
                )

    def _set_absorber_design_inputs(self):
        # Fix column design variables
        # Absorber diameter
        self.absorber.diameter_column.fix(12) 

        # Absorber length
        if not self.config.has_absorber_intercooling:
            self.absorber.length_column.fix(20)  # meter

        # Fix operating conditions

        # Flue gas inlet
        self.absorber.vapor_inlet.flow_mol.fix(12000)  # mol/sec
        self.absorber.vapor_inlet.temperature.fix(313.15)  # K
        self.absorber.vapor_inlet.pressure.fix(105000)  # Pa
        self.absorber.vapor_inlet.mole_frac_comp[0, "CO2"].fix(0.042)
        self.absorber.vapor_inlet.mole_frac_comp[0, "H2O"].fix(0.054)
        self.absorber.vapor_inlet.mole_frac_comp[0, "N2"].fix(0.770)
        self.absorber.vapor_inlet.mole_frac_comp[0, "O2"].fix(0.13)

        # Makeup H2O
        self.makeup_mixer.h2o_makeup.flow_mol.fix(2500)  # mol/sec
        self.makeup_mixer.h2o_makeup.temperature.fix(313.15)  # K
        self.makeup_mixer.h2o_makeup.pressure.fix(183700)  # Pa
        # Pressure determined by pressure equality with lean solvent stream
        self.makeup_mixer.h2o_makeup.mole_frac_comp[0, "CO2"].fix(1e-5)
        self.makeup_mixer.h2o_makeup.mole_frac_comp[0, "H2O"].fix(1.0)
        self.makeup_mixer.h2o_makeup.mole_frac_comp[0, "MEA"].fix(1e-5)

        # Set operating conditions
        self.rich_solvent_pump.deltaP.fix(103700)
        self.rich_solvent_pump.efficiency_pump.fix(0.9)

        self.lean_solvent_pump.deltaP.fix(100000)
        self.lean_solvent_pump.efficiency_pump.fix(0.9)

        # self.lean_rich_heat_exchanger.area.fix(4630)
        self.lean_rich_heat_exchanger.area.fix(1500)
        self.lean_rich_heat_exchanger.overall_heat_transfer_coefficient.fix(850) # Used in Aspen flowsheet
        self.lean_rich_heat_exchanger.hot_side.deltaP.fix(
            -25000
        )  # These are very low pressure drops
        self.lean_rich_heat_exchanger.cold_side.deltaP.fix(-25000)
        
    def _scaling(self):
        gsf = iscale.get_scaling_factor
        ssf = iscale.set_scaling_factor

        def cst(con, s):
            iscale.constraint_scaling_transform(con, s, overwrite=False)

        scale_mea_vapor_params(self.vapor_properties, scaling_factor_flow_mol=3e-4)
        scale_mea_liquid_params(self.liquid_properties, ions=True, scaling_factor_flow_mol=3e-4)
        scale_mea_liquid_params(self.liquid_properties_no_ions, ions=False, scaling_factor_flow_mol=3e-4)

        # Stripper column
        if not self.config.has_absorber_intercooling:
            column_list = [self.absorber]
        else:
            column_list = [getattr(self.absorber, f'absorber_bed_{i}') for i in range(1,4)]
        
        if not self.config.has_absorber_intercooling:
            column = self.absorber
            for t in self.time:
                for x in column.liquid_phase.length_domain:
                    ssf(column.velocity_liq[t, x], 20)
                    ssf(column.interphase_mass_transfer[t, x, "CO2"], 1 / 20)
                    ssf(column.interphase_mass_transfer[t, x, "H2O"], 1 / 100)
                    # ssf(column.mass_transfer_driving_force[t, x, "CO2"], 1 / 200 * 20)
                    # ssf(column.mass_transfer_driving_force[t, x, "H2O"], 1e-4 * 20)
                    # if self.config.surrogate_enhancement_factor_model is None:
                    #     ssf(column.conc_CO2_bulk[t, x], 3)
                        # if x != column.liquid_phase.length_domain.last():
                        #     cst(column.conc_CO2_equil_bulk_eqn[t, x], 10)
                    ssf(column.liquid_phase.heat[t, x], 1e-4)



            for x in column.vapor_phase.length_domain:
                ssf(column.heat_transfer_coeff[t, x], 1 / 3e6)
                ssf(column.vapor_phase.heat[t, x], 1e-4)

        # Have to do these scaling factors first in order for the column to scale right
        iscale.calculate_scaling_factors(self.liquid_properties)
        iscale.calculate_scaling_factors(self.vapor_properties)
        for column in column_list:
            iscale.calculate_scaling_factors(column)
        
        column = self.absorber
        for t in self.time:
            sf_mol = gsf(column.vapor_inlet.flow_mol[t])
            sf_CO2 = gsf(column.vapor_inlet.mole_frac_comp[t, "CO2"])
            ssf(self.absorber.co2_capture[t], 1/100)
            cst(self.absorber.co2_capture_eqn[t], sf_mol*sf_CO2/100)

        for t in self.time:
            ssf(self.lean_solvent_pump.control_volume.work[t], 1e-4)
            ssf(self.lean_solvent_pump.control_volume.deltaP[t], 1e-5)
            ssf(self.lean_rich_heat_exchanger.hot_side.heat[t], 1e-8)
            ssf(self.lean_rich_heat_exchanger.cold_side.heat[t], 1e-8)
            ssf(
                self.lean_rich_heat_exchanger.overall_heat_transfer_coefficient[t],
                1/self.lean_rich_heat_exchanger.overall_heat_transfer_coefficient[t].value
            )
            ssf(self.lean_rich_heat_exchanger.area, 1/4000)


            
    def print_column_design_parameters(self):
        print("\n ******* Printing some results *******")
        print("\nColumn diameter: ", pyo.value(self.absorber.diameter_column), "m")
        print("Column height: ", pyo.value(self.absorber.length_column), "m")
        print("Column volume: ", pyo.value(self.volume_column), "m3")
        print("Column volume with heads: ", pyo.value(self.volume_column_withheads), "m3")
        print("Column height to diameter ratio: ", pyo.value(self.HDratio))
        print("\nSolvent inlet molar flowrate: ", pyo.value(self.absorber.liquid_inlet.flow_mol[0]), "mol/s")
        print("L/G ratio: ", pyo.value(self.LGratio))
        inlet = self.absorber.liquid_inlet
        print("Inlet CO2 Loading: ", pyo.value(inlet.mole_frac_comp[0, "CO2"]/inlet.mole_frac_comp[0, "MEA"]))
        print("Inlet H2O Loading: ", pyo.value(inlet.mole_frac_comp[0, "H2O"]/inlet.mole_frac_comp[0, "MEA"]))
        outlet = self.absorber.liquid_outlet
        print("Outlet CO2 Loading: ", pyo.value(outlet.mole_frac_comp[0, "CO2"]/outlet.mole_frac_comp[0, "MEA"]))
        print("Outlet H2O Loading: ", pyo.value(outlet.mole_frac_comp[0, "H2O"]/outlet.mole_frac_comp[0, "MEA"]))

        print("\nCO2 capture: ", pyo.value(self.absorber.co2_capture[0]), "%")

    def strip_statevar_bounds_for_absorber_initialization(self):
        def strip_var_bounds(var):
            for j in var.keys():
                var[j].setlb(None)
                var[j].setub(None)
        def strip_property_blk_bounds(props, ions=True):
            for k in props.keys():
                if ions:
                    strip_var_bounds(props[k].flow_mol_phase_comp_true)
                    strip_var_bounds(props[k].mole_frac_phase_comp_true)
                strip_var_bounds(props[k].flow_mol)
                strip_var_bounds(props[k].temperature)
                strip_var_bounds(props[k].pressure)
                strip_var_bounds(props[k].mole_frac_phase_comp)
                strip_var_bounds(props[k].mole_frac_comp)
        # Strip variable bounds
        
        if not self.config.has_absorber_intercooling:
            prop_list = [
                self.absorber.liquid_phase.properties,
            ]
        else:
            prop_list = [
                getattr(self.absorber, f'absorber_bed_{i}.liquid_phase.properties') for i in range(1,4)
            ]
        for control_volume in [
            self.lean_rich_heat_exchanger.hot_side,
            self.lean_rich_heat_exchanger.cold_side,
            self.lean_solvent_pump.control_volume,
        ]:
            prop_list.append(control_volume.properties_in)
            prop_list.append(control_volume.properties_out)
        for props in prop_list:
            strip_property_blk_bounds(props)
        
        if not self.config.has_absorber_intercooling:
            prop_list_no_ions = [
                self.makeup_mixer.h2o_makeup_state,
                self.absorber.vapor_phase.properties,
            ]
        else:
            prop_list_no_ions = [
                self.makeup_mixer.h2o_makeup_state,
                self.absorber.absorber_bed_1.vapor_phase.properties,
                self.absorber.absorber_bed_2.vapor_phase.properties,
                self.absorber.absorber_bed_3.vapor_phase.properties,
            ]
        for props in prop_list_no_ions:
            strip_property_blk_bounds(props, ions=False)

        
    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg=None,
        liquid_temperature_guess=313.15,
        liquid_pressure_guess=101325,
    ):
        solver_obj = get_solver(solver, optarg)
        init_log = idaeslog.getInitLogger(self.name, outlvl)
        solve_log = idaeslog.getSolveLogger(self.name, outlvl)

        def safe_solve(blk):
            assert degrees_of_freedom(blk) == 0
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = solver_obj.solve(blk, tee=slc.tee)
            pyo.assert_optimal_termination(results)

        # First, initialize lean solvent pump
        self.lean_solvent_pump.initialize(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
        )

        # Next, initialize absorber column
        self.absorber_liquid_feed.expanded_block.deactivate()
        propagate_state(
            source=self.lean_solvent_pump.inlet, destination=self.absorber.liquid_inlet
        )
        # TODO: Bug in column model does not check if inlets are fixed
        # Fix liquid feed flowrate and composition
        self.absorber.liquid_inlet.flow_mol.fix()
        self.absorber.liquid_inlet.mole_frac_comp[...].fix()
        # Fix liquid feed temperature and pressure with initial guesses
        # TODO: Add initial guesses
        self.absorber.liquid_inlet.temperature.fix(liquid_temperature_guess)  # K
        self.absorber.liquid_inlet.pressure.fix(liquid_pressure_guess)  # Pa
        # Initialize
        self.absorber.initialize(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
        )

        # Rich solvent pump
        propagate_state(self.absorber_liquid_product)
        self.rich_solvent_pump.initialize(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
        )

        # HX unit
        propagate_state(self.hx_cold_feed)
        propagate_state(self.hx_hot_feed)
        self.lean_rich_heat_exchanger.initialize(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
            duty=(30.0, pyo.units.MW),
        )

        # Make-up mixer
        propagate_state(self.solvent_to_makeup)
        self.makeup_mixer.initialize(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
        )

        # Solve sub-flowsheet with absorber liquid feed uncoupled
        # init_log.info("Solving with absorber liquid feed fixed")
        # safe_solve(self)

        # Unfix liquid feed and activate coupling stream
        self.absorber_liquid_feed.expanded_block.activate()
        self.absorber.liquid_inlet.flow_mol.unfix()
        self.absorber.liquid_inlet.mole_frac_comp[...].unfix()
        self.absorber.liquid_inlet.temperature.unfix()
        self.absorber.liquid_inlet.pressure.unfix()

        # Finally, solve full sub-flowsheet
        init_log.info("Solving with coupled sub-flowsheet")
        safe_solve(self)
