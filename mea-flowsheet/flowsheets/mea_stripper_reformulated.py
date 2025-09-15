"""
Stripper sub-flowsheet for Monoethanolamine solvent carbon capture system
"""
# Python imports
import copy
import math
import logging

# Pyomo imports
import pyomo.environ as pyo
from pyomo.network import Arc, Port
from pyomo.common.config import ConfigValue, ConfigDict, In, Bool, ListOf
from pyomo.common.collections import ComponentSet

# IDAES imports
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlockData, declare_process_block_class

from idaes.models.unit_models import Mixer, MomentumMixingType

from idaes.models_extra.column_models.solvent_reboiler import SolventReboiler
from idaes.models_extra.column_models.solvent_condenser import SolventCondenser
from idaes.models_extra.column_models.MEAsolvent_column import MEAColumn

from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale

# External IDAES imports (workspace)

from mea_properties import (
    MEALiquidParameterBlock,
    MEAVaporParameterBlock,
    scale_mea_liquid_params,
    scale_mea_vapor_params,
    switch_liquid_to_parmest_params,
    state_bounds_default
)

logging.getLogger('pyomo.repn.plugins.nl_writer').setLevel(logging.ERROR)

@declare_process_block_class("MEAStripperFlowsheet")
class MEAStripperFlowsheetData(FlowsheetBlockData):
    CONFIG = FlowsheetBlockData.CONFIG()
    CONFIG.declare(
        "column_finite_element_set",
        ConfigValue(
            domain=ListOf(float),
            description="List containing coordinates finite element faces "
            "of stripper column. Coordinates must start with zero, be "
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
    def build(self):
        super().build()
        self._add_properties()
        self._add_units()
        self._add_arcs()
        self._add_performance_math()
        self._set_stripper_design_inputs()
        self._scaling()

    def _add_properties(self):
        if len(self.config.property_package_bounds) == 0:
            bounds = None
        else:
            bounds = self.config.property_package_bounds
        self.vapor_properties = MEAVaporParameterBlock(state_bounds=bounds)
        self.liquid_properties = MEALiquidParameterBlock(ions=True, state_bounds=bounds)
        self.liquid_properties_no_ions = MEALiquidParameterBlock(ions=False, state_bounds=bounds)

        self.liquid_properties.CO2.dh_abs_co2.fix(-97000)
        self.liquid_properties_no_ions.CO2.dh_abs_co2.fix(-97000)

    def _add_units(self):
        self.reflux_mixer = Mixer(
            property_package=self.liquid_properties_no_ions,
            inlet_list=["rich_solvent", "reflux"],
            momentum_mixing_type=MomentumMixingType.none
        )

        @self.reflux_mixer.Constraint(self.time)
        def pressure_equality_eqn(b, t):
            return b.mixed_state[t].pressure == b.rich_solvent_state[t].pressure
        
        fe_set = self.config.column_finite_element_set
        self.stripper = MEAColumn(
            finite_elements=len(fe_set) - 1,
            length_domain_set=fe_set,
            vapor_phase={"property_package": self.vapor_properties},
            liquid_phase={"property_package": self.liquid_properties},
            # corrected_hx_coeff_eqn=False,
            # surrogate_enhancement_factor_model=self.config.surrogate_enhancement_factor_model,
        )
        
        self.condenser = SolventCondenser(
            liquid_property_package=self.liquid_properties_no_ions,
            vapor_property_package=self.vapor_properties,
        )
        
        self.reboiler = SolventReboiler(
            liquid_property_package=self.liquid_properties,
            vapor_property_package=self.vapor_properties,
            has_pressure_change=True,
        )
        
    def _add_arcs(self):
        # Define Arcs (streams)
        self.column_liquid_feed = Arc(
            source=self.reflux_mixer.outlet,
            destination=self.stripper.liquid_inlet,
        )
        self.column_vapor_outlet = Arc(
            source=self.stripper.vapor_outlet,
            destination=self.condenser.inlet,
        )
        self.column_reflux = Arc(
            source=self.condenser.reflux,
            destination=self.reflux_mixer.reflux,
        )
        self.column_liquid_outlet = Arc(
            source=self.stripper.liquid_outlet,
            destination=self.reboiler.inlet,
        )
        self.column_boilup = Arc(
            source=self.reboiler.vapor_reboil,
            destination=self.stripper.vapor_inlet,
        )

        # Transform Arcs
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)

        # Add sub-flowsheet level ports for easy connectivity
        self.rich_solvent = Port(
            extends=self.reflux_mixer.rich_solvent
        )
        self.lean_solvent = Port(
            extends=self.reboiler.bottoms
        )
        self.distillate = Port(
            extends=self.condenser.vapor_outlet
        )
        
    def _add_performance_math(self):
        @self.Expression(doc="Column height to column diameter ratio (-)")
        def HDratio(b):
            return b.stripper.length_column / b.stripper.diameter_column

        @self.Expression(self.time, doc="Inlet liquid (solvent) flowrate to inlet gas flowrate ratio (-)")
        def LGratio(b, t):
            return b.stripper.liquid_inlet.flow_mol[t] / b.stripper.vapor_inlet.flow_mol[t]

        @self.Expression(doc="Volume of column (cubic m)")
        def volume_column(b):
            return math.pi * b.stripper.length_column * (b.stripper.diameter_column / 2) ** 2

        @self.Expression(doc="Volume of column with heads (cubic m)")
        def volume_column_withheads(b):
            return (
                math.pi/4 * b.stripper.diameter_column ** 2 * b.stripper.length_column
                + math.pi / 3 * b.stripper.diameter_column ** 3
            )
        @self.Expression(self.time, doc="Captured CO2 mass flow rate")
        def distillate_flow_mass_CO2(b, t):
            return(
                b.condenser.vapor_phase.properties_out[t].flow_mass_comp["CO2"]
            )
        
        @self.Expression(self.time, doc="Lean solvent mass flow rate")
        def lean_solvent_flow_mass(b, t):
            return(
                b.reboiler.liquid_phase.properties_out[t].flow_mass
            )
        @self.Expression(self.time, doc="Lean solvent volumetric flow rate")
        def lean_solvent_flow_vol(b, t):
            return(
                b.reboiler.liquid_phase.properties_out[t].flow_vol
            )

    def _set_stripper_design_inputs(self):
        # Stripper geometry
        self.stripper.diameter_column.fix(4)
        self.stripper.length_column.fix(12)

        # Condenser
        self.condenser.heat_duty.fix(-4e7)

        # Reboiler
        self.reboiler.heat_duty.fix(1e8)
        self.reboiler.inlet.pressure.fix(183700)
        
    def _scaling(self):
        gsf = iscale.get_scaling_factor
        ssf = iscale.set_scaling_factor

        def cst(con, s):
            iscale.constraint_scaling_transform(con, s, overwrite=False)

        scale_mea_vapor_params(self.vapor_properties, scaling_factor_flow_mol=3e-4)
        scale_mea_liquid_params(self.liquid_properties, ions=True, scaling_factor_flow_mol=3e-4)
        scale_mea_liquid_params(self.liquid_properties_no_ions, ions=False, scaling_factor_flow_mol=3e-4)

        # Stripper column
        column = self.stripper
        for t in self.time:
            for x in column.liquid_phase.length_domain:
                ssf(column.velocity_liq[t, x], 20)
                ssf(column.interphase_mass_transfer[t, x, "CO2"], 1 / 20)
                ssf(column.interphase_mass_transfer[t, x, "H2O"], 1 / 100)
                ssf(column.mass_transfer_driving_force[t, x, "CO2"], 1 / 25000 * 20)
                ssf(column.mass_transfer_driving_force[t, x, "H2O"], 1e-5 * 20)
                if self.config.surrogate_enhancement_factor_model is None:
                    ssf(column.conc_CO2_bulk[t, x], 1)
                ssf(column.liquid_phase.heat[t, x], 1e-4)

                # if x != column.liquid_phase.length_domain.last():
                    # ssf(column.omega[t, x], 1e-1)
                    # cst(column.enhancement_factor_eqn1[t, x], 1e-1)
                    # cst(column.enhancement_factor_eqn2[t, x], 1e-1)

            for x in column.vapor_phase.length_domain:
                ssf(column.heat_transfer_coeff[t, x], 1 / 500000)
                ssf(column.vapor_phase.heat[t, x], 1e-4)

        for t in self.time:
            ssf(self.reboiler.liquid_phase.mass_transfer_term[t, "Liq", "CO2"], 0.001)
            ssf(self.reboiler.liquid_phase.mass_transfer_term[t, "Liq", "H2O"], 0.001)
            ssf(self.reboiler.liquid_phase.properties_out[t].fug_phase_comp["Liq", "CO2"], 0.0001)
            ssf(self.reboiler.liquid_phase.properties_out[t].fug_phase_comp["Liq", "H2O"], 0.0001)
            ssf(self.reboiler.liquid_phase.properties_out[t].fug_phase_comp["Liq", "MEA"], 0.001)
            ssf(self.reboiler.liquid_phase.enthalpy_transfer[t], 1e-7)
            ssf(self.reboiler.liquid_phase.properties_out[t].pressure, 1e-5)

            cst(self.reflux_mixer.pressure_equality_eqn[t], 1e-5)


            
    def print_column_design_parameters(self):
        print("\n ******* Printing some results *******")
        print("\nColumn diameter: ", pyo.value(self.stripper.diameter_column), "m")
        print("Column height: ", pyo.value(self.stripper.length_column), "m")
        print("Column volume: ", pyo.value(self.volume_column), "m3")
        print("Column volume with heads: ", pyo.value(self.volume_column_withheads), "m3")
        print("Column height to diameter ratio: ", pyo.value(self.HDratio))
        print("\nSolvent inlet molar flowrate: ", pyo.value(self.stripper.liquid_inlet.flow_mol[0]), "mol/s")
        print("L/G ratio: ", pyo.value(self.LGratio[0]))

        inlet = self.stripper.liquid_inlet
        print("Inlet CO2 Loading: ", pyo.value(inlet.mole_frac_comp[0, "CO2"]/inlet.mole_frac_comp[0, "MEA"]))
        print("Inlet H2O Loading: ", pyo.value(inlet.mole_frac_comp[0, "H2O"]/inlet.mole_frac_comp[0, "MEA"]))
        outlet = self.stripper.liquid_outlet
        print("Outlet CO2 Loading: ", pyo.value(outlet.mole_frac_comp[0, "CO2"]/outlet.mole_frac_comp[0, "MEA"]))
        print("Outlet H2O Loading: ", pyo.value(outlet.mole_frac_comp[0, "H2O"]/outlet.mole_frac_comp[0, "MEA"]))

    def strip_statevar_bounds_for_stripper_initialization(self):
        # Strip variable bounds 
        for t in self.time:
            for props in [
                self.condenser.vapor_phase.properties_in[t],
                self.condenser.vapor_phase.properties_out[t],
                self.condenser.liquid_phase[t],
                self.reboiler.liquid_phase.properties_in[t],
                self.reboiler.liquid_phase.properties_out[t],
                self.reboiler.vapor_phase[t]
            ]:
                for idx, var in props.mole_frac_comp.items():
                    var.domain = pyo.Reals
                    var.bounds = (None, None)
                for idx, var in props.mole_frac_phase_comp.items():
                    var.domain = pyo.Reals
                    var.bounds = (None, None)
                # for idx, var in props.flow_mol_phase_comp.items():
                #     var.domain = pyo.Reals
                #     var.bounds = (None, None)

            for props in [
                self.reboiler.liquid_phase.properties_in[t],
                self.reboiler.liquid_phase.properties_out[t],
            ]:
                for idx, var in props.mole_frac_phase_comp_true.items():
                    var.domain = pyo.Reals
                    var.bounds = (None, None)
                for idx, var in props.flow_mol_phase_comp_true.items():
                    var.domain = pyo.Reals
                    var.bounds = (None, None)
                
            for x in self.stripper.liquid_phase.length_domain:
                for j in self.stripper.config.liquid_phase.property_package.true_phase_component_set:
                    self.stripper.liquid_phase.properties[t, x].flow_mol_phase_comp_true[j].domain = pyo.Reals
                    self.stripper.liquid_phase.properties[t, x].flow_mol_phase_comp_true[j].setlb(None)

                    self.stripper.liquid_phase.properties[t, x].mole_frac_phase_comp_true[j].domain = pyo.Reals
                    self.stripper.liquid_phase.properties[t, x].mole_frac_phase_comp_true[j].setlb(None)

                for j in self.stripper.config.liquid_phase.property_package.apparent_species_set:
                    self.stripper.liquid_phase.properties[t, x].mole_frac_comp[j].domain = pyo.Reals
                    self.stripper.liquid_phase.properties[t, x].mole_frac_comp[j].setlb(None)

                for j in self.stripper.config.liquid_phase.property_package.apparent_phase_component_set:
                    self.stripper.liquid_phase.properties[t, x].mole_frac_phase_comp.domain = pyo.Reals
                    self.stripper.liquid_phase.properties[t, x].mole_frac_phase_comp.setlb(None)
                if self.config.surrogate_enhancement_factor_model is None:
                    self.stripper.conc_interface_MEA[t, x].setub(None)
                    self.stripper.log_conc_interface_MEA[t, x].setub(None)

            for x in self.stripper.vapor_phase.length_domain:
                for j in self.stripper.config.vapor_phase.property_package.component_list:
                    self.stripper.vapor_phase.properties[t, x].mole_frac_comp[j].domain = pyo.Reals
                    self.stripper.vapor_phase.properties[t, x].mole_frac_comp[j].setlb(None)

                for j in self.stripper.config.vapor_phase.property_package._phase_component_set:
                    self.stripper.vapor_phase.properties[t, x].mole_frac_phase_comp[j].domain = pyo.Reals
                    self.stripper.vapor_phase.properties[t, x].mole_frac_phase_comp[j].setlb(None)
    
        
    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg=None,
        boilup_guess=None
    ):
        solver_obj = get_solver(solver, optarg)
        init_log = idaeslog.getInitLogger(self.name, outlvl)
        solve_log = idaeslog.getSolveLogger(self.name, outlvl)

        def safe_solve(blk):
            assert degrees_of_freedom(blk) == 0
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = solver_obj.solve(blk, tee=slc.tee, symbolic_solver_labels=True)
            pyo.assert_optimal_termination(results)


        # Start with initializing stripper column
        # Get guess for liquid feed conditions from mixer rich solvent
        propagate_state(
            source=self.reflux_mixer.rich_solvent, destination=self.stripper.liquid_inlet
        )

        # Set initial guesses for reboil stream
        if boilup_guess is None:
            boilup_guess = {
                "flow_mol": 2500,
                "temperature": 392.75,
                "mole_frac_comp": {
                    "CO2": 0.06,
                    "H2O": 0.94,
                },
            }

        self.column_boilup.expanded_block.deactivate()
        for t in self.time:
            self.stripper.vapor_inlet.flow_mol[t].fix(boilup_guess["flow_mol"])  # mol/sec
            self.stripper.vapor_inlet.temperature[t].fix(boilup_guess["temperature"])  # K
            self.stripper.vapor_inlet.pressure[t].fix(self.reboiler.inlet.pressure[t])  # Fix to reboiler pressure
            self.stripper.vapor_inlet.mole_frac_comp[t, "CO2"].fix(
                boilup_guess["mole_frac_comp"]["CO2"]
            )
            self.stripper.vapor_inlet.mole_frac_comp[t, "H2O"].fix(
                boilup_guess["mole_frac_comp"]["H2O"]
            )

        # Initialize stripper column
        self.stripper.initialize(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
            mode="stripper",
        )

        # Initialize condenser
        propagate_state(self.column_vapor_outlet)
        # TODO: Condenser requires 0 DOF in order to initialize at the moment, should be fixed
        self.condenser.inlet.flow_mol.fix()
        self.condenser.inlet.temperature.fix()
        self.condenser.inlet.pressure.fix()
        self.condenser.inlet.mole_frac_comp[...].fix()
        self.condenser.initialize(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
        )

        # TODO: Manually unfix streams, condenser initialization should be updated
        self.condenser.inlet.flow_mol.unfix()
        self.condenser.inlet.temperature.unfix()
        self.condenser.inlet.pressure.unfix()
        self.condenser.inlet.mole_frac_comp[...].unfix()

        # Next, initialize reflux mixer
        propagate_state(self.column_reflux)
        self.reflux_mixer.initialize(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
        )
        propagate_state(self.column_liquid_feed)

        # Solve the reflux loop for the stripper assembly
        init_log.info("Solving column with condenser")
        self.reboiler.deactivate()
        self.column_liquid_outlet.expanded_block.deactivate()
        safe_solve(self)

        # Next, initialize the reboiler
        self.reboiler.activate()
        self.column_liquid_outlet.expanded_block.activate()
        propagate_state(self.column_liquid_outlet)
        # TODO: Bug in reboiler initialization too
        self.reboiler.inlet.flow_mol.fix()
        self.reboiler.inlet.temperature.fix()
        self.reboiler.inlet.pressure.fix()
        self.reboiler.inlet.mole_frac_comp[...].fix()
        # For initialization, need to fix deltaP to zero
        self.reboiler.deltaP.fix(0)
        self.reboiler.initialize(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
        )
        # TODO: Bug in reboiler initialization too
        self.reboiler.inlet.flow_mol.unfix()
        self.reboiler.inlet.temperature.unfix()
        # Don't unfix pressure - that is an operating condition
        self.reboiler.inlet.mole_frac_comp[...].unfix()
        # Do need to unfix column inlet vapor pressure at this point however
        self.stripper.vapor_inlet.pressure.unfix()

        # Finally, try to close the boil-up recycle solve the full flowsheet
        init_log.info("Solving full flowsheet")

        self.column_boilup.expanded_block.activate()
        self.stripper.vapor_inlet.flow_mol.unfix()
        self.stripper.vapor_inlet.temperature.unfix()

        self.stripper.vapor_inlet.mole_frac_comp[...].unfix()
        self.reboiler.deltaP.unfix()

        # xfrm = pyo.TransformationFactory("contrib.strip_var_bounds")
        # xfrm.apply_to(self, reversible=True)
        safe_solve(self)
