#!/usr/bin/env python
# coding: utf-8





from pyomo.environ import (Constraint,
                           Var,
                           ConcreteModel,
                           Expression,
                           Objective,
                           TransformationFactory,
                           value)




# Todo: Import the above mentioned tools from pyomo.network



# Todo: Import the above mentioned tools from pyomo.network
from pyomo.network import Arc, SequentialDecomposition




from idaes.core import FlowsheetBlock



from idaes.generic_models.unit_models import (PressureChanger,
                                              Mixer,
                                              Separator as Splitter,
                                              Heater,
                                              CSTR,
                                              Flash,
                                              Translator)

from idaes.generic_models.unit_models.column_models import TrayColumn
from idaes.generic_models.unit_models.column_models.condenser     import CondenserType, TemperatureSpec 




# Utility tools to put together the flowsheet and calculate the degrees of freedom
from idaes.generic_models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from idaes.core.util.misc import get_solver

# Import idaes logger to set output levels
import idaes.logger as idaeslog




import hda_reaction_kinetic as reaction_props
from idaes.generic_models.properties.activity_coeff_models.    BTX_activity_coeff_VLE import BTXParameterBlock

from hda_ideal_VLE import HDAParameterBlock


def model():

    # Create a Pyomo Concrete Model to contain the problem
    m = ConcreteModel()

    # Add a steady state flowsheet block to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})




    # Property package for benzene, toluene, hydrogen, methane mixture
    m.fs.BTHM_params = HDAParameterBlock()

    # Property package for the benzene-toluene mixture
    m.fs.BT_params = BTXParameterBlock(default={
            "valid_phase": ('Liq', 'Vap'),
            "activity_coeff_model": "Ideal"
    })

    # Reaction package for the HDA reaction
    m.fs.reaction_params = reaction_props.HDAReactionParameterBlock(
            default={"property_package": m.fs.BTHM_params})




    # Adding the mixer M101 to the flowsheet
    m.fs.M101 = Mixer(default={"property_package": m.fs.BTHM_params,
                               "inlet_list": ["toluene_feed", "hydrogen_feed", "vapor_recycle"]})

    # Adding the heater H101 to the flowsheet
    m.fs.H101 = Heater(default={"property_package": m.fs.BTHM_params,
                                "has_phase_equilibrium": True})




    # Todo: Add reactor with the specifications above



    # Todo: Add reactor with the specifications above
    m.fs.R101 = CSTR(
                default={"property_package": m.fs.BTHM_params,
                         "reaction_package": m.fs.reaction_params,
                         "has_heat_of_reaction": True,
                         "has_heat_transfer": True})




    # Adding the flash tank F101 to the flowsheet
    m.fs.F101 = Flash(default={"property_package": m.fs.BTHM_params,
                               "has_heat_transfer": True,
                               "has_pressure_change": True})

    # Adding the splitter S101 to the flowsheet
    m.fs.S101 = Splitter(default={"property_package": m.fs.BTHM_params,
                                  "outlet_list": ["purge", "recycle"]})

    # Adding the compressor C101 to the flowsheet
    m.fs.C101 = PressureChanger(default={
                "property_package": m.fs.BTHM_params,
                "compressor": True,
                "thermodynamic_assumption": ThermodynamicAssumption.isothermal})





    # Add translator block to convert between property packages
    m.fs.translator = Translator(default={
            "inlet_property_package": m.fs.BTHM_params,
            "outlet_property_package": m.fs.BT_params
    })




    # Add constraint: Total flow = benzene flow + toluene flow (molar)
    m.fs.translator.eq_total_flow = Constraint(
        expr=m.fs.translator.outlet.flow_mol[0] ==
        m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] +
        m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"])

    # Add constraint: Outlet temperature = Inlet temperature
    m.fs.translator.eq_temperature = Constraint(
        expr=m.fs.translator.outlet.temperature[0] ==
        m.fs.translator.inlet.temperature[0])





    # Todo: Add constraint: Outlet pressure = Inlet pressure



    # Todo: Add constraint: Outlet pressure = Inlet pressure
    m.fs.translator.eq_pressure = Constraint(
        expr=m.fs.translator.outlet.pressure[0] ==
        m.fs.translator.inlet.pressure[0])



    # Remaining constraints on the translator block

    # Add constraint: Benzene mole fraction definition
    m.fs.translator.eq_mole_frac_benzene = Constraint(
        expr=m.fs.translator.outlet.mole_frac_comp[0, "benzene"] ==
        m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] /
        (m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] +
        m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"]))

    # Add constraint: Toluene mole fraction definition
    m.fs.translator.eq_mole_frac_toluene = Constraint(
        expr=m.fs.translator.outlet.mole_frac_comp[0, "toluene"] ==
        m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"] /
        (m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] +
        m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"]))




    # Todo: Add the Heater H102 to the flowsheet




    # Todo: Add the Heater H102 to the flowsheet
    m.fs.H102 = Heater(default={"property_package": m.fs.BT_params,
                                "has_pressure_change": True,
                                "has_phase_equilibrium": True})




    m.fs.s03 = Arc(source=m.fs.M101.outlet, destination=m.fs.H101.inlet)




    # Todo: Connect the H101 outlet to R101 inlet



    # Todo: Connect the H101 outlet to R101 inlet
    m.fs.s04 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)




    m.fs.s05 = Arc(source=m.fs.R101.outlet, destination=m.fs.F101.inlet)
    m.fs.s06 = Arc(source=m.fs.F101.vap_outlet, destination=m.fs.S101.inlet)
    m.fs.s08 = Arc(source=m.fs.S101.recycle, destination=m.fs.C101.inlet)
    m.fs.s09 = Arc(source=m.fs.C101.outlet,
                       destination=m.fs.M101.vapor_recycle)
    m.fs.s10a = Arc(source=m.fs.F101.liq_outlet,
                       destination=m.fs.translator.inlet)
    m.fs.s10b = Arc(source=m.fs.translator.outlet,
                       destination=m.fs.H102.inlet)




    TransformationFactory("network.expand_arcs").apply_to(m)




    # Define the conversion variables using 'Var'
    m.fs.R101.conversion = Var(initialize=0.75, bounds=(0, 1))

    # Append the constraint to the model
    m.fs.R101.conv_constraint = Constraint(
        expr=m.fs.R101.conversion*m.fs.R101.inlet.
        flow_mol_phase_comp[0, "Vap", "toluene"] ==
        (m.fs.R101.inlet.flow_mol_phase_comp[0, "Vap", "toluene"] -
        m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "toluene"]))




    print(degrees_of_freedom(m))



    # Check the degrees of freedom
    assert degrees_of_freedom(m) == 29




    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(0.30)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5)
    m.fs.M101.toluene_feed.temperature.fix(303.2)
    m.fs.M101.toluene_feed.pressure.fix(350000)




    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(0.30)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(0.02)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5)
    m.fs.M101.hydrogen_feed.temperature.fix(303.2)
    m.fs.M101.hydrogen_feed.pressure.fix(350000)




    # Fix the temperature of the outlet from the heater H101
    m.fs.H101.outlet.temperature.fix(600)




    # Todo: Fix the 'conversion' of the reactor R101


    # Todo: Fix the 'heat_duty' of the reactor R101



    # Todo: Fix the 'conversion' of the reactor R101
    m.fs.R101.conversion.fix(0.75)

    # Todo: Fix the 'heat_duty' of the reactor R101
    m.fs.R101.heat_duty.fix(0)




    # Fix the temperature of the vapor outlet from F101
    m.fs.F101.vap_outlet.temperature.fix(325.0)

    # Fix the pressure drop in the flash F101
    m.fs.F101.deltaP.fix(0)


    # Fix the split fraction of the 'purge' stream from S101
    m.fs.S101.split_fraction[0, "purge"].fix(0.2)

    # Fix the pressure of the outlet from the compressor C101
    m.fs.C101.outlet.pressure.fix(350000)




    # Fix the temperature of the outlet from the heater H102
    m.fs.H102.outlet.temperature.fix(375)

    # Fix the pressure drop in the heater H102
    m.fs.H102.deltaP.fix(-200000)




    # Todo: Check the degrees of freedom



    # Todo: Check the degrees of freedom
    print(degrees_of_freedom(m))



    # Check the degrees of freedom
    assert degrees_of_freedom(m) == 0




    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 3

    # Using the SD tool
    G = seq.create_graph(m)
    heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    order = seq.calculation_order(G)




    for o in heuristic_tear_set:
        print(o.name)




    for o in order:
        print(o[0].name)




    tear_guesses = {
            "flow_mol_phase_comp": {
                    (0, "Vap", "benzene"): 1e-5,
                    (0, "Vap", "toluene"): 1e-5,
                    (0, "Vap", "hydrogen"): 0.30,
                    (0, "Vap", "methane"): 0.02,
                    (0, "Liq", "benzene"): 1e-5,
                    (0, "Liq", "toluene"): 0.30,
                    (0, "Liq", "hydrogen"): 1e-5,
                    (0, "Liq", "methane"): 1e-5},
            "temperature": {0: 303},
            "pressure": {0: 350000}}

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.H101.inlet, tear_guesses)




    def function(unit):
            unit.initialize(outlvl=idaeslog.INFO)




    seq.run(m, function)




    # Create the solver object
    solver = get_solver()

    # Solve the model
    results = solver.solve(m, tee=True)



    # Check solver solve status
    from pyomo.environ import TerminationCondition
    assert results.solver.termination_condition == TerminationCondition.optimal




    # Add distillation column to the flowsheet
    m.fs.D101 = TrayColumn(default={
                            "number_of_trays": 10,
                            "feed_tray_location": 5,
                            "condenser_type":
                                CondenserType.totalCondenser,
                            "condenser_temperature_spec":
                                TemperatureSpec.atBubblePoint,
                            "property_package": m.fs.BT_params})

    # Connect the outlet from the heater H102 to the distillation column
    m.fs.s11 = Arc(source=m.fs.H102.outlet,
                   destination=m.fs.D101.feed)

    # Add the necessary equality constraints
    TransformationFactory("network.expand_arcs").apply_to(m)

    # Propagate the state
    propagate_state(m.fs.s11)

    # Fix the reflux ratio, boilup ratio, and the condenser pressure
    m.fs.D101.condenser.reflux_ratio.fix(0.5)
    m.fs.D101.reboiler.boilup_ratio.fix(0.5)
    m.fs.D101.condenser.condenser_pressure.fix(150000)

    # Initialize the distillation column
    m.fs.D101.initialize(outlvl=idaeslog.INFO)




    # Expression to compute the total cooling cost
    m.fs.cooling_cost = Expression(expr=0.25e-7 * (-m.fs.F101.heat_duty[0]) +
                                    0.2e-7 * (-m.fs.D101.condenser.heat_duty[0]))

    # Expression to compute the total heating cost
    m.fs.heating_cost = Expression(expr=2.2e-7 * m.fs.H101.heat_duty[0] +
                                    1.2e-7 * m.fs.H102.heat_duty[0] +
                                    1.9e-7 * m.fs.D101.reboiler.heat_duty[0])

    # Expression to compute the total operating cost
    m.fs.operating_cost = Expression(expr=(3600 * 24 * 365 *
                                    (m.fs.heating_cost + m.fs.cooling_cost)))

    # Expression to compute the total capital cost
    m.fs.capital_cost = Expression(expr=1e5*m.fs.R101.volume[0])




    # Check that the degrees of freedom is zero
    assert degrees_of_freedom(m) == 0



    solver.solve(m, tee=True)



    # Check solver solve status
    from pyomo.environ import TerminationCondition
    assert results.solver.termination_condition == TerminationCondition.optimal




    print('total cost = $', value(m.fs.capital_cost) + value(m.fs.operating_cost))
    print('operating cost = $', value(m.fs.operating_cost))
    print('capital cost = $', value(m.fs.capital_cost))
    print()
    print('Distillate flowrate = ', value(m.fs.D101.condenser.distillate.flow_mol[0]()), 'mol/s')
    print('Benzene purity = ', 100 * value(m.fs.D101.    condenser.distillate.mole_frac_comp[0, "benzene"]), '%')
    print('Residue flowrate = ', value(m.fs.D101.reboiler.bottoms.flow_mol[0]()), 'mol/s')
    print('Toluene purity = ', 100 * value(m.fs.D101.    reboiler.bottoms.mole_frac_comp[0, "toluene"]), '%')
    print()
    print('Conversion = ', 100 * value(m.fs.R101.conversion), '%')
    print()
    print('Overhead benzene loss in F101 = ',     100 * value(m.fs.F101.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"]) /         value(m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "benzene"]), '%' )



    import pytest
    assert value(m.fs.operating_cost) == pytest.approx(427593.007, abs=100)
    assert value(m.fs.capital_cost) == pytest.approx(14704.740, abs=100)




    m.fs.R101.report()




    m.fs.F101.report()




    from idaes.core.util.tables import create_stream_table_dataframe, stream_table_dataframe_to_string

    st = create_stream_table_dataframe({"Reactor": m.fs.s05, "Light Gases": m.fs.s06})
    print(stream_table_dataframe_to_string(st))






    m.fs.objective = Objective(expr=m.fs.operating_cost + m.fs.capital_cost)




    m.fs.H101.outlet.temperature.unfix()
    m.fs.R101.conversion.unfix()
    m.fs.F101.vap_outlet.temperature.unfix()
    m.fs.D101.condenser.condenser_pressure.unfix()
    m.fs.D101.condenser.reflux_ratio.unfix()
    m.fs.D101.reboiler.boilup_ratio.unfix()




    # Todo: Unfix the temperature of the outlet from H102



    # Todo: Unfix the temperature of the outlet from H102
    m.fs.H102.outlet.temperature.unfix()




    # Set bounds on the temperature of the outlet from H101
    m.fs.H101.outlet.temperature[0].setlb(500)
    m.fs.H101.outlet.temperature[0].setub(600)

    # Set bounds on the temperature of the outlet from R101
    m.fs.R101.outlet.temperature[0].setlb(600)
    m.fs.R101.outlet.temperature[0].setub(900)

    # Set bounds on the volume of the reactor R101
    m.fs.R101.volume[0].setlb(0)

    # Set bounds on the temperature of the vapor outlet from F101
    m.fs.F101.vap_outlet.temperature[0].setlb(298)
    m.fs.F101.vap_outlet.temperature[0].setub(450.0)

    # Set bounds on the temperature of the outlet from H102
    m.fs.H102.outlet.temperature[0].setlb(350)
    m.fs.H102.outlet.temperature[0].setub(400)

    # Set bounds on the pressure inside the condenser
    m.fs.D101.condenser.condenser_pressure.setlb(101325)
    m.fs.D101.condenser.condenser_pressure.setub(150000)




    # Todo: Set bounds on the reflux ratio



    # Todo: Set bounds on the boilup ratio



    # Todo: Set bounds on the reflux ratio
    m.fs.D101.condenser.reflux_ratio.setlb(0.1)
    m.fs.D101.condenser.reflux_ratio.setub(5)

    # Todo: Set bounds on the boilup ratio
    m.fs.D101.reboiler.boilup_ratio.setlb(0.1)
    m.fs.D101.reboiler.boilup_ratio.setub(5)




    # Ensure that the overhead loss of benzene from F101 <= 20%
    m.fs.overhead_loss = Constraint(
        expr=m.fs.F101.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"] <=
        0.20 * m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "benzene"])




    # Todo: Add minimum product flow constraint



    # Todo: Add minimum product flow constraint
    m.fs.product_flow = Constraint(
        expr=m.fs.D101.condenser.distillate.flow_mol[0] >= 0.18)




    m.fs.product_purity = Constraint(
        expr=m.fs.D101.condenser.
        distillate.mole_frac_comp[0, "benzene"] >= 0.99)



    return m


if __name__ == "__main__":
    import sys
    from util import fsvis_main
    sys.exit(fsvis_main(model(), name="HDA"))

