##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
# =============================================================================
"""
Created on Tue Jun  9 18:17:01 2020
Goal: Demonstrate Data reconciliation framework for a SCPC Boiler subflowsheet

Inputs:
    BFW - boiler feed water (from Feed water heaters)
    Coal from pulverizers

Main Assumptions:
    Coal flowrate as a function of load, coal HHV is fixed and heat dutty
    splitt from fire side to water wall and platen superheater is fixed.

    Boiler heat exchanger network:
        Water Flow:
            BFW -> ECONOMIZER -> Water Wall -> Primary SH -> Platen SH -> Finishing Superheate -> HP Turbine -> Reheater -> IP Turbine
        Flue Gas Flow:
            Fire Ball -> Platen SH -> Finishing SH -> Reheater  -> o -> Economizer -> Air Preheater
                                                   -> Primary SH --^

        *HP Turbine, IP Turbine, Air Preheater ==> not included in this release

    Models used:
        - Mixers: Attemperator, Flue gas mix
        - Heater: Platen SH, Fire/Water side (simplified model)
        - BoilerHeatExchanger: Economizer, Primary SH, Finishing SH, Reheater
            + Shell and tube heat exchanger
                - tube side: Steam (side 1 holdup)
                - shell side: flue gas (side 2 holdup)

    Property packages used:
        - IAPWS: Water/steam side
        - IDEAL GAS: Flue Gas side
"""
__author__ = "Miguel Zamarripa"
# Import Pyomo libraries
from pyomo.environ import value
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES core
from idaes.generic_models.properties import iapws95
from idaes.core.util import model_serializer as ms
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale

def deactivate_performance_constraints(m):
    for blk in m.fs.component_objects(pyo.Block, descend_into=False):
        try:
            Uval = value(blk.overall_heat_transfer_coefficient[0])
            Dps = value(blk.deltaP_shell[0])
            Dpt = value(blk.deltaP_tube[0])
            blk.overall_heat_transfer_coefficient[0].setlb(Uval*0.25)
            blk.overall_heat_transfer_coefficient[0].setub(Uval*2.5)
            blk.deltaP_shell.setlb(Dps*10)
            blk.deltaP_shell.setub(0)
            blk.deltaP_tube.setlb(Dpt*10)
            blk.deltaP_tube.setub(0)

            blk.overall_heat_transfer_coefficient_eqn.deactivate()
            blk.rcond_wall_eqn.deactivate()
            blk.hconv_shell_total_eqn.deactivate()
            blk.hconv_shell_conv_eqn.deactivate()
            blk.N_Nu_shell_eqn.deactivate()
            blk.N_Pr_shell_eqn.deactivate()
            blk.deltaP_shell_eqn.deactivate()
            blk.friction_factor_shell_eqn.deactivate()
            blk.N_Re_shell_eqn.deactivate()
            blk.v_shell_eqn.deactivate()
            blk.hconv_tube_eqn.deactivate()
            blk.N_Nu_tube_eqn.deactivate()
            blk.N_Pr_tube_eqn.deactivate()
            blk.deltaP_tube_eqn.deactivate()
            blk.deltaP_tube_uturn_eqn.deactivate()
            blk.deltaP_tube_friction_eqn.deactivate()
            blk.friction_factor_tube_eqn.deactivate()
            blk.N_Re_tube_eqn.deactivate()
            blk.v_tube_eqn.deactivate()

            if blk.config.has_radiation is True:
                blk.hconv_shell_rad_eqn.deactivate()
                blk.frad_gas_shell_eqn.deactivate()
                blk.gas_gray_fraction_eqn.deactivate()
                blk.gas_emissivity_mul2_eqn.deactivate()
                blk.gas_emissivity_div2_eqn.deactivate()
                blk.gas_emissivity_eqn.deactivate()
        except AttributeError:
            pass
    print('degrees of freedom = ' + str(degrees_of_freedom(m)))


def boiler_flowsheet():
    """
    Make the flowsheet object, fix some variables, and solve the problem
    """
    import idaes.power_generation.flowsheets.\
        supercritical_power_plant.boiler_subflowsheet_build as blr
    # First, we import the boler model from boiler_subflowsheet_build
    m, solver = blr.main()
    # initialize boiler flowsheet
    blr.initialize(m)
    # unfix inlet values to each unit and build the flowsheet connectivity
    blr.unfix_inlets(m)
    results = solver.solve(m, tee=True)
    print(results)

    # adding boiler block
    # The model solved in line 116 considers heat to platen superheater,
    # heat to water wall, and flue gas inlet to finishing superheater as fixed
    # variables. This code uses simplified surrogate models for these inputs.
    # Therefore, in the new model, these variables will change with plant load,
    # More details regarding this methodology can be found in Zamarripa et al.,
    # FOCAPD 2 page paper included in this directory as FOCAPD_paper.pdf
    m.fs.SR = pyo.Var(initialize=1.15, doc='stoichiometric ratio')
    m.fs.SR.setlb(1.0)
    m.fs.SR.setub(1.25)
    m.fs.SR.fix(1.15)
    m.fs.coal_flow = pyo.Var(initialize=54.01, doc='coal flowrate kg.s')
    m.fs.FGET = pyo.Var(initialize=1300.00, doc='flue gas exit temperature K')
    m.fs.FG_flow = pyo.Var(m.fs.prop_fluegas.component_list,
                           initialize=2000,
                           domain=pyo.Reals, doc='component molar flow rate')
    init_fg = {"H2O": (8.69/100), "CO2": (14.49/100),
               "N2": (0.6999), "O2": (2.47/100), "NO": 0.0006, "SO2": (0.002)}
    m.fs.fg_frac = pyo.Var(m.fs.prop_fluegas.component_list,
                           initialize=init_fg,
                           doc='flue gas molar fraction')

    # surrogate models valid for coal flowrate (30 to 70 kg/s) and SR (1-2.5)
    def fg_mol_frac_rule(b, i):
        if i == 'CO2':
            return m.fs.fg_frac[i] == 0.23235491508307534735955\
                * m.fs.SR**-0.5
        elif i == 'H2O':
            return m.fs.fg_frac[i] == 0.45009703293861530459807E-001 \
                * m.fs.SR
        elif i == 'N2':
            return m.fs.fg_frac[i] == 1 - (m.fs.fg_frac['CO2']
                                           + m.fs.fg_frac['H2O']
                                           + m.fs.fg_frac['NO']
                                           + m.fs.fg_frac['O2']
                                           + m.fs.fg_frac['SO2'])
        elif i == 'NO':
            return m.fs.fg_frac[i] == 0.38148020722713599076938E-005 \
                * m.fs.coal_flow
        elif i == 'O2':
            return m.fs.fg_frac[i] == 0.60149266916011525155317E-003 \
                * m.fs.coal_flow
        elif i == 'SO2':
            return m.fs.fg_frac[i] == 0.21991717807021016347323E-004 \
                * m.fs.coal_flow
    m.fs.fg_mol_cn = pyo.Constraint(m.fs.prop_fluegas.component_list,
                                    rule=fg_mol_frac_rule)

    def fg_flow_rule(b, i):
        mass_flow = 0.87757407893140026988732 * m.fs.coal_flow \
            - 0.68240933480416252066014E-001 * m.fs.SR + 8.6386912890637752582\
            * m.fs.coal_flow*m.fs.SR - 0.11737247790640543564002E-003 \
            * (m.fs.coal_flow/m.fs.SR)**2  # ~ 28.3876e3 kg/s

        # flow mol component in mol/s = kg/s/kg/mol
        return m.fs.FSH.side_2.properties_in[0].flow_mol_comp[i] == \
            (mass_flow / sum(m.fs.fg_frac[c]*m.fs.prop_fluegas.mw_comp[c]
                             for c in m.fs.prop_fluegas.component_list))\
            * m.fs.fg_frac[i]
    m.fs.fg_flow_cn = pyo.Constraint(m.fs.prop_fluegas.component_list,
                                     rule=fg_flow_rule)

    def fg_temp_rule(b):
        return m.fs.FSH.side_2.properties_in[0].temperature == \
            59836.381548557488713413 * m.fs.coal_flow**-0.5 \
            + 791.74814907302368283126 * m.fs.coal_flow**0.5 \
            + 0.60200443235342349090899 * m.fs.coal_flow**1.5 \
            + 285.74858049626226375040 * m.fs.SR**3 - 29865.27413896147845662\
            * (m.fs.coal_flow*m.fs.SR)**-.5 - 518.58090213915738786454 \
            * (m.fs.coal_flow*m.fs.SR)**0.5 - 0.86351166781748790735040E-002 \
            * (m.fs.coal_flow*m.fs.SR)**2 - 30141.308801694875000976 \
            * (m.fs.SR/m.fs.coal_flow)**0.5 - 18.109911868067683826666 \
            * m.fs.coal_flow/m.fs.SR + 1017.0807559525446777116 \
            * m.fs.SR/m.fs.coal_flow
    m.fs.fg_temp_cn = pyo.Constraint(rule=fg_temp_rule)

    def heat_2_PLSH_rule(b):
        return m.fs.PlSH.heat_duty[0] == - 42149808.046329699456692 \
            * pyo.log(m.fs.coal_flow) + 13125913.817196270450950 \
            * pyo.exp(m.fs.SR) - 82168941.403612509369850 \
            * m.fs.coal_flow**-.5 - 20751.165176131220505340 \
            * (m.fs.coal_flow*m.fs.SR)**1.5 + 35303843.583323523402214 \
            * (m.fs.coal_flow/m.fs.SR)**0.5
    m.fs.heat_2_PLSH_cn = pyo.Constraint(rule=heat_2_PLSH_rule)

    def heat_2_ww_rule(b):
        heat_total = 15383517.068522246554494 * m.fs.coal_flow \
            - 145195958.70188459753990 * m.fs.SR**0.5 + 91548063.4268338829278\
            * (m.fs.coal_flow*m.fs.SR)**0.5 - 11732787.822234204038978 \
            * m.fs.coal_flow*m.fs.SR - 19639.552666366322227987 \
            * (m.fs.coal_flow/m.fs.SR)**2
        return m.fs.Water_wall.heat_duty[0] == heat_total \
            - m.fs.PlSH.heat_duty[0]
    m.fs.heat_2_ww_cn = pyo.Constraint(rule=heat_2_ww_rule)

    # unfix variables to build flowsheet connections
    m.fs.PlSH.heat_duty.unfix()
    m.fs.Water_wall.heat_duty.unfix()
    m.fs.FSH.side_2.properties_in[:].flow_mol_comp.unfix()
    m.fs.FSH.side_2.properties_in[:].temperature.unfix()
    m.fs.coal_flow.fix(50.15)
    print('degrees of freedom = ' + str(degrees_of_freedom(m)))
    # initialize surrogate models (better initial point before solve)
    calculate_variable_from_constraint(
            m.fs.PlSH.heat_duty[0],
            m.fs.heat_2_PLSH_cn)
    calculate_variable_from_constraint(
        m.fs.Water_wall.heat_duty[0],
        m.fs.heat_2_ww_cn)
    # set scaling parameters
    iscale.calculate_scaling_factors(m)
    # final solve
    results = solver.solve(m, tee=False)
    return m


if __name__ == "__main__":
    m = boiler_flowsheet()
    deactivate_performance_constraints(m)
