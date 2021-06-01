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
"""
This file builds the kriging SOFC reduced order model (developed by Pacific
Northwest National Laboratory) in a Pyomo block. A data file of kriging
coefficients must be provided. The ROM is specifically designed for use with
the NGFC flowsheet. The code for generating SOFC ROMS can be found on PNNL's
Github page: https://github.com/NGFC-Lib/NGFC-Lib.
"""

import os
from pyomo.environ import Block, Constraint, Param, Var, exp, value, units
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.fileutils import this_file_dir


# creates a dictionary from a list of indices and values
def build_dict(index, values):
    d = {}
    for i, v in zip(index, values):
        d[i] = v
    return d


# creates a dictionary with double indices
def build_matrix(index1, index2, values):
    d = {}
    count = 0
    for i in index1:
        for j in index2:
            d[(i, j)] = values[count]
            count += 1
    return d


def build_SOFC_ROM(m):
    m.SOFC = b = Block()

    # load kriging coefficients
    fname = 'kriging_coefficients.dat'
    path = os.path.join(this_file_dir(), fname)

    with open(path) as file:
        text = file.readlines()

    kriging = [float(line) for line in text]

    n_inputs = 9
    n_outputs = 48
    n_samples = 13424

    # create indecies for vars and params
    input_index = list(range(n_inputs))
    input_plus_index = list(range(n_inputs+1))
    output_index = list(range(n_outputs))
    samples_index = list(range(n_samples))

    # read through data file and create params
    start_count = 0
    end_count = n_inputs
    values = kriging[start_count:end_count]
    b.mean_input = Param(input_index,
                         initialize=build_dict(input_index, values),
                         mutable=False)

    start_count += n_inputs
    end_count += n_inputs
    values = kriging[start_count:end_count]
    b.sigma_input = Param(input_index,
                          initialize=build_dict(input_index, values),
                          mutable=False)

    start_count += n_inputs
    end_count += n_outputs
    values = kriging[start_count:end_count]
    b.mean_output = Param(output_index,
                          initialize=build_dict(output_index, values),
                          mutable=False)

    start_count += n_outputs
    end_count += n_outputs
    values = kriging[start_count:end_count]
    b.sigma_output = Param(output_index,
                           initialize=build_dict(output_index, values),
                           mutable=False)

    start_count += n_outputs
    end_count += n_inputs*n_samples
    values = kriging[start_count:end_count]
    b.ds_input = Param(samples_index, input_index,
                       initialize=build_matrix(
                           samples_index, input_index, values),
                       mutable=False)

    start_count += n_inputs*n_samples
    end_count += n_inputs
    values = kriging[start_count:end_count]
    b.theta = Param(input_index,
                    initialize=build_dict(input_index, values),
                    mutable=False)

    start_count += n_inputs
    end_count += (n_inputs+1)*n_outputs
    values = kriging[start_count:end_count]
    b.beta = Param(input_plus_index, output_index,
                   initialize=build_matrix(
                       input_plus_index, output_index, values),
                   mutable=False)

    start_count += (n_inputs+1)*n_outputs
    end_count += n_samples*n_outputs
    values = kriging[start_count:end_count]
    b.gamma = Param(samples_index, output_index,
                    initialize=build_matrix(
                        samples_index, output_index, values),
                    mutable=False)

    # create input vars for the user to interface with
    b.current_density = Var(initialize=4000,
                            units=units.A/units.m**2,
                            bounds=(2000, 6000))
    # units for T should be degC but pyomo doesn't support conversion from K
    b.fuel_temperature = Var(initialize=500,
                             units=None,
                             bounds=(15, 600))

    b.internal_reforming = Var(initialize=0.4,
                               units=None,
                               bounds=(0, 1))

    b.air_temperature = Var(initialize=700,
                            units=None,
                            bounds=(550, 800))

    b.air_recirculation = Var(initialize=0.5,
                              units=None,
                              bounds=(0, 0.8))

    b.OTC = Var(initialize=2.1,
                units=None,
                bounds=(1.5, 3))

    b.fuel_util = Var(initialize=0.85,
                      units=None,
                      bounds=(0.4, 0.95))

    b.air_util = Var(initialize=0.5,
                     units=None,
                     bounds=(0.125, 0.833))

    b.pressure = Var(initialize=1,
                     units=units.atm,
                     bounds=(1, 2.5))

    # create vars for intermediate calculations
    ROM_initialize_values = [4000, 500, 0.4, 700, 0.5, 2.1, .85, 0.5, 1]
    b.ROM_input = Var(input_index,
                      initialize=build_dict(input_index,
                                            ROM_initialize_values))

    b.norm_input = Var(input_index, initialize=0)

    b.F = Var(input_plus_index, initialize=0)
    b.F[0].fix(1)

    b.R = Var(samples_index, initialize=0)

    b.norm_output = Var(output_index, initialize=0)

    b.ROM_output = Var(output_index)

    # create kriging regression constraints

    # this dict maps the index values to the input vars
    input_map = {0: b.current_density,
                 1: b.fuel_temperature,
                 2: b.internal_reforming,
                 3: b.air_temperature,
                 4: b.air_recirculation,
                 5: b.OTC,
                 6: b.fuel_util,
                 7: b.air_util,
                 8: b.pressure}

    def input_rule(b, i):
        if units.get_units(input_map[i]) is None:
            return b.ROM_input[i] == input_map[i]
        else:
            unit_conversion = units.get_units(input_map[i])
            return b.ROM_input[i] == input_map[i]/unit_conversion

    b.input_mapping_eqs = Constraint(input_index, rule=input_rule)

    def norm_input_rule(b, i):
        return (b.norm_input[i] ==
                (b.ROM_input[i] - b.mean_input[i])/b.sigma_input[i])

    b.norm_input_eqs = Constraint(input_index, rule=norm_input_rule)

    def F_rule(b, i):
        return b.F[i+1] == b.norm_input[i]

    b.F_eqs = Constraint(input_index, rule=F_rule)

    def R_rule(b, i):
        return (b.R[i] == exp(-1*sum(b.theta[j] *
                                     (b.ds_input[i, j] - b.norm_input[j])**2
                                     for j in input_index)))

    b.R_eqs = Constraint(samples_index, rule=R_rule)

    def norm_output_rule(b, i):
        return (b.norm_output[i] ==
                sum(b.F[j]*b.beta[j, i] for j in input_plus_index)
                + sum(b.R[k]*b.gamma[k, i] for k in samples_index))

    b.norm_output_eqs = Constraint(output_index, rule=norm_output_rule)

    def ROM_output_rule(b, i):
        return (b.ROM_output[i] ==
                b.mean_output[i] + b.norm_output[i]*b.sigma_output[i])

    b.ROM_output_eqs = Constraint(output_index, rule=ROM_output_rule)

    # create output variables and constraints
    b.anode_outlet_temperature = Var(initialize=600, units=None)
    b.cathode_outlet_temperature = Var(initialize=600, units=None)
    b.stack_voltage = Var(initialize=1, units=units.V)
    b.max_cell_temperature = Var(initialize=750, units=None)
    b.deltaT_cell = Var(initialize=100, units=None)

    def anode_outlet_rule(b):
        return b.anode_outlet_temperature == b.ROM_output[11]
    b.anode_outlet_eq = Constraint(rule=anode_outlet_rule)

    def cathode_outlet_rule(b):
        return b.cathode_outlet_temperature == b.ROM_output[13]
    b.cathode_outlet_eq = Constraint(rule=cathode_outlet_rule)

    def stack_voltage_rule(b):
        return b.stack_voltage == b.ROM_output[1]*units.V
    b.stack_voltage_eq = Constraint(rule=stack_voltage_rule)

    def max_cell_temp_rule(b):
        return b.max_cell_temperature == b.ROM_output[8]
    b.max_cell_temp_eq = Constraint(rule=max_cell_temp_rule)

    def deltaT_cell_rule(b):
        return b.deltaT_cell == b.ROM_output[10]
    b.deltaT_cell_eq = Constraint(rule=deltaT_cell_rule)


# initialization procedure for the ROM
def initialize_SOFC_ROM(b):

    def cvfc_indexed(variable, constraint):
        for i in constraint.keys():
            calculate_variable_from_constraint(variable[i],
                                               constraint[i])

    # F needs to be done manually
    def calculate_F():
        for i in b.norm_input.keys():
            b.F[i+1].fix(value(b.norm_input[i]))
            b.F[i+1].unfix()

    print('Starting ROM initialization')

    cvfc_indexed(b.ROM_input, b.input_mapping_eqs)
    cvfc_indexed(b.norm_input, b.norm_input_eqs)
    calculate_F()
    cvfc_indexed(b.R, b.R_eqs)
    cvfc_indexed(b.norm_output, b.norm_output_eqs)
    cvfc_indexed(b.ROM_output, b.ROM_output_eqs)
    cvfc_indexed(b.anode_outlet_temperature, b.anode_outlet_eq)
    cvfc_indexed(b.cathode_outlet_temperature, b.cathode_outlet_eq)
    cvfc_indexed(b.stack_voltage, b.stack_voltage_eq)
    cvfc_indexed(b.max_cell_temperature, b.max_cell_temp_eq)
    cvfc_indexed(b.deltaT_cell, b.deltaT_cell_eq)

    print('ROM initialization completed')
