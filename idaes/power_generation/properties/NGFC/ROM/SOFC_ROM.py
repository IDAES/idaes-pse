##############################################################################
# The development of this flowsheet/code is funded by the ARPA-E DIFFERENTIATE
# project: “Machine Learning for Natural Gas to Electric Power System Design”
# Project number: DE-FOA-0002107-1625.
# This project is a collaborative effort between the Pacific Northwest National
# Laboratory, the National Energy Technology Laboratory, and the University of
# Washington to design NGFC systems with high efficiencies and low CO2
# emissions.
##############################################################################
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

from pyomo.environ import Constraint, Param, Var, exp, value, units
from pyomo.util.calc_var_value import calculate_variable_from_constraint


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


def make_SOFC_ROM(self):
    # load kriging coefficients
    ROM_filename = 'ROM/kriging_coefficients.dat'

    with open(ROM_filename) as file:
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
    self.mean_input = Param(input_index,
                            initialize=build_dict(input_index, values),
                            mutable=False)

    start_count += n_inputs
    end_count += n_inputs
    values = kriging[start_count:end_count]
    self.sigma_input = Param(input_index,
                             initialize=build_dict(input_index, values),
                             mutable=False)

    start_count += n_inputs
    end_count += n_outputs
    values = kriging[start_count:end_count]
    self.mean_output = Param(output_index,
                             initialize=build_dict(output_index, values),
                             mutable=False)

    start_count += n_outputs
    end_count += n_outputs
    values = kriging[start_count:end_count]
    self.sigma_output = Param(output_index,
                              initialize=build_dict(output_index, values),
                              mutable=False)

    start_count += n_outputs
    end_count += n_inputs*n_samples
    values = kriging[start_count:end_count]
    self.ds_input = Param(samples_index, input_index,
                          initialize=build_matrix(
                              samples_index, input_index, values),
                          mutable=False)

    start_count += n_inputs*n_samples
    end_count += n_inputs
    values = kriging[start_count:end_count]
    self.theta = Param(input_index,
                       initialize=build_dict(input_index, values),
                       mutable=False)

    start_count += n_inputs
    end_count += (n_inputs+1)*n_outputs
    values = kriging[start_count:end_count]
    self.beta = Param(input_plus_index, output_index,
                      initialize=build_matrix(
                          input_plus_index, output_index, values),
                      mutable=False)

    start_count += (n_inputs+1)*n_outputs
    end_count += n_samples*n_outputs
    values = kriging[start_count:end_count]
    self.gamma = Param(samples_index, output_index,
                       initialize=build_matrix(
                           samples_index, output_index, values),
                       mutable=False)

    # create input vars for the user to interface with
    self.current_density = Var(initialize=4000,
                               units=units.A/units.m**2,
                               bounds=(2000, 6000))

    self.fuel_temperature = Var(initialize=500,
                                units=units.C,
                                bounds=(15, 600))

    self.internal_reforming = Var(initialize=0.4,
                                  units=None,
                                  bounds=(0, 1))

    self.air_temperature = Var(initialize=700,
                               units=units.C,
                               bounds=(550, 800))

    self.air_recirculation = Var(initialize=0.5,
                                 units=None,
                                 bounds=(0, 0.8))

    self.OTC = Var(initialize=2.1,
                   units=None,
                   bounds=(1.5, 3))

    self.fuel_util = Var(initialize=0.85,
                         units=None,
                         bounds=(0.4, 0.95))

    self.air_util = Var(initialize=0.5,
                        units=None,
                        bounds=(0.125, 0.833))

    self.pressure = Var(initialize=1,
                        units=units.atm,
                        bounds=(1, 2.5))

    # create vars for intermediate calculations
    ROM_initialize_values = [4000, 500, 0.4, 700, 0.5, 2.1, .85, 0.5, 1]
    self.ROM_input = Var(input_index,
                         initialize=build_dict(input_index,
                                               ROM_initialize_values))

    self.norm_input = Var(input_index, initialize=0)

    self.F = Var(input_plus_index, initialize=0)
    self.F[0].fix(1)

    self.R = Var(samples_index, initialize=0)

    self.norm_output = Var(output_index, initialize=0)

    self.ROM_output = Var(output_index)

    # create kriging regression constraints

    # this dict maps the index values to the input vars
    input_map = {0: self.current_density,
                 1: self.fuel_temperature,
                 2: self.internal_reforming,
                 3: self.air_temperature,
                 4: self.air_recirculation,
                 5: self.OTC,
                 6: self.fuel_util,
                 7: self.air_util,
                 8: self.pressure}

    def input_rule(b, i):
        return b.ROM_input[i] == input_map[i]

    self.input_mapping_eqs = Constraint(input_index, rule=input_rule)

    def norm_input_rule(b, i):
        return (b.norm_input[i] ==
                (b.ROM_input[i] - b.mean_input[i])/b.sigma_input[i])

    self.norm_input_eqs = Constraint(input_index, rule=norm_input_rule)

    def F_rule(b, i):
        return b.F[i+1] == b.norm_input[i]

    self.F_eqs = Constraint(input_index, rule=F_rule)

    def R_rule(b, i):
        return (b.R[i] == exp(-1*sum(b.theta[j] *
                                     (b.ds_input[i, j] - b.norm_input[j])**2
                                     for j in input_index)))

    self.R_eqs = Constraint(samples_index, rule=R_rule)

    def norm_output_rule(b, i):
        return (b.norm_output[i] ==
                sum(b.F[j]*b.beta[j, i] for j in input_plus_index)
                + sum(b.R[k]*b.gamma[k, i] for k in samples_index))

    self.norm_output_eqs = Constraint(output_index, rule=norm_output_rule)

    def ROM_output_rule(b, i):
        return (b.ROM_output[i] ==
                b.mean_output[i] + b.norm_output[i]*b.sigma_output[i])

    self.ROM_output_eqs = Constraint(output_index, rule=ROM_output_rule)

    # create output variables and constraints
    self.anode_outlet_temperature = Var(initialize=600, units=units.C)
    self.cathode_outlet_temperature = Var(initialize=600, units=units.C)
    self.stack_voltage = Var(initialize=1, units=units.V)
    self.max_cell_temperature = Var(initialize=750, units=units.C)
    self.deltaT_cell = Var(initialize=100, units=units.C)

    def anode_outlet_rule(b):
        return b.anode_outlet_temperature == b.ROM_output[11]*units.C
    self.anode_outlet_eq = Constraint(rule=anode_outlet_rule)

    def cathode_outlet_rule(b):
        return b.cathode_outlet_temperature == b.ROM_output[13]*units.C
    self.cathode_outlet_eq = Constraint(rule=cathode_outlet_rule)

    def stack_voltage_rule(b):
        return b.stack_voltage == self.ROM_output[1]*units.V
    self.stack_voltage_eq = Constraint(rule=stack_voltage_rule)

    def max_cell_temp_rule(b):
        return b.max_cell_temperature == self.ROM_output[8]*units.C
    self.max_cell_temp_eq = Constraint(rule=max_cell_temp_rule)

    def deltaT_cell_rule(b):
        return b.deltaT_cell == self.ROM_output[10]*units.C
    self.deltaT_cell_eq = Constraint(rule=deltaT_cell_rule)


# initialization procedure for the ROM
def initialize_SOFC_ROM(self):

    def cvfc_indexed(variable, constraint):
        for i in constraint.keys():
            calculate_variable_from_constraint(variable[i],
                                               constraint[i])

    # F needs to be done manually
    def calculate_F():
        for i in self.norm_input.keys():
            self.F[i+1].fix(value(self.norm_input[i]))
            self.F[i+1].unfix()

    print('Starting ROM initialization')

    cvfc_indexed(self.ROM_input, self.input_mapping_eqs)
    cvfc_indexed(self.norm_input, self.norm_input_eqs)
    calculate_F()
    cvfc_indexed(self.R, self.R_eqs)
    cvfc_indexed(self.norm_output, self.norm_output_eqs)
    cvfc_indexed(self.ROM_output, self.ROM_output_eqs)
    cvfc_indexed(self.anode_outlet_temperature, self.anode_outlet_eq)
    cvfc_indexed(self.cathode_outlet_temperature, self.cathode_outlet_eq)
    cvfc_indexed(self.stack_voltage, self.stack_voltage_eq)
    cvfc_indexed(self.max_cell_temperature, self.max_cell_temp_eq)
    cvfc_indexed(self.deltaT_cell, self.deltaT_cell_eq)

    print('ROM initialization completed')
