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

import os
import numpy as np

from pyomo.environ import SolverFactory, Block, Constraint, Param, Var, Set,\
    exp, value

path = os.path.dirname(os.path.abspath(__file__))


def build_dict(index, values):
    d = {}
    for i, v in zip(index, values):
        d[i] = v
    return d


def build_matrix(index1, index2, values):
    d = {}
    for i in index1:
        for j in index2:
            d[(i, j)] = values[i, j]
    return d


def build_SOFC_ROM(m):
    m.SOFC = self = Block()

    ROM_name = "ADV50_CCS"
    folder = "/sofc_rom_data"

    w1 = np.loadtxt(path+folder+'/rom_w1.dat')
    w2 = np.loadtxt(path+folder+'/rom_w2.dat')
    w3 = np.loadtxt(path+folder+'/rom_w3.dat')
    w4 = np.loadtxt(path+folder+'/rom_w4.dat')
    w5 = np.loadtxt(path+folder+'/rom_w5.dat')

    b1 = np.loadtxt(path+folder+'/rom_b1.dat')
    b2 = np.loadtxt(path+folder+'/rom_b2.dat')
    b3 = np.loadtxt(path+folder+'/rom_b3.dat')
    b4 = np.loadtxt(path+folder+'/rom_b4.dat')
    b5 = np.loadtxt(path+folder+'/rom_b5.dat')

    meanX = np.loadtxt(path+folder+'/rom_meanx.dat')
    meanY = np.loadtxt(path+folder+'/rom_meany.dat')
    stdX = np.loadtxt(path+folder+'/rom_stdx.dat')
    stdY = np.loadtxt(path+folder+'/rom_stdy.dat')

    all_weights = [w1, w2, w3, w4, w5]
    all_biases = [b1, b2, b3, b4, b5]

    layers = len(all_weights)

    if layers == 1:
        raise Exception('''Networks with no hidden layers are currently
                         not supported.''')

    # prepare IndexedSets
    n_X = np.shape(all_weights[0])[0]
    self.X_set = Set(initialize=np.arange(n_X))

    all_sets = [self.X_set]

    # get the list of all_sets from the list of all weights
    for i, w in enumerate(all_weights):
        if i == layers - 1:
            set_name = 'Y_set'
        else:
            set_name = 'hl%d_set' % (i+1)

        n_nodes = np.shape(w)[1]
        set_value = Set(initialize=np.arange(n_nodes))

        setattr(self, set_name, set_value)

        all_sets.append(getattr(self, set_name))

    # Standardization of inputs and outputs

    self.mean_X = Param(self.X_set,
                        mutable=False,
                        initialize=build_dict(self.X_set, meanX))

    self.std_X = Param(self.X_set,
                       mutable=False,
                       initialize=build_dict(self.X_set, stdX))

    self.mean_Y = Param(self.Y_set,
                        mutable=False,
                        initialize=build_dict(self.Y_set, meanY))

    self.std_Y = Param(self.Y_set,
                       mutable=False,
                       initialize=build_dict(self.Y_set, stdY))

    self.input = Var(self.X_set, initialize=0)

    self.norm_input = Var(self.X_set, initialize=0)

    self.output = Var(self.Y_set, initialize=0)

    self.norm_output = Var(self.Y_set, initialize=0)

    def norm_input_rule(b, i):
        return b.norm_input[i] == (b.input[i] - b.mean_X[i])/b.std_X[i]
    self.norm_input_constrant = Constraint(self.X_set, rule=norm_input_rule)

    def output_rule(b, i):
        return b.output[i] == b.norm_output[i]*b.std_Y[i] + b.mean_Y[i]
    self.output_constrant = Constraint(self.Y_set, rule=output_rule)

    # this loop creates parameter and variables
    for lr in range(layers):
        # in this section we create the parameters for the weights and biases
        previous_set = all_sets[lr]
        current_set = all_sets[lr+1]

        w = Param(previous_set, current_set,
                  initialize=build_matrix(previous_set, current_set,
                                          all_weights[lr]),
                  mutable=False)
        setattr(self, 'w%d' % (lr+1), w)

        b = Param(current_set,
                  initialize=build_dict(current_set, all_biases[lr]),
                  mutable=False)
        setattr(self, 'b%d' % (lr+1), b)

        # this section creates the input and output vars for each hidden layer
        if lr != 0:  # not the output layer

            hl_input = Var(previous_set, initialize=0)
            setattr(self, 'hl%d_input' % (lr), hl_input)

            hl_output = Var(previous_set, initialize=0)
            setattr(self, 'hl%d_output' % (lr), hl_output)

    # Creating Hidden Layer Constraints
    # current layer input, previous layer output, weights, biases
    def build_layer_input_constraint(previous_output, current_input,
                                     weights, biases):
        i_set = current_input.index_set()
        j_set = previous_output.index_set()

        def layer_input_rule(b, i):
            return (current_input[i] ==
                    sum(previous_output[j]*weights[j, i] for j in j_set)
                    + biases[i])

        return Constraint(i_set, rule=layer_input_rule)

    def build_layer_output_constraint(current_input, current_output):
        i_set = current_input.index_set()

        def sigmoid(b, i):
            return current_output[i] == 1/(1 + exp(-1*current_input[i]))

        return Constraint(i_set, rule=sigmoid)

    # this loop creates constraints
    for lr in range(layers):
        if lr == 0:  # case for norm_input to hl1_input
            self.hl1_input_constraint = \
                build_layer_input_constraint(self.norm_input, self.hl1_input,
                                             self.w1, self.b1)

            self.hl1_output_constraint = \
                build_layer_output_constraint(self.hl1_input, self.hl1_output)

        elif lr == layers-1:  # case for last hidden layer to output layer
            previous_output = getattr(self, 'hl%d_output' % (lr))
            w = getattr(self, 'w%d' % (lr+1))
            b = getattr(self, 'b%d' % (lr+1))

            self.norm_output_constraint = \
                build_layer_input_constraint(previous_output, self.norm_output,
                                             w, b)

        else:  # case for one hidden layer to the next
            previous_output = getattr(self, 'hl%d_output' % (lr))
            current_input = getattr(self, 'hl%d_input' % (lr+1))
            current_output = getattr(self, 'hl%d_output' % (lr+1))
            w = getattr(self, 'w%d' % (lr+1))
            b = getattr(self, 'b%d' % (lr+1))

            hl_input_constraint = \
                build_layer_input_constraint(previous_output, current_input,
                                             w, b)

            setattr(self, 'hl%d_input_constraint' % (lr+1),
                    hl_input_constraint)

            hl_output_constraint = \
                build_layer_output_constraint(current_input, current_output)

            setattr(self, 'hl%d_output_constraint' % (lr+1),
                    hl_output_constraint)

    # constraints for input vars
    self.current_density = Var(initialize=4000, bounds=(2000, 6000))
    self.fuel_temperature = Var(initialize=500, bounds=(15, 600))
    self.internal_reforming = Var(initialize=0.4, bounds=(0, 1))
    self.air_temperature = Var(initialize=700, bounds=(550, 800))
    self.air_recirculation = Var(initialize=0.5, bounds=(0, 0.8))
    self.OTC = Var(initialize=2.1, bounds=(1.5, 3))
    self.fuel_util = Var(initialize=0.85, bounds=(0.4, 0.95))
    self.air_util = Var(initialize=0.5, bounds=(0.125, 0.833))

    if ROM_name == "ADV50_VGR" or ROM_name == "ADV50_VGR_NOWGS":
        self.VGR_temperature = Var(initialize=300)
        self.VGR_rate = Var(initialize=0.94)

        input_dict = {0: self.current_density,
                      1: self.fuel_temperature,
                      2: self.internal_reforming,
                      3: self.air_temperature,
                      4: self.air_recirculation,
                      5: self.OTC,
                      6: self.fuel_util,
                      7: self.air_util,
                      8: self.VGR_temperature,
                      9: self.VGR_rate}
    else:
        self.pressure = Var(initialize=1, bounds=(1, 5))

        input_dict = {0: self.current_density,
                      1: self.fuel_temperature,
                      2: self.internal_reforming,
                      3: self.air_temperature,
                      4: self.air_recirculation,
                      5: self.OTC,
                      6: self.fuel_util,
                      7: self.air_util,
                      8: self.pressure}

    def input_rule(b, i):
        return b.input[i] == input_dict[i]

    self.input_mapping_eqs = Constraint(self.X_set, rule=input_rule)

    # constraints for output vars
    self.anode_outlet_temperature = Var(initialize=600)
    self.cathode_outlet_temperature = Var(initialize=600)
    self.stack_voltage = Var(initialize=1)
    self.max_cell_temperature = Var(initialize=750)
    self.deltaT_cell = Var(initialize=100)

    def anode_outlet_rule(b):
        return b.anode_outlet_temperature == b.output[10]
    self.anode_outlet_eq = Constraint(rule=anode_outlet_rule)

    def cathode_outlet_rule(b):
        return b.cathode_outlet_temperature == b.output[12]
    self.cathode_outlet_eq = Constraint(rule=cathode_outlet_rule)

    def stack_voltage_rule(b):
        return b.stack_voltage == b.output[0]
    self.stack_voltage_eq = Constraint(rule=stack_voltage_rule)

    def max_cell_temp_rule(b):
        return b.max_cell_temperature == b.output[7]
    self.max_cell_temp_eq = Constraint(rule=max_cell_temp_rule)

    def deltaT_cell_rule(b):
        return b.deltaT_cell == b.output[9]
    self.deltaT_cell_eq = Constraint(rule=deltaT_cell_rule)


def initialize_SOFC_ROM(b):
    # fix all input vars
    input_vars = [b.current_density, b.fuel_temperature, b.internal_reforming,
                  b.air_temperature, b.air_recirculation, b.OTC,
                  b.fuel_util, b.air_util, b.pressure]

    vars_to_unfix = []
    for var in input_vars:
        if var.fixed is False:
            vars_to_unfix.append(var)
        var.fix()

    # unfix all output vars
    output_vars = [b.anode_outlet_temperature, b.cathode_outlet_temperature,
                   b.stack_voltage, b.max_cell_temperature,
                   b.deltaT_cell]

    vars_to_refix = []
    var_values = []
    for var in output_vars:
        if var.fixed is True:
            vars_to_refix.append(var)
            var_values.append(value(var))
        var.unfix()

    # temporary check
    from idaes.core.util.model_statistics import degrees_of_freedom
    assert degrees_of_freedom(b) == 0

    # solve
    solver = SolverFactory("ipopt")
    solver.solve(b, tee=False)

    # correct variables
    for var in vars_to_unfix:
        var.unfix()

    for i in range(len(vars_to_refix)):
        vars_to_refix[i].fix(var_values[i])

