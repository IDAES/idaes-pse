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
Power Plant costing library
This method leverages NETL costing capabilities. Two main methods have been
developed to calculate the capital cost of power generation plants:

1.- Fossil fueled power plants (from SCPC to IGCC) (get_PP_costing)
2.- supercritical CO2 power cycles (direct and indirect) (get_sCO2_unit_cost)

other methods:
* get_ASU_cost() to cost air separation units
* costing_initialization() to initialize costing blocks
* display_total_plant_costs() to display total plant cost
* display_bare_erected_costs() to display BEC costs
* build_flowsheet_cost_constraint() to display the total cost of the entire
flowsheet
* display_flowsheet_cost() to display flowsheet cost
* check_sCO2_costing_bounds() to display a warnning if costing model have been
used outside the range that where designed for
"""
__author__ = "Costing Team (A. Noring and M. Zamarripa)"
__version__ = "1.0.0"

from pyomo.environ import Param, Var, Block, Constraint, Expression, value, \
    Expr_if
import idaes.core.util.scaling as iscale
from idaes.power_generation.costing.costing_dictionaries import \
    BB_costing_exponents, BB_costing_params, sCO2_costing_params
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# -----------------------------------------------------------------------------
# Power Plant Costing Library
# -----------------------------------------------------------------------------


def get_PP_costing(self, cost_accounts,
                   scaled_param, units, tech, ccs='B'):
    '''
    Power Plant Costing Method
    This method relies on the capital cost scaling methodologies developed
    by NETL. Report #DOE/NETL-341/013113
    Multiple vendors quotes have been used to determine the cost of several
    plant equipments (i.e. boiler, pumps, heat exchangers, etc.), other cost
    incurred during the plant operation (i.e. solids handling, etc.)

    Scaling approach uses one main equation:
        SC = RC*(SP/RP)^Exp
    where:
        SC is the scaled cost
        RC is the reference cost
        SP is the scaled operational parameter
        RP is the reference operational parameter
        Exp is the scaling exponent

    The scaled cost is computed using ref values for different technoligies.
    Categories:
    1 - Supercritical PC, air-fired, with and without CO2 capture,
    Illinois No. 6 coal
    2 - Subcritical PC, air-fired, with and without CO2 capture,
    Illinois No. 6 coal
    3 - Two-stage, slurry-feed, oxygen-blown gasifier with and without
    CO2 capture, Illinois No. 6 coal
    4 - Single-stage, slurry-feed, oxygen-blown gasifier with and without
    CO2 capture, Illinois No. 6 coal
    5 - Single-stage, dry-feed, oxygen-blown, up-flow gasifier with
    and without CO2 capture, Illinois No. 6 coal
    6 - Natural gas, air-fired, with and without CO2 capture
    7 - Advanced Ultrasupercritical PC

    This method computes the capital cost of units and main components of the
    power plant, and requires a few arguments to build a constraint as part of
    your main model.

    Args:
    * self: A block or unit model where costing constraints can be added to
    * accounts: A list of accounts to be included in the total cost,
    they should all use the same reference parameter
    * scaled_param: the process parameter for the system(s) being costed
    * units: the units of the scaled_param, used for verification
    * tech: int 1-7 representing the above catagories
    * ccs: 'A' or 'B' representing no CCS or CCS

    The appropriate scaling parameters for various cost accounts can be found
    in the QGESS on capital cost scaling (Report #DOE/NETL-341/013113).
    The correct units for the reference parameters are found in the BBR4 COE
    spreadsheet.

    '''
    # ------------------------ Power Plant Cost ------------------------

    # check to see if a costing block already exists
    if hasattr(self, 'costing'):
        raise AttributeError("{} already has an attribute costing. "
                             "Check that you are not calling get_costing"
                             " twice on the same model".format(self.name))

    # create a costing Block
    self.costing = Block()
    self.costing.library = 'PP'

    # find flowsheet block to create global costing parameters
    try:
        fs = self.flowsheet()
    except AttributeError:
        fs = self.parent_block()

    # build flowsheet level parameters CE_index = year 
    if not hasattr(fs, 'costing'):
        fs.get_costing(year='2018')

    CE_index = fs.costing.CE_index

    # define preloaded accounts
    PC_preloaded_accounts = {'Coal Handling': ['1.1', '1.2',
                                               '1.3', '1.4', '1.9a'],
                             'Sorbent Handling': ['1.5', '1.6',
                                                  '1.7', '1.8', '1.9b'],
                             'Coal Feed': ['2.1', '2.2', '2.9a'],
                             'Sorbent Feed': ['2.5', '2.6', '2.9b'],
                             'Feedwater System': ['3.1', '3.3'],
                             'PC Boiler': ['4.9'],
                             'Steam Turbine': ['8.1'],
                             'Condenser': ['8.3'],
                             'Cooling Tower': ['9.1'],
                             'Circulating Water System': ['9.2', '9.3',
                                                          '9.4', '9.6', '9.7'],
                             'Ash Handling': ['10.6', '10.7', '10.9']}

    IGCC_preloaded_accounts = {'Coal Handling': ['1.1', '1.2',
                                                 '1.3', '1.4', '1.9'],
                               'Coal Feed': ['2.1', '2.2',
                                             '2.3', '2.4', '2.9'],
                               'Feedwater System': ['3.1', '3.3'],
                               'Gasifier': ['4.1'],
                               'Syngas Cooler': ['4.2'],
                               'ASU': ['4.3a'],
                               'ASU Oxidant Compression': ['4.3b'],
                               'Combustion Turbine': ['6.1', '6.3'],
                               'Syngas Expander': ['6.2'],
                               'HRSG': ['7.1', '7.2'],
                               'Steam Turbine': ['8.1'],
                               'Condenser': ['8.3'],
                               'Cooling Tower': ['9.1'],
                               'Circulating Water System': ['9.2', '9.3',
                                                            '9.4', '9.6',
                                                            '9.7'],
                               'Slag Handling': ['10.1', '10.2',
                                                 '10.3', '10.6',
                                                 '10.7', '10.8',
                                                 '10.9']}

    NGCC_preloaded_accounts = {'Feedwater System': ['3.1', '3.3'],
                               'Combustion Turbine': ['6.1', '6.3'],
                               'HRSG': ['7.1', '7.2'],
                               'Steam Turbine': ['8.1'],
                               'Condenser': ['8.3'],
                               'Cooling Tower': ['9.1'],
                               'Circulating Water System': ['9.2', '9.3',
                                                            '9.4', '9.6',
                                                            '9.7']}

    AUSC_preloaded_accounts = {'PC Boiler': ['4.9'],
                               'Steam Turbine': ['8.1'],
                               'Steam Piping': ['8.4']}

    # preloaded account handling
    if type(cost_accounts) == str:
        if tech in [1, 2]:
            cost_accounts = PC_preloaded_accounts[cost_accounts]
        elif tech in [3, 4, 5]:
            cost_accounts = IGCC_preloaded_accounts[cost_accounts]
        elif tech == 6:
            cost_accounts = NGCC_preloaded_accounts[cost_accounts]
        elif tech == 7:
            cost_accounts = AUSC_preloaded_accounts[cost_accounts]
        else:
            AttributeError("{} technology not supported".format(self.name))

    # check that all accounts use the same process parameter
    param_check = None
    for account in cost_accounts:
        param = BB_costing_exponents[str(tech)][account]['Process Parameter']
        if param_check is None:
            param_check = param
        elif param != param_check:
            raise ValueError("{} cost accounts selected do not use "
                             " the same process parameter".format(self.name))

    # check that the user passed the correct units
    ref_units = BB_costing_params[str(tech)][ccs][cost_accounts[0]]['Units']
    if units != ref_units:
        raise ValueError('Account %s uses units of %s. '
                         'Units of %s were passed.'
                         % (cost_accounts[0], ref_units, units))

    # construct dictionaries
    account_names = {}
    exponents = {}
    reference_costs = {}
    reference_params = {}
    engineering_fees = {}
    process_contingencies = {}
    project_contingencies = {}

    for account in cost_accounts:
        account_names[account] = BB_costing_exponents[str(
            tech)][account]['Account Name']
        exponents[account] = float(
            BB_costing_exponents[str(tech)][account]['Exponent'])
        reference_costs[account] = BB_costing_params[str(
            tech)][ccs][account]['BEC']
        reference_params[account] = BB_costing_params[str(
            tech)][ccs][account]['RP Value']
        engineering_fees[account] = BB_costing_params[str(
            tech)][ccs][account]['Eng Fee']
        process_contingencies[account] = BB_costing_params[str(
            tech)][ccs][account]['Process Contingency']
        project_contingencies[account] = BB_costing_params[str(
            tech)][ccs][account]['Project Contingency']

    # Used by other functions for reporting results
    self.costing.account_names = account_names

    # define parameters
    self.costing.exp = Param(cost_accounts,
                             mutable=True,
                             initialize=exponents,
                             doc='exponential parameter for account')

    self.costing.ref_cost = Param(cost_accounts,
                                  mutable=True,
                                  initialize=reference_costs,
                                  doc='reference cost for account')

    self.costing.ref_param = Param(cost_accounts,
                                   mutable=True,
                                   initialize=reference_params,
                                   doc='reference parameter for account')

    self.costing.eng_fee = Param(cost_accounts,
                                 mutable=True,
                                 initialize=engineering_fees,
                                 doc='engineering fee percentage')

    self.costing.process_conting = Param(cost_accounts,
                                         mutable=True,
                                         initialize=process_contingencies,
                                         doc='process contingency percentage')

    self.costing.project_conting = Param(cost_accounts,
                                         mutable=True,
                                         initialize=project_contingencies,
                                         doc='project contingency percentage')

    # define variables
    self.costing.bare_erected_cost = Var(cost_accounts,
                                         initialize=reference_costs,
                                         bounds=(0, 1e4),
                                         doc='scaled bare erected cost in $MM')

    self.costing.total_plant_cost = Var(cost_accounts,
                                        initialize=reference_costs,
                                        bounds=(0, 1e4),
                                        doc='total plant cost in $MM')

    # rule for scaling BEC
    # reference cost is in 2018 dollars, 671.1 is CE index for 2018
    def bare_erected_cost_rule(costing, i):
        return (costing.bare_erected_cost[i]*1e3 ==
                (CE_index/671.1)*costing.ref_cost[i] *
                (scaled_param/costing.ref_param[i])**costing.exp[i])

    self.costing.bare_erected_cost_eq = Constraint(
        cost_accounts, rule=bare_erected_cost_rule)

    # rule for calculating TPC
    def total_plant_cost_rule(costing, i):
        return (costing.total_plant_cost[i] == costing.bare_erected_cost[i] *
                (1 + costing.eng_fee[i] + costing.process_conting[i]) *
                (1 + costing.project_conting[i]))

    self.costing.total_plant_cost_eq = Constraint(
        cost_accounts, rule=total_plant_cost_rule)

    # rule for sum of BEC
    def BEC_sum_rule(costing):
        return sum(costing.bare_erected_cost[i] for i in cost_accounts)
    self.costing.bare_erected_cost_sum = Expression(rule=BEC_sum_rule)

    # rule for sum of TPC
    def TPC_sum_rule(costing):
        return sum(costing.total_plant_cost[i] for i in cost_accounts)
    self.costing.total_plant_cost_sum = Expression(rule=TPC_sum_rule)

    #  # add variable and constraint scaling
    for i in cost_accounts:
        iscale.set_scaling_factor(self.costing.bare_erected_cost[i], 1)
        iscale.set_scaling_factor(self.costing.total_plant_cost[i], 1)
        iscale.constraint_scaling_transform(self.
                                            costing.bare_erected_cost_eq[i],
                                            1e-3,
                                            overwrite=False)
        iscale.constraint_scaling_transform(self.
                                            costing.total_plant_cost_eq[i],
                                            1,
                                            overwrite=False)


# -----------------------------------------------------------------------------
# Supercritical CO2 Costing Library
# -----------------------------------------------------------------------------
def get_sCO2_unit_cost(self, equipment, scaled_param, temp_C=None, n_equip=1):
    '''
    Args:
    self - pyomo Block where constraints will be made
    unit_name - the name of the SCO2 equipment to cost
    scaling_param - the scaling parameter (in appropriate units) for the
                    selected equipment
    temp - the maximum temperature of the equipment. Not all types of
           equipment use a temperature correction factor, so it is optional
    n_equip - the number of pieces of equipment to cost

    Cost is in M$
    '''
    # check to see if a costing block already exists
    if hasattr(self, 'costing'):
        raise AttributeError("{} already has an attribute costing. "
                             "Check that you are not calling get_costing"
                             " twice on the same model".format(self.name))

    # create a costing Block
    self.costing = Block()
    self.costing.library = 'sCO2'
    self.costing.equipment = equipment

    # find flowsheet block to create global costing parameters
    try:
        fs = self.flowsheet()
    except AttributeError:
        fs = self.parent_block()

    # build flowsheet level parameters CE_index = year 
    if not hasattr(fs, 'costing'):
        fs.get_costing(year='2017')

    CE_index = fs.costing.CE_index

    param_dict = sCO2_costing_params[equipment]

    # define parameters
    self.costing.ref_cost = Param(mutable=True,
                                  initialize=param_dict['a'],
                                  doc='Reference cost')

    self.costing.exp = Param(mutable=True,
                             initialize=param_dict['b'],
                             doc='Scaling exponent')

    self.costing.c = Param(mutable=True,
                           initialize=param_dict['c'],
                           doc='coefficient for temperature correction')

    self.costing.d = Param(mutable=True,
                           initialize=param_dict['d'],
                           doc='coefficient for temperature correction')

    self.costing.material_cost = Param(mutable=True,
                                       doc='material installation cost',
                                       initialize=param_dict['Material Cost'])

    self.costing.labor_cost = Param(mutable=True,
                                    initialize=param_dict['Labor Cost'],
                                    doc='labor installation cost')

    # estimates for the percentages of TPC will be added later

    self.costing.eng_fee = Param(mutable=True,
                                 initialize=0,
                                 doc='engineering fee percentage')

    self.costing.process_conting = Param(mutable=True,
                                         initialize=0,
                                         doc='process contingency percentage')

    self.costing.project_conting = Param(mutable=True,
                                         initialize=0,
                                         doc='project contingency percentage')

    # define variables
    # n_equip is left as a fixed variable to support MINLP optimization
    self.costing.n_equip = Var(initialize=n_equip,
                               doc='number of pieces of equipment')

    self.costing.n_equip.fix(n_equip)

    self.costing.scaled_param = Var(initialize=scaled_param,
                                    bounds=(0, 1e12),
                                    doc='scaled parameter')

    self.costing.temp_factor = Var(initialize=1,
                                   bounds=(0.9, 100),
                                   doc='temperature correction factor')

    self.costing.equipment_cost = Var(initialize=self.costing.ref_cost,
                                      bounds=(0, 1e4),
                                      doc='equipment cost of sCO2 unit in $MM')

    self.costing.bare_erected_cost = Var(initialize=self.costing.ref_cost,
                                         bounds=(0, 1e4),
                                         doc='bare erected cost of sCO2 unit'
                                         'in $MM')

    self.costing.total_plant_cost = Var(initialize=self.costing.ref_cost,
                                        bounds=(0, 1e4),
                                        doc='total plant cost of sCO2 unit'
                                        'in $MM')

    # divides the scaled parameter by the number of pieces of equipment
    def scaled_param_rule(costing):
        return costing.scaled_param*costing.n_equip == scaled_param

    self.costing.scaled_param_eq = Constraint(rule=scaled_param_rule)

    # check if equipment requires a temperature correction factor
    if equipment in ['Axial turbine', 'Radial turbine', 'Coal-fired heater',
                     'Natural gas-fired heater', 'Recuperator']:

        if temp_C is None:
            raise ValueError('Temperature argument is '
                             'required to cost %s equipment' % equipment)

        else:
            self.costing.temperature = Var(initialize=500,
                                           bounds=(0, 1e6),
                                           doc='dummy var for temperature')

            self.costing.temp_eq = Constraint(expr=(self.costing.temperature
                                                    == temp_C))

            def temp_correction_rule(costing):  # rule for temp correction
                return (Expr_if(costing.temperature < 550,
                                1e-6*costing.temperature + 1,
                                1 + costing.c*(costing.temperature - 550)
                                + costing.d*(costing.temperature - 550)**2) ==
                        costing.temp_factor)

            self.costing.temp_correction_eq = Constraint(
                rule=temp_correction_rule)
    else:
        self.costing.temp_factor.fix(1)

    # rule for equipment cost
    def equipment_cost_rule(costing):
        return (costing.equipment_cost*1e6 ==
                (CE_index/567.5) * costing.n_equip * costing.ref_cost *
                (costing.scaled_param**costing.exp) * costing.temp_factor)

    self.costing.equipment_cost_eq = Constraint(rule=equipment_cost_rule)

    # rule for bare erected cost
    def bare_erected_cost_rule(costing):
        return (costing.bare_erected_cost == costing.equipment_cost *
                (1 + costing.material_cost + costing.labor_cost))

    self.costing.bare_erected_cost_eq = Constraint(rule=bare_erected_cost_rule)

    # rule for calculating total plant cost
    def total_plant_cost_rule(costing):
        return (costing.total_plant_cost == costing.bare_erected_cost *
                (1 + costing.eng_fee + costing.process_conting
                 + costing.project_conting))
    self.costing.total_plant_cost_eq = Constraint(rule=total_plant_cost_rule)

    # add variable and constraint scaling
    if equipment in ["Recuperator", "Direct air cooler"]:
        iscale.set_scaling_factor(self.costing.scaled_param, 1e-5)
    else:
        iscale.set_scaling_factor(self.costing.scaled_param, 1)

    iscale.set_scaling_factor(self.costing.equipment_cost, 1e3)
    iscale.set_scaling_factor(self.costing.bare_erected_cost, 1e3)
    iscale.set_scaling_factor(self.costing.total_plant_cost, 1e3)
    iscale.constraint_scaling_transform(
        self.costing.equipment_cost_eq, 1e-6, overwrite=False)
    iscale.constraint_scaling_transform(
        self.costing.bare_erected_cost_eq, 1e3, overwrite=False)
    iscale.constraint_scaling_transform(
        self.costing.bare_erected_cost_eq, 1e3, overwrite=False)


# -----------------------------------------------------------------------------
# Air Separation Unit Costing Library
# -----------------------------------------------------------------------------
def get_ASU_cost(self, scaled_param):
    # scaled parameter is O2 flowrate in TPD

    params = {'Reference Cost': 3.26e6,
              'Reference Parameter': 13078,
              'Exponent': 0.7,
              'Eng Fee': 0.097,
              'Process': 0,
              'Project': 0.110}

    # check to see if a costing block already exists
    if hasattr(self, 'costing'):
        raise AttributeError("{} already has an attribute costing. "
                             "Check that you are not calling get_costing"
                             " twice on the same model".format(self.name))

    # create a costing Block
    self.costing = Block()
    self.costing.library = 'ASU'

    # find flowsheet block to create global costing parameters
    try:
        fs = self.flowsheet()
    except AttributeError:
        fs = self.parent_block()

    # build flowsheet level parameters CE_index = year 
    if not hasattr(fs, 'costing'):
        fs.get_costing(year='2017')

    CE_index = fs.costing.CE_index

    # define parameters
    self.costing.ref_cost = Param(initialize=params['Reference Cost'],
                                  mutable=True,
                                  doc='ASU reference cost')

    self.costing.ref_param = Param(initialize=params['Reference Parameter'],
                                   mutable=True,
                                   doc='ASU reference parameter value')

    self.costing.exp = Param(initialize=params['Exponent'],
                             mutable=True,
                             doc='ASU scaling exponent')

    self.costing.eng_fee = Param(mutable=True,
                                 initialize=params['Eng Fee'],
                                 doc='engineering fee percentage')

    self.costing.process_conting = Param(mutable=True,
                                         initialize=params['Process'],
                                         doc='process contingency percentage')

    self.costing.project_conting = Param(mutable=True,
                                         initialize=params['Project'],
                                         doc='project contingency percentage')

    # define variables
    self.costing.bare_erected_cost = Var(initialize=params['Reference Cost'],
                                         bounds=(0, 1e4),
                                         doc='scaled bare erected cost in $MM')

    self.costing.total_plant_cost = Var(initialize=params['Reference Cost'],
                                        bounds=(0, 1e4),
                                        doc='total plant cost in $MM')

    # rule for scaling BEC
    # reference cost is in 2008 dollars, 566.2 is CE index for Nov 2008

    def bare_erected_cost_rule(costing):
        return (costing.bare_erected_cost*1e3 ==
                (CE_index/566.2)*costing.ref_cost *
                (scaled_param/costing.ref_param)**costing.exp)

    self.costing.bare_erected_cost_eq = Constraint(rule=bare_erected_cost_rule)

    # rule for calculating TPC
    def total_plant_cost_rule(costing):
        return (costing.total_plant_cost == costing.bare_erected_cost *
                (1 + costing.eng_fee + costing.process_conting
                 + costing.project_conting))
    self.costing.total_plant_cost_eq = Constraint(rule=total_plant_cost_rule)

    # add variable and constraint scaling
    iscale.set_scaling_factor(self.costing.bare_erected_cost, 1)
    iscale.set_scaling_factor(self.costing.total_plant_cost, 1)

    iscale.constraint_scaling_transform(
        self.costing.bare_erected_cost_eq, 1e-3, overwrite=False)
    iscale.constraint_scaling_transform(
        self.costing.total_plant_cost_eq, 1, overwrite=False)


# -----------------------------------------------------------------------------
# Costing Library Utility Functions
# -----------------------------------------------------------------------------
def costing_initialization(fs):
    for o in fs.component_objects(descend_into=False):
        # look for costing blocks
        if hasattr(o, 'costing'):
            if o.costing.library == 'sCO2':

                if o.costing.equipment in ['Axial turbine',
                                           'Radial turbine',
                                           'Coal-fired heater',
                                           'Natural gas-fired heater',
                                           'Recouperator']:
                    calculate_variable_from_constraint(o.costing.temperature,
                                                       o.costing.temp_eq)
                    calculate_variable_from_constraint(o.costing.temp_factor,
                                                       o.costing.
                                                       temp_correction_eq)

                calculate_variable_from_constraint(o.costing.scaled_param,
                                                   o.costing.scaled_param_eq)
                calculate_variable_from_constraint(o.costing.equipment_cost,
                                                   o.costing.equipment_cost_eq)
                calculate_variable_from_constraint(o.costing.bare_erected_cost,
                                                   o.costing.
                                                   bare_erected_cost_eq)
                calculate_variable_from_constraint(o.costing.total_plant_cost,
                                                   o.costing.
                                                   total_plant_cost_eq)

            elif o.costing.library in ['PP', 'ASU']:
                for key in o.costing.bare_erected_cost.keys():
                    calculate_variable_from_constraint(o.costing.
                                                       bare_erected_cost[key],
                                                       o.costing.
                                                       bare_erected_cost_eq[
                                                           key])
                    calculate_variable_from_constraint(o.costing.
                                                       total_plant_cost[key],
                                                       o.costing.
                                                       total_plant_cost_eq[
                                                           key])


def display_total_plant_costs(fs):
    print('-----Total Plant Costs-----')
    for o in fs.component_objects(descend_into=False):
        # look for costing blocks
        if hasattr(o, 'costing') and hasattr(o.costing, 'total_plant_cost'):
            print('%s: $%.2f Million' % (value(o.name),
                                         value(o.costing.total_cost)))


def display_bare_erected_costs(fs):
    print('-----Bare Erected Costs-----')
    for o in fs.component_objects(descend_into=False):
        # look for costing blocks
        if hasattr(o, 'costing') and hasattr(o.costing, 'bare_erected_cost'):
            print('%s: $%.2f Million' % (value(o.name),
                                         value(o.costing.bare_erected_cost)))


def display_equipment_costs(fs):
    print('-----Equipment Costs-----')
    for o in fs.component_objects(descend_into=False):
        # look for costing blocks
        if hasattr(o, 'costing') and hasattr(o.costing, 'equipment_cost'):
            print('%s: $%.2f Million' % (value(o.name),
                                         value(o.costing.equipment_cost)))


def build_flowsheet_cost_constraint(m):
    total_cost_list = []
    for o in m.fs.component_objects(descend_into=False):
        # look for costing blocks
        if hasattr(o, 'costing'):
            for key in o.costing.total_plant_cost.keys():
                total_cost_list.append(o.costing.total_plant_cost[key])

    m.fs.flowsheet_cost = Var(initialize=0,
                              bounds=(0, 1e12),
                              doc='cost of entire process')

    def flowsheet_cost_rule(fs):
        return fs.flowsheet_cost == sum(total_cost_list)

    m.fs.flowsheet_cost_eq = Constraint(rule=flowsheet_cost_rule)


def display_flowsheet_cost(m):
    print('\n')
    print('Total flowsheet cost: $%.3f Million' %
          value(m.fs.flowsheet_cost_exp))


def check_sCO2_costing_bounds(fs):
    # interate through the children of the flowsheet
    for o in fs.component_objects(descend_into=False):
        # look for costing blocks
        if hasattr(o, 'costing'):
            costing = o.costing
            if costing.library == 'sCO2':
                equipment = costing.equipment
                lower_bound = sCO2_costing_params[equipment]['Lower Bound']
                upper_bound = sCO2_costing_params[equipment]['Upper Bound']
                if value(costing.scaled_param) < lower_bound:
                    print('''%s: The scaled parameter (%f) is below the lower
                          bound (%f).''' % (value(o.name),
                                            value(costing.scaled_param),
                                            lower_bound))
                elif value(costing.scaled_param) > upper_bound:
                    print('''%s: The scaled parameter (%f) is above the upper
                          bound (%f).''' % (value(o.name),
                                            value(costing.scaled_param),
                                            upper_bound))
                else:
                    print('''%s: The scaled parameter is within the bounds.'''
                          % value(o.name))
