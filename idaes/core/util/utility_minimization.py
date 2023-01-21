# -*- coding: UTF-8 -*-
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
General purpose heat integration block for IDAES models
"""
# Import Plotting and Numpy libraries for composite curves
import matplotlib.pyplot as plt
import numpy as np

# Import Pyomo libraries
from pyomo.environ import Constraint, sqrt, Expression, value, Var

from pyomo.environ import units as pyunits

__author__ = "Alejandro Garciadiego, Alexander Dowling"

# Temperature Epsilon to add or subsract 1 degree to avoid division over
# 0 in equipment not changing temperature
EpsT = 1


def min_utility(blk, heating, cooling, DTmin, eps=1e-6, DG_units=pyunits.Mwatt):
    """
    Function for Duran-Grossmann for minimization of utilities
    This function adds variables and constraints to the flowsheet to
    minimize utilities

    Args:
        blk: flowsheet
        heating: Equipment that requires Heat
        cooling: Equipment from wich heat is being removed
        DTmin: HRAT (Heat Recovery Approximation Temperature)
        eps: Epsilon for smoothing operation

    Returns:
        Constraint and variable object representing the calculation
        for hot utility and cold utility
    """

    # Generate lists out of strings from Args
    exchanger_list = heating + cooling
    pinch_streams = exchanger_list

    # Generate dictionaries out of strings to convert to pyomo objects
    pinch_streamsdict = {}
    coolingdict = {}
    heatingdict = {}
    exchangerdict = {}

    for i, el in enumerate(exchanger_list):
        pinch_streamsdict[str(el)] = el

    for i, el in enumerate(cooling):
        coolingdict[str(el)] = el

    for i, el in enumerate(heating):
        heatingdict[str(el)] = el

    for i, el in enumerate(exchanger_list):
        exchangerdict[str(el)] = el

    Q = {}

    for i in exchangerdict:
        Q[i] = pyunits.convert(pinch_streamsdict[i].heat_duty[0], to_units=DG_units)

    # Run function to fill exchanger data for pinch calculations for
    # initialization of variables
    exchangeData = heat_data(blk, heating, cooling, DG_units)

    # Call pinch calculation for initialization of variables in PD class
    PD = pinch_calc(heating, cooling, exchangeData, DTmin, eps)

    # Define dictionary for addition of DTmin to heating equipment
    dT = {}

    for i, el in enumerate(cooling):
        dT[str(el)] = 0

    for i, el in enumerate(heating):
        dT[str(el)] = DTmin

    def T_in(blk, i):
        return pinch_streamsdict[i].control_volume.properties_in[0].temperature

    blk.Tin = Expression(
        pinch_streamsdict.keys(), rule=T_in, doc="Inlet temperature in exchangers"
    )

    def T_out(blk, i):
        return pinch_streamsdict[i].control_volume.properties_out[0].temperature

    blk.Tout = Expression(
        pinch_streamsdict.keys(), rule=T_out, doc="Outlet temperature in exchangers"
    )

    # Expression for cp of equimpent with heat exchange
    def Theta_(blk, i):
        if i in heatingdict.keys():
            return Q[i] / (
                0.5
                * (
                    (blk.Tout[i] - blk.Tin[i] + (EpsT))
                    + sqrt((blk.Tout[i] - blk.Tin[i] - (EpsT)) ** 2 + eps)
                )
            )
        else:
            return Q[i] / (
                -0.5
                * (
                    (-(blk.Tout[i] - blk.Tin[i]) - (-EpsT))
                    + sqrt(((-EpsT) - (blk.Tout[i] - blk.Tin[i])) ** 2 + eps)
                )
            )

    blk.Theta = Expression(
        pinch_streamsdict.keys(), rule=Theta_, doc="FCp in exchangers"
    )

    # Define expression for pinch candidate temperature
    def T_(blk, i):
        return pinch_streamsdict[i].control_volume.properties_in[0].temperature + dT[i]

    blk.T_ = Expression(
        pinch_streamsdict.keys(), rule=T_, doc="Pinch candidate temperature"
    )

    # Define variable for heat content above the pinch point
    blk.QAh = Var(
        pinch_streamsdict.keys(),
        initialize=PD.initQAh,
        bounds=(1e-8, None),
        doc="Heat content above pinch",
        units=pyunits.get_units(Q[str(pinch_streams[0])]),
    )

    # Define a cpnstraint to calculate the varaible QAh
    def rule_heat_above_pinch(blk, p):
        return blk.QAh[p] == sum(
            blk.Theta[i]
            * (
                0.5
                * ((blk.Tin[i] - blk.T_[p]) + sqrt((blk.Tin[i] - blk.T_[p]) ** 2 + eps))
                - 0.5
                * (
                    (blk.Tout[i] - blk.T_[p])
                    + sqrt((blk.Tout[i] - blk.T_[p]) ** 2 + eps)
                )
            )
            for i in coolingdict.keys()
        )

    blk.heat_above_pinch = Constraint(
        pinch_streamsdict.keys(), rule=rule_heat_above_pinch
    )

    # Define variable for heat content below the pinch point
    blk.QAc = Var(
        pinch_streamsdict.keys(),
        initialize=PD.initQAc,
        bounds=(1e-8, None),
        doc="Heat content bellow pinch",
        units=pyunits.get_units(Q[str(pinch_streams[0])]),
    )

    # Define a constraint to calculate the varaible QAc
    def rule_heat_below_pinch(blk, p):
        return blk.QAc[p] == sum(
            blk.Theta[i]
            * (
                0.5
                * (
                    (blk.Tout[i] - blk.T_[p] + DTmin)
                    + sqrt((blk.Tout[i] - blk.T_[p] + DTmin) ** 2 + eps)
                )
                - 0.5
                * (
                    (blk.Tin[i] - blk.T_[p] + DTmin)
                    + sqrt((blk.Tin[i] - blk.T_[p] + DTmin) ** 2 + eps)
                )
            )
            for i in heatingdict.keys()
        )

    blk.heat_below_pinch = Constraint(
        pinch_streamsdict.keys(), rule=rule_heat_below_pinch
    )

    # Define variable for Heat of hot utility
    blk.Qs = Var(
        initialize=PD.initQs,
        bounds=(1e-8, None),
        doc="Heating utilities",
        units=pyunits.get_units(Q[str(pinch_streams[0])]),
    )

    # Define a constraint to solve for Qs
    # Where 1E-6 is added to both sides of the constraint as a scaling factor
    def rule_heating_utility(blk, p):
        return blk.Qs >= (blk.QAc[p] - blk.QAh[p])

    blk.heating_utility = Constraint(
        pinch_streamsdict.keys(), rule=rule_heating_utility
    )

    # Define variable for Heat of cold utility
    blk.Qw = Var(
        initialize=PD.initQw,
        bounds=(1e-8, None),
        doc="Cooling utilities",
        units=pyunits.get_units(Q[str(pinch_streams[0])]),
    )

    # Define a constraint to solve for Qw
    # Where 1E-6 is added to both sides of the constraint as a scaling factor
    def rule_cooling_utility(blk):
        return blk.Qw == -sum(Q[i] for i in exchangerdict.keys()) + blk.Qs

    blk.cooling_utility = Constraint(rule=rule_cooling_utility)


def heat_data(blk, heating, cooling, DG_units=pyunits.Mwatt):
    """
    Function for generating necesary data for Duran-Grossmann initialization
    Allows the calculation of reactors
    Args:
        blk: flowsheet
        heating: Equipment that requires Heat
        cooling: Equipment from wich heat is being removed
    Returns:
        Dictionary for the equipment exchanging heat containing:
            Tin
            Tout
            FCp
            Q
    """
    # Generate lists out of strings from Args
    exchanger_list = heating + cooling

    # Generate dictionaries out of strings to convert to pyomo objects
    pinch_streamsdict = {}
    coolingdict = {}
    heatingdict = {}
    exchangerdict = {}

    for i, el in enumerate(exchanger_list):
        pinch_streamsdict[str(el)] = el

    for i, el in enumerate(cooling):
        coolingdict[str(el)] = el

    for i, el in enumerate(heating):
        heatingdict[str(el)] = el

    for i, el in enumerate(exchanger_list):
        exchangerdict[str(el)] = el

    # Generate dictionaries for information for heat data
    Q = {}
    T_in = {}
    T_out = {}
    FCp_ = {}

    # Defining inlet temperature for inlet from equipment control volume
    for i in heatingdict.keys():
        T_in[i] = value(
            pinch_streamsdict[i].control_volume.properties_in[0].temperature
        )
        T_out[i] = (
            value(pinch_streamsdict[i].control_volume.properties_out[0].temperature)
            + EpsT
        )

    # Defining inlet temperature for utlet from equipment control volume
    for i in coolingdict.keys():
        T_in[i] = (
            value(pinch_streamsdict[i].control_volume.properties_in[0].temperature)
            + EpsT
        )
        T_out[i] = value(
            pinch_streamsdict[i].control_volume.properties_out[0].temperature
        )

    # Calculating FCp out of heat in control volume
    # Obtaines Q from equipment's control volume heat
    for i in pinch_streamsdict.keys():
        Q[i] = value(
            pyunits.convert(pinch_streamsdict[i].heat_duty[0], to_units=DG_units)
        )
        FCp_[i] = Q[i] / (T_out[i] - T_in[i])

    # Generate a large dictioary containing all the data obtained
    # from the equipment
    exchangeData = {}
    for i in pinch_streamsdict.keys():
        exchangeData[i] = {
            "T_in": T_in[i],
            "T_out": T_out[i],
            "FCp_": FCp_[i],
            "Q": Q[i],
        }

    return exchangeData


def pinch_calc(heating, cooling, exchangeData, DTmin, eps):
    """
    Function for calculating heat data for Duran-Grossmann initialization

    Args:
        heating: Equipment that requires Heat
        cooling: Equipment from wich heat is being removed
        exchangeData: Dictionary containing Tin, Tout, FCp and Q
        for each equipment in lists
        DTmin: HRAT (Heat Recovery Approximation Temperature)
        eps: Epsilon for smoothing operation
    Returns:
        PD (Dictionary containint initialized values for QAh, QAc, Qs and Qw)
    """
    # Generate lists out of strings from Args
    exchanger_list = heating + cooling

    # Generate dictionaries out of strings to convert to pyomo objects
    pinch_streamsdict = {}
    coolingdict = {}
    heatingdict = {}
    exchangerdict = {}

    for i, el in enumerate(exchanger_list):
        pinch_streamsdict[str(el)] = el

    for i, el in enumerate(cooling):
        coolingdict[str(el)] = el

    for i, el in enumerate(heating):
        heatingdict[str(el)] = el

    for i, el in enumerate(exchanger_list):
        exchangerdict[str(el)] = el

    # Generate dictionaries to contain initialized data
    T_ = {}
    initQAh = {}
    initQAc = {}
    b = []

    # Define dictionary for addition of DTmin to heating equipment
    dT = {}
    for i, el in enumerate(cooling):
        dT[str(el)] = 0

    for i, el in enumerate(heating):
        dT[str(el)] = DTmin

    # Calculate pinch temperature candidate
    # Calculate QAh and QAc
    for i in pinch_streamsdict.keys():
        T_[i] = exchangeData[i]["T_in"] + dT[i]
        initQAh[i] = sum(
            exchangeData[j]["FCp_"]
            * (
                0.5
                * (
                    (exchangeData[j]["T_in"] - T_[i])
                    + sqrt((exchangeData[j]["T_in"] - T_[i]) ** 2 + eps)
                )
                - 0.5
                * (
                    (exchangeData[j]["T_out"] - T_[i])
                    + sqrt((exchangeData[j]["T_out"] - T_[i]) ** 2 + eps)
                )
            )
            for j in coolingdict.keys()
        )
        initQAc[i] = sum(
            exchangeData[j]["FCp_"]
            * (
                0.5
                * (
                    (exchangeData[j]["T_out"] - T_[i] + DTmin)
                    + sqrt((exchangeData[j]["T_out"] - T_[i] + DTmin) ** 2 + eps)
                )
                - 0.5
                * (
                    (exchangeData[j]["T_in"] - T_[i] + DTmin)
                    + sqrt((exchangeData[j]["T_in"] - T_[i] + DTmin) ** 2 + eps)
                )
            )
            for j in heatingdict.keys()
        )

    # Generate array with all possible QS
    for i in exchangerdict.keys():
        b.append(initQAc[i] - initQAh[i])

    # Define largest value of QS
    c = max([max(b), 0.0])
    initQs = c

    # Calculate Qw from largest value of Qs
    initQw = -sum(value(exchangeData[i]["Q"]) for i in exchangerdict.keys()) + initQs
    initQw = max([initQw, 0.0])

    # Fill Class with all the data to initialioze Duran-Grossmann variables
    PD = PinchDataClass(initQs, initQw)
    PD.initQAh = initQAh
    PD.initQAc = initQAc

    return PD


def generate_curves(CD):
    """
    Function for plotting composite curves

    Args:
        CD: Class containing curve data
    Returns:
        Composite curves
    """
    # Transform information from CD class into arrays for plotting function

    # Transforms cooling equipment data
    CTin = CD.Cooling_Tin
    CTout = CD.Cooling_Tout
    CQ = CD.Cooling_Q
    CTin, CTout, CQ = np.array((CTin, CTout, CQ))

    # Transforms heating equipment data
    HTin = CD.Heating_Tin
    HTout = CD.Heating_Tout
    HQ = CD.Heating_Q
    HTin, HTout, HQ = np.array((HTin, HTout, HQ))

    Thot, Qhot = gen_curves(CTin, CTout, CQ)
    Tcold, Qcold = gen_curves(HTin, HTout, -HQ)
    # Brings Qw scalar into scope out of the class
    Qw = CD.Qw

    plt.figure(figsize=(10, 10))
    ax = plt.subplot(111)

    # Plot values for Hot streams
    # Negative sign corrects for Q < 0 for hot streams
    plt.plot(-Qhot, Thot, color="r", label="Hot Streams / Streams that are Cooled")

    # Need to shift curve to account for the cooling utility
    Qcold = Qcold + sum(HQ)
    # Plot values for cold streams
    plt.plot(
        (Qcold + value(Qw)),
        Tcold,
        color="b",
        label="Cold Streams / Streams that are Heated",
    )
    plt.xlabel(
        "Cumulative Process-Wide Heat Exchange [" + str(CD.Q_unit) + "]", fontsize=18
    )
    plt.ylabel("Temperature [" + str(CD.T_unit) + "]", fontsize=18)
    plt.title("Composite Curves for minimum utilities", size=20)
    plt.legend(loc="upper left", fontsize=14)
    plt.grid()
    plt.show()
    return


def heat_ex_data(blk, heating, cooling):
    """
    Function for turning IDAES heat exchanging equipment into a class for
    use in plotting

    Args:
        blk: Flowsheet
        heating: Equipment that heats streams
        cooling: Equipment that cools streams
    Returns:
        CD: Class with heating and cooling equipment as arrays
    """
    # Generate dictionaries out of strings to convert to pyomo objects
    exchanger_list = heating + cooling
    pinch_streams = exchanger_list

    pinch_streamsdict = {}
    coolingdict = {}
    heatingdict = {}
    exchangerdict = {}

    T_units = pyunits.get_units(
        pinch_streams[0].control_volume.properties_in[0].temperature
    )
    DG_units = pyunits.get_units(blk.Qw)

    for i, el in enumerate(exchanger_list):
        pinch_streamsdict[str(el)] = el

    for i, el in enumerate(cooling):
        coolingdict[str(el)] = el

    for i, el in enumerate(heating):
        heatingdict[str(el)] = el

    for i, el in enumerate(exchanger_list):
        exchangerdict[str(el)] = el

    # Generate zero arrays from length of the cooling list
    CHX = len(coolingdict)
    CTin = np.zeros(CHX)
    CTout = np.zeros(CHX)
    CQ = np.zeros(CHX)

    # Convert Pyomo model values into arrays
    j = 0
    for i in coolingdict.keys():
        CTin[j] = value(coolingdict[i].control_volume.properties_in[0].temperature) + 1
        CTout[j] = value(coolingdict[i].control_volume.properties_out[0].temperature)
        CQ[j] = value(
            pyunits.convert(coolingdict[i].control_volume.heat[0], to_units=DG_units)
        )
        j += 1

    # Generate zero arrays from length of the heating list
    HHX = len(heatingdict)
    HTin = np.zeros(HHX)
    HTout = np.zeros(HHX)
    HQ = np.zeros(HHX)

    # Convert Pyomo model values into arrays
    j = 0
    for i in heatingdict.keys():
        HTin[j] = value(heatingdict[i].control_volume.properties_in[0].temperature)
        HTout[j] = (
            value(heatingdict[i].control_volume.properties_out[0].temperature) + 1
        )
        HQ[j] = value(
            pyunits.convert(heatingdict[i].control_volume.heat[0], to_units=DG_units)
        )
        j += 1

        # Fill class with values and arrays
    Qw = value(blk.Qw)
    CD = CurveData(Qw, T_units, DG_units)
    CD.Cooling_Tin = CTin
    CD.Cooling_Tout = CTout
    CD.Cooling_Q = CQ
    CD.Heating_Tin = HTin
    CD.Heating_Tout = HTout
    CD.Heating_Q = HQ
    return CD


def gen_curves(Tin, Tout, Q):
    """
    Function to do add cumulative heat arrays
    Args:
        Tin: Inlet temperatures of cooling/heating equipment
        Toout: Outlet temperatures of cooling/heating equipment
        Q: Heat of cooling/heating equipment
    Returns:
        Tstar: Cumulative Temperature array
        NQstar: Cumulative heat array
    """
    # Q < 0 for hot streams = heat removing = cooling units
    # Q > 0 for cold streams = heat added = heating units

    # Ignoring edge cases for phase changes
    # Shaping Temperature array
    # Calls unique function to avoid  repeating
    Tstart = unique(np.vstack([Tin, Tout]).reshape((1, -1)))
    Tstar = np.sort(Tstart[0])

    # Generate vector of Tin size and Tunique
    nTin = len(Tin)
    nTstar = len(Tstar)
    Qmat = np.zeros((nTin, nTstar))

    # Generate cumulative arays for temperature and heat
    for i in range(nTin):
        for j in range(nTstar):
            Qmat[i, j] = linear_interpolation([Tin[i], Tout[i]], [0, Q[i]], Tstar[j])
    Qstar = sum(Qmat, 0)
    NQstar = Qstar.reshape(nTstar)
    return Tstar, NQstar


def linear_interpolation(x, y, t):
    """
    Function to do Linear interpolation with nearest neighbor extrapolation
    Args:
        x: Inlet and Outlet Temperature values
        y: 0 and Heat value
        t: Unique inlet and outlet temperatures
    Returns:
        Qstar: New array of Q Values
    """
    # Set the upper or lower value
    if x[0] < x[1]:
        lo = 0
        hi = 1
    else:
        lo = 1
        hi = 0
    # do Linear interpolation with nearest neighbor extrapolation
    alpha = (x[hi] - t) / (x[hi] - x[lo])
    alpha = max([0, alpha])
    alpha = min([alpha, 1])
    Qsta = alpha * (y[hi] - y[lo]) + y[lo]
    return Qsta


def print_HX_results(blk, exchanger_list):
    """
    Function to print results of Heat Exchangers
    Args:
        blk: flowsheet of the equipment
        exchanger_list: List of equipment to print data
    Returns:
        Printed List
    """
    # Initialize dictionary and fill with exchanger list
    exchangerdict = {}

    for i, el in enumerate(exchanger_list):
        exchangerdict[str(el)] = el

    # initialize null dictionaries for data to be printed
    Tin_ = {}
    Tout_ = {}
    f_ = {}
    Q_ = {}

    # Loop over heat exchangers
    for i in exchangerdict.keys():
        Tin_[i] = value(exchangerdict[i].control_volume.properties_in[0].temperature)
        Tout_[i] = value(exchangerdict[i].control_volume.properties_out[0].temperature)
        f_[i] = value(exchangerdict[i].control_volume.properties_out[0].flow_mol)
        Q_[i] = value(exchangerdict[i].control_volume.heat[0])
        T_units = pyunits.get_units(
            exchangerdict[i].control_volume.properties_in[0].temperature
        )
        DG_units = pyunits.get_units(exchangerdict[i].heat_duty[0])

    # Print the header
    print("Heat Exchanger Summary: ")

    # Print Inlet Temperature, Outlet Temperature and Heat
    for i in exchangerdict.keys():
        print("Heat exchanger: ", exchangerdict[i])
        print(f'Inlet T: {" "*3} {Tin_[i] : 0.3f} {T_units}')
        print(f'Outlet T: {" "*2} {Tout_[i] : 0.3f} {T_units}')
        print(f'Q : {" "*9} {Q_[i]: 0.3f} {DG_units}')


def unique(list1):
    """
    Function to remove not unique elements of a list
    Args:
        list1: List from where to obtain only unique values
    Returns:
        unique_list
    """
    # intilize a null list
    unique_list = []

    # traverse for all elements
    for i in list1:
        # check if exists in unique_list or not
        if i not in unique_list:
            unique_list.append(i)

    return unique_list


class PinchDataClass:
    """
    Class containing all the initialized values for all heat Variables
    """

    def __init__(self, initQs, initQw):
        """
        Args:
            initQs: Initial value of heat of hot utility
            initQw: Initial value of heat to be removed by cold utility
        Returns:
            None
        """
        self.initQs = initQs
        self.initQw = initQw
        self.initQAh = {}
        self.initQAc = {}

    def HeatAbove(self, data_dic):
        """
        Args:
            data_dic: Initial value of heat above pinch
        Returns:
            None
        """
        self.initQAh = data_dic

    def HeatBellow(self, data_dic2):
        """
        Args:
            data_dic2: Initial value of heat below pinch
        Returns:
            None
        """
        self.initQAc = data_dic2


class CurveData:
    """
    Class containing necessary data to generate composite curves
    """

    def __init__(self, Qw, T_unit, Q_unit):
        self.Qw = Qw
        self.T_unit = T_unit
        self.Q_unit = Q_unit

    def Cooling_Tin(self, list):
        Cooling_Tin = list

    def Cooling_Tout(self, list):
        Cooling_Tout = list

    def Cooling_Q(self, list):
        Cooling_Q = list

    def Heating_Tin(self, list):
        Heating_Tin = list

    def Heating_Tout(self, list):
        Heating_Tout = list

    def Heating_Q(self, list):
        Heating_Q = list
