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

__author__ = "John Eslick, Douglas Allan"

import enum

from pyomo.common.config import ConfigValue, In
import pyomo.environ as pyo


from idaes.core import useDefault
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.core.util.model_statistics as mstat
from idaes.core.util.constants import Constants
import idaes.logger as idaeslog


class CV_Bound(enum.Enum):
    EXTRAPOLATE = 1
    NODE_VALUE = 2


def _set_and_get_attr(obj, name, val):
    setattr(obj, name, val)
    return getattr(obj, name)


def _set_default_factor(c, s):
    """Iterate over an indexed component, and set individual scaling factors
    if there is not an existing scaling factor

    Args:
        c: indexed component to be scaled
        s: scaling factor
    """
    for i in c:
        if iscale.get_scaling_factor(c[i]) is None:
            iscale.set_scaling_factor(c[i], s)


def _set_if_unfixed(v, val):
    """Set a variable's value so long as it isn't fixed

    Args:
        v: scalar variable
        val: value to be assigned to v
    """
    if not v.fixed:
        v.set_value(pyo.value(val))


def _create_if_none(blk, var_name, idx_set, units):
    """Looks in a block's config to see whether or not a variable has been
    supplied. If it has been, creates a Reference to that variable. If it
    has not been supplied, creates a local variable.

    Args:
        blk: Pyomo block
        var_name: variable name as string
        idx_set: index set of variable var. If scalar, pass None
        units: variable units
    """
    attr = getattr(blk.config, var_name)
    if attr is None:
        if idx_set is None:
            setattr(
                blk,
                var_name,
                pyo.Var(
                    initialize=0,
                    units=units,
                ),
            )
        else:
            setattr(
                blk,
                var_name,
                pyo.Var(
                    *idx_set,
                    initialize=0,
                    units=units,
                ),
            )
    else:
        setattr(blk, var_name, pyo.Reference(attr))


def _init_solve_block(blk, solver, log):
    """Checks whether solving a block is a square problem. If not, raise
    InitializationError. If it is, solve block while redirecting output to
    IDAES log. Then check whether solve was successful. If not successful,
    raise InitializationError

    Args:
        blk: Pyomo block to solve
        solver: solver object (NOT name string)
        log: IDAES solver log
    """
    if not mstat.degrees_of_freedom(blk) == 0:
        raise InitializationError(
            f"{blk.name} encountered a nonsquare problem during "
            "initialization. Check whether all parameter values are fixed and "
            "keep in mind that the SOC submodel initialization methods assume "
            "that certain degrees of freedom have been fixed by the cell "
            "initialization method."
        )
    with idaeslog.solver_log(log, idaeslog.DEBUG) as slc:
        results = solver.solve(blk, tee=slc.tee)
    if not pyo.check_optimal_termination(results):
        raise InitializationError(
            f"{blk.name} failed to initialize successfully. Please check "
            f"the output logs for more information."
        )


def _interpolate_channel(iz, ifaces, nodes, faces, phi_func, phi_inlet, opposite_flow):
    """PRIVATE Function: Interpolate faces of control volumes in 1D

    Args:
        iz: index of the face
        ifaces: set of face indexes
        nodes: set of node z locations
        faces: set of face z locations
        phi_func: function that returns an expression for the quantity to be
            interpolated at the node as a function of node index
        phi_inlet: expression for the value of the quantity to be interpolated
            at the channel inlet face (not the inlet to the CV)
        method: interpolation method
        opposite_flow: if True assume velocity is negative

    Returns:
        expression for phi at face[iz]
    """
    # I don't always need these, but it doesn't take long to calculate them
    if not opposite_flow:
        izu = iz - 1  # adjacent node upstream of the face
        izuu = iz - 2  # node upstream adacjent to node upstream adjacent to face
        izd = iz  # downstream node adjacent to face
    else:
        izu = iz  # adjacent node upstream of the face
        izuu = iz + 1  # node upstream adacjent to node upstream adjacent to face
        izd = iz - 1  # downstream node adjacent to face
    if iz == ifaces.first() and not opposite_flow:
        return phi_inlet
    if iz == ifaces.last() and opposite_flow:
        return phi_inlet
    return phi_func(izu)


def _interpolate_2D(
    ic,
    ifaces,
    nodes,
    faces,
    phi_func,
    phi_bound_0,
    phi_bound_1,
    derivative=False,
):
    """PRIVATE Function: Interpolate faces of control volumes in 2D

    Args:
        ic: face index in x or z direction matching specified direction
        ifaces: set of face indexes
        nodes: set of node locations in specified direction
        faces: set of face locations in specified direction
        phi_func: function that returns an expression for the quantity to be
            interpolated at a node center as a function of node index
        phi_bound_0: expression for the value of the quantity to be interpolated
            at the 0 bound
        phi_bound_1: expression for the value of the quantity to be interpolated
            at the 1 bound
        derivative: If True estimate derivative
        method: interpolation method currently only CDS is supported

    Returns:
        expression for phi at face
    """
    icu = ic - 1
    icd = ic
    if ic == ifaces.first():
        if not isinstance(phi_bound_0, CV_Bound):
            return phi_bound_0
        icu = ic + 1
    if ic == ifaces.last():
        if not isinstance(phi_bound_1, CV_Bound):
            return phi_bound_1
        icd = ic - 2
    # For now CDS is the only option
    cu = nodes.at(icu)
    cd = nodes.at(icd)
    if not derivative:
        cf = faces.at(ic)
        lambf = (cd - cf) / (cd - cu)
        return (1 - lambf) * phi_func(icu) + lambf * phi_func(icd)
    else:
        # Since we are doing linear interpolation derivative is the slope
        # between node centers even if they are not evenly spaced
        return (phi_func(icd) - phi_func(icu)) / (cd - cu)


class _SubsetOf(object):
    # Copied and modified from Pyomo's implementation of In
    #  ___________________________________________________________________________
    #
    #  Pyomo: Python Optimization Modeling Objects
    #  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
    #  Under the terms of Contract DE-NA0003525 with National Technology and
    #  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain
    #  rights in this software.
    #  This software is distributed under the 3-clause BSD License.
    #  ___________________________________________________________________________
    #
    #  This module was originally developed as part of the PyUtilib project
    #  Copyright (c) 2008 Sandia Corporation.
    #  This software is distributed under the BSD License.
    #  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
    #  the U.S. Government retains certain rights in this software.
    #  ___________________________________________________________________________
    """_SubsetOf(domain, cast=None)
    Domain validation class admitting a Container of possible values

    This will admit any iterable containing values that are in the
     `domain` Container (i.e., Container.__contains__() returns True).
    Most common domains are list, set, and dict objects.

    Parameters
    ----------
    domain: Container
        The container that specifies the allowable values.  Incoming
        values are passed to ``domain.__contains__()``, and if ``True``
        is returned, the value is accepted and returned.
    """
    #  TODO Need to determine what to do about repeated entries before this can become a Pyomo PR
    def __new__(cls, domain=None, cast=None):
        return super(_SubsetOf, cls).__new__(cls)

    def __init__(self, domain):
        self._domain = domain

    def __call__(self, possible_subset):
        # We often like to use None to represent the empty set in Config blocks because using a mutable object
        # for a default function value results in the same object being used every time the function is called
        # Not sure the same is true for Config blocks, but it probably is.
        if possible_subset is None:
            return None
        for value in possible_subset:
            if value not in self._domain:
                raise ValueError("value %s not in domain %s" % (value, self._domain))
        return possible_subset


_element_list = ["Ar", "H", "O", "C", "N", "S"]
_gas_species_list = [
    "Ar",
    "CH4",
    "CO",
    "CO2",
    "C2H4",
    "C2H6",
    "C3H8",
    "H2",
    "H2O",
    "H2S",
    "NO",
    "N2",
    "N2O",
    "O2",
    "SO2",
]
_all_species_list = _gas_species_list.copy()
_all_species_list.extend(["Vac", "O^2-", "e^-"])
_element_dict = {
    "Ar": {"Ar": 1},
    "H": {"CH4": 4, "C2H4": 4, "C2H6": 6, "C3H8": 8, "H2": 2, "H2O": 2, "H2S": 2},
    "O": {"CO": 1, "CO2": 2, "H2O": 1, "NO": 1, "N2O": 1, "O2": 2, "SO2": 2},
    "C": {"CH4": 1, "CO": 1, "CO2": 1, "C2H4": 4, "C2H6": 2, "C3H8": 3},
    "N": {"NO": 1, "N2": 2, "N2O": 2},
    "S": {"H2S": 1, "SO2": 1},
}

# If a molecule is not listed under an element, it has zero atoms of that element
for value in _element_dict.values():
    for species in _gas_species_list:
        if species not in value.keys():
            value[species] = 0

# Parameters for binary diffusion coefficients from:
#  "Properties of Gases and Liquids" 5th Ed.
bin_diff_sigma = {
    "Ar": 3.542,
    "Air": 3.711,
    "CH4": 3.758,
    "CO": 3.690,
    "CO2": 3.941,
    "C2H4": 4.163,
    "C2H6": 4.443,
    "C3H8": 5.118,
    "H2": 2.827,
    "H2O": 2.641,
    "H2S": 3.623,
    "NO": 3.492,
    "N2": 3.798,
    "N2O": 3.826,
    "O2": 3.467,
    "SO2": 4.112,
}
bin_diff_epsok = {
    "Ar": 93.3,
    "Air": 78.6,
    "CH4": 148.6,
    "CO": 91.7,
    "CO2": 195.2,
    "C2H4": 224.7,
    "C2H6": 215.7,
    "C3H8": 237.1,
    "H2": 59.7,
    "H2O": 809.1,
    "H2S": 301.1,
    "NO": 116.7,
    "N2": 71.4,
    "N2O": 232.4,
    "O2": 106.7,
    "SO2": 335.4,
}
bin_diff_M = {
    "Ar": 39.948,
    "Air": 28.96,
    "CH4": 16.0425,
    "CO": 28.0101,
    "CO2": 44.0095,
    "C2H4": 28.0532,
    "C2H6": 30.0690,
    "C3H8": 44.0956,
    "H2": 2.01588,
    "H2O": 18.0153,
    "H2S": 34.081,
    "NO": 30.0061,
    "N2": 28.0134,
    "N2O": 44.0128,
    "O2": 31.9988,
    "SO2": 64.064,
}

# Shomate equation parameters
# NIST Webbook
h_params = {
    "Ar": {
        "A": 20.78600,
        "B": 2.825911e-7,
        "C": -1.464191e-7,
        "D": 1.092131e-8,
        "E": -3.661371e-8,
        "F": -6.197350,
        "G": 179.9990,
        "H": 0.000000,
    },
    "CH4": {
        "A": -0.703029,
        "B": 108.4773,
        "C": -42.52157,
        "D": 5.862788,
        "E": 0.678565,
        "F": -76.84376,
        "G": 158.7163,
        "H": -74.87310,
    },
    "CO": {
        "A": 25.56759,
        "B": 6.096130,
        "C": 4.054656,
        "D": -2.671301,
        "E": 0.131021,
        "F": -118.0089,
        "G": 227.3665,
        "H": -110.5271,
    },
    "CO2": {
        "A": 24.99735,
        "B": 55.18696,
        "C": -33.69137,
        "D": 7.948387,
        "E": -0.136638,
        "F": -403.6075,
        "G": 228.2431,
        "H": -393.5224,
    },
    "C2H4": {
        "A": -6.387880,
        "B": 184.4019,
        "C": -112.9718,
        "D": 28.49593,
        "E": 0.315540,
        "F": 48.17332,
        "G": 163.1568,
        "H": 52.46694,
    },
    "H2": {  # 1000K-2500K
        "A": 33.066178,
        "B": -11.363417,
        "C": 11.432816,
        "D": -2.772874,
        "E": -0.158558,
        "F": -9.980797,
        "G": 172.707974,
        "H": 0.0,
    },
    "H2O": {  # 500K-1700K
        "A": 30.09200,
        "B": 6.832514,
        "C": 6.793435,
        "D": -2.534480,
        "E": 0.082139,
        "F": -250.8810,
        "G": 223.3967,
        "H": -241.8264,
    },
    "H2S": {
        "A": 26.88412,
        "B": 18.67809,
        "C": 3.434203,
        "D": -3.378702,
        "E": 0.135882,
        "F": -28.91211,
        "G": 233.3747,
        "H": -20.50202,
    },
    "N2": {
        "A": 19.50583,
        "B": 19.88705,
        "C": -8.598535,
        "D": 1.369784,
        "E": 0.527601,
        "F": -4.935202,
        "G": 212.3900,
        "H": 0.0,
    },
    "O2": {  # 700 to 2000
        "A": 30.03235,
        "B": 8.772972,
        "C": -3.988133,
        "D": 0.788313,
        "E": -0.741599,
        "F": -11.32468,
        "G": 236.1663,
        "H": 0.0,
    },
    "SO2": {
        "A": 21.43049,
        "B": 74.35094,
        "C": -57.75217,
        "D": 16.35534,
        "E": 0.086731,
        "F": -305.7688,
        "G": 254.8872,
        "H": -296.8422,
    },
    "Vac": {
        "A": 0.0,
        "B": 0.0,
        "C": 0.0,
        "D": 0.0,
        "E": 0.0,
        "F": 0.0,
        "G": 0.0,
        "H": 0.0,
    },
    "O^2-": {  # Chosen to match O2 coeffs besides enthalpy and entropy of formation
        "A": 0.5 * 30.03235,
        "B": 0.5 * 8.772972,
        "C": 0.5 * -3.988133,
        "D": 0.5 * 0.788313,
        "E": 0.5 * -0.741599,
        "F": 0.5 * -11.32468,
        "G": 0.5 * 236.1663,  # 0,
        "H": 0,  # -236.0,
    },
    "e^-": {
        "A": 0.0,
        "B": 0.0,
        "C": 0.0,
        "D": 0.0,
        "E": 0.0,
        "F": 0.0,
        "G": 0.0,
        "H": 0.0,
    },
}


def _binary_diffusion_coefficient_expr(temperature, p, c1, c2):
    mab = 2 * (1.0 / bin_diff_M[c1] + 1.0 / bin_diff_M[c2]) ** (-1)
    sab = (bin_diff_sigma[c1] + bin_diff_sigma[c2]) / 2.0
    epsok = (bin_diff_epsok[c1] * bin_diff_epsok[c2]) ** 0.5
    tr = temperature / (epsok * pyo.units.K)
    a = 1.06036
    b = 0.15610
    c = 0.19300
    d = 0.47635
    e = 1.03587
    f = 1.52996
    g = 1.76474
    h = 3.89411
    omega = (
        a / tr**b + c / pyo.exp(d * tr) + e / pyo.exp(f * tr) + g / pyo.exp(h * tr)
    )
    cm2_to_m2 = 0.01 * 0.01
    Pa_to_bar = 1e-5
    return (
        (
            0.00266
            * cm2_to_m2
            * (temperature / pyo.units.K) ** (3 / 2)
            / (p / pyo.units.Pa)
            / Pa_to_bar
            / mab**0.5
            / sab**2
            / omega
        )
        * pyo.units.m**2
        / pyo.units.s
    )


def _pure_component_cp(temperature, comp):
    d = h_params[comp]
    t = temperature / 1000.0 / pyo.units.K
    return (
        (d["A"] + d["B"] * t + d["C"] * t**2 + d["D"] * t**3 + d["E"] / t**2)
        * pyo.units.J
        / pyo.units.mol
        / pyo.units.K
    )


def _pure_component_cv(temperature, comp):
    return _pure_component_cp(temperature, comp) - Constants.gas_constant


def _pure_component_visc(temperature, comp):
    # The properties of gases and liquids 5th ed.  eqn 9-4.6
    # low pressure gas from theory
    return (
        (
            1e-7
            * 16.64
            * bin_diff_M[comp] ** 0.5
            * temperature
            / (bin_diff_epsok[comp] ** 0.5 * bin_diff_sigma[comp] ** 2)
        )
        * pyo.units.Pa
        * pyo.units.s
    )


def _pure_component_thermal_comductivity(temperature, c):
    # The properties of gases and liquids 5th ed.  eqn 10-3.6
    # low pressure gas Stiel and Thodos 1964.
    cv = _pure_component_cv(temperature, c)
    visc = _pure_component_visc(temperature, c)
    m = bin_diff_M[c] / 1000
    units = pyo.units.W / pyo.units.m / pyo.units.K
    return (1.15 + Constants.gas_constant * 2.03 / cv) * visc * cv / m * units


def _comp_enthalpy_expr(temperature, comp):
    # ideal gas enthalpy
    d = h_params[comp]
    t = temperature / 1000.0 / pyo.units.K
    return (
        1000
        * (
            d["A"] * t
            + d["B"] * t**2 / 2.0
            + d["C"] * t**3 / 3.0
            + d["D"] * t**4 / 4.0
            - d["E"] / t
            + d["F"]
        )
        * pyo.units.J
        / pyo.units.mol
    )


_monotomic_gas_standard_state = ["He", "Ne", "Ar", "Kr", "Xe", "Ra"]
_diatomic_gas_standard_state = ["F", "Cl", "H", "N", "O"]


def _comp_int_energy_expr(temperature, comp):
    # ideal gas internal energy
    # NIST has 298.15 K as a reference state, so adjust internal energy expression for that
    T_ref = 298.15
    dn_form = 1
    for element, molecule_dict in _element_dict.items():
        if element in _monotomic_gas_standard_state:
            dn_form -= molecule_dict[comp]
        elif element in _diatomic_gas_standard_state:
            dn_form -= 0.5 * molecule_dict[comp]
    return (
        _comp_enthalpy_expr(temperature, comp)
        - Constants.gas_constant * (temperature - T_ref)
        + dn_form * Constants.gas_constant * T_ref
    )


def _comp_entropy_expr(temperature, comp):
    # ideal gas entropy
    d = h_params[comp]
    t = temperature / 1000.0 / pyo.units.K
    return (
        (
            d["A"] * pyo.log(t)
            + d["B"] * t
            + d["C"] * t**2 / 2.0
            + d["D"] * t**3 / 3.0
            - d["E"] / 2.0 / t**2
            + d["G"]
        )
        * pyo.units.J
        / pyo.units.mol
        / pyo.units.K
    )


def _submodel_boilerplate_config(CONFIG):
    # Add common fields to submodel CONFIG dict
    CONFIG.declare(
        "control_volume_zfaces",
        ConfigValue(
            description="List containing coordinates of control volume faces "
            "in z direction. Coordinates must start with zero, be strictly "
            "increasing, and end with one"
        ),
    )
    CONFIG.declare(
        "length_z",
        ConfigValue(
            default=None, description="Length in the direction of flow (z-direction)"
        ),
    )
    CONFIG.declare(
        "length_y", ConfigValue(default=None, description="Width of cell (y-direction)")
    )
    CONFIG.declare(
        "current_density",
        ConfigValue(default=None, description="Optional current_density variable"),
    )
    CONFIG.declare(
        "include_temperature_x_thermo",
        ConfigValue(
            domain=In([useDefault, True, False]),
            default=True,
            description="Whether to consider temperature variations in "
            "x direction in thermodynamic equations",
        ),
    )


def _submodel_boilerplate_create_if_none(unit):
    # Check for certain common variables on unit block. If they were provided
    # in CONFIG, create references, otherwise create them new on unit block
    tset = unit.flowsheet().config.time
    iznodes = unit.iznodes
    _create_if_none(unit, "length_z", idx_set=None, units=pyo.units.m)
    _create_if_none(unit, "length_y", idx_set=None, units=pyo.units.m)
    _create_if_none(
        unit,
        "current_density",
        idx_set=(tset, iznodes),
        units=pyo.units.A / pyo.units.m**2,
    )


def _thermal_boundary_conditions_config(CONFIG, thin):
    # Add fieldnames relating to thermal boundary conditions to submodel
    # CONFIG dict. If submodel is thin, there is only a single temperature,
    # otherwise there is one at both ends
    CONFIG.declare(
        "temperature_z",
        ConfigValue(
            default=None, description="Temperature as indexed by time " "and z"
        ),
    )
    if thin:
        CONFIG.declare(
            "temperature_deviation_x",
            ConfigValue(
                default=None,
                description="Deviation of temperature " "from temperature_z",
            ),
        )
    else:
        CONFIG.declare(
            "temperature_deviation_x0",
            ConfigValue(
                default=None,
                description="Deviation of temperature at x=0 " "from temperature_z",
            ),
        )
        CONFIG.declare(
            "temperature_deviation_x1",
            ConfigValue(
                default=None,
                description="Deviation of temperature at x=1 " "from temperature_z",
            ),
        )
    CONFIG.declare(
        "heat_flux_x0",
        ConfigValue(
            default=None, description="Heat flux through x=0 " "(positive is in)"
        ),
    )
    CONFIG.declare(
        "heat_flux_x1",
        ConfigValue(
            default=None, description="Heat flux through x=1 " "(positive is out)"
        ),
    )


def _create_thermal_boundary_conditions_if_none(unit, thin):
    # Check for variables relating to thermal boundary conditions on unit
    # block. If they were provided in CONFIG, create references, otherwise
    # create them new on unit block
    tset = unit.flowsheet().config.time
    include_temp_x_thermo = unit.config.include_temperature_x_thermo
    iznodes = unit.iznodes

    _create_if_none(unit, "temperature_z", idx_set=(tset, iznodes), units=pyo.units.K)

    if thin:
        _create_if_none(
            unit, "temperature_deviation_x", idx_set=(tset, iznodes), units=pyo.units.K
        )

        @unit.Expression(tset, iznodes)
        def temperature(b, t, iz):
            if include_temp_x_thermo:
                return b.temperature_z[t, iz] + b.temperature_deviation_x[t, iz]
            else:
                return b.temperature_z[t, iz]

    else:
        _create_if_none(
            unit, "temperature_deviation_x0", idx_set=(tset, iznodes), units=pyo.units.K
        )

        @unit.Expression(tset, iznodes)
        def temperature_x0(b, t, iz):
            if include_temp_x_thermo:
                return b.temperature_z[t, iz] + b.temperature_deviation_x0[t, iz]
            else:
                return b.temperature_z[t, iz]

        _create_if_none(
            unit, "temperature_deviation_x1", idx_set=(tset, iznodes), units=pyo.units.K
        )

        @unit.Expression(tset, iznodes)
        def temperature_x1(b, t, iz):
            if include_temp_x_thermo:
                return b.temperature_z[t, iz] + b.temperature_deviation_x1[t, iz]
            else:
                return b.temperature_z[t, iz]

    _create_if_none(
        unit,
        "heat_flux_x0",
        idx_set=(tset, iznodes),
        units=pyo.units.W / pyo.units.m**2,
    )

    _create_if_none(
        unit,
        "heat_flux_x1",
        idx_set=(tset, iznodes),
        units=pyo.units.W / pyo.units.m**2,
    )


def _material_boundary_conditions_config(CONFIG, thin):
    # Add fieldnames relating to material boundary conditions to submodel
    # CONFIG dict. If submodel is thin, there is only a single concentration,
    # otherwise either end of unit have different concentrations
    if thin:
        CONFIG.declare(
            "material_flux_x",
            ConfigValue(default=None, description="Variable for material flux"),
        )
        CONFIG.declare(
            "conc_mol_comp_deviation_x",
            ConfigValue(
                default=None,
                description="Deviation of concentration from channel bulk "
                "concentration",
            ),
        )
    else:
        CONFIG.declare(
            "conc_mol_comp_deviation_x0",
            ConfigValue(
                default=None,
                description="Deviation of concentration at x= 0 from "
                "channel bulk concentration",
            ),
        )
        CONFIG.declare(
            "conc_mol_comp_deviation_x1",
            ConfigValue(
                default=None,
                description="Deviation of concentration at x= 0 from "
                "channel bulk concentration",
            ),
        )
        CONFIG.declare(
            "material_flux_x0",
            ConfigValue(default=None, description="Variable for material flux at x=0"),
        )
        CONFIG.declare(
            "material_flux_x1",
            ConfigValue(default=None, description="Variable for material flux at x=1"),
        )


def _create_material_boundary_conditions_if_none(unit, thin):
    # Check for variables relating to material boundary conditions on unit
    # block. If they were provided in CONFIG, create references, otherwise
    # create them new on unit block
    tset = unit.flowsheet().config.time
    iznodes = unit.iznodes
    comps = unit.component_list
    if thin:
        _create_if_none(
            unit,
            "conc_mol_comp_deviation_x",
            idx_set=(tset, iznodes, comps),
            units=pyo.units.mol / pyo.units.m**3,
        )
        _create_if_none(
            unit,
            "material_flux_x",
            idx_set=(tset, iznodes, comps),
            units=pyo.units.mol / (pyo.units.s * pyo.units.m**2),
        )
    else:
        _create_if_none(
            unit,
            "conc_mol_comp_deviation_x0",
            idx_set=(tset, iznodes, comps),
            units=pyo.units.mol / pyo.units.m**3,
        )
        _create_if_none(
            unit,
            "conc_mol_comp_deviation_x1",
            idx_set=(tset, iznodes, comps),
            units=pyo.units.mol / pyo.units.m**3,
        )
        _create_if_none(
            unit,
            "material_flux_x0",
            idx_set=(tset, iznodes, comps),
            units=pyo.units.mol / (pyo.units.s * pyo.units.m**2),
        )
        _create_if_none(
            unit,
            "material_flux_x1",
            idx_set=(tset, iznodes, comps),
            units=pyo.units.mol / (pyo.units.s * pyo.units.m**2),
        )


def _face_initializer(blk, faces, direction):
    dfaces = direction + "faces"
    dnodes = direction + "nodes"
    if faces[0] != 0:
        raise ConfigurationError(
            f"Smallest control volume face provided in "
            f"{dfaces} to block {blk.name} is not zero."
        )
    if faces[-1] != 1:
        raise ConfigurationError(
            f"Largest control volume face provided in "
            f"{dfaces} to block {blk.name} is not one."
        )
    for i in range(len(faces) - 1):
        if not faces[i + 1] - faces[i] > 0:
            raise ConfigurationError(
                f"Sequence of control volume face distances in {dfaces} "
                "is not strictly increasing."
            )

    face_set = _set_and_get_attr(
        blk,
        dfaces,
        pyo.Set(
            initialize=faces,
            ordered=True,
            doc=f"{direction} coordinates for control volume faces",
        ),
    )
    node_set = _set_and_get_attr(
        blk,
        dnodes,
        pyo.Set(
            initialize=[
                (face_set.at(i) + face_set.at(i + 1)) / 2.0
                for i in range(1, len(faces))
            ],
            ordered=True,
            doc=f"{direction} coordinates for control volume centers",
        ),
    )
    # This sets provide an integer index for nodes and faces
    iface_set = _set_and_get_attr(
        blk,
        "i" + dfaces,
        pyo.Set(
            initialize=range(1, len(face_set) + 1),
            ordered=True,
            doc=f"Integer index set for {dfaces}",
        ),
    )
    inode_set = _set_and_get_attr(
        blk,
        "i" + dnodes,
        pyo.Set(
            initialize=range(1, len(node_set) + 1),
            ordered=True,
            doc=f"Integer index set for {dnodes}",
        ),
    )
    return iface_set, inode_set
