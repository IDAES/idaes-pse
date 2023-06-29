#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""This provides functions to create and check Helmholtz equation of state
parameter and expression files.
"""

import logging
import json
import numpy
import pyomo.environ as pyo

from idaes.models.properties.general_helmholtz.expressions import (
    phi_residual_types,
    phi_ideal_types,
    delta_sat_types,
    surface_tension_types,
)

_log = logging.getLogger("idaes.helmholtz_parameters")


def _parse_int_key(pairs):
    """Since json can only store dictionary keys as strings this converts the string
    keys to integers or integer tuples if possible
    """
    d = {}
    for x in pairs:
        try:
            if x[0][0] == "(" and x[0][-1] == ")":
                d[tuple(map(int, map(str.strip, x[0][1:-1].split(","))))] = x[1]
            else:
                d[int(x[0])] = x[1]
        except ValueError:
            d[x[0]] = x[1]
    return d


class WriteParameters(object):
    """This class allows users to set parameters and expressions to generate parameter
    and NL expression files that define a Helmholtz equation of state and optionally
    transport properties for a substance.
    """

    # Dictionary used to track the location of variables in NL files.
    variables = {
        "delta": 0,
        "tau": 1,
        "p": 2,
        "T": 3,
    }

    # Dictionary used to track the location of expressions in NL files.
    expressions = {
        # phi is dimensionless Helmholtz free energy
        "phii": 0,  # ideal part of phi(delta, tau)
        "phii_d": 1,  # ideal part of phi partial wrt delta
        "phii_dd": 2,  # ideal part of phi partial wrt delta and delta
        "phii_t": 3,  # ideal part of phi partial wrt tau
        "phii_tt": 4,  # ideal part of phi partial tau delta and tau
        "phii_dt": 5,  # ideal part of phi partial wrt delta and tau
        "phir": 6,  # residual part of phi(delta, tau)
        "phir_d": 7,  # residual part of phi partial wrt delta
        "phir_dd": 8,  # residual part of phi partial wrt delta and delta
        "phir_t": 9,  # residual part of phi partial wrt tau
        "phir_tt": 10,  # residual part of phi partial wrt tau and tau
        "phir_dt": 11,  # residual part of phi partial wrt delta and tau
        # Initial guess on saturated liquid and vapor reduced density
        "delta_v_sat_approx": 12,  # approximate delta_sat_vapor(tau)
        "delta_l_sat_approx": 13,  # approximate delta_sat_liquid(tau)
    }

    # A dictionary of optional expressions for transport properties, the value
    # in the dictionary is a short name used in the NL file names.
    optional_expressions = {
        "thermal_conductivity": "tcx",
        "surface_tension": "st",
        "viscosity": "visc",
    }

    def __init__(
        self,
        parameters,
    ):
        """Create a parameter writer object.  Parameters can either be specified in
        dictionary or read from a json file.

        Args:
            parameters (str|dict):  Either a dictionary of parameters or a path to
                a json file.

        Returns:
            WriteParameters: Object used to write Helmholtz EoS parameter and expression files
        """
        # A list that will store the names of defined expressions
        self.has_expression = []

        # if parameters is a string read from json
        if isinstance(parameters, str):
            with open(parameters, "r") as fp:
                parameters = json.load(fp, object_pairs_hook=_parse_int_key)
        # we may need to access parameters that are provided but not stored
        # directly as attributes, so store a link to the parameter dictionary
        self.parameters = parameters
        # Store the basic parameters as attributes.
        self.comp = parameters["comp"]
        self.R = parameters["basic"]["R"]
        self.MW = parameters["basic"]["MW"]
        self.T_star = parameters["basic"]["T_star"]
        self.rho_star = parameters["basic"]["rho_star"]
        self.Tc = parameters["basic"]["Tc"]
        self.rhoc = parameters["basic"]["rhoc"]
        self.Pc = parameters["basic"]["Pc"]
        self.Tt = parameters["basic"]["Tt"]
        self.Pt = parameters["basic"]["Pt"]
        self.rhot_l = parameters["basic"]["rhot_l"]
        self.rhot_v = parameters["basic"]["rhot_v"]
        self.P_min = parameters["basic"]["P_min"]
        self.P_max = parameters["basic"]["P_max"]
        self.rho_max = parameters["basic"]["rho_max"]
        self.T_min = parameters["basic"]["T_min"]
        self.T_max = parameters["basic"]["T_max"]
        # If provided, set the reference state offset, otherwise set to (0, 0).
        try:
            self.reference_state_offset = parameters["eos"]["reference_state_offset"]
        except KeyError:
            self.reference_state_offset = [0.0, 0.0]

        # Create the main Pyomo model to write equation of state expressions
        self.model = self.make_model("delta", "tau")

        # Optional models (Thermal Conductivity, Viscosity, and Surface Tension)
        self.model_tcx = self.make_model("delta", "tau")
        self.model_visc = self.make_model("delta", "tau")

        # Since surface tension is only valid for the two phase region, surface tension
        # is a function of temperature only
        self.model_st = self.make_model("tau")

        # Check if a predefined form of the ideal part of Helmholtz free energy is used
        try:
            phi_ideal_type = parameters["eos"]["phi_ideal_type"]
            if phi_ideal_type:
                self.add(
                    phi_ideal_types[phi_ideal_type](
                        model=self.model, parameters=parameters
                    )
                )
        except KeyError:
            phi_ideal_type = 0

        # Check if a predefined form of the residual part of Helmholtz free energy is used
        try:
            phi_residual_type = parameters["eos"]["phi_residual_type"]
            if phi_residual_type:
                self.add(
                    phi_residual_types[phi_residual_type](
                        model=self.model, parameters=parameters
                    )
                )
        except KeyError:
            phi_residual_type = 0

        # Check if predefined approximate liquid and vapor saturated density curves are used
        aux_parameters = parameters.get("aux", None)
        if aux_parameters is not None:
            delta_l_sat_parameters = parameters["aux"].get("delta_l_sat_approx", None)
            delta_v_sat_parameters = parameters["aux"].get("delta_v_sat_approx", None)
            if delta_l_sat_parameters is not None:
                etype = delta_l_sat_parameters.get("type", 0)
                if etype:
                    self.add(
                        {
                            "delta_l_sat_approx": delta_sat_types[etype](
                                model=self.model,
                                name="delta_l_sat_approx",
                                parameters=parameters,
                            ),
                        }
                    )
            if delta_v_sat_parameters is not None:
                etype = delta_v_sat_parameters.get("type", 0)
                if etype:
                    self.add(
                        {
                            "delta_v_sat_approx": delta_sat_types[etype](
                                model=self.model,
                                name="delta_v_sat_approx",
                                parameters=parameters,
                            ),
                        }
                    )

        # Check if using predefined surface tension expression
        try:
            etype = parameters["transport"]["surface_tension"]["type"]
            if etype:
                self.add(
                    {
                        "surface_tension": surface_tension_types[etype](
                            model=self.model_st, parameters=parameters
                        ),
                    }
                )
        except KeyError:  # No surface tension to add
            pass

    def calculate_pressure(self, rho, T):
        """From the expressions provided, calculate pressure from density and
        temperature. This can be used for testing and calculating the critical
        pressure based on the critical temperature and density.

        Args:
            rho (float): density in kg/m3
            T (float): temperature in K

        Returns:
            float: pressure in kPa
        """
        self.model.delta = rho / self.rho_star
        self.model.tau = self.T_star / T
        return pyo.value(rho * self.R * T * (1 + self.model.delta * self.model.phir_d))

    def calculate_enthalpy(self, rho, T):
        """From the expressions provided, calculate enthalpy from density and
        temperature. This can be used for testing.

        Args:
            rho (float): density in kg/m3
            T (float): temperature in K

        Returns:
            float: enthalpy in kJ/kg
        """
        self.model.delta = rho / self.rho_star
        self.model.tau = self.T_star / T
        return pyo.value(
            self.R
            * T
            * (
                1
                + self.model.tau * (self.model.phii_t + self.model.phir_t)
                + self.model.delta * self.model.phir_d
            )
        )

    def calculate_entropy(self, rho, T):
        """From the expressions provided, calculate entropy from density and
        temperature. This can be used for testing.

        Args:
            rho (float): density in kg/m3
            T (float): temperature in K

        Returns:
            float: entropy in kJ/kg/K
        """
        self.model.delta = rho / self.rho_star
        self.model.tau = self.T_star / T
        return pyo.value(
            self.R
            * (
                self.model.tau * (self.model.phii_t + self.model.phir_t)
                - self.model.phii
                - self.model.phir
            )
        )

    def make_model(self, *args):
        """Make a Pyomo model used to define expression NL files.  The arguments are strings
        for variables to create.  The basic parameters are also added as parameters in the
        model.

        Args:
            (str): Variables to add to the model

        Returns:
            ConcreteModel: Pyomo model with variables from args
        """
        m = pyo.ConcreteModel()
        for a in args:
            setattr(m, a, pyo.Var())
        m.R = pyo.Param(initialize=self.R)
        m.MW = pyo.Param(initialize=self.MW)
        m.T_star = pyo.Param(initialize=self.T_star)
        m.rho_star = pyo.Param(initialize=self.rho_star)
        m.Tc = pyo.Param(initialize=self.Tc)
        m.rhoc = pyo.Param(initialize=self.rhoc)
        m.Pc = pyo.Param(initialize=self.Pc)
        m.Tt = pyo.Param(initialize=self.Tt)
        m.Pt = pyo.Param(initialize=self.Pt)
        m.rhot_l = pyo.Param(initialize=self.rhot_l)
        m.rhot_v = pyo.Param(initialize=self.rhot_v)
        m.P_min = pyo.Param(initialize=self.P_min)
        m.P_max = pyo.Param(initialize=self.P_max)
        m.rho_max = pyo.Param(initialize=self.rho_max)
        m.T_min = pyo.Param(initialize=self.T_min)
        m.T_max = pyo.Param(initialize=self.T_max)
        return m

    def add(self, expressions):
        """This adds expressions to the to the object.  These expressions are written
        to NL files to be used by external functions.

        Args:
            expressions (dict): Dictionary where the key is an expression name and the
                value is a Pyomo expression.

        Returns:
            None
        """
        for name, expr in expressions.items():
            if name in self.expressions:  # check if in thermo model
                m = self.model
            elif name == "thermal_conductivity":
                m = self.model_tcx
            elif name == "surface_tension":
                m = self.model_st
            elif name == "viscosity":
                m = self.model_visc
            else:
                raise RuntimeError(f"Unknown expression {name}")
            if callable(expr):
                setattr(m, name, pyo.Objective(rule=expr))
            else:
                setattr(m, name, pyo.Objective(expr=expr))
            self.has_expression.append(name)

    def approx_sat_curves(self, trange):
        """_log.infos a table to verify that the approximate saturated density curves
        are correct.  Since the approximate curves are used as an initial guess to
        the phase equilibrium problems, they are not directly testable.

        Args:
            trange (iterable): temperature points in K

        Returns:
            tuple: list of liquid densities and list of vapor densities
        """
        _log.info(
            "\n====================================================================="
        )
        _log.info(" Check approx sat delta curves")
        _log.info(
            "====================================================================="
        )
        _log.info(
            f"{'T [K]':7s}  {'T [C]':8s}  {'~rho_l [kg/m3]':14s}  {'~rho_v [kg/m3]':14s}  {'~P [kPa]':9s}"
        )
        _log.info(
            "---------------------------------------------------------------------"
        )
        rhol_list = []
        rhov_list = []
        for T in trange:
            self.model.tau = self.model.T_star / T
            delta_l = pyo.value(self.model.delta_l_sat_approx)
            delta_v = pyo.value(self.model.delta_v_sat_approx)
            rho_l = pyo.value(delta_l * self.model.rho_star)
            rho_v = pyo.value(delta_v * self.model.rho_star)
            rhol_list.append(rho_l)
            rhov_list.append(rho_v)
            self.model.delta = delta_v
            P_v = pyo.value(
                rho_v * self.R * T * (1 + self.model.delta * self.model.phir_d)
            )
            _log.info(
                f"{T:7.3f}, {T - 273.15: 8.3f}, {rho_l:14.4f}, {rho_v:14.4f}, {P_v:9.4f}"
            )
        _log.info(
            "=====================================================================\n"
        )
        return (rhol_list, rhov_list)

    def calculate_reference_offset(self, delta, tau, s0, h0):
        """Given delta and tau for a reference state and the reference entropy and enthalpy, calculate
        the reference state offset parameters. Since temperature and density are independent of reference
        state, they can be calculated at a new reference state.

        Args:
            delta (float): reduced density at new reference state
            tau (float): 1/reduced temperature at new reference state
            s0 (float): Entropy at new reference state [kJ/kg/K]
            h0 (float): Enthalpy at new reference state [kJ/kg]

        Returns:
            (tuple): Offset for given reference state
        """
        self.model.tau = tau
        self.model.delta = delta
        s1 = pyo.value(
            self.model.tau * (self.model.phii_t + self.model.phir_t)
            - self.model.phii
            - self.model.phir
        )
        n1_off = s1 - s0
        h1 = pyo.value(
            1.0
            / self.model.tau
            * (
                1
                + self.model.tau * (self.model.phii_t + self.model.phir_t)
                + self.model.delta * self.model.phir_d
            )
        )
        n2_off = h0 - h1
        return (n1_off, n2_off)

    def write_model(self, model, model_name, expressions=None):
        """Write an NL file and create an expression and variable map for the EoS model or just a
        variable map for transport property expressions where there is only one expression per file.

        Args:
            model (ConcreteModel): a Pyomo model
            model_name (str): the model name used in the NL file

        Returns:
            tuple: NL file, expression map, variable map for EOS
        """
        nl_file, smap_id = model.write(f"{self.comp}_expressions_{model_name}.nl")
        for v in model.component_data_objects(pyo.Var):
            v.unfix()
        smap = model.solutions.symbol_map[smap_id]
        var_map = [1000] * 4
        for s, c in smap.bySymbol.items():
            if s.startswith("v"):
                j = int(s[1:])
                var_map[j] = self.variables[c.name]
        if expressions is not None:
            expr_map = [0] * len(expressions)
            for s, c in smap.bySymbol.items():
                if s.startswith("o"):
                    i = expressions[c.name]
                    j = int(s[1:])
                    expr_map[i] = j
            return nl_file, expr_map, var_map
        return nl_file, None, var_map

    def write(self, dry_run=False):
        """Write the parameter and expression files, and _log.info some diagnostics."""
        if dry_run:
            # TODO<jce> this isn't important to IDAES function, so I need to come back later
            #    and figure out how to test this.
            return

        _log.info(
            "\n======================================================================="
        )
        _log.info(f" Writing expressions for {self.comp}")
        _log.info(
            "======================================================================="
        )
        pc = self.calculate_pressure(self.rhoc, self.Tc)
        _log.info(f"Calculated Pc = {pc}, Pc given {self.Pc}")
        _log.info(f"Tc = {self.Tc}, T^* = {self.T_star}")
        _log.info(f"rhoc = {self.rhoc}, rho^* = {self.rho_star}")
        self.Pc = pc
        for name in self.expressions:
            if name not in self.has_expression:
                raise RuntimeError(f"Required expression {name} not provided.")
        nl_file, expr_map, var_map = self.write_model(
            self.model, "eos", self.expressions
        )
        param_dict = {
            "nl_file": nl_file,
            "expr_map": expr_map,
            "var_map": var_map,
            "param": {
                "R": self.R,
                "MW": self.MW,
                "T_star": self.T_star,
                "rho_star": self.rho_star,
                "Tc": self.Tc,
                "rhoc": self.rhoc,
                "Pc": self.Pc,
                "Tt": self.Tt,
                "Pt": self.Pt,
                "rhot_l": self.rhot_l,
                "rhot_v": self.rhot_v,
                "P_min": self.P_min,
                "P_max": self.P_max,
                "rho_max": self.rho_max,
                "T_min": self.T_min,
                "T_max": self.T_max,
                "reference_state_offset": self.reference_state_offset,
            },
        }

        # Add optional models
        for name, short_name in self.optional_expressions.items():
            if name in self.has_expression:
                model = getattr(self, f"model_{short_name}")
                nl_file, _, var_map = self.write_model(model, short_name)
                param_dict[f"nl_file_{short_name}"] = nl_file
                param_dict[f"var_map_{short_name}"] = var_map
                param_dict[f"have_{short_name}"] = True
            else:
                _log.warning(f"Missing optional expression {name}")
                param_dict[f"have_{short_name}"] = False

            with open(f"{self.comp}_parameters.json", "w") as f:
                json.dump(param_dict, f, indent=4)

        trng = []
        linspc = numpy.linspace(self.Tt, self.T_star, 15)
        for i, t in enumerate(linspc):
            if i != 0 and i != len(linspc) - 1:
                t = round(t, 3)
            trng.append(t)
        self.approx_sat_curves(trng)
