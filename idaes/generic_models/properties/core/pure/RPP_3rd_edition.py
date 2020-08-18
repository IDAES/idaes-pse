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
Methods for calculating pure component properties from:

The Properties of Gases & Liquids, 3rd Edition
Reid, Prausnitz and Polling, 1987, McGraw-Hill

All parameter indicies based on conventions used by the source
"""

from pyomo.environ import exp, value, log, Var, units as pyunits

from idaes.generic_models.properties.core.generic.utility import \
    set_param_value

# -----------------------------------------------------------------------------
# Saturation pressure
# Note that this equation in not valid beyond the critical temperature
class pressure_sat_comp():
    def build_parameters(cobj):
        cobj.pressure_sat_comp_coeff_A = Var(
                doc="Antoine A coefficient for calculating Psat",
                units=None)
        set_param_value(cobj,
                        param="pressure_sat_comp_coeff",
                        units=None,
                        index="A")

        cobj.pressure_sat_comp_coeff_B = Var(
                doc="Antoine B coefficient for calculating Psat",
                units=pyunits.K)
        set_param_value(cobj,
                        param="pressure_sat_comp_coeff",
                        units=pyunits.K,
                        index="B")

        cobj.pressure_sat_comp_coeff_C = Var(
                doc="Antoine C coefficient for calculating Psat",
                units=pyunits.K)
        set_param_value(cobj,
                        param="pressure_sat_comp_coeff",
                        units=pyunits.K,
                        index="C")

    def return_expression(b, cobj, T, dT=False):
        if dT:
            return pressure_sat_comp.dT_expression(b, cobj, T)

        psat = (exp(cobj.pressure_sat_comp_coeff_A -
                    cobj.pressure_sat_comp_coeff_B /
                    (pyunits.convert(T, to_units=pyunits.K) +
                     cobj.pressure_sat_comp_coeff_C)))*pyunits.mmHg

        base_units = b.params.get_metadata().default_units
        p_units = (base_units["mass"] *
                   base_units["length"]**-1 *
                   base_units["time"]**-2)
        return pyunits.convert(psat, to_units=p_units)

    def dT_expression(b, cobj, T):
        p_sat_dT = (pressure_sat_comp.return_expression(b, cobj, T) *
                    cobj.pressure_sat_comp_coeff_B /
                    (pyunits.convert(T, to_units=pyunits.K) +
                             cobj.pressure_sat_comp_coeff_C)**2)

        base_units = b.params.get_metadata().default_units
        dp_units = (base_units["mass"] *
                    base_units["length"]**-1 *
                    base_units["time"]**-2 *
                    base_units["temperature"]**-1)
        return pyunits.convert(p_sat_dT, to_units=dp_units)
