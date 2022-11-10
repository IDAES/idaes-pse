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
NRTL property model for a benzene-toluene mixture.
The example model is from the IDAES tutorial,
https://github.com/IDAES/examples-pse/blob/main/src/Tutorials/Advanced/ParamEst/
"""
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Flash
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
import idaes.logger as idaeslog
from pyomo.environ import *


def NRTL_model(data):
    """This function generates an instance of the NRTL Pyomo model using 'data' as the input argument

    Parameters
    ----------
    data: pandas DataFrame, list of dictionaries, or list of json file names
        Data that is used to build an instance of the Pyomo model

    Returns
    -------
    m: an instance of the Pyomo model
        for estimating parameters and covariance
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BTXParameterBlock(
        valid_phase=("Liq", "Vap"), activity_coeff_model="NRTL"
    )
    m.fs.flash = Flash(property_package=m.fs.properties)

    # Initialize at a certain inlet condition
    m.fs.flash.inlet.flow_mol.fix(1)
    m.fs.flash.inlet.temperature.fix(368)
    m.fs.flash.inlet.pressure.fix(101325)
    m.fs.flash.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.flash.inlet.mole_frac_comp[0, "toluene"].fix(0.5)

    # Set Flash unit specifications
    m.fs.flash.heat_duty.fix(0)
    m.fs.flash.deltaP.fix(0)

    # Fix NRTL specific variables
    # alpha values (set at 0.3)
    m.fs.properties.alpha["benzene", "benzene"].fix(0)
    m.fs.properties.alpha["benzene", "toluene"].fix(0.3)
    m.fs.properties.alpha["toluene", "toluene"].fix(0)
    m.fs.properties.alpha["toluene", "benzene"].fix(0.3)

    # initial tau values
    m.fs.properties.tau["benzene", "benzene"].fix(0)
    m.fs.properties.tau["benzene", "toluene"].fix(0.1690)
    m.fs.properties.tau["toluene", "toluene"].fix(0)
    m.fs.properties.tau["toluene", "benzene"].fix(-0.1559)

    # Initialize the flash unit
    m.fs.flash.initialize(outlvl=idaeslog.INFO_LOW)

    # Fix at actual temperature
    m.fs.flash.inlet.temperature.fix(float(data["temperature"]))

    # Set bounds on variables to be estimated
    m.fs.properties.tau["benzene", "toluene"].setlb(-5)
    m.fs.properties.tau["benzene", "toluene"].setub(5)

    m.fs.properties.tau["toluene", "benzene"].setlb(-5)
    m.fs.properties.tau["toluene", "benzene"].setub(5)

    # Return initialized flash model
    return m


def NRTL_model_opt():
    """This function generates an instance of the NRTL Pyomo model

    Returns
    -------
    m: an instance of the Pyomo model
        for uncertainty propagation
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BTXParameterBlock(
        valid_phase=("Liq", "Vap"), activity_coeff_model="NRTL"
    )
    m.fs.flash = Flash(property_package=m.fs.properties)

    # Initialize at a certain inlet condition
    m.fs.flash.inlet.flow_mol.fix(1)
    m.fs.flash.inlet.temperature.fix(368)
    m.fs.flash.inlet.pressure.fix(101325)
    m.fs.flash.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.flash.inlet.mole_frac_comp[0, "toluene"].fix(0.5)

    # Set Flash unit specifications
    m.fs.flash.heat_duty.fix(0)
    m.fs.flash.deltaP.fix(0)

    # Fix NRTL specific variables
    # alpha values (set at 0.3)
    m.fs.properties.alpha["benzene", "benzene"].fix(0)
    m.fs.properties.alpha["benzene", "toluene"].fix(0.3)
    m.fs.properties.alpha["toluene", "toluene"].fix(0)
    m.fs.properties.alpha["toluene", "benzene"].fix(0.3)

    # initial tau values
    m.fs.properties.tau["benzene", "benzene"].fix(0)
    m.fs.properties.tau["benzene", "toluene"].fix(0.1690)
    m.fs.properties.tau["toluene", "toluene"].fix(0)
    m.fs.properties.tau["toluene", "benzene"].fix(-0.1559)

    # Initialize the flash unit
    m.fs.flash.initialize(outlvl=idaeslog.INFO_LOW)

    # Fix at actual temperature
    m.fs.flash.inlet.temperature.fix(float(368))

    # Set bounds on variables to be estimated
    m.fs.properties.tau["benzene", "toluene"].setlb(-5)
    m.fs.properties.tau["benzene", "toluene"].setub(5)

    m.fs.properties.tau["toluene", "benzene"].setlb(-5)
    m.fs.properties.tau["toluene", "benzene"].setub(5)

    # To use kaug
    # objective function required
    # need to unfix the variables
    m.obj = Objective(
        expr=0 * m.fs.properties.tau["benzene", "toluene"]
        + exp(
            -m.fs.properties.alpha["toluene", "benzene"].value
            * m.fs.properties.tau["toluene", "benzene"]
        ),
        sense=minimize,
    )
    m.fs.properties.tau["benzene", "toluene"].fixed = False  # To use kaug
    m.fs.properties.tau["toluene", "benzene"].fixed = False
    return m


def NRTL_model_opt_infeasible():
    """This function generates an instance of the NRTL Pyomo model

    Returns
    -------
    m: an instance of the Pyomo model
        for uncertainty propagation
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BTXParameterBlock(
        valid_phase=("Liq", "Vap"), activity_coeff_model="NRTL"
    )
    m.fs.flash = Flash(property_package=m.fs.properties)

    # Initialize at a certain inlet condition
    m.fs.flash.inlet.flow_mol.fix(1)
    m.fs.flash.inlet.temperature.fix(368)
    m.fs.flash.inlet.pressure.fix(101325)
    m.fs.flash.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.flash.inlet.mole_frac_comp[0, "toluene"].fix(0.5)

    # Set Flash unit specifications
    m.fs.flash.heat_duty.fix(0)
    m.fs.flash.deltaP.fix(0)

    # Fix NRTL specific variables
    # alpha values (set at 0.3)
    m.fs.properties.alpha["benzene", "benzene"].fix(0)
    m.fs.properties.alpha["benzene", "toluene"].fix(0.3)
    m.fs.properties.alpha["toluene", "toluene"].fix(0)
    m.fs.properties.alpha["toluene", "benzene"].fix(0.3)

    # initial tau values
    m.fs.properties.tau["benzene", "benzene"].fix(0)
    m.fs.properties.tau["benzene", "toluene"].fix(0.1690)
    m.fs.properties.tau["toluene", "toluene"].fix(0)
    m.fs.properties.tau["toluene", "benzene"].fix(-0.1559)

    # Set bounds on variables to be estimated
    m.fs.properties.tau["benzene", "toluene"].setlb(-5)
    m.fs.properties.tau["benzene", "toluene"].setub(5)

    m.fs.properties.tau["toluene", "benzene"].setlb(-5)
    m.fs.properties.tau["toluene", "benzene"].setub(5)

    # To use kaug
    # objective function required
    # need to unfix the variables
    m.obj = Objective(
        expr=0 * m.fs.properties.tau["benzene", "toluene"]
        + exp(
            -m.fs.properties.alpha["toluene", "benzene"].value
            * m.fs.properties.tau["toluene", "benzene"]
        ),
        sense=minimize,
    )
    m.fs.properties.tau["benzene", "toluene"].fixed = False
    m.fs.properties.tau["toluene", "benzene"].fixed = False
    return m
