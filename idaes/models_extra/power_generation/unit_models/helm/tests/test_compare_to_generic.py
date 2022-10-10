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
"""Tests that helmholtz specific models match generic models"""
import pytest

import pyomo.environ as pyo
import idaes.core
import idaes.models.unit_models as cmodels
import idaes.models_extra.power_generation.unit_models.helm as hmodels

from idaes.models.properties import iapws95


@pytest.mark.component
def test_pump():
    m = pyo.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit1 = cmodels.Pump(property_package=m.fs.properties)
    m.fs.unit2 = hmodels.HelmPump(property_package=m.fs.properties)
    # set inputs

    Fin = 1e4  # mol/s
    hin = 4000  # J/mol
    Pin = 101325  # Pa
    Pout = 2 * Pin  # Pa
    eff = 0.7
    m.fs.unit1.inlet.flow_mol[0].fix(Fin)
    m.fs.unit2.inlet.flow_mol[0].fix(Fin)
    m.fs.unit1.inlet.enth_mol[0].fix(hin)
    m.fs.unit2.inlet.enth_mol[0].fix(hin)
    m.fs.unit1.inlet.pressure[0].fix(Pin)
    m.fs.unit2.inlet.pressure[0].fix(Pin)
    m.fs.unit1.outlet.pressure[0].fix(Pout)
    m.fs.unit2.outlet.pressure[0].fix(Pout)
    m.fs.unit1.efficiency_pump.fix(eff)
    m.fs.unit2.efficiency_pump.fix(eff)
    m.fs.unit1.initialize()
    m.fs.unit2.initialize()

    assert pyo.value(m.fs.unit1.control_volume.work[0]) == pytest.approx(
        pyo.value(m.fs.unit2.control_volume.work[0]), rel=1e-7
    )
    assert pyo.value(
        m.fs.unit1.control_volume.properties_out[0].temperature
    ) == pytest.approx(
        pyo.value(m.fs.unit2.control_volume.properties_out[0].temperature), rel=1e-7
    )


@pytest.mark.component
def test_turbine():
    m = pyo.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit1 = cmodels.Turbine(property_package=m.fs.properties)
    m.fs.unit2 = hmodels.HelmIsentropicTurbine(property_package=m.fs.properties)

    # set inputs
    Fin = 1000  # mol/s
    Tin = 500  # K
    Pin = 1e6  # Pa
    Pout = 7e5  # Pa
    hin = pyo.value(iapws95.htpx(Tin * pyo.units.K, Pin * pyo.units.Pa))  # J/mol

    m.fs.unit1.inlet.flow_mol[0].fix(Fin)
    m.fs.unit2.inlet.flow_mol[0].fix(Fin)
    m.fs.unit1.inlet.enth_mol[0].fix(hin)
    m.fs.unit2.inlet.enth_mol[0].fix(hin)
    m.fs.unit1.inlet.pressure[0].fix(Pin)
    m.fs.unit2.inlet.pressure[0].fix(Pin)
    m.fs.unit1.outlet.pressure[0].fix(Pout)
    m.fs.unit2.outlet.pressure[0].fix(Pout)
    m.fs.unit1.efficiency_isentropic.fix(0.9)
    m.fs.unit2.efficiency_isentropic.fix(0.9)
    m.fs.unit1.initialize()
    m.fs.unit2.initialize()

    assert pyo.value(
        m.fs.unit1.control_volume.properties_out[0].temperature
    ) == pytest.approx(
        pyo.value(m.fs.unit2.control_volume.properties_out[0].temperature), rel=1e-7
    )
    assert pyo.value(m.fs.unit1.control_volume.work[0]) == pytest.approx(
        pyo.value(m.fs.unit2.control_volume.work[0]), rel=1e-7
    )


@pytest.mark.component
def test_compressor():
    m = pyo.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit1 = cmodels.Compressor(property_package=m.fs.properties)
    m.fs.unit2 = hmodels.HelmIsentropicCompressor(property_package=m.fs.properties)

    # set inputs
    Fin = 1000  # mol/s
    Tin = 500  # K
    Pin = 2222782.4  # Pa
    Pout = 2.7 * Pin  # Pa
    hin = pyo.value(iapws95.htpx(Tin * pyo.units.K, Pin * pyo.units.Pa))  # J/mol
    eff = 0.324

    m.fs.unit1.inlet.flow_mol[0].fix(Fin)
    m.fs.unit2.inlet.flow_mol[0].fix(Fin)
    m.fs.unit1.inlet.enth_mol[0].fix(hin)
    m.fs.unit2.inlet.enth_mol[0].fix(hin)
    m.fs.unit1.inlet.pressure[0].fix(Pin)
    m.fs.unit2.inlet.pressure[0].fix(Pin)
    m.fs.unit1.outlet.pressure[0].fix(Pout)
    m.fs.unit2.outlet.pressure[0].fix(Pout)
    m.fs.unit1.efficiency_isentropic.fix(eff)
    m.fs.unit2.efficiency_isentropic.fix(eff)
    m.fs.unit1.initialize()
    m.fs.unit2.initialize()

    assert pyo.value(
        m.fs.unit1.control_volume.properties_out[0].temperature
    ) == pytest.approx(
        pyo.value(m.fs.unit2.control_volume.properties_out[0].temperature), rel=1e-7
    )
    assert pyo.value(m.fs.unit1.control_volume.work[0]) == pytest.approx(
        pyo.value(m.fs.unit2.control_volume.work[0]), rel=1e-7
    )


@pytest.mark.component
def test_compressor_pump_compare():
    m = pyo.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit1 = hmodels.HelmPump(property_package=m.fs.properties)
    m.fs.unit2 = hmodels.HelmIsentropicCompressor(property_package=m.fs.properties)

    Fin = 1e4  # mol/s
    hin = 4000  # J/mol
    Pin = 101325  # Pa
    Pout = 2 * Pin  # Pa
    eff = 0.7
    m.fs.unit1.inlet.flow_mol[0].fix(Fin)
    m.fs.unit1.inlet.enth_mol[0].fix(hin)
    m.fs.unit1.inlet.pressure[0].fix(Pin)
    m.fs.unit1.outlet.pressure[0].fix(Pout)
    m.fs.unit1.efficiency_pump.fix(eff)

    m.fs.unit2.inlet.flow_mol[0].fix(Fin)
    m.fs.unit2.inlet.enth_mol[0].fix(hin)
    m.fs.unit2.inlet.pressure[0].fix(Pin)
    m.fs.unit2.outlet.pressure[0].fix(Pout)
    m.fs.unit2.efficiency_isentropic.fix(eff)

    m.fs.unit1.initialize()
    m.fs.unit2.initialize()

    # The pump calculations are a bit more approximate assuming incompressible
    # fluid entropy is independent of pressure, so the results here should be
    # close, but not exactly the same.
    assert pyo.value(
        m.fs.unit1.control_volume.properties_out[0].temperature
    ) == pytest.approx(
        pyo.value(m.fs.unit2.control_volume.properties_out[0].temperature), rel=1e-3
    )
    assert pyo.value(m.fs.unit1.control_volume.work[0]) == pytest.approx(
        pyo.value(m.fs.unit2.control_volume.work[0]), rel=1e-3
    )
