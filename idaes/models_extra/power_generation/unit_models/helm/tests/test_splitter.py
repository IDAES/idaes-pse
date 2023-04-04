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
import pytest
import pyomo.environ as pyo
import idaes.core
from idaes.models_extra.power_generation.unit_models.helm import HelmSplitter
from idaes.models.properties import iapws95
from idaes.models.properties.general_helmholtz import helmholtz_available


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
@pytest.mark.component
def test_splitter():
    m = pyo.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = HelmSplitter(
        property_package=m.fs.properties, outlet_list=["o1", "o2", "o3"]
    )

    Fin = 1e4  # mol/s
    hin = 4000  # J/mol
    Pin = 101325  # Pa

    m.fs.unit.inlet.flow_mol[0].fix(Fin)
    m.fs.unit.inlet.enth_mol[0].fix(hin)
    m.fs.unit.inlet.pressure[0].fix(Pin)
    m.fs.unit.split_fraction[0, "o2"] = 0.3
    m.fs.unit.split_fraction[0, "o3"] = 0.2

    m.fs.unit.initialize()

    assert pyo.value(m.fs.unit.o1.flow_mol[0]) == pytest.approx(1e4 * 0.5, rel=1e-7)
    assert pyo.value(m.fs.unit.o2.flow_mol[0]) == pytest.approx(1e4 * 0.3, rel=1e-7)
    assert pyo.value(m.fs.unit.o3.flow_mol[0]) == pytest.approx(1e4 * 0.2, rel=1e-7)
    assert pyo.value(m.fs.unit.o1.pressure[0]) == pytest.approx(101325, rel=1e-7)
    assert pyo.value(m.fs.unit.o2.pressure[0]) == pytest.approx(101325, rel=1e-7)
    assert pyo.value(m.fs.unit.o3.pressure[0]) == pytest.approx(101325, rel=1e-7)
    assert pyo.value(m.fs.unit.o1.enth_mol[0]) == pytest.approx(4000, rel=1e-7)
    assert pyo.value(m.fs.unit.o2.enth_mol[0]) == pytest.approx(4000, rel=1e-7)
    assert pyo.value(m.fs.unit.o3.enth_mol[0]) == pytest.approx(4000, rel=1e-7)


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
@pytest.mark.unit
def test_get_stream_table_contents():
    m = pyo.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = HelmSplitter(
        property_package=m.fs.properties, outlet_list=["o1", "o2", "o3"]
    )

    stable = m.fs.unit._get_stream_table_contents()

    expected = {
        "Units": {
            "Mass Flow": getattr(pyo.units.pint_registry, "kg/s"),
            "Molar Flow": getattr(pyo.units.pint_registry, "mol/s"),
            "Molar Enthalpy": getattr(pyo.units.pint_registry, "J/mol"),
            "P": getattr(pyo.units.pint_registry, "Pa"),
            "T": getattr(pyo.units.pint_registry, "K"),
            "Vapor Fraction": getattr(pyo.units.pint_registry, "dimensionless"),
        },
        "inlet": {
            "Mass Flow": pytest.approx(0.01801527, rel=1e-5),
            "Molar Flow": pytest.approx(1.0, rel=1e-5),
            "Molar Enthalpy": pytest.approx(0.01102139, rel=1e-5),
            "P": pytest.approx(11032300, rel=1e-5),
            "T": pytest.approx(270.4877, rel=1e-5),
            "Vapor Fraction": pytest.approx(0.0, abs=1e-5),
        },
        "o1": {
            "Mass Flow": pytest.approx(0.01801527, rel=1e-5),
            "Molar Flow": pytest.approx(1.0, rel=1e-5),
            "Molar Enthalpy": pytest.approx(0.01102139, rel=1e-5),
            "P": pytest.approx(11032300, rel=1e-5),
            "T": pytest.approx(270.4877, rel=1e-5),
            "Vapor Fraction": pytest.approx(0.0, abs=1e-5),
        },
        "o2": {
            "Mass Flow": pytest.approx(0.01801527, rel=1e-5),
            "Molar Flow": pytest.approx(1.0, rel=1e-5),
            "Molar Enthalpy": pytest.approx(0.01102139, rel=1e-5),
            "P": pytest.approx(11032300, rel=1e-5),
            "T": pytest.approx(270.4877, rel=1e-5),
            "Vapor Fraction": pytest.approx(0.0, abs=1e-5),
        },
        "o3": {
            "Mass Flow": pytest.approx(0.01801527, rel=1e-5),
            "Molar Flow": pytest.approx(1.0, rel=1e-5),
            "Molar Enthalpy": pytest.approx(0.01102139, rel=1e-5),
            "P": pytest.approx(11032300, rel=1e-5),
            "T": pytest.approx(270.4877, rel=1e-5),
            "Vapor Fraction": pytest.approx(0.0, abs=1e-5),
        },
    }

    assert stable.to_dict() == expected
