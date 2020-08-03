import pytest
import pyomo.environ as pyo
import idaes.core
from idaes.power_generation.unit_models.helm import HelmMixer, MomentumMixingType
from idaes.generic_models.properties import iapws95


@pytest.mark.component
def test_mixer():
    m = pyo.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = HelmMixer(
        default={
            "property_package": m.fs.properties,
            "inlet_list": ["i1", "i2", "i3"]
        }
    )

    Fin1 = 1.1e4 # mol/s
    hin1 = 4000 # J/mol
    Pin1 = 1.2e5 # Pa
    Fin2 = 1e4 # mol/s
    hin2 = 5000 # J/mol
    Pin2 = 2e5 # Pa
    Fin3 = 1.3e4 # mol/s
    hin3 = 6000 # J/mol
    Pin3 = 3e5 # Pa
    Pout = 1.5e5 # Pa

    m.fs.unit.i1.flow_mol[0].fix(Fin1)
    m.fs.unit.i1.enth_mol[0].fix(hin1)
    m.fs.unit.i1.pressure[0].fix(Pin1)
    m.fs.unit.i2.flow_mol[0].fix(Fin2)
    m.fs.unit.i2.enth_mol[0].fix(hin2)
    m.fs.unit.i2.pressure[0].fix(Pin2)
    m.fs.unit.i3.flow_mol[0].fix(Fin3)
    m.fs.unit.i3.enth_mol[0].fix(hin3)
    m.fs.unit.i3.pressure[0].fix(Pin3)

    m.fs.unit.initialize()
    Fout = Fin1 + Fin2 + Fin3
    hout = (hin1*Fin1 + hin2*Fin2 + hin3*Fin3)/Fout

    assert pyo.value(m.fs.unit.outlet.flow_mol[0]) == pytest.approx(Fout, rel=1e-7)
    assert pyo.value(m.fs.unit.outlet.enth_mol[0]) == pytest.approx(hout, rel=1e-7)
    assert pyo.value(m.fs.unit.outlet.pressure[0]) == pytest.approx(Pin1, rel=1e-7)


@pytest.mark.component
def test_mixer2():
    m = pyo.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = HelmMixer(
        default={
            "momentum_mixing_type": MomentumMixingType.equality,
            "property_package": m.fs.properties,
            "inlet_list": ["i1", "i2", "i3"]
        }
    )

    Fin1 = 1.1e4 # mol/s
    hin1 = 4000 # J/mol
    Pin1 = 1.2e5 # Pa
    Fin2 = 1e4 # mol/s
    hin2 = 5000 # J/mol
    Pin2 = 2e5 # Pa
    Fin3 = 1.3e4 # mol/s
    hin3 = 6000 # J/mol
    Pin3 = 3e5 # Pa
    Pout = 1.5e5 # Pa

    m.fs.unit.i1.flow_mol[0].fix(Fin1)
    m.fs.unit.i1.enth_mol[0].fix(hin1)
    m.fs.unit.i1.pressure[0] = Pin1
    m.fs.unit.i2.flow_mol[0].fix(Fin2)
    m.fs.unit.i2.enth_mol[0].fix(hin2)
    m.fs.unit.i2.pressure[0] = Pin2
    m.fs.unit.i3.flow_mol[0].fix(Fin3)
    m.fs.unit.i3.enth_mol[0].fix(hin3)
    m.fs.unit.i3.pressure[0] = Pin3
    m.fs.unit.outlet.pressure[0].fix(Pout)

    m.fs.unit.initialize()
    Fout = Fin1 + Fin2 + Fin3
    hout = (hin1*Fin1 + hin2*Fin2 + hin3*Fin3)/Fout

    assert pyo.value(m.fs.unit.outlet.flow_mol[0]) == pytest.approx(Fout, rel=1e-7)
    assert pyo.value(m.fs.unit.outlet.enth_mol[0]) == pytest.approx(hout, rel=1e-7)
    assert pyo.value(m.fs.unit.outlet.pressure[0]) == pytest.approx(Pout, rel=1e-7)
    assert pyo.value(m.fs.unit.i1.pressure[0]) == pytest.approx(Pout, rel=1e-7)
    assert pyo.value(m.fs.unit.i2.pressure[0]) == pytest.approx(Pout, rel=1e-7)
    assert pyo.value(m.fs.unit.i3.pressure[0]) == pytest.approx(Pout, rel=1e-7)
