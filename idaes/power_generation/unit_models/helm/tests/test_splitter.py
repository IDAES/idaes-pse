import pytest
import pyomo.environ as pyo
import idaes.core
from idaes.power_generation.unit_models.helm import HelmSplitter
from idaes.generic_models.properties import iapws95


@pytest.mark.component
def test_splitter():
    m = pyo.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = HelmSplitter(
        default={
            "property_package": m.fs.properties,
            "outlet_list": ["o1", "o2", "o3"]
        }
    )

    Fin = 1e4 # mol/s
    hin = 4000 # J/mol
    Pin = 101325 # Pa

    m.fs.unit.inlet.flow_mol[0].fix(Fin)
    m.fs.unit.inlet.enth_mol[0].fix(hin)
    m.fs.unit.inlet.pressure[0].fix(Pin)
    m.fs.unit.split_fraction[0, "o2"] = 0.3
    m.fs.unit.split_fraction[0, "o3"] = 0.2

    m.fs.unit.initialize()

    assert pyo.value(m.fs.unit.o1.flow_mol[0]) == pytest.approx(1e4*0.5, rel=1e-7)
    assert pyo.value(m.fs.unit.o2.flow_mol[0]) == pytest.approx(1e4*0.3, rel=1e-7)
    assert pyo.value(m.fs.unit.o3.flow_mol[0]) == pytest.approx(1e4*0.2, rel=1e-7)
    assert pyo.value(m.fs.unit.o1.pressure[0]) == pytest.approx(101325, rel=1e-7)
    assert pyo.value(m.fs.unit.o2.pressure[0]) == pytest.approx(101325, rel=1e-7)
    assert pyo.value(m.fs.unit.o3.pressure[0]) == pytest.approx(101325, rel=1e-7)
    assert pyo.value(m.fs.unit.o1.enth_mol[0]) == pytest.approx(4000, rel=1e-7)
    assert pyo.value(m.fs.unit.o2.enth_mol[0]) == pytest.approx(4000, rel=1e-7)
    assert pyo.value(m.fs.unit.o3.enth_mol[0]) == pytest.approx(4000, rel=1e-7)
