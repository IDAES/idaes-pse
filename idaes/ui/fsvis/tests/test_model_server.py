"""
Tests for model_server module
"""
# stdlib
# ext
import pytest
from pyomo.environ import ConcreteModel
# pkg
from idaes.ui.fsvis import model_server, errors, persist
from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.generic_models.unit_models import Flash


@pytest.mark.unit
def test_flowsheet_server_class():
    srv = model_server.FlowsheetServer()
    assert srv.port is not None


@pytest.mark.unit
def test_save_flowsheet(flash_model):
    srv = model_server.FlowsheetServer()
    with pytest.raises(errors.ProcessingError):
        srv.save_flowsheet("oscar", {})
    fs = flash_model.fs
    srv.add_flowsheet("oscar", fs, persist.MemoryDataStore())
    with pytest.raises(errors.ProcessingError):
        srv.save_flowsheet("oscar", -1)


@pytest.fixture(scope="module")
def flash_model():
    """Flash unit model. Use '.fs' attribute to get the flowsheet.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    # Flash properties
    m.fs.properties = BTXParameterBlock(
        default={
            "valid_phase": ("Liq", "Vap"),
            "activity_coeff_model": "Ideal",
            "state_vars": "FTPz",
        }
    )
    # Flash unit
    m.fs.flash = Flash(default={"property_package": m.fs.properties})
    m.fs.flash.inlet.flow_mol.fix(1)
    m.fs.flash.inlet.temperature.fix(368)
    m.fs.flash.inlet.pressure.fix(101325)
    m.fs.flash.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.flash.inlet.mole_frac_comp[0, "toluene"].fix(0.5)
    m.fs.flash.heat_duty.fix(0)
    m.fs.flash.deltaP.fix(0)
    return m
