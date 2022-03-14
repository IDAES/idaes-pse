import sys
import time
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.generic_models.unit_models import Flash
from idaes.ui.fsvis import visualize

def model():
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


def main():
    m = model()
    result = visualize(m.fs, name="Flash", browser=False)
    print(f"[scraper] http://localhost:{result.port}/app?id=Flash\n")
    sys.stdout.flush()
    time.sleep(60)
    return 0


if __name__ == "__main__":
    sys.exit(main())