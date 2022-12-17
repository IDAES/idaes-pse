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
Tests for Reboiler unit model. Tests 2 sets of state vars using the
ideal property package (FTPz and FcTP).

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import check_optimal_termination, ConcreteModel, value
from pyomo.environ import units as pyunits
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    declare_process_block_class,
)
from idaes.models_extra.column_models import Tray
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)

from idaes.models.properties.modular_properties.state_definitions.FPhx import (
    FPhx,
    define_state,
    set_metadata,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.solvers import get_solver

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.state_definitions import FPhx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity
from idaes.models.properties.modular_properties.pure.RPP4 import RPP4
from idaes.models.properties.modular_properties.pure import Perrys

from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# Mark module as all unit tests
pytestmark = pytest.mark.unit


class TestenthIdeal:
    @pytest.fixture(scope="class")
    def btx(self):

        config_dict = {
            "components": {
                "benzene": {
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP4,
                    "pressure_sat_comp": RPP4,
                    "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                    "parameter_data": {
                        "mw": (78.1136e-3, pyunits.kg / pyunits.mol),  # [1]
                        "pressure_crit": (48.9e5, pyunits.Pa),  # [1]
                        "temperature_crit": (562.2, pyunits.K),  # [1]
                        "dens_mol_liq_comp_coeff": {
                            "eqn_type": 1,
                            "1": (
                                1.0162,
                                pyunits.kmol * pyunits.m**-3,
                            ),  # [2] pg. 2-98
                            "2": (0.2655, None),
                            "3": (562.16, pyunits.K),
                            "4": (0.28212, None),
                        },
                        "cp_mol_ig_comp_coeff": {
                            "A": (-3.392e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                            "B": (4.739e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                            "C": (-3.017e-4, pyunits.J / pyunits.mol / pyunits.K**3),
                            "D": (7.130e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                        },
                        "cp_mol_liq_comp_coeff": {
                            "1": (1.29e2, pyunits.J / pyunits.kmol / pyunits.K),  # [2]
                            "2": (-1.7e-1, pyunits.J / pyunits.kmol / pyunits.K**2),
                            "3": (6.48e-4, pyunits.J / pyunits.kmol / pyunits.K**3),
                            "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                            "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                        },
                        "enth_mol_form_liq_comp_ref": (
                            49.0e3,
                            pyunits.J / pyunits.mol,
                        ),  # [3]
                        "enth_mol_form_vap_comp_ref": (
                            82.9e3,
                            pyunits.J / pyunits.mol,
                        ),  # [3]
                        "pressure_sat_comp_coeff": {
                            "A": (-6.98273, None),  # [1]
                            "B": (1.33213, None),
                            "C": (-2.62863, None),
                            "D": (-3.33399, None),
                        },
                    },
                },
                "toluene": {
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP4,
                    "pressure_sat_comp": RPP4,
                    "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                    "parameter_data": {
                        "mw": (92.1405e-3, pyunits.kg / pyunits.mol),  # [1]
                        "pressure_crit": (41e5, pyunits.Pa),  # [1]
                        "temperature_crit": (591.8, pyunits.K),  # [1]
                        "dens_mol_liq_comp_coeff": {
                            "eqn_type": 1,
                            "1": (
                                0.8488,
                                pyunits.kmol * pyunits.m**-3,
                            ),  # [2] pg. 2-98
                            "2": (0.26655, None),
                            "3": (591.8, pyunits.K),
                            "4": (0.2878, None),
                        },
                        "cp_mol_ig_comp_coeff": {
                            "A": (-2.435e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                            "B": (5.125e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                            "C": (-2.765e-4, pyunits.J / pyunits.mol / pyunits.K**3),
                            "D": (4.911e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                        },
                        "cp_mol_liq_comp_coeff": {
                            "1": (1.40e2, pyunits.J / pyunits.kmol / pyunits.K),  # [2]
                            "2": (-1.52e-1, pyunits.J / pyunits.kmol / pyunits.K**2),
                            "3": (6.95e-4, pyunits.J / pyunits.kmol / pyunits.K**3),
                            "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                            "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                        },
                        "enth_mol_form_liq_comp_ref": (
                            12.0e3,
                            pyunits.J / pyunits.mol,
                        ),  # [3]
                        "enth_mol_form_vap_comp_ref": (
                            50.1e3,
                            pyunits.J / pyunits.mol,
                        ),  # [3]
                        "pressure_sat_comp_coeff": {
                            "A": (-7.28607, None),  # [1]
                            "B": (1.38091, None),
                            "C": (-2.83433, None),
                            "D": (-2.79168, None),
                        },
                    },
                },
            },
            "phases": {
                "Vap": {"equation_of_state": Ideal},
                "Liq": {"equation_of_state": Ideal},
            },
            "state_definition": FPhx,
            "pressure_ref": 100000.0,
            "temperature_ref": 300,
            "state_bounds": {
                "flow_mol": (0, 100, 200),
                "temperature": (290, 345, 400),
                "pressure": (100000.0, 300000.0, 500000.0),
                "enth_mol": (0, 500, 1000),
            },
            "base_units": {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
        }

        m = ConcreteModel()

        m.props = GenericParameterBlock(**config_dict)

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.unit = Tray(
            property_package=m.props,
            has_liquid_side_draw=False,
            has_heat_transfer=True,
            has_pressure_change=True,
        )

        assert len(m.fs.unit.config) == 9

        assert not m.fs.unit.config.is_feed_tray
        assert not m.fs.unit.config.has_liquid_side_draw
        assert not m.fs.unit.config.has_vapor_side_draw

        # Set inputs
        m.fs.unit.liq_in.flow_mol.fix(1)
        m.fs.unit.liq_in.enth_mol.fix(369)
        m.fs.unit.liq_in.pressure.fix(101325)
        m.fs.unit.liq_in.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.liq_in.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.vap_in.flow_mol.fix(1)
        m.fs.unit.vap_in.enth_mol.fix(372)
        m.fs.unit.vap_in.pressure.fix(101325)
        m.fs.unit.vap_in.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.vap_in.mole_frac_comp[0, "toluene"].fix(0.5)

        m.display()

        return m

    @pytest.mark.unit
    def test_build(self, btx):

        assert hasattr(btx.fs.unit, "material_mixing_equations")
        assert hasattr(btx.fs.unit, "enthalpy_mixing_equations")
        assert hasattr(btx.fs.unit, "pressure_drop_equation")

        assert btx.fs.unit.config.has_heat_transfer
        assert hasattr(btx.fs.unit, "heat_duty")

        assert btx.fs.unit.config.has_pressure_change
        assert hasattr(btx.fs.unit, "deltaP")

        # State blocks
        assert hasattr(btx.fs.unit, "properties_in_liq")
        assert hasattr(btx.fs.unit, "properties_in_vap")
        assert hasattr(btx.fs.unit, "properties_out")

        # Ports
        assert hasattr(btx.fs.unit, "liq_in")
        assert hasattr(btx.fs.unit.liq_in, "flow_mol")
        assert hasattr(btx.fs.unit.liq_in, "mole_frac_comp")
        assert hasattr(btx.fs.unit.liq_in, "enth_mol")
        assert hasattr(btx.fs.unit.liq_in, "pressure")

        assert hasattr(btx.fs.unit, "vap_in")
        assert hasattr(btx.fs.unit.vap_in, "flow_mol")
        assert hasattr(btx.fs.unit.vap_in, "mole_frac_comp")
        assert hasattr(btx.fs.unit.vap_in, "enth_mol")
        assert hasattr(btx.fs.unit.vap_in, "pressure")

        assert hasattr(btx.fs.unit, "liq_out")
        assert hasattr(btx.fs.unit.liq_out, "flow_mol")
        assert hasattr(btx.fs.unit.liq_out, "mole_frac_comp")
        assert hasattr(btx.fs.unit.liq_out, "enth_mol")
        assert hasattr(btx.fs.unit.liq_out, "pressure")

        assert hasattr(btx.fs.unit, "vap_out")
        assert hasattr(btx.fs.unit.vap_out, "flow_mol")
        assert hasattr(btx.fs.unit.vap_out, "mole_frac_comp")
        assert hasattr(btx.fs.unit.vap_out, "enth_mol")
        assert hasattr(btx.fs.unit.vap_out, "pressure")

        assert not hasattr(btx.fs.unit, "liq_side_draw")
        assert not hasattr(btx.fs.unit, "liq_side_sf")

        assert not hasattr(btx.fs.unit, "vap_side_draw")
        assert not hasattr(btx.fs.unit, "vap_side_sf")


class TestenthIdeal:
    @pytest.fixture(scope="class")
    def btx(self):

        config_dict = {
            "components": {
                "benzene": {
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP4,
                    "pressure_sat_comp": RPP4,
                    "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                    "parameter_data": {
                        "mw": (78.1136e-3, pyunits.kg / pyunits.mol),  # [1]
                        "pressure_crit": (48.9e5, pyunits.Pa),  # [1]
                        "temperature_crit": (562.2, pyunits.K),  # [1]
                        "dens_mol_liq_comp_coeff": {
                            "eqn_type": 1,
                            "1": (
                                1.0162,
                                pyunits.kmol * pyunits.m**-3,
                            ),  # [2] pg. 2-98
                            "2": (0.2655, None),
                            "3": (562.16, pyunits.K),
                            "4": (0.28212, None),
                        },
                        "cp_mol_ig_comp_coeff": {
                            "A": (-3.392e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                            "B": (4.739e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                            "C": (-3.017e-4, pyunits.J / pyunits.mol / pyunits.K**3),
                            "D": (7.130e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                        },
                        "cp_mol_liq_comp_coeff": {
                            "1": (1.29e2, pyunits.J / pyunits.kmol / pyunits.K),  # [2]
                            "2": (-1.7e-1, pyunits.J / pyunits.kmol / pyunits.K**2),
                            "3": (6.48e-4, pyunits.J / pyunits.kmol / pyunits.K**3),
                            "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                            "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                        },
                        "enth_mol_form_liq_comp_ref": (
                            49.0e3,
                            pyunits.J / pyunits.mol,
                        ),  # [3]
                        "enth_mol_form_vap_comp_ref": (
                            82.9e3,
                            pyunits.J / pyunits.mol,
                        ),  # [3]
                        "pressure_sat_comp_coeff": {
                            "A": (-6.98273, None),  # [1]
                            "B": (1.33213, None),
                            "C": (-2.62863, None),
                            "D": (-3.33399, None),
                        },
                    },
                },
                "toluene": {
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP4,
                    "pressure_sat_comp": RPP4,
                    "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                    "parameter_data": {
                        "mw": (92.1405e-3, pyunits.kg / pyunits.mol),  # [1]
                        "pressure_crit": (41e5, pyunits.Pa),  # [1]
                        "temperature_crit": (591.8, pyunits.K),  # [1]
                        "dens_mol_liq_comp_coeff": {
                            "eqn_type": 1,
                            "1": (
                                0.8488,
                                pyunits.kmol * pyunits.m**-3,
                            ),  # [2] pg. 2-98
                            "2": (0.26655, None),
                            "3": (591.8, pyunits.K),
                            "4": (0.2878, None),
                        },
                        "cp_mol_ig_comp_coeff": {
                            "A": (-2.435e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                            "B": (5.125e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                            "C": (-2.765e-4, pyunits.J / pyunits.mol / pyunits.K**3),
                            "D": (4.911e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                        },
                        "cp_mol_liq_comp_coeff": {
                            "1": (1.40e2, pyunits.J / pyunits.kmol / pyunits.K),  # [2]
                            "2": (-1.52e-1, pyunits.J / pyunits.kmol / pyunits.K**2),
                            "3": (6.95e-4, pyunits.J / pyunits.kmol / pyunits.K**3),
                            "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                            "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                        },
                        "enth_mol_form_liq_comp_ref": (
                            12.0e3,
                            pyunits.J / pyunits.mol,
                        ),  # [3]
                        "enth_mol_form_vap_comp_ref": (
                            50.1e3,
                            pyunits.J / pyunits.mol,
                        ),  # [3]
                        "pressure_sat_comp_coeff": {
                            "A": (-7.28607, None),  # [1]
                            "B": (1.38091, None),
                            "C": (-2.83433, None),
                            "D": (-2.79168, None),
                        },
                    },
                },
            },
            "phases": {
                "Vap": {"equation_of_state": Ideal},
                "Liq": {"equation_of_state": Ideal},
            },
            "state_definition": FPhx,
            "pressure_ref": 100000.0,
            "temperature_ref": 300,
            "state_bounds": {
                "flow_mol": (0, 100, 200),
                "temperature": (290, 345, 400),
                "pressure": (100000.0, 300000.0, 500000.0),
                "enth_mol": (0, 500, 1000),
            },
            "base_units": {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
        }

        m = ConcreteModel()

        m.props = GenericParameterBlock(**config_dict)

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.unit = Tray(
            property_package=m.props,
            has_liquid_side_draw=False,
            has_heat_transfer=True,
            has_pressure_change=True,
        )

        assert len(m.fs.unit.config) == 9

        assert not m.fs.unit.config.is_feed_tray
        assert not m.fs.unit.config.has_liquid_side_draw
        assert not m.fs.unit.config.has_vapor_side_draw

        # Set inputs
        m.fs.unit.liq_in.flow_mol.fix(1)
        m.fs.unit.liq_in.enth_mol.fix(369)
        m.fs.unit.liq_in.pressure.fix(101325)
        m.fs.unit.liq_in.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.liq_in.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.vap_in.flow_mol.fix(1)
        m.fs.unit.vap_in.enth_mol.fix(372)
        m.fs.unit.vap_in.pressure.fix(101325)
        m.fs.unit.vap_in.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.vap_in.mole_frac_comp[0, "toluene"].fix(0.5)

        m.display()

        return m

    @pytest.mark.unit
    def test_build(self, btx):

        assert hasattr(btx.fs.unit, "material_mixing_equations")
        assert hasattr(btx.fs.unit, "enthalpy_mixing_equations")
        assert hasattr(btx.fs.unit, "pressure_drop_equation")

        assert btx.fs.unit.config.has_heat_transfer
        assert hasattr(btx.fs.unit, "heat_duty")

        assert btx.fs.unit.config.has_pressure_change
        assert hasattr(btx.fs.unit, "deltaP")

        # State blocks
        assert hasattr(btx.fs.unit, "properties_in_liq")
        assert hasattr(btx.fs.unit, "properties_in_vap")
        assert hasattr(btx.fs.unit, "properties_out")

        # Ports
        assert hasattr(btx.fs.unit, "liq_in")
        assert hasattr(btx.fs.unit.liq_in, "flow_mol")
        assert hasattr(btx.fs.unit.liq_in, "mole_frac_comp")
        assert hasattr(btx.fs.unit.liq_in, "enth_mol")
        assert hasattr(btx.fs.unit.liq_in, "pressure")

        assert hasattr(btx.fs.unit, "vap_in")
        assert hasattr(btx.fs.unit.vap_in, "flow_mol")
        assert hasattr(btx.fs.unit.vap_in, "mole_frac_comp")
        assert hasattr(btx.fs.unit.vap_in, "enth_mol")
        assert hasattr(btx.fs.unit.vap_in, "pressure")

        assert hasattr(btx.fs.unit, "liq_out")
        assert hasattr(btx.fs.unit.liq_out, "flow_mol")
        assert hasattr(btx.fs.unit.liq_out, "mole_frac_comp")
        assert hasattr(btx.fs.unit.liq_out, "enth_mol")
        assert hasattr(btx.fs.unit.liq_out, "pressure")

        assert hasattr(btx.fs.unit, "vap_out")
        assert hasattr(btx.fs.unit.vap_out, "flow_mol")
        assert hasattr(btx.fs.unit.vap_out, "mole_frac_comp")
        assert hasattr(btx.fs.unit.vap_out, "enth_mol")
        assert hasattr(btx.fs.unit.vap_out, "pressure")

        assert not hasattr(btx.fs.unit, "liq_side_draw")
        assert not hasattr(btx.fs.unit, "liq_side_sf")

        assert not hasattr(btx.fs.unit, "vap_side_draw")
        assert not hasattr(btx.fs.unit, "vap_side_sf")


class TestmassIdeal:
    @pytest.fixture(scope="class")
    def btx_mass(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(
            valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal", state_vars="FTPz"
        )
        m.fs.unit = Tray(
            property_package=m.fs.properties,
            has_liquid_side_draw=True,
            has_heat_transfer=True,
            has_pressure_change=True,
        )

        # Set inputs
        m.fs.unit.liq_in.flow_mol.fix(1)
        m.fs.unit.liq_in.temperature.fix(369)
        m.fs.unit.liq_in.pressure.fix(101325)
        m.fs.unit.liq_in.mass_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.liq_in.mass_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.vap_in.flow_mol.fix(1)
        m.fs.unit.vap_in.temperature.fix(372)
        m.fs.unit.vap_in.pressure.fix(101325)
        m.fs.unit.vap_in.mass_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.vap_in.mass_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.deltaP.fix(0)
        m.fs.unit.heat_duty.fix(0)

        m.fs.unit.liq_side_sf.fix(0.5)

        return m
