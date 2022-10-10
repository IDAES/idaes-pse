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
import pyomo.common.unittest as unittest
import pytest
import pyomo.environ as pyo
from pyomo.util.check_units import (
    assert_units_equivalent,
    assert_units_consistent,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.models_extra.gas_distribution.properties.natural_gas import (
    NaturalGasParameterBlock,
)

"""
Test for simple ideal methane property package created with generic
property framework.
"""


@pytest.mark.unit
class TestNaturalGasPropertyPackage(unittest.TestCase):
    def test_build(self):
        m = pyo.ConcreteModel()
        m.properties = NaturalGasParameterBlock()
        m.state = m.properties.build_state_block(defined_state=False)
        state = m.state

        # The following tests are very specific to this property package.
        # They test units and values, which may be different for other
        # property packages that are valid for the pipeline models.

        assert_units_equivalent(state.pressure.get_units(), pyo.units.bar)
        assert_units_equivalent(
            state.flow_mol.get_units(), pyo.units.kmol / pyo.units.hr
        )
        assert_units_equivalent(state.temperature.get_units(), pyo.units.K)
        component_list = m.state.config.parameters.component_list
        j = next(iter(component_list))
        self.assertIs(state.mole_frac_comp[j].get_units(), None)

        mw = m.state.mw
        assert_units_equivalent(pyo.units.get_units(mw), pyo.units.kg / pyo.units.kmol)
        self.assertEqual(pyo.value(mw), 18.0)

        flow_mass = m.state.flow_mass
        assert_units_equivalent(
            pyo.units.get_units(flow_mass), pyo.units.kg / pyo.units.hr
        )
        m.state.flow_mol.set_value(100.0 * pyo.units.kmol / pyo.units.hr)
        self.assertAlmostEqual(pyo.value(flow_mass), 1800.0)

        kJkmolK = pyo.units.kJ / pyo.units.kmol / pyo.units.K
        cp_mol_comp = m.state.cp_mol_comp[j]
        assert_units_equivalent(pyo.units.get_units(cp_mol_comp), kJkmolK)
        self.assertAlmostEqual(pyo.value(cp_mol_comp), 2.34 * 18.0)

        cp_mol = m.state.cp_mol
        assert_units_equivalent(pyo.units.get_units(cp_mol), kJkmolK)
        self.assertAlmostEqual(pyo.value(cp_mol), 2.34 * 18.0)

        kJkgK = pyo.units.kJ / pyo.units.kg / pyo.units.K
        cp_mass = m.state.cp_mass
        assert_units_equivalent(cp_mass, kJkgK)
        self.assertEqual(pyo.value(cp_mass), 2.34)

        cv_mol_comp = m.state.cv_mol_comp[j]
        assert_units_equivalent(pyo.units.get_units(cv_mol_comp), kJkmolK)
        self.assertAlmostEqual(pyo.value(cv_mol_comp), 1.85 * 18.0)

        cv_mol = m.state.cv_mol
        assert_units_equivalent(pyo.units.get_units(cv_mol_comp), kJkmolK)
        self.assertAlmostEqual(pyo.value(cv_mol), 1.85 * 18.0)

        cv_mass = m.state.cv_mass
        assert_units_equivalent(pyo.units.get_units(cv_mass), kJkgK)
        self.assertAlmostEqual(pyo.value(cv_mass), 1.85)

        heat_capacity_ratio = m.state.heat_capacity_ratio
        assert_units_equivalent(
            pyo.units.get_units(heat_capacity_ratio), pyo.units.dimensionless
        )
        self.assertAlmostEqual(pyo.value(heat_capacity_ratio), 2.34 / 1.85)

        speed_of_sound = m.state.speed_of_sound
        speed_of_sound_eq = m.state.speed_of_sound_eq
        assert_units_equivalent(speed_of_sound, pyo.units.m / pyo.units.s)
        assert_units_consistent(speed_of_sound_eq)
        m.state.temperature.set_value(293.15 * pyo.units.K)
        calculate_variable_from_constraint(speed_of_sound, speed_of_sound_eq)
        # 370.2 = (0.8*2.34/1.85*8314*293.15/18.0)**0.5
        self.assertAlmostEqual(pyo.value(speed_of_sound), 370.1628645)

        dens_mol = m.state.dens_mol
        dens_mol_eq = m.state.dens_mol_eq
        assert_units_equivalent(dens_mol, pyo.units.kmol / pyo.units.m**3)
        assert_units_consistent(dens_mol_eq)
        m.state.temperature.set_value(293.15 * pyo.units.K)
        m.state.pressure.set_value(50 * pyo.units.bar)
        calculate_variable_from_constraint(dens_mol, dens_mol_eq)
        # 2.654 = 50*1e2/0.8/8.314/293.15
        self.assertAlmostEqual(pyo.value(dens_mol), 2.5642238411)

        component_list = m.state.config.parameters.component_list
        j = next(iter(component_list))
        dens_mol_comp = m.state.dens_mol_comp[j]
        assert_units_equivalent(
            pyo.units.get_units(dens_mol_comp), pyo.units.kmol / pyo.units.m**3
        )
        self.assertAlmostEqual(pyo.value(dens_mol_comp), 2.5642238411)

        # flow_mol was set to 100 kmol/hr earlier.
        m.state.mole_frac_comp[j].set_value(1.0)
        self.assertEqual(
            pyo.value(m.state.flow_mol_comp[j]),
            100.0,
        )
        assert_units_equivalent(
            pyo.units.get_units(m.state.flow_mol_comp[j]),
            pyo.units.kmol / pyo.units.hr,
        )

        enth_flow = m.state.get_enthalpy_flow_terms("Vap")
        pred_enth_flow = pyo.value(
            (293.15 - 298.15) * m.state.cp_mol * m.state.flow_mol
        )
        self.assertAlmostEqual(
            pyo.value(enth_flow),
            pred_enth_flow,
        )
        kJhr = pyo.units.kJ / pyo.units.hr
        assert_units_equivalent(
            pyo.units.get_units(enth_flow),
            kJhr,
        )

    def test_nominal_density(self):
        m = pyo.ConcreteModel()
        m.properties = NaturalGasParameterBlock()
        self.assertEqual(m.properties.dens_nominal.value, 0.72)
        assert_units_equivalent(
            m.properties.dens_nominal.get_units(),
            pyo.units.kg / pyo.units.m**3,
        )

    def test_compressibility(self):
        m = pyo.ConcreteModel()
        m.properties = NaturalGasParameterBlock()
        m.state = m.properties.build_state_block()
        self.assertEqual(m.state.compressibility.value, 0.80)

    def test_temperature_ref(self):
        m = pyo.ConcreteModel()
        m.properties = NaturalGasParameterBlock()
        self.assertEqual(pyo.value(m.properties.temperature_ref), 298.15)
        assert_units_equivalent(
            m.properties.temperature_ref,
            pyo.units.K,
        )

    def test_pprint_noerror(self):
        m = pyo.ConcreteModel()
        m.properties = NaturalGasParameterBlock()
        m.state = m.properties.build_state_block()
        m.state.dens_mol
        m.state.pprint()


if __name__ == "__main__":
    unittest.main()
