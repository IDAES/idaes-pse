import pyomo.common.unittest as unittest
import pytest
import pyomo.environ as pyo
from pyomo.util.check_units import (
    assert_units_equivalent,
    assert_units_consistent,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.gas_distribution.properties.natural_gas import (
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
        m.state = m.properties.build_state_block(
            # As usual, if defined_state is True, we need to specify
            # every mole fraction
            default={"defined_state": False},
        )
        state = m.state

        # The following tests are very specific to this property package.
        # They test units and values, which may be different for other
        # property packages that are valid for the pipeline models.

        assert_units_equivalent(state.pressure.get_units(), pyo.units.bar)
        assert_units_equivalent(
            state.flow_mol.get_units(), pyo.units.kmol/pyo.units.hr
        )
        assert_units_equivalent(state.temperature.get_units(), pyo.units.K)
        component_list = m.state.config.parameters.component_list
        j = next(iter(component_list))
        self.assertIs(state.mole_frac_comp[j].get_units(), None)
        # TODO: How do I go from component_list to the actual component
        # object?

        mw = m.state.mw
        assert_units_equivalent(
            pyo.units.get_units(mw), pyo.units.kg/pyo.units.kmol
        )
        self.assertEqual(pyo.value(mw), 18.0)

        flow_mass = m.state.flow_mass
        assert_units_equivalent(
            pyo.units.get_units(flow_mass), pyo.units.kg/pyo.units.hr
        )
        m.state.flow_mol.set_value(100.0*pyo.units.kmol/pyo.units.hr)
        self.assertAlmostEqual(pyo.value(flow_mass), 1800.0)

        kJkmolK = pyo.units.kJ/pyo.units.kmol/pyo.units.K
        cp_mol_comp = m.state.cp_mol_comp[j]
        assert_units_equivalent(
            pyo.units.get_units(cp_mol_comp), kJkmolK
        )
        self.assertAlmostEqual(pyo.value(cp_mol_comp), 2.34*18.0)

        cp_mol = m.state.cp_mol
        assert_units_equivalent(pyo.units.get_units(cp_mol), kJkmolK)
        self.assertAlmostEqual(pyo.value(cp_mol), 2.34*18.0)

        kJkgK = pyo.units.kJ/pyo.units.kg/pyo.units.K
        cp_mass = m.state.cp_mass
        assert_units_equivalent(cp_mass, kJkgK)
        self.assertEqual(pyo.value(cp_mass), 2.34)

        cv_mol_comp = m.state.cv_mol_comp[j]
        assert_units_equivalent(pyo.units.get_units(cv_mol_comp), kJkmolK)
        self.assertAlmostEqual(pyo.value(cv_mol_comp), 1.85*18.0)

        cv_mol = m.state.cv_mol
        assert_units_equivalent(pyo.units.get_units(cv_mol_comp), kJkmolK)
        self.assertAlmostEqual(pyo.value(cv_mol), 1.85*18.0)

        cv_mass = m.state.cv_mass
        assert_units_equivalent(pyo.units.get_units(cv_mass), kJkgK)
        self.assertAlmostEqual(pyo.value(cv_mass), 1.85)

        heat_capacity_ratio = m.state.heat_capacity_ratio
        assert_units_equivalent(
            pyo.units.get_units(heat_capacity_ratio), pyo.units.dimensionless
        )
        self.assertAlmostEqual(pyo.value(heat_capacity_ratio), 2.34/1.85)

        speed_of_sound = m.state.speed_of_sound
        speed_of_sound_eq = m.state.speed_of_sound_eq
        assert_units_equivalent(speed_of_sound, pyo.units.m/pyo.units.s)
        assert_units_consistent(speed_of_sound_eq)
        m.state.temperature.set_value(293.15*pyo.units.K)
        calculate_variable_from_constraint(speed_of_sound, speed_of_sound_eq)
        # 370.2 = (0.8*2.34/1.85*8314*293.15/18.0)**0.5
        self.assertAlmostEqual(pyo.value(speed_of_sound), 370.1628645)

        dens_mol = m.state.dens_mol
        dens_mol_eq = m.state.dens_mol_eq
        assert_units_equivalent(dens_mol, pyo.units.kmol/pyo.units.m**3)
        assert_units_consistent(dens_mol_eq)
        m.state.temperature.set_value(293.15*pyo.units.K)
        m.state.pressure.set_value(50*pyo.units.bar)
        calculate_variable_from_constraint(dens_mol, dens_mol_eq)
        # 2.654 = 50*1e2/0.8/8.314/293.15
        self.assertAlmostEqual(pyo.value(dens_mol), 2.5642238411)

        component_list = m.state.config.parameters.component_list
        j = next(iter(component_list))
        dens_mol_comp = m.state.dens_mol_comp[j]
        assert_units_equivalent(
            pyo.units.get_units(dens_mol_comp),
            pyo.units.kmol/pyo.units.m**3
        )
        self.assertAlmostEqual(pyo.value(dens_mol_comp), 2.5642238411)


if __name__ == "__main__":
    unittest.main()
