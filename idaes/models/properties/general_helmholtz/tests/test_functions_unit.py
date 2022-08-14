import pytest

import pyomo.environ as pyo

from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    HelmholtzThermoExpressions,
    PhaseType,
    StateVars,
    AmountBasis,
    add_helmholtz_external_functions,
    available,
)


@pytest.mark.unit
def test_available():
    assert available()


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_available2():
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE
    )
    assert m.hparam.available()


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_add_funcion():
    """Test mixed phase form with P-H state vars and phase mass balances"""
    m = pyo.ConcreteModel()

    # Test add one
    add_helmholtz_external_functions(m, "p_func")
    assert isinstance(m.p_func, pyo.ExternalFunction)

    # Test add list
    add_helmholtz_external_functions(m, "h_func")
    add_helmholtz_external_functions(m, "s_func")
    assert isinstance(m.h_func, pyo.ExternalFunction)
    assert isinstance(m.s_func, pyo.ExternalFunction)

    # Test add all
    add_helmholtz_external_functions(m)
    assert isinstance(m.phi0_func, pyo.ExternalFunction)


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_htpx_mass():
    """Test mixed phase form with P-H state vars and phase mass balances"""
    param = HelmholtzParameterBlock(pure_component="h2o", amount_basis=AmountBasis.MASS)
    param.construct()

    T = 300 * pyo.units.K
    units = pyo.units.kJ / pyo.units.kg

    h = param.htpx(T=T, x=0.0, with_units=False, units=units)
    assert h == pytest.approx(112.56, rel=1e-4)
    h = param.htpx(T=T, x=1.0, with_units=False, units=units)
    assert h == pytest.approx(2549.9, rel=1e-4)
    h = param.htpx(T=T, x=0.0, with_units=False, units=None)
    assert h == pytest.approx(112.56*1000, rel=1e-4)
    h = param.htpx(T=T, x=1.0, with_units=False, units=None)
    assert h == pytest.approx(2549.9*1000, rel=1e-4)
    h = param.htpx(T=T, x=0.0, with_units=True, units=None)
    assert pyo.value(pyo.units.convert(h, units)) == pytest.approx(112.56, rel=1e-4)
    h = param.htpx(T=T, x=1.0, with_units=True, units=None)
    assert pyo.value(pyo.units.convert(h, units)) == pytest.approx(2549.9, rel=1e-4)

    u = param.utpx(T=T, x=0.0, units=units)
    assert u == pytest.approx(112.56, rel=1e-4)
    u = param.utpx(T=T, x=1.0, units=units)
    assert u == pytest.approx(2411.6, rel=1e-4)

    u = param.utpx(T=T, p=3.537*pyo.units.kPa, units=units)
    assert u == pytest.approx(112.56, rel=1e-4)
    u = param.utpx(T=T, p=3.536*pyo.units.kPa, units=units)
    assert u == pytest.approx(2411.64, rel=1e-4)

    units = pyo.units.kJ / pyo.units.kg / pyo.units.K
    s = param.stpx(T=T, x=0.0, units=units)
    assert s == pytest.approx(0.39309, rel=1e-4)
    s = param.stpx(T=T, x=1.0, units=units)
    assert s == pytest.approx(8.5174, rel=1e-4)

@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_htpx_mole():
    """Test mixed phase form with P-H state vars and phase mass balances"""
    # Make 2 param blocks here just to check the default amount basis
    param = HelmholtzParameterBlock(pure_component="h2o")
    param2 = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE
    )
    param.construct()
    param2.construct()

    T = 300 * pyo.units.K
    units = pyo.units.kJ / pyo.units.mol

    h = param.htpx(T=T, x=0.0, with_units=False, units=units)
    assert h == pytest.approx(2.0279, rel=1e-4)
    h = param.htpx(T=T, x=1.0, with_units=False, units=units)
    assert h == pytest.approx(45.936, rel=1e-4)
    h = param.htpx(T=T, x=0.0, with_units=False, units=None)
    assert h == pytest.approx(2.0279*1000, rel=1e-4)
    h = param.htpx(T=T, x=1.0, with_units=False, units=None)
    assert h == pytest.approx(45.936*1000, rel=1e-4)
    h = param.htpx(T=T, x=0.0, with_units=True, units=None)
    assert pyo.value(pyo.units.convert(h, units)) == pytest.approx(2.0279, rel=1e-4)
    h = param.htpx(T=T, x=1.0, with_units=True, units=None)
    assert pyo.value(pyo.units.convert(h, units)) == pytest.approx(45.936, rel=1e-4)

    u = param.utpx(T=T, x=0.0, units=units)
    assert u == pytest.approx(2.0278, rel=1e-4)
    u = param.utpx(T=T, x=1.0, units=units)
    assert u == pytest.approx(43.446, rel=1e-4)
    u = param.utpx(T=T, x=0.0, units=None)
    assert u == pytest.approx(2.0278*1000, rel=1e-4)
    u = param.utpx(T=T, x=1.0, units=None)
    assert u == pytest.approx(43.446*1000, rel=1e-4)

    u = param.utpx(T=T, p=3.537*pyo.units.kPa, units=units)
    assert u == pytest.approx(2.0278, rel=1e-4)
    u = param.utpx(T=T, p=3.536*pyo.units.kPa, units=units)
    assert u == pytest.approx(43.446, rel=1e-4)

    units = pyo.units.kJ / pyo.units.kmol / pyo.units.K
    s = param2.stpx(T=T, x=0.0, units=units)
    assert s == pytest.approx(7.0816, rel=1e-4)
    s = param2.stpx(T=T, x=1.0, units=units)
    assert s == pytest.approx(153.44, rel=1e-4)
    s = param2.stpx(T=T, x=0.0)
    assert s == pytest.approx(7.0816, rel=1e-4)
    s = param2.stpx(T=T, x=1.0)
    assert s == pytest.approx(153.44, rel=1e-4)

@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_expression_writter_mass():
    """Test mixed phase form with P-H state vars and phase mass balances"""
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MASS
    )
    te = HelmholtzThermoExpressions(m, m.hparam)
