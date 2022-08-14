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
    param2 = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE
    )
    param.construct()
    param2.construct()

    T = 300 * pyo.units.K
    units = pyo.units.kJ / pyo.units.kg

    h = param.htpx(T=T, x=0.0, with_units=False, units=units)
    assert h == pytest.approx(112.56, rel=1e-4)
    h = param.htpx(T=T, x=1.0, with_units=False, units=units)
    assert h == pytest.approx(2549.9, rel=1e-4)

    units = pyo.units.kJ / pyo.units.kg / pyo.units.K
    s = param.stpx(T=T, x=0.0, with_units=False, units=units)
    assert s == pytest.approx(0.39309, rel=1e-4)
    s = param.stpx(T=T, x=1.0, with_units=False, units=units)
    assert s == pytest.approx(8.5174, rel=1e-4)


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_expression_writter_mass():
    """Test mixed phase form with P-H state vars and phase mass balances"""
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MASS
    )
    te = HelmholtzThermoExpressions(m, m.hparam)
