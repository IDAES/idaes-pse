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
    assert h == pytest.approx(112.56 * 1000, rel=1e-4)
    h = param.htpx(T=T, x=1.0, with_units=False, units=None)
    assert h == pytest.approx(2549.9 * 1000, rel=1e-4)
    h = param.htpx(T=T, x=0.0, with_units=True, units=None)
    assert pyo.value(pyo.units.convert(h, units)) == pytest.approx(112.56, rel=1e-4)
    h = param.htpx(T=T, x=1.0, with_units=True, units=None)
    assert pyo.value(pyo.units.convert(h, units)) == pytest.approx(2549.9, rel=1e-4)

    u = param.utpx(T=T, x=0.0, units=units)
    assert u == pytest.approx(112.56, rel=1e-4)
    u = param.utpx(T=T, x=1.0, units=units)
    assert u == pytest.approx(2411.6, rel=1e-4)

    u = param.utpx(T=T, p=3.537 * pyo.units.kPa, units=units)
    assert u == pytest.approx(112.56, rel=1e-4)
    u = param.utpx(T=T, p=3.536 * pyo.units.kPa, units=units)
    assert u == pytest.approx(2411.64, rel=1e-4)

    units = pyo.units.kJ / pyo.units.kg / pyo.units.K
    s = param.stpx(T=T, x=0.0, units=units)
    assert s == pytest.approx(0.39309, rel=1e-4)
    s = param.stpx(T=T, x=1.0, units=None)
    assert s == pytest.approx(8.5174 * 1000, rel=1e-4)


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
    assert h == pytest.approx(2.0279 * 1000, rel=1e-4)
    h = param.htpx(T=T, x=1.0, with_units=False, units=None)
    assert h == pytest.approx(45.936 * 1000, rel=1e-4)
    h = param.htpx(T=T, x=0.0, with_units=True, units=None)
    assert pyo.value(pyo.units.convert(h, units)) == pytest.approx(2.0279, rel=1e-4)
    h = param.htpx(T=T, x=1.0, with_units=True, units=None)
    assert pyo.value(pyo.units.convert(h, units)) == pytest.approx(45.936, rel=1e-4)

    u = param.utpx(T=T, x=0.0, units=units)
    assert u == pytest.approx(2.0278, rel=1e-4)
    u = param.utpx(T=T, x=1.0, units=units)
    assert u == pytest.approx(43.446, rel=1e-4)
    u = param.utpx(T=T, x=0.0, units=None, with_units=True)
    assert pyo.value(u) == pytest.approx(2.0278 * 1000, rel=1e-4)
    u = param.utpx(T=T, x=1.0, units=None)
    assert u == pytest.approx(43.446 * 1000, rel=1e-4)

    u = param.utpx(T=T, p=3.537 * pyo.units.kPa, units=units)
    assert u == pytest.approx(2.0278, rel=1e-4)
    u = param.utpx(T=T, p=3.536 * pyo.units.kPa, units=None)
    assert u == pytest.approx(43.446 * 1000, rel=1e-4)

    units = pyo.units.kJ / pyo.units.kmol / pyo.units.K
    s = param2.stpx(T=T, x=0.0, units=units)
    assert s == pytest.approx(7.0816, rel=1e-4)
    s = param2.stpx(T=T, x=1.0, units=units)
    assert s == pytest.approx(153.44, rel=1e-4)
    s = param2.stpx(T=T, x=1.0, units=units, with_units=True)
    assert pyo.value(s) == pytest.approx(153.44, rel=1e-4)
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

    assert pytest.approx(0.39309, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.s(
                h=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
        )
    )

    assert pytest.approx(0.39309, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.s_liq(
                h=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
        )
    )

    assert pytest.approx(8.5174, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.s_vap(
                h=2549.9 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
        )
    )

    assert pytest.approx(112.56, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u(
                h=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(112.56, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u_liq(
                h=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(2411.6, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u_vap(
                h=2549.9 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(112.56, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h(
                u=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(112.56, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_liq(
                u=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(2549.9, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_vap(
                u=2411.6 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(112.56, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h(
                s=0.39309 * 1000 * pyo.units.J / pyo.units.kg / pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(112.56, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_liq(
                s=0.39309 * 1000 * pyo.units.J / pyo.units.kg / pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(2549.9, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_vap(
                s=8.5174 * 1000 * pyo.units.J / pyo.units.kg / pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(112.56, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h(x=0.0, p=3.5368 * 1000 * pyo.units.Pa),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(112.56, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_liq(x=0.0, p=3.5368 * 1000 * pyo.units.Pa),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(2549.9, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_vap(x=1.0, p=3.5368 * 1000 * pyo.units.Pa),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(112.56, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u(
                h=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(112.56, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u_liq(
                h=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(2411.6, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u_vap(
                h=2549.9 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(3.5368, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.p(T=300.0 * pyo.units.K, x=0),
            pyo.units.kPa,
        )
    )

    assert pytest.approx(-5.361849, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.g(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(-5.361849, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.g_liq(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(-5.361849, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.g_vap(T=300 * pyo.units.K, x=1),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(-5.3654, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.f(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(-5.3654, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.f_liq(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(-143.57411, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.f_vap(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(1.8078e-02, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.v_mol(T=300 * pyo.units.K, x=0),
            pyo.units.m**3 / pyo.units.kmol,
        )
    )

    assert pytest.approx(1.8078e-02, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.v_mol_liq(T=300 * pyo.units.K, x=0),
            pyo.units.m**3 / pyo.units.kmol,
        )
    )

    assert pytest.approx(704.01, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.v_mol_vap(T=300 * pyo.units.K, x=0),
            pyo.units.m**3 / pyo.units.kmol,
        )
    )

    assert pytest.approx(1.0, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.x(
                h=2549.9 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.dimensionless,
        )
    )

    assert pytest.approx(300.0, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.T(
                h=2549.9 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.K,
        )
    )

    assert pytest.approx(647.096 / 300.0, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.tau(
                h=2549.9 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.dimensionless,
        )
    )

    assert pytest.approx(3.09476, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.delta_liq(T=300 * pyo.units.K, x=0),
            pyo.units.dimensionless,
        )
    )

    assert pytest.approx(7.9471e-05, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.delta_vap(T=300 * pyo.units.K, x=1),
            pyo.units.dimensionless,
        )
    )

    assert pytest.approx(1 / 1.8078e-02, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho_mol_liq(T=300 * pyo.units.K, x=0),
            pyo.units.kmol / pyo.units.m**3,
        )
    )

    assert pytest.approx(1 / 704.01, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho_mol_vap(T=300 * pyo.units.K, x=1),
            pyo.units.kmol / pyo.units.m**3,
        )
    )

    assert pytest.approx(1 / 1.8078e-02 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho_liq(T=300 * pyo.units.K, x=0),
            pyo.units.kg / pyo.units.m**3,
        )
    )

    assert pytest.approx(1 / 704.01 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho_vap(T=300 * pyo.units.K, x=1),
            pyo.units.kg / pyo.units.m**3,
        )
    )

    assert pytest.approx(4.1305, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.cv_liq(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
        )
    )

    assert pytest.approx(1.4422, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.cv_vap(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
        )
    )

    assert pytest.approx(4.1809, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.cp_liq(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
        )
    )

    assert pytest.approx(1.9141, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.cp_vap(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
        )
    )

    assert pytest.approx(1501.4, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.w(T=300 * pyo.units.K, x=0),
            pyo.units.m / pyo.units.s,
        )
    )

    assert pytest.approx(1501.4, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.w_liq(T=300 * pyo.units.K, x=0),
            pyo.units.m / pyo.units.s,
        )
    )

    assert pytest.approx(427.89, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.w_vap(T=300 * pyo.units.K, x=0),
            pyo.units.m / pyo.units.s,
        )
    )


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_expression_writter_mole():
    """Test mixed phase form with P-H state vars and phase mass balances"""
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE
    )
    te = HelmholtzThermoExpressions(m, m.hparam)

    assert pytest.approx(0.39309 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.s(
                h=112.56 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol / pyo.units.K,
        )
    )

    assert pytest.approx(0.39309 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.s_liq(
                h=112.56 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol / pyo.units.K,
        )
    )

    assert pytest.approx(8.5174 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.s_vap(
                h=2549.9 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol / pyo.units.K,
        )
    )

    assert pytest.approx(112.56 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u(
                h=112.56 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(112.56 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u_liq(
                h=112.56 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(2411.6 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u_vap(
                h=2549.9 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(112.56 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h(
                u=112.56 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(112.56 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_liq(
                u=112.56 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(2549.9 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_vap(
                u=2411.6 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(112.56 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h(
                s=0.39309 * 18.015268 * pyo.units.J / pyo.units.mol / pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(112.56 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_liq(
                s=0.39309 * 18.015268 * pyo.units.J / pyo.units.mol / pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(2549.9 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_vap(
                s=8.5174 * 18.015268 * pyo.units.J / pyo.units.mol / pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(112.56 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h(x=0.0, p=3.5368 * 1000 * pyo.units.Pa),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(112.56 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_liq(x=0.0, p=3.5368 * 1000 * pyo.units.Pa),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(2549.9 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_vap(x=1.0, p=3.5368 * 1000 * pyo.units.Pa),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(112.56 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u(
                h=112.56 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(112.56 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u_liq(
                h=112.56 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(2411.6 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u_vap(
                h=2549.9 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(3.5368, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.p(T=300.0 * pyo.units.K, x=0),
            pyo.units.kPa,
        )
    )

    assert pytest.approx(-5.361849 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.g(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(-5.361849 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.g_liq(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(-5.361849 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.g_vap(T=300 * pyo.units.K, x=1),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(-5.3654 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.f(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(-5.3654 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.f_liq(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(-143.57411 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.f_vap(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kmol,
        )
    )

    assert pytest.approx(1.8078e-02, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.v_mol(T=300 * pyo.units.K, x=0),
            pyo.units.m**3 / pyo.units.kmol,
        )
    )

    assert pytest.approx(1.8078e-02, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.v_mol_liq(T=300 * pyo.units.K, x=0),
            pyo.units.m**3 / pyo.units.kmol,
        )
    )

    assert pytest.approx(704.01, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.v_mol_vap(T=300 * pyo.units.K, x=0),
            pyo.units.m**3 / pyo.units.kmol,
        )
    )

    assert pytest.approx(1.0, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.x(
                h=2549.9 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.dimensionless,
        )
    )

    assert pytest.approx(300.0, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.T(
                h=2549.9 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.K,
        )
    )

    assert pytest.approx(647.096 / 300.0, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.tau(
                h=2549.9 * 18.015268 * pyo.units.J / pyo.units.mol,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.dimensionless,
        )
    )

    assert pytest.approx(3.09476, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.delta_liq(T=300 * pyo.units.K, x=0),
            pyo.units.dimensionless,
        )
    )

    assert pytest.approx(7.9471e-05, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.delta_vap(T=300 * pyo.units.K, x=1),
            pyo.units.dimensionless,
        )
    )

    assert pytest.approx(1 / 1.8078e-02, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho_mol_liq(T=300 * pyo.units.K, x=0),
            pyo.units.kmol / pyo.units.m**3,
        )
    )

    assert pytest.approx(1 / 704.01, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho_mol_vap(T=300 * pyo.units.K, x=1),
            pyo.units.kmol / pyo.units.m**3,
        )
    )

    assert pytest.approx(1 / 1.8078e-02 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho_liq(T=300 * pyo.units.K, x=0),
            pyo.units.kg / pyo.units.m**3,
        )
    )

    assert pytest.approx(1 / 704.01 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho_vap(T=300 * pyo.units.K, x=1),
            pyo.units.kg / pyo.units.m**3,
        )
    )

    assert pytest.approx(4.1305 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.cv_liq(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kmol / pyo.units.K,
        )
    )

    assert pytest.approx(1.4422 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.cv_vap(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kmol / pyo.units.K,
        )
    )

    assert pytest.approx(4.1809 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.cp_liq(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kmol / pyo.units.K,
        )
    )

    assert pytest.approx(1.9141 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.cp_vap(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kmol / pyo.units.K,
        )
    )


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_expression_writter_sat():
    """Test mixed phase form with P-H state vars and phase mass balances"""
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MASS
    )
    te = HelmholtzThermoExpressions(m, m.hparam)

    assert pytest.approx(0.0035368, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.p_sat(T=300 * pyo.units.K),
            pyo.units.MPa,
        )
    )

    assert pytest.approx(0.0035368, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.p_sat(tau=m.hparam.temperature_star / 300 / pyo.units.K),
            pyo.units.MPa,
        )
    )

    assert pytest.approx(996.51, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.delta_liq_sat(T=300 * pyo.units.K) * m.hparam.dens_mass_star,
            pyo.units.kg / pyo.units.m**3,
        )
    )

    assert pytest.approx(996.51, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.delta_liq_sat(tau=m.hparam.temperature_star / 300 / pyo.units.K)
            * m.hparam.dens_mass_star,
            pyo.units.kg / pyo.units.m**3,
        )
    )

    assert pytest.approx(0.025590, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.delta_vap_sat(T=300 * pyo.units.K) * m.hparam.dens_mass_star,
            pyo.units.kg / pyo.units.m**3,
        )
    )

    assert pytest.approx(0.025590, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.delta_vap_sat(tau=m.hparam.temperature_star / 300 / pyo.units.K)
            * m.hparam.dens_mass_star,
            pyo.units.kg / pyo.units.m**3,
        )
    )

    assert pytest.approx(300, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.T_sat(p=3536.8 * pyo.units.Pa),
            pyo.units.K,
        )
    )

    assert pytest.approx(300, rel=1e-4) == pyo.value(
        pyo.units.convert(
            m.hparam.temperature_star / te.tau_sat(p=3536.8 * pyo.units.Pa),
            pyo.units.K,
        )
    )


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_h2o_transport():
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MASS
    )
    te = HelmholtzThermoExpressions(m, m.hparam)

    assert pytest.approx(0.00085375, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.viscosity_liq(T=300 * pyo.units.K, x=0),
            pyo.units.Pa * pyo.units.s,
        )
    )

    assert pytest.approx(5.5307e-05, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.viscosity_liq(T=640 * pyo.units.K, x=0),
            pyo.units.Pa * pyo.units.s,
        )
    )

    assert pytest.approx(9.7596e-06, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.viscosity_vap(T=300 * pyo.units.K, x=1),
            pyo.units.Pa * pyo.units.s,
        )
    )

    assert pytest.approx(2.7860e-05, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.viscosity_vap(T=640 * pyo.units.K, x=1),
            pyo.units.Pa * pyo.units.s,
        )
    )


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_co2_transport():
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="co2", amount_basis=AmountBasis.MASS
    )
    te = HelmholtzThermoExpressions(m, m.hparam)

    assert pytest.approx(0.00017249, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.viscosity_liq(T=240 * pyo.units.K, x=0),
            pyo.units.Pa * pyo.units.s,
        )
    )

    assert pytest.approx(5.3192e-05, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.viscosity_liq(T=300 * pyo.units.K, x=0),
            pyo.units.Pa * pyo.units.s,
        )
    )

    assert pytest.approx(1.1062e-05, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.viscosity_vap(T=220 * pyo.units.K, x=1),
            pyo.units.Pa * pyo.units.s,
        )
    )

    assert pytest.approx(2.0811e-05, rel=1e-1) == pyo.value(
        pyo.units.convert(
            te.viscosity_vap(T=300 * pyo.units.K, x=1),
            pyo.units.Pa * pyo.units.s,
        )
    )


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_r134a_transport():
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="r134a", amount_basis=AmountBasis.MASS
    )
    te = HelmholtzThermoExpressions(m, m.hparam)

    assert pytest.approx(1590.7, rel=1e-3) == pyo.value(
        pyo.units.convert(
            te.rho_liq(T=170 * pyo.units.K, x=0),
            pyo.units.kg / pyo.units.m**3,
        )
    )

    assert pytest.approx(0.0021397, rel=3e-1) == pyo.value(
        pyo.units.convert(
            te.viscosity_liq(T=170 * pyo.units.K, x=0),
            pyo.units.Pa * pyo.units.s,
        )
    )

    assert pytest.approx(5.7956e-05, rel=1e-1) == pyo.value(
        pyo.units.convert(
            te.viscosity_liq(T=370 * pyo.units.K, x=0),
            pyo.units.Pa * pyo.units.s,
        )
    )
    assert pytest.approx(6.8353e-06, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.viscosity_vap(T=170 * pyo.units.K, x=1),
            pyo.units.Pa * pyo.units.s,
        )
    )

    assert pytest.approx(2.1336e-05, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.viscosity_vap(T=370 * pyo.units.K, x=1),
            pyo.units.Pa * pyo.units.s,
        )
    )
