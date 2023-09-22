#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

import pytest

import pyomo.environ as pyo

from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    HelmholtzThermoExpressions,
    AmountBasis,
    add_helmholtz_external_functions,
    helmholtz_available as available,
    HelmholtzEoSInitializer,
)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.initialization.initializer_base import (
    InitializationStatus,
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
def test_add_function():
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
def test_expression_writer_mass():
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
            te.s(
                u=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
        )
    )

    assert pytest.approx(0.39309, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.s(
                s=0.39309 * 1000 * pyo.units.J / pyo.units.kg / pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
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
            te.u(
                s=0.39309 * 1000 * pyo.units.J / pyo.units.kg / pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(112.56, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u(
                u=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(112.56, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h(
                h=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
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

    assert pytest.approx(-5.361836, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.g(
                u=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(-5.361836, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.g(
                h=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(-5.361836, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.g(
                s=0.39309 * 1000 * pyo.units.J / pyo.units.kg / pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(-5.3652750407, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.f(
                u=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(3.5368, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.p(
                u=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kPa,
        )
    )

    assert pytest.approx(-5.3652750407, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.f(
                h=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(-5.3652750407, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.f(
                u=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(1.81e-05, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.v(
                u=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
                result_basis=AmountBasis.MOLE,
            ),
            pyo.units.m**3 / pyo.units.mol,
        )
    )

    assert pytest.approx(1.81e-05, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.v(
                h=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
                result_basis=AmountBasis.MOLE,
            ),
            pyo.units.m**3 / pyo.units.mol,
        )
    )

    assert pytest.approx(1.81e-05, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.v(
                s=0.39309 * 1000 * pyo.units.J / pyo.units.kg / pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
                result_basis=AmountBasis.MOLE,
            ),
            pyo.units.m**3 / pyo.units.mol,
        )
    )

    assert pytest.approx(0.001003499, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.v(
                u=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.m**3 / pyo.units.kg,
        )
    )

    assert pytest.approx(0.001003499, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.v(
                h=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.m**3 / pyo.units.kg,
        )
    )

    assert pytest.approx(0.001003499, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.v(
                s=0.39309 * 1000 * pyo.units.J / pyo.units.kg / pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.m**3 / pyo.units.kg,
        )
    )

    assert pytest.approx(0.001003499, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.v(
                T=300 * pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
                x=0,
            ),
            pyo.units.m**3 / pyo.units.kg,
        )
    )

    assert pytest.approx(0.001003499, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.v_liq(
                T=300 * pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
                x=0,
            ),
            pyo.units.m**3 / pyo.units.kg,
        )
    )

    assert pytest.approx(39.07833781228, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.v_vap(
                T=300 * pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
                x=1.0,
            ),
            pyo.units.m**3 / pyo.units.kg,
        )
    )

    assert pytest.approx(300, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.T(
                T=300 * pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
                x=1.0,
            ),
            pyo.units.K,
        )
    )

    assert pytest.approx(300, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.T(
                s=0.39309 * 1000 * pyo.units.J / pyo.units.kg / pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.K,
        )
    )

    assert pytest.approx(2.1569869, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.tau(
                T=300 * pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
                x=1.0,
            ),
            pyo.units.dimensionless,
        )
    )

    assert pytest.approx(300, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.T(
                u=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.K,
        )
    )

    assert pytest.approx(2.1569869, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.tau(
                s=0.39309 * 1000 * pyo.units.J / pyo.units.kg / pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.dimensionless,
        )
    )

    assert pytest.approx(2.1569869, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.tau(
                u=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.dimensionless,
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
            te.v(T=300 * pyo.units.K, x=0, result_basis=AmountBasis.MOLE),
            pyo.units.m**3 / pyo.units.kmol,
        )
    )

    assert pytest.approx(1.8078e-02, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.v_liq(T=300 * pyo.units.K, x=0, result_basis=AmountBasis.MOLE),
            pyo.units.m**3 / pyo.units.kmol,
        )
    )

    assert pytest.approx(704.01, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.v_vap(T=300 * pyo.units.K, x=0, result_basis=AmountBasis.MOLE),
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

    assert pytest.approx(1.0, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.x(
                u=2549.9 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.dimensionless,
        )
    )

    assert pytest.approx(0.19022335, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.x(
                s=549.9 * 1000 * pyo.units.J / pyo.units.kg,
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

    assert pytest.approx(1 / 1.8078e-02, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho(T=300 * pyo.units.K, x=0, result_basis=AmountBasis.MOLE),
            pyo.units.kmol / pyo.units.m**3,
        )
    )

    assert pytest.approx(1 / 1.8078e-02, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho_liq(T=300 * pyo.units.K, x=0, result_basis=AmountBasis.MOLE),
            pyo.units.kmol / pyo.units.m**3,
        )
    )

    assert pytest.approx(1 / 704.01, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho_vap(T=300 * pyo.units.K, x=1, result_basis=AmountBasis.MOLE),
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
            te.cv(T=300 * pyo.units.K, x=0),
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
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

    assert pytest.approx(0.07441200836, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.cv_liq(T=300 * pyo.units.K, x=0, result_basis=AmountBasis.MOLE),
            pyo.units.kJ / pyo.units.mol / pyo.units.K,
        )
    )

    assert pytest.approx(0.025982313, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.cv_vap(T=300 * pyo.units.K, x=0, result_basis=AmountBasis.MOLE),
            pyo.units.kJ / pyo.units.mol / pyo.units.K,
        )
    )

    assert pytest.approx(0.0753202555, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.cp_liq(T=300 * pyo.units.K, x=0, result_basis=AmountBasis.MOLE),
            pyo.units.kJ / pyo.units.mol / pyo.units.K,
        )
    )

    assert pytest.approx(0.034482648393, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.cp_vap(T=300 * pyo.units.K, x=0, result_basis=AmountBasis.MOLE),
            pyo.units.kJ / pyo.units.mol / pyo.units.K,
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

    assert pytest.approx(1501.3549, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.w(
                h=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.m / pyo.units.s,
        )
    )

    assert pytest.approx(1501.3549, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.w(
                u=112.56 * 1000 * pyo.units.J / pyo.units.kg,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.m / pyo.units.s,
        )
    )

    assert pytest.approx(1501.3549, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.w(
                s=0.39309 * 1000 * pyo.units.J / pyo.units.kg / pyo.units.K,
                p=3.5368 * 1000 * pyo.units.Pa,
            ),
            pyo.units.m / pyo.units.s,
        )
    )


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_expression_writer_mole():
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
            te.v(T=300 * pyo.units.K, x=0),
            pyo.units.m**3 / pyo.units.kmol,
        )
    )

    assert pytest.approx(1.8078e-02, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.v_liq(T=300 * pyo.units.K, x=0),
            pyo.units.m**3 / pyo.units.kmol,
        )
    )

    assert pytest.approx(704.01, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.v_vap(T=300 * pyo.units.K, x=0),
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

    assert pytest.approx(1 / 1.8078e-02, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho_liq(T=300 * pyo.units.K, x=0),
            pyo.units.kmol / pyo.units.m**3,
        )
    )

    assert pytest.approx(1 / 704.01, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho_vap(T=300 * pyo.units.K, x=1),
            pyo.units.kmol / pyo.units.m**3,
        )
    )

    assert pytest.approx(1 / 1.8078e-02 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho_liq(T=300 * pyo.units.K, x=0, result_basis=AmountBasis.MASS),
            pyo.units.kg / pyo.units.m**3,
        )
    )

    assert pytest.approx(1 / 704.01 * 18.015268, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.rho_vap(T=300 * pyo.units.K, x=1, result_basis=AmountBasis.MASS),
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
def test_expression_writer_sat():
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

    assert pytest.approx(300, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.T_sat(p=3536.8 * pyo.units.Pa),
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
            te.viscosity(T=300 * pyo.units.K, x=0),
            pyo.units.Pa * pyo.units.s,
        )
    )

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

    assert pytest.approx(0.60944, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity(T=300 * pyo.units.K, x=0),
            pyo.units.W / pyo.units.m / pyo.units.K,
        )
    )

    assert pytest.approx(0.60944, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity_liq(T=300 * pyo.units.K, x=0),
            pyo.units.W / pyo.units.m / pyo.units.K,
        )
    )

    assert pytest.approx(0.018563, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity_vap(T=300 * pyo.units.K, x=1),
            pyo.units.W / pyo.units.m / pyo.units.K,
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

    assert pytest.approx(0.011424, rel=1e-1) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity_vap(T=220 * pyo.units.K, x=1),
            pyo.units.W / pyo.units.m / pyo.units.K,
        )
    )


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_r134a_thermo():
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="r134a", amount_basis=AmountBasis.MASS
    )
    te = HelmholtzThermoExpressions(m, m.hparam)

    assert pytest.approx(209.27, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_liq_sat(T=280 * pyo.units.K),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(402.55, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_vap_sat(T=280 * pyo.units.K),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(209.27, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_liq_sat(p=0.37271 * 1e6 * pyo.units.Pa),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(402.55, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.h_vap_sat(p=0.37271 * 1e6 * pyo.units.Pa),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(1.0333, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.s_liq_sat(T=280 * pyo.units.K),
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
        )
    )

    assert pytest.approx(1.7236, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.s_vap_sat(T=280 * pyo.units.K),
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
        )
    )

    assert pytest.approx(1.0333, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.s_liq_sat(p=0.37271 * 1e6 * pyo.units.Pa),
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
        )
    )

    assert pytest.approx(1.7236, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.s_vap_sat(p=0.37271 * 1e6 * pyo.units.Pa),
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
        )
    )

    assert pytest.approx(208.98, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u_liq_sat(T=280 * pyo.units.K),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(382.10, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u_vap_sat(T=280 * pyo.units.K),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(208.98, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u_liq_sat(p=0.37271 * 1e6 * pyo.units.Pa),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(382.10, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.u_vap_sat(p=0.37271 * 1e6 * pyo.units.Pa),
            pyo.units.kJ / pyo.units.kg,
        )
    )

    assert pytest.approx(0.00078629, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.v_liq_sat(T=280 * pyo.units.K),
            pyo.units.m**3 / pyo.units.kg,
        )
    )

    assert pytest.approx(0.054861, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.v_vap_sat(T=280 * pyo.units.K),
            pyo.units.m**3 / pyo.units.kg,
        )
    )

    assert pytest.approx(0.00078629, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.v_liq_sat(p=0.37271 * 1e6 * pyo.units.Pa),
            pyo.units.m**3 / pyo.units.kg,
        )
    )

    assert pytest.approx(0.054861, rel=1e-4) == pyo.value(
        pyo.units.convert(
            te.v_vap_sat(p=0.37271 * 1e6 * pyo.units.Pa),
            pyo.units.m**3 / pyo.units.kg,
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

    assert pytest.approx(0.0021397, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.viscosity_liq(T=170 * pyo.units.K, x=0),
            pyo.units.Pa * pyo.units.s,
        )
    )

    assert pytest.approx(5.7956e-05, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.viscosity_liq(T=370 * pyo.units.K, x=0),
            pyo.units.Pa * pyo.units.s,
        )
    )

    assert pytest.approx(0.14516, rel=1e-3) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity_liq(T=170 * pyo.units.K, x=0),
            pyo.units.W / pyo.units.m / pyo.units.K,
        )
    )

    assert pytest.approx(0.093414, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity_liq(T=270 * pyo.units.K, x=0),
            pyo.units.W / pyo.units.m / pyo.units.K,
        )
    )

    assert pytest.approx(0.011871, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.surface_tension(T=270 * pyo.units.K, x=0),
            pyo.units.N / pyo.units.m,
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

    assert pytest.approx(0.0030921, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity_vap(T=170 * pyo.units.K, x=0),
            pyo.units.W / pyo.units.m / pyo.units.K,
        )
    )

    assert pytest.approx(0.011241, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity_vap(T=270 * pyo.units.K, x=0),
            pyo.units.W / pyo.units.m / pyo.units.K,
        )
    )


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_r1234ze_transport():
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="r1234ze", amount_basis=AmountBasis.MASS
    )
    te = HelmholtzThermoExpressions(m, m.hparam)

    assert pytest.approx(0.0098503, rel=1e-3) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity_vap(
                T=250 * pyo.units.K, p=0.05e6 * pyo.units.Pa, x=1
            ),
            pyo.units.W / pyo.units.m / pyo.units.K,
        )
    )

    assert pytest.approx(0.013933, rel=1e-3) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity_vap(
                T=300 * pyo.units.K, p=0.1e6 * pyo.units.Pa, x=1
            ),
            pyo.units.W / pyo.units.m / pyo.units.K,
        )
    )

    assert pytest.approx(0.10066, rel=1e-3) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity_liq(
                T=250 * pyo.units.K, p=20e6 * pyo.units.Pa, x=0
            ),
            pyo.units.W / pyo.units.m / pyo.units.K,
        )
    )

    assert pytest.approx(0.085389, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity_liq(
                T=300 * pyo.units.K, p=20e6 * pyo.units.Pa, x=0
            ),
            pyo.units.W / pyo.units.m / pyo.units.K,
        )
    )

    # roughly 0 pressure
    assert pytest.approx(11.777, rel=1e-3) == pyo.value(
        pyo.units.convert(
            te.viscosity_vap(T=300 * pyo.units.K, p=1e-3 * pyo.units.Pa, x=1),
            pyo.units.microPa * pyo.units.s,
        )
    )

    assert pytest.approx(12.041, rel=1e-3) == pyo.value(
        pyo.units.convert(
            te.viscosity_vap(T=300 * pyo.units.K, p=1e5 * pyo.units.Pa, x=1),
            pyo.units.microPa * pyo.units.s,
        )
    )

    assert pytest.approx(10.522, rel=1e-3) == pyo.value(
        pyo.units.convert(
            te.rho_liq(
                T=300 * pyo.units.K,
                p=10e6 * pyo.units.Pa,
                x=0,
                result_basis=AmountBasis.MOLE,
            ),
            pyo.units.mol / pyo.units.l,
        )
    )

    assert pytest.approx(217.89, rel=1e-3) == pyo.value(
        pyo.units.convert(
            te.viscosity_liq(T=300 * pyo.units.K, p=10e6 * pyo.units.Pa, x=0),
            pyo.units.microPa * pyo.units.s,
        )
    )


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_propane_transport():
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="propane", amount_basis=AmountBasis.MASS
    )
    te = HelmholtzThermoExpressions(m, m.hparam)

    # Test viscosity with data from
    #
    # Vogel E, C Küchenmeister, E Bich, A Laesecke, Reference Correlation of the
    #     Viscosity of Propane. Journal of Physical and Chemical Reference Data
    #     27, 947–970 (1998)
    assert pytest.approx(8.731, rel=1e-3) == pyo.value(
        pyo.units.convert(
            te.viscosity(
                T=320.0 * pyo.units.K,
                p=1e5 * pyo.units.Pa,
                x=1,
            ),
            pyo.units.uPa * pyo.units.s,
        )
    )

    assert pytest.approx(78.97826, rel=1e-3) == pyo.value(
        pyo.units.convert(
            te.viscosity(
                T=320.0 * pyo.units.K,
                p=2.5e6 * pyo.units.Pa,
                x=0,
            ),
            pyo.units.uPa * pyo.units.s,
        )
    )

    assert pytest.approx(82.16, rel=1e-2) == pyo.value(
        pyo.units.convert(
            te.viscosity(
                T=320.0 * pyo.units.K,
                p=4.0e6 * pyo.units.Pa,
                x=0,
            ),
            pyo.units.uPa * pyo.units.s,
        )
    )

    # Test thermal conductivity with data from
    #
    # NIST Chemistry WebBook
    assert pytest.approx(20.78, rel=1e-3) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity(
                T=320.0 * pyo.units.K,
                p=1e5 * pyo.units.Pa,
                x=1,
            ),
            pyo.units.mW / pyo.units.m / pyo.units.K,
        )
    )

    assert pytest.approx(85.215, rel=1e-3) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity(
                T=320.0 * pyo.units.K,
                p=2.5e6 * pyo.units.Pa,
                x=0,
            ),
            pyo.units.mW / pyo.units.m / pyo.units.K,
        )
    )

    assert pytest.approx(87.159, rel=1e-3) == pyo.value(
        pyo.units.convert(
            te.thermal_conductivity(
                T=320.0 * pyo.units.K,
                p=4.0e6 * pyo.units.Pa,
                x=0,
            ),
            pyo.units.mW / pyo.units.m / pyo.units.K,
        )
    )


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_initialize_param_block():
    # this should do absolutely nothing, so just make sure there is no exception
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="r1234ze", amount_basis=AmountBasis.MASS
    )
    m.hparam.initialize()


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_errors():
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="r1234ze", amount_basis=AmountBasis.MASS
    )
    te = HelmholtzThermoExpressions(m, m.hparam)
    with pytest.raises(TypeError):
        te.p_sat()
    with pytest.raises(RuntimeError):
        te.x(T=100 * pyo.units.K, p=3 * pyo.units.Pa)
    with pytest.raises(RuntimeError):
        te.h(T=100 * pyo.units.K)
    with pytest.raises(RuntimeError):
        te.h(T=100 * pyo.units.K, tacos=3)
    with pytest.raises(RuntimeError):
        te.h(T=100 * pyo.units.K, p=3 * pyo.units.Pa)
    with pytest.raises(RuntimeError):
        te.u(T=100 * pyo.units.K, p=3 * pyo.units.Pa)
    with pytest.raises(RuntimeError):
        te.s(T=100 * pyo.units.K, p=3 * pyo.units.Pa)
    with pytest.raises(RuntimeError):
        te.cp(T=100 * pyo.units.K, p=3 * pyo.units.Pa)
    with pytest.raises(ConfigurationError):
        m.err_param = HelmholtzParameterBlock(
            pure_component="not a real thing", amount_basis=AmountBasis.MASS
        )
    with pytest.raises(RuntimeError):
        # Can't specify all tree T, P, x
        h = m.hparam.htpx(T=300 * pyo.units.K, p=10 * pyo.units.Pa, x=0.0)
    with pytest.raises(RuntimeError):
        # Temperature way too high
        h = m.hparam.htpx(T=300e7 * pyo.units.K, x=0.0)
    with pytest.raises(RuntimeError):
        # Vapor fraction > 1.0
        h = m.hparam.htpx(T=300 * pyo.units.K, x=10.0)
    with pytest.raises(RuntimeError):
        # Pressure way too high
        h = m.hparam.htpx(p=300e17 * pyo.units.Pa, x=0.0)


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_plot_no_exception():
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="r1234ze", amount_basis=AmountBasis.MASS
    )
    m.hparam.ph_diagram(isotherms=True)
    m.hparam.st_diagram()
    m.hparam.pt_diagram()


@pytest.mark.unit
@pytest.mark.skipif(not available(), reason="General Helmholtz not available")
def test_default_initializer():
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="r1234ze", amount_basis=AmountBasis.MASS
    )
    m.state = m.hparam.build_state_block([0])

    assert m.state.default_initializer is HelmholtzEoSInitializer


@pytest.mark.unit
def test_HelmholtzEoSInitializer():
    m = pyo.ConcreteModel()
    m.b = pyo.Block()

    init = HelmholtzEoSInitializer()

    assert init.initialize(m.b) == InitializationStatus.Ok
    assert init.summary[m.b]["status"] == InitializationStatus.Ok
