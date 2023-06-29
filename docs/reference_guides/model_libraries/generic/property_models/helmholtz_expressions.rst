Helmholtz EoS Expression Writer
===============================

.. index::
  pair: idaes.models.properties.general_helmholtz.helmholtz_functions; HelmholtzThermoExpressions

.. module:: idaes.models.properties.general_helmholtz.helmholtz_functions

In addition to the usual IDAES property block, property functions can be accessed directly
as Pyomo expressions using the HelmholtzThermoExpressions class.

Class
-----

.. autoclass:: HelmholtzThermoExpressions
    :members:


All the property expression methods that take indeterminate ``**kwarg`` take
the same set of arguments. The state variables should be provided as arguments and include
units. The state variables are passed to external functions to calculate properties.  If
``convert_units=False`` is passed the state variable expressions are assumed to include 
units native to the external functions and no unit conversion expression is needed;
otherwise, units are assumed to be in SI. While you do need to provide units to check for 
consistency, arbitrary unit conversion is not currently supported, due to potentially slow
model building.  The state variables are on a mass or mole basis as determined by the 
``amount_basis`` argument given when creating the HelmholtzThermoExpressions object. An 
optional ``result_basis`` argument can be provided as an ``AmountBasis`` Enum to set whether 
the resulting property expression is given on a mass or mole basis. If not provided, the 
`amount_basis` argument is used for the result. 

State variables should be in one of the sets below.

* {h, p} 
* {u, p} 
* {s, p}
* {s, T}
* {T, x}
* {p, x}
* {T, P, x}

Saturation properties are given as a function of either temperature or pressure, so only the temperature 
or pressure state variable are required.  The ``convert_units`` and ``result_basis`` arguments are also
taken by the saturated property methods.

Example
-------

The sample code below gives an example of how to use the expression writer in constructing a Pyomo model.

.. testcode::

    import pytest
    import pyomo.environ as pyo
    from idaes.models.properties.general_helmholtz import (
        HelmholtzParameterBlock,
        HelmholtzThermoExpressions,
        AmountBasis,
    )
    m = pyo.ConcreteModel()
    m.hparam = HelmholtzParameterBlock(
        pure_component="r134a", amount_basis=AmountBasis.MASS
    )
    te = HelmholtzThermoExpressions(m, m.hparam)
    m.density = pyo.Expression(expr=te.rho_liq(T=170 * pyo.units.K, x=0))

    assert pytest.approx(1.5907, rel=1e-3) == pyo.value(
        pyo.units.convert(m.density,
            pyo.units.g / pyo.units.cm**3,
        )
    )
