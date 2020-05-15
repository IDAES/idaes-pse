Cubic Equations of State
========================

.. contents:: Contents 
    :depth: 2

Introduction
------------

This module implements a general form of a cubic equation of state which can be used for most cubic-type equations of state. The following forms are currently supported:

* Peng-Robinson
* Soave-Redlich-Kwong

General Cubic Equation of State
-------------------------------
All equations come from "The Properties of Gases and Liquids, 4th Edition" by Reid, Prausnitz and Poling. The general cubic equation of state is represented by the following equations:

.. math:: 0 = Z^3 - (1+B-uB)Z^2 + (A-uB-(u-w)B^2)Z - AB-wB^2-wB^3
.. math:: A = \frac{a_mP}{R^2T^2}
.. math:: B = \frac{b_mP}{RT}

where :math:`Z` is the compressibility factor of the mixture, :math:`a_m` and :math:`b_m` are properties of the mixture and :math:`u` and :math:`w` are parameters which depend on the specific equation of state being used as show in the table below.

.. csv-table::
   :header: "Equation", ":math:`u`", ":math:`w`", ":math:`\Omega_A`", ":math:`\Omega_B`", ":math:`\kappa_j`"

   "Peng-Robinson", "2", "-1", "0.45724", "0.07780", ":math:`(1+(1-T_r^2)(0.37464+1.54226\omega_j-0.26992\omega_j^2))^2`"
   "Soave-Redlich-Kwong", "1", "0", "0.42748", "0.08664", ":math:`(1+(1-T_r^2)(0.48+1.574\omega_j-0.176\omega_j^2))^2`"

The properties :math:`a_m` and :math:`b_m` are calculated from component specific properties :math:`a_j` and :math:`b_j` as shown below:

.. math:: a_j = \frac{\Omega_AR^2T_{c,j}^2}{P_{c, j}}\kappa_j
.. math:: b_j = \frac{\Omega_BRT_{c,j}}{P_{c,j}}
.. math:: a_m = \sum_i{\sum_j{y_iy_j(a_ia_j)^{1/2}(1-k_{ij})}}
.. math:: b_m = \sum_i{y_ib_i}

where :math:`P_{c,j}` and :math:`T_{c,j}` are the component critical pressures and temperatures, :math:`y_j` is the mole fraction of component :math`j`, :math:`k_{ij}` are a set of binary interaction parameters which are specific to the equation of state and :math:`\Omega_A`, :math:`\Omega_B` and :math:`\kappa_j` are taken from the table above. :math:`\omega_j` is the Pitzer acentric factor of each component.

The cubic equation of state is solved for each phase via a call to an external function which automatically identifies the correct root of the cubic and returns the value of :math:`Z` as a function of :math:`A` and :math:`B` along with the first and second partial derivatives.

Property Package Options
------------------------

When using the general cubic equation of state module, users must specify the type of cubic to use. This is done using by providing a `type` option in the `equation_of_state_options` argument in the `Phase` definition, as shown in the example below.

.. code-block:: python

    from idaes.generic_models.properties.core.eos.ceos import Cubic, CubicType

    configuration = {
        "phases": {
            "Liquid": {
                "type": LiquidPhase,
                "equation_of_state": Cubic,
                "equation_of_state_options": {
                    "type": CubicType.PR}}}

Required Parameters
-------------------

Cubic equations of state require the following parameters to be defined:

1. `omega` (Pitzer acentricity factor) needs to be defined for each component (in the `parameter_data` for each component).
2. `kappa` (binary interaction parameters) need to be defined for each component pair in the system. This parameter needs to be defined in the general `parameter_data` argument for the overall property package (as it can be used in multiple phases).

Calculation of Properties
-------------------------

Many thermophysical properties are calculated using an ideal and residual term, such that:

.. math:: p = p^0 + p^r

The residual term is derived from the partial derivatives of the cubic equation of state, whilst the ideal term is determined using pure component properties for the ideal gas phase defined for each component.

Mass Density by Phase
---------------------

The following equation is used for both liquid and vapor phases, where :math:`p` indicates a given phase:

.. math:: \rho_{mass, p} = \rho_{mol, p} \times MW_p

where :math:`MW_p` is the mixture molecular weight of phase :math:`p`.

Molar Density by Phase
----------------------

Molar density is calculated using the following equation

.. math:: \rho_{mol, Vap} = \frac{P}{ZRT}

Molar Enthalpy by Phase
-----------------------

The residual enthalpy term is given by:

.. math:: h_{i}^{r}b_m\sqrt{u^2-4w} = \left(T\frac{da}{dT}-a_m\right)\ln{\left(\frac{2Z+B(u+\sqrt{u^2-4w})}{2Z+B(u-\sqrt{u^2-4w})}\right)} +RT(Z-1)b_m\sqrt{u^2-4w}

.. math:: \frac{da}{dT}\sqrt{T} = -\frac{R}{2}\sqrt{\Omega_A}\sum_i{\sum_j{y_iy_j(1-k_{ij})\left(f_{w,j}\sqrt{a_i\frac{T_{c,j}}{P_{c,j}}}+f_{w,i}\sqrt{a_j\frac{T_{c,i}}{P_{c,i}}}\right)}}

The ideal component is calculated from the weighted sum of the (ideal) component molar enthalpies.

Component Molar Enthalpy by Phase
---------------------------------

Component molar enthalpies by phase are calculated using the pure component method provided by the users in the property package configuration arguments.

Molar Entropy by Phase
-----------------------

The residual entropy term is given by:

.. math:: s_{i}^{r}b_m\sqrt{u^2-4w} = R\ln{\frac{Z-B}{Z}}b_m\sqrt{u^2-4w} + R\ln{\frac{ZP^{ref}}{P}}b_m\sqrt{u^2-4w} + \frac{da}{dT}\ln{\left(\frac{2Z+B(u+\sqrt{u^2-4w})}{2Z+B(u-\sqrt{u^2-4w})}\right)}

The ideal component is calculated from the weighted sum of the (ideal) components molar enthalpies.

Component Molar Entropy by Phase
--------------------------------

Component molar entropies by phase are calculated using the pure component methods provided by the users in the property package configuration arguments.

Component Fugacity by Phase
---------------------------

Fugacity is calculated from the system pressure and fugacity coefficients as follows:

.. math :: f_{i, p} = \phi_{i, p} P

Component Fugacity Coefficient by Phase
---------------------------------------

The fugacity coefficient is calculated from the departure function of the cubic equation of state as shown below:

.. math:: \ln{\phi_{i}} = \frac{b_i}{b_m}(Z-1) - \ln{(Z-B)} + \frac{A}{B\sqrt{u^2-4w}}\left(\frac{b_i}{b_m}-\delta_i\right)\ln{\left(\frac{2Z+B(u+\sqrt{u^2-4w})}{2Z+B(u-\sqrt{u^2-4w})}\right)}

.. math:: \delta_i = \frac{2a_i^{1/2}}{a_m} \sum_j{x_j a_j^{1/2}(1-k_{ij})}

Molar Gibbs Energy by Phase
---------------------------

For both liquid and vapor phases, the molar Gibbs energy is calculated as the weighted sum of the component molar Gibbs energies for the given phase:

.. math:: g_{mol, p} = \sum_j{x_{p, j} \times g_{mol, p, j}}

where :math:`x_{p, j}` is the mole fraction of component :math:`j` in the phase :math:`p`.

Component Gibbs Energy by Phase
-------------------------------

Component molar Gibbs energies are calculated using the definition of Gibbs energy:

.. math:: g_{mol, p, j} = h_{mol, p, j} - s_{mol, p, j} \times T
