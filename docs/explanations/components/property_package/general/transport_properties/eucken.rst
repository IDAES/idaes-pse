Eucken Thermal Conductivity Model
=================================

.. index::
   pair: idaes.models.properties.modular_properties.pure.Eucken; Eucken

.. currentmodule:: idaes.models.properties.modular_properties.pure.Eucken

The :code:`Eucken` pure component gas thermal conductivity model can be imported from :code:`idaes.models.properties.modular_properties.pure`.
In order to use it, methods to compute `visc_d_phase_comp` and `cp_mol_ig_comp` are also required.

Formulation
-----------

The version of the Eucken model used here is derived from that in :ref:`The Properties of Gases and Liquids, Section 10-3-1 <poling-2001>`.
Pure component gas thermal conductivity is given by the equation

.. math::

  k = \frac{\mu}{M}\left(f_\text{int}C_p \left(15/4 - 5/2 f_\text{int}\right) \right)

in which :math:`k` is the thermal conductivity, :math:`\mu` is the dynamic viscosity, :math:`M` is the molar mass,
:math:`C_p` is the pure gas constant pressure heat capacity, :math:`R` is the universal gas constant,
and :math:`f_\text{int}` is a dimensionless constant that accounts for thermal conductivity due to a molecule's
internal degrees of freedom. Eucken originally chose :math:`f_\text{int}=1`, but subsequent authors have variously
chosen :math:`f_\text{int}=1.32` or :math:`f_\text{int}=1.15`. No value of :math:`f_\text{int}` is clearly
superior at predicting thermal conductivity.

This equation does not appear in :ref:`The Properties of Gases and Liquids<poling-2001>`, but was derived
from Equation 10-3.2 by the approximation, exact for ideal gases, that :math:`C_p = C_v + R`, in which 
:math:`C_v` is the constant volume heat capacity.

List of Parameters
------------------
.. csv-table::
   :header: "Parameter Name", "Description", "Units"

   "``mw``", "Molecular weight :math:`M`", "Mass/Amount"
   "``f_int_eucken``", "Eucken internal DOF factor :math:`f_\text{int}`", "Dimensionless"
   "``visc_d_phase_comp``", "Method to compute pure gas viscosity :math:`\mu`", "Pressure\*Time"
   "``cp_mol_ig_comp``", "Method to compute pure ideal gas constant pressure heat capacity", "Energy/(Amount\*Temperature)"

Example
-------

The code snippet below demonstrates how to specify use of the :code:`Eucken` model for pure component gas thermal conductivity 
as part of the modular property framework. Note that if you specify :code:`therm_cond_phase_comp` for one phase, you must specify 
it for all phases, even if only to pass :code:`None` as the method.

.. code-block:: python

  from idaes.models.properties.modular_properties.pure import Eucken, ChapmanEnskogLennardJones, NIST

  configuration = {
    "components":{
      "H2O": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase, PhaseType.liquidPhase],
          ...
        "cp_mol_ig_comp": NIST,
        "visc_d_phase_comp": {"Vap": ChapmanEnskogLennardJones, "Liq": None},
        "therm_cond_phase_comp": {"Vap": Eucken, "Liq": None}
        "parameter_data": {
          "mw": (0.01801528, pyunits.kg / pyunits.mol),
          ...
          "f_int_eucken": 1,
        }
        ...
      }
      ...
    }
    ...
  }

References
----------

.. _poling-2001:

Poling, Bruce,  E. et al. *The Properties of Gases and Liquids*. 5th ed. New York: NcGraw-Hill, 2001. 
