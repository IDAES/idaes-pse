Wilke Method for Low-pressure Gas Mixture Viscosity
===================================================

.. index::
   pair: idaes.models.properties.modular_properties.transport_properties.viscosity_wilke; ViscosityWilke

.. currentmodule:: idaes.models.properties.modular_properties.transport_properties.viscosity_wilke

The :code:`ViscosityWilke` mixture gas dynamic viscosity model can be imported from :code:`idaes.models.properties.modular_properties.transport_properties`,
while additional rules and utility functions can be imported from ``idaes.models.properties.modular_properties.transport_properties.viscosity_wilke``.

Formulation
-----------

This description of the Wilke method to calculate gas mixture viscosity is derived from that in :ref:`The Properties of Gases and Liquids, Section 9-5 <poling-2001>`.
Viscosity is given by the formula

.. math::

  \mu_m = \sum_{i=1}^n \frac{y_i \mu_i}{\sum_{j=1}^n y_j \phi_{ij}}

in which :math:`\mu_m` is the mixture dynamic viscosity, :math:`\mu_i` is the pure component viscosity of component :math:`i`, 
:math:`y_i` is the gas mole fraction of component :math:`i`, and :math:`\phi_{ij}` is a mixing function defined below. 


Callbacks for Mixing Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These callbacks provide expressions for the mixing function :math:`\phi_{ij}`.

.. autofunction:: wilke_phi_ij_callback


.. autofunction:: herring_zimmer_phi_ij_callback


List of Parameters
------------------
.. csv-table::
   :header: "Parameter Name", "Description", "Units"

   "``mw``", "Molecular weight :math:`M_i`", "Mass/Amount"
   "``visc_d_phase_comp``", "Pure component viscosity :math:`\mu_i`", "Pressure\*Time"
   "``viscosity_phi_ij_callback``", "Callback to use for :math:`\phi_{ij}`", "n/a"

Example
-------

The code snippet below demonstrates how to specify use of the :code:`ViscosityWilke` model for pure component vapor viscosity as part of the modular property 
framework.

.. code-block:: python

  from idaes.models.properties.modular_properties.pure import ChapmanEnskogLennardJones
  from idaes.models.properties.modular_properties.transport_properties import ViscosityWilke
  from idaes.models.properties.modular_properties.transport_properties.viscosity_wilke import wilke_phi_ij_callback, herring_zimmer_phi_ij_callback

  configuration = {
    "components":{
      "H2O": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        ...
        "visc_d_phase_comp": {"Vap": ChapmanEnskogLennardJones}
        ...
      }
      "H2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        ...
        "visc_d_phase_comp": {"Vap": ChapmanEnskogLennardJones}
        ...
      }
    }
    "phases":{
      "Vap": {
        "type": VaporPhase,
        ...
        "visc_d_phase": ViscosityWilke,
        "transport_property_options": {
          "viscosity_phi_ij_callback": wilke_phi_ij_callback,
        }
      },
      ...
    }
    ...
  }

References
----------

.. _poling-2001:

Poling, Bruce,  E. et al. *The Properties of Gases and Liquids*. 5th ed. New York: NcGraw-Hill, 2001. 
