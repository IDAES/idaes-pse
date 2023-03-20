Wassiljewa-Mason-Saxena Method for Low-pressure Gas Mixture Thermal Conductivity
================================================================================

.. index::
   pair: idaes.models.properties.modular_properties.transport_properties.thermal_conductivity_wms; ThermalConductivityWMS

.. currentmodule:: idaes.models.properties.modular_properties.transport_properties.thermal_conductivity_wms

The :code:`ThermalConductivityWMS` mixture gas thermal conductivity model can be imported from :code:`idaes.models.properties.modular_properties.transport_properties`.

Formulation
-----------

This description of the Wassiljewa-Mason-Saxena method to calculate gas mixture thermal conductivity is derived from that in :ref:`The Properties of Gases and Liquids, Section 10-6 <poling-2001>`.
Thermal conductivity is given by the formula

.. math::

  k_m = \sum_{i=1}^n \frac{y_i k_i}{\sum_{j=1}^n y_j \phi_{ij}}

in which :math:`k_m` is the mixture thermal conductivity, :math:`k_i` is the pure component thermal conductivity of component :math:`i`, 
:math:`y_i` is the gas mole fraction of component :math:`i`, and :math:`\phi_{ij}` is a mixing function shared with the 
:ref:`Wilke method for calculating mixture viscosity<explanations/components/property_package/general/transport_properties/viscosity_wilke:Callbacks for Mixing Function>`. 

List of Parameters
------------------
.. csv-table::
   :header: "Parameter Name", "Description", "Units"

   "``therm_cond_phase_comp``", "Pure component thermal conductivity :math:`k_i`", "Power / Length / Time"
   "``viscosity_phi_ij_callback``", "Callback to use for :math:`\phi_{ij}`", "n/a"

Example
-------

The code snippet below demonstrates how to specify use of the :code:`ThermalConductivityWMS` model for pure component vapor viscosity as part of the modular property 
framework.

.. code-block:: python

  from idaes.models.properties.modular_properties.pure import Eucken
  from idaes.models.properties.modular_properties.transport_properties import ThermalConductivityWMS
  from idaes.models.properties.modular_properties.transport_properties.viscosity_wilke import wilke_phi_ij_callback

  configuration = {
    "components":{
      "H2O": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        ...
        "therm_cond_phase_comp": {"Vap": Eucken}
        ...
      }
      "H2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        ...
        "therm_cond_phase_comp": {"Vap": Eucken}
        ...
      }
    }
    "phases":{
      "Vap": {
        "type": VaporPhase,
        ...
        "therm_cond_phase": ThermalConductivityWMS,
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
