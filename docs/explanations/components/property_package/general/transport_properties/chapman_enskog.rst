Chapman Enskog Viscosity Model
==============================

.. index::
   pair: idaes.models.properties.modular_properties.pure.ChapmanEnskog; ChapmanEnskogLennardJones

.. currentmodule:: idaes.models.properties.modular_properties.pure.ChapmanEnskog

The :code:`ChapmanEnskogLennardJones` pure component gas dynamic viscosity model can be imported from :code:`idaes.models.properties.modular_properties.pure`,
while additional rules and utility functions can be imported from ``idaes.models.property_models.modular_properties.pure.ChapmanEnskog``.

Formulation
-----------

The description of Chapman Enskog theory here is derived from that in :ref:`The Properties of Gases and Liquids, Section 9-3 <poling-2001>`.
Pure component gas viscosity is given by the formula

.. math::

  \mu = \frac{C \sqrt{MT}}{\sigma^2 \Omega_v(T)}

in which :math:`\mu` is the dynamic viscosity, :math:`M` is the molar mass, :math:`\sigma` is the collision diameter of the molecule 
(equal to the molecule's diameter for hard spheres and typically taken to be the Lennard-Jones parameter :math:`\sigma` for real molecules), 
and :math:`\Omega_v(T)` is a dimensionless quantity known as the collision integral (equal to 1 for hard spheres but a function of temperature 
for real molecules). :math:`C` is a constant with value :math:`26.69\;\mu\text{PÅ}^2/\sqrt{(\text{g}/\text{mol})\text{K}}`.

The collision integral is specific for the property of viscosity (e.g., there is a different collision integral for diffusivity calculations)
and is given in terms of the dimensionless temperature

.. math::

  T^* = \frac{k_BT}{\varepsilon}

in which :math:`k_B` is the Boltzmann constant and :math:`\varepsilon` is the Lennard-Jones well depth. For good viscosity predictions, 
**it is vital that** :math:`\sigma` **and** :math:`\varepsilon` **come from the same source**, preferably from fitting experimental
viscosity data (see :ref:`Reichenberg (1973) <reichenberg-1973>` for more details). Several different forms of the collision integral, 
of varying complexity, are available to use. If no callback is specified, the form proposed by Neufeld is used.

Defined Collision Integral Callbacks 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These callbacks provide expressions for the viscosity collision integral.

.. autofunction:: collision_integral_neufeld_callback


.. autofunction:: collision_integral_kim_ross_callback


List of Parameters
------------------
.. csv-table::
   :header: "Parameter Name", "Description", "Units"

   "``mw``", "Molecular weight :math:`M`", "Mass/Amount"
   "``lennard_jones_sigma``", "Lennard-Jones 'particle size' :math:`\sigma`", "Length"
   "``lennard_jones_epsilon_reduced``", "Reduced Lennard-Jones well depth :math:`\varepsilon /k_B`", "Temperature"
   "``viscosity_collision_integral_callback``", "Callback to use for viscosity integral", "n/a"

Example
-------

The code snippet below demonstrates how to specify use of the :code:`ChapmanEnskogLennardJones` model for pure component vapor viscosity as part of the modular property 
framework. Note that if you specify :code:`visc_d_phase_comp` for one phase, you must specify it for all phases, even if only to pass :code:`None` as the
method.

.. code-block:: python

  from idaes.models.properties.modular_properties.pure import ChapmanEnskogLennardJones
  from idaes.models.properties.modular_properties.pure.ChapmanEnskog import collision_integral_neufeld_callback, collision_integral_kim_ross_callback

  configuration = {
    "components":{
      "H2O": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase, PhaseType.liquidPhase],
        ...
        "visc_d_phase_comp": {"Vap": ChapmanEnskogLennardJones, "Liq": None},
        "viscosity_collision_integral_callback": collision_integral_neufeld_callback,
        "parameter_data": {
          "mw": (0.01801528, pyunits.kg / pyunits.mol),
          "lennard_jones_sigma": (2.641, pyunits.angstrom),
          "lennard_jones_epsilon_reduced": (809.1, pyunits.K),
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

.. _reichenberg-1973:

Reichenberg, Daniel. "The indeterminacy of the values of potential parameters as derived from transport and virial coefficients." AIChE Journal 19.4 (1973): 854-856.
