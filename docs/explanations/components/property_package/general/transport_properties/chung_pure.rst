Chung et al. Viscosity Model
============================

.. index::
   pair: idaes.models.properties.modular_properties.pure.ChungPure; ChungViscosityPure

.. currentmodule:: idaes.models.properties.modular_properties.pure.ChungPure

The :code:`ChungViscosityPure` pure component gas dynamic viscosity model can be imported from :code:`idaes.models.properties.modular_properties.pure`.

Formulation
-----------

The description of the Chung et al. model for viscosity here is derived from that in :ref:`The Properties of Gases and Liquids, Section 9-4-2 <poling-2001>`.
Pure component gas viscosity is given by the formula:

.. math::

  \mu = C\frac{F_c \sqrt{MT}}{V_c^{2/3} \Omega_v(T)}

in which :math:`\mu` is the dynamic viscosity, :math:`M` is the molar mass, :math:`V_c` is the critical molar volume, 
and :math:`\Omega_v(T)` is a dimensionless quantity known as the collision integral. :math:`C` is a constant with value 
:math:`40.785\;\mu\text{P}(\text{mL}/\text{mol})^{2/3}/\sqrt{(\text{g}/\text{mol})\text{K}}`, and 

.. math::

    F_c = 1 - 0.2756 \omega + 0.059035 p_r^4 + \kappa

in which :math:`\omega` is the acentric factor, :math:`p_r` is the reduced molecular dipole moment, and :math:`\kappa` 
is a correction factor, called the association factor, given for certain polar molecules. Values for :math:`\kappa` can be found in 
:ref:`The Properties of Gases and Liquids, Table 9-1 <poling-2001>`, which in turn are taken from :ref:`Chung et al. (1988) <chung-1988>`.
The reduced dipole moment is given by:

.. math::

    p_r = K \frac{p}{\sqrt{V_cT_c}}

in which :math:`p` is the (dimensionful) dipole moment and :math:`K` is a constant with value :math:`131.3\;\sqrt{(\text{g}/\text{mol})\text{K}}/\text{debye}`.

The collision integral is specific for the property of viscosity (e.g., there is a different collision integral for diffusivity calculations)
and is given in terms of dimensionless temperature. For the purpose of this correlation, the dimensionless temperature used is a scaled
version of the reduced temperature:

.. math::

  T^* = 1.2593\frac{T}{T_c}

in which :math:`T_c` is the critical temperature. Callbacks for the collision integral are described 
:ref:`here<explanations/components/property_package/general/transport_properties/chapman_enskog:Defined Collision Integral Callbacks>`.

List of Parameters
------------------
.. csv-table::
   :header: "Parameter Name", "Description", "Units"

   "``mw``", "Molecular weight :math:`M`", "Mass/Amount"
   "``temperature_crit``", "Critical Temperature :math:`T_c`", "Temperature"
   "``dens_mol_crit``", "Inverse critical molar volume :math:`1/V_c`", "Amount/Volume"
   "``omega``", "Acentric factor :math:`\omega`", "Dimensionless"
   "``dipole_moment``", "Molecular dipole moment :math:`p`", "Charge/Length"
   "``association_factor_chung``", Association factor :math:`\kappa`, "Dimensionless"
   "``viscosity_collision_integral_callback``", "Callback to use for viscosity integral", "n/a"

Example
-------

The code snippet below demonstrates how to specify use of the :code:`ChungViscosityPure` model for pure component vapor viscosity as part of the modular property 
framework. Note that if you specify :code:`visc_d_phase_comp` for one phase, you must specify it for all phases, even if only to pass :code:`None` as the
method.

.. code-block:: python

  from idaes.models.properties.modular_properties.pure import ChungViscosityPure
  from idaes.models.properties.modular_properties.pure.ChapmanEnskog import collision_integral_neufeld_callback, collision_integral_kim_ross_callback

  configuration = {
    "components":{
      "H2O": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase, PhaseType.liquidPhase],
        ...
        "visc_d_phase_comp": {"Vap": ChungViscosityPure, "Liq": None},
        "viscosity_collision_integral_callback": collision_integral_neufeld_callback,
        "parameter_data": {
          "mw": (0.01801528, pyunits.kg / pyunits.mol),
          "temperature_crit": (647.3, pyunits.K),
          "dens_mol_crit": (0.01787, pyunits.mol/pyunits.mL),
          "omega": 0.344,
          "dipole_moment": (1.8546, pyunits.debye),
          "association_factor_chung": 0.076
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

.. _chung-1988:

Chung, Ting Horng, et al. "Generalized multiparameter correlation for nonpolar and polar fluid transport properties." Industrial & Engineering Chemistry Research 27.4 (1988): 671-679.
