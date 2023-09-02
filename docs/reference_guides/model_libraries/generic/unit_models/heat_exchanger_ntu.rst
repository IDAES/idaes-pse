Heat Exchanger using the NTU Method
===================================

.. index::
   pair: idaes.models.unit_models.heat_exchanger_ntu;HeatExchangerNTU

.. currentmodule:: idaes.models.unit_models.heat_exchanger_ntu

The HeatExchangerNTU model can be imported from :code:`idaes.models.unit_models`, and models a heat exchanger using the effectiveness-NTU method.
This model derived from the

Degrees of Freedom
------------------

Aside from the inlet conditions, an NTU heat exchanger model has three degrees of freedom which must be fixed for it to be fully specified. Additionally, users may include a pressure change in the unit which adds two additional degrees of freedom.

Standard design variables are:

    * heat transfer area,
    * heat transfer coefficient, and
    * effectiveness.

If pressure change is included, additional degrees of freedom are:

    * hot and cold side pressure changes.

Model Structure
---------------

The ``HeatExchanger`` model contains two ``ControlVolume0DBlock`` blocks named ``hot_side`` and the ``cold side``.
These names are configurable using the ``hot_side_name`` and ``cold_side_name`` configuration arguments, in which case
aliases are assigned to the control volumes and associated Ports using the names provided (note that ``hot_side`` and
``cold_side`` will always work).
The sign convention is that duty is positive for heat flowing from the hot side to the cold
side.  Aside from the sign convention there is no requirement that the hot side be hotter
than the cold side.

The ``HeatExchanger`` has two inlet ports and two outlet ports. By default these are
``hot_side_inlet``, ``cold_side_inlet``, ``hot_side_outlet``, and ``cold_side_outlet``. If the user
supplies different hot and cold side names the inlet and outlets are named accordingly.

Variables
---------

=========================== ================== =========== =============================================================================
Variable                    Symbol             Index Sets  Doc
=========================== ================== =========== =============================================================================
heat_duty                   :math:`Q`          t           Heat transferred from hot side to the cold side
area                        :math:`A`          None        Heat transfer area
heat_transfer_coefficient   :math:`U`          t           Heat transfer coefficient
effectiveness               :math:`\epsilon`   t           Effectiveness factor
=========================== ================== =========== =============================================================================

Expressions
-----------

The following Expressions are constructed by the model and can be used in correlations to determine the effectiveness factor.

Minimum heat capacitance:

.. math::
  C_{min} = min((F_{mol, hot} \times c_{p, mol, hot}), (F_{mol, hot} \times c_{p, mol, hot}))

Maximum heat capacitance:

.. math::
  C_{max} = max((F_{mol, hot} \times c_{p, mol, hot}), (F_{mol, hot} \times c_{p, mol, hot}))

Min and max operators are implemented using smooth approximation using the :math:`\epsilon_{cmin}` parameter.

Heat capacitance ratio:

.. math::
  C_{ratio} = \frac{C_{min}}{C_{max}}

Number of theoretical heat transfer units:

.. math::
  NTU = \frac{U \times A}{C_{min}}

Constraints
-----------

The effectiveness-NTU method is a method to approximate the heat transferred in a heat exchanger using the following calculation:

.. math::
  Q_{cold} = \epsilon \times C_{min} \times (T_{hot} - T_{cold})

where :math:`Q_{cold}` is the heat transferred from the hot side to the cold side, :math:`\epsilon` is the effectiveness factor for the heat exchanger, :math:`C_{min}` is the minimum heat capacitance between the hot and cold inlet streams and :math:`T_{hot}` and :math:`T_{cold}` are the temperatures of the hot and cold inlet streams respectively.

Additionally, and overall energy balance constraint is written:

.. math::
  Q_{hot} = -Q_{cold}
  
Initialization
--------------

.. autoclass:: HXNTUInitializer
   :members: initialization_routine

Class Documentation
-------------------

.. autoclass:: HeatExchangerNTU
   :members:

.. autoclass:: HeatExchangerNTUData
   :members:
