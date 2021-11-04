Plate Heat Exchanger
====================

The thermal model of IDAES Plate Heat Exchanger (PHE) as part of the MEA
scrubbing process for post-combustion carbon capture (PCC) is based on the
Effectiveness-Number of Transfer Units (e-NTU) approach. In amine-based PCC,
the rich solvent leaving the absorber is pre-heated in the PHE using heat recovered
from the lean solvent leaving the stripper to reduce the regeneration energy
requirement. In the PHE unit, the series of plates stacked together form channels
where hot and cold fluids flow alternatively as shown in Figure 1(A). Divider plates
enable the partitioning of PHEs into different operating zones. The main dimensions
of a gasket plate are shown in Figure 1(B). The PHE is a viable alternative to the
conventional Shell and Tube Heat Exchanger specifically because of its lower
approach temperature difference capability. For more information on the PHE model
see `Akula et al. (2019) <https://doi.org/10.1016/B978-0-12-818597-1.50008-4>`_.

.. figure:: ./phe_a.png
  :alt: Z-configuration Plate Heat Exchanger with P passes
  :width: 90%
  :align: center

  **Figure 1(A). Z-configuration Plate Heat Exchanger with P passes**

.. figure:: ./phe_b.png
  :alt: Basic details of a Chevron Plate
  :width: 100%
  :align: center

  **Figure 1(B). Basic details of a Chevron Plate**


Degrees of Freedom
------------------

Once the configuration parameters (construction arguments of the PHE Class)
have been specified, the PHE unit model has 12 degrees of freedom which are the
operating parameters as listed in the Specification Table below. The indexed
components for ``mole_frac_comp`` are ``H2O, MEA and CO2``.


Specification
^^^^^^^^^^^^^

.. csv-table::
   :header: "Variable Name", "Description", "Units"
   :widths: 16,29,8

   "**Configuration parameters**", " "," "
   "``passes``","Number of passes of the fluids", "None "
   "``channel_list``","Number of channels in each pass as a list", "None "
   "``divider_plate_number``","Number of divider plates", "None "
   "``port_diameter``","Diameter of the plate ports (Dp)", ":math:`m` "
   "``plate_thermal_cond``","Thermal conductivity of the plate material", ":math:`W/m.K`"
   "``total_area``","Total heat transfer area as specified by the manufacturer", ":math:`m^{2}`"
   "``plate_thickness``","Plate thickness", ":math:`m`"
   "``plate_vertical_dist``","Vertical distance between centers of ports (Lv)", ":math:`m`"
   "``plate_horizontal_dist``","Horizontal distance between centers of ports (Lh)", ":math:`m`"
   "``plate_pact_length``","Compressed plate pact length (optional)", ":math:`m`"
   "``surface_enlargement_factor``","Ratio of single plate area obtained from the total area
   to the projected plate area (optional)", "None"
   "``plate_gap``","Distance between two adjacent plates that forms a flow channel", ":math:`m` "
   "**Operating parameters**", " ", " "
   "``hot_inlet.flow_mol``", "Hot fluid inlet total molar flowrate", ":math:`mol/s`"
   "``hot_inlet.temperature``", "Hot fluid inlet temperature", ":math:`K`"
   "``hot_inlet.pressure``", "Hot fluid inlet pressure", ":math:`Pa`"
   "``hot_inlet.mole_frac_comp``", "Hot fluid inlet mole fraction indexed by component", "None"
   "``cold_inlet.flow_mol``", "Cold fluid inlet total molar flowrate", ":math:`mol/s`"
   "``cold_inlet.temperature``", "Cold fluid inlet temperature", ":math:`K`"
   "``cold_inlet.pressure``", "Cold fluid inlet pressure", ":math:`Pa`"
   "``cold_inlet.mole_frac_comp``", "Cold fluid inlet mole fraction indexed by component", "None"


Model Structure
---------------

The PHE unit model consists of two
:ref:`ControlVolume0D Blocks <technical_specs/core/control_volume_0d:0D Control Volume Class>`
(named ``hot_side`` and ``cold_side``), each with one Inlet Port (named ``hot_inlet``
and ``cold_inlet``) and one Outlet Port (named ``hot_outlet`` and ``cold_outlet``).
The ``hot_side`` and ``cold_side`` ControlVolume0D Blocks use the
:ref:`Liquid Phase Property Methods <technical_specs/model_libraries/power_generation/carbon_capture/mea_solvent_system/properties/liquid_prop:Liquid Phase Property Methods>` which is built off of the
:ref:`Physical Property Package Class <technical_specs/core/physical_property_class:Physical Property Package Classes>`.
The Energy balance is based on the Effectiveness Number of Transfer Units
(e-NTU method) and is included as performance equations (Additional Constraints).
Hence, the control volume energy balances are not added.


Additional Constraints
----------------------

The PHE unit model writes additional ``Constraints`` beyond those written by the
:ref:`ControlVolume0D Blocks <technical_specs/core/control_volume_0d:0D Control Volume Class>`
to describe the heat exchange between the rich and lean solvent for
post-combustion carbon capture using MEA solvent.



PHE Class
----------

.. module:: idaes.power_generation.carbon_capture.mea_solvent_system.unit_models.phe

.. autoclass:: PHE
  :members:

PHEData Class
--------------

.. autoclass:: PHEData
  :members:

References
------------

1. Akula, P., Eslick, J., Bhattacharyya, D., & Miller, D. C. (2019).
   "Modelling and Parameter Estimation of a Plate Heat Exchanger as Part of a
   Solvent-Based Post-Combustion CO2 Capture System", In Computer Aided Chemical
   Engineering (Vol. 47, pp. 47-52). Elsevier. https://doi.org/10.1016/B978-0-12-818597-1.50008-4
2. Kakac, S., Liu, H., & Pramuanjaroenkij, A. (2012). Heat exchangers:
   selection, rating, and thermal design. CRC press.
