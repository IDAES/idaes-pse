Defining State Variables
========================

An important part of defining a set of property calculations is choosing the set of variables which will describe the state of the material. The set of state variables needs to include information on the extensive flow, composition and thermodynamic state of the material. However, there are many ways in which this information can be described, and the best choice of state variables depends on many factors.

Within the IDAES Generic Property Package Framework, the definition of state variables is done using sub-modules which create the necessary variables supporting information for the property package. A state definition sub-module may define any set of state variables the user feel appropriate, but must define the following components as either state variables or function of the state variables:

* `temperature`
* `pressure`
* `mole_frac_phase_comp`
* `phase_frac`

The IDAES Generic Property Package Framework has a library of prebuilt state definition sub-modules for users to use which are listed below.

State Definition Libraries
--------------------------

.. toctree::
    :maxdepth: 1

    state/FTPx
