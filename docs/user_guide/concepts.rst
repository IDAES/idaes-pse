Concepts
========

The IDAES computational platform integrates an extensible, equation-oriented Process Modelling
Framework with advanced solver and computer architectures to enable advanced process systems
engineering capabilities. The platform is based on the Python-based algebraic modeling language
Pyomo, and while not necessary, users may benefit from a basic familiarity with Pyomo
(`link to documentation <https://pyomo.readthedocs.io/en/stable/index.html>`_).

The IDAES computational platform was designed to be modular, and is based on the
block-hierarchical structure shown below:

.. image:: ../_images/IDAES_structure.png

An IDAES process model begins with a Process Flowsheet, which is the canvas on which the
representation of the user’s process will be constructed. Each process consists of a network of
Unit Operations which represent different pieces of equipment within the process (such as
reactors, heater and pumps) and are connected together to form the overall process. Each Unit
Operation in turn is made up of modular components – a Unit model which describes the behavior
and performance of the piece of equipment, a Thermophysical Property Package which represents
the material being processed by the unit operation, and a Reaction Package (if applicable) which
represents any chemical reactions that may occur within the unit. Each of these components can
be further broken down into sub-modules:

* Unit Models consist of a set of material, energy and momentum balance equations which describe how
  material flows through the system, coupled with a set of performance equation which describe
  phenomena such as heat and mass transfer.
* Thermophysical Property Packages (generally) consist of a set of ideal, pure component properties
  for each component, a set of mixing rules and departure functions which describe how the mixture
  properties depend on the ideal properties, and a set of equations describing phase-equilibrium
  phenomena.

At the other end of the spectrum, IDAES process models are designed to be general purpose and
to be applicable to a wide range of modeling activities. By providing access to a wide range of
different numerical solvers and modeling tools, IDAES process models can be applied to a wide
range of different problems, such as:

* process optimization and simulation of both steady-state and dynamic applications,
* data reconciliation,
* parameter estimation and uncertainty quantification,
* optimization under uncertainty, and
* conceptual design (superstructure problems).

Modeling Components
-------------------
The IDAES Process Modeling Framework represents each level within the hierarchy above using
“modeling components”. Each of these components represents a part of the overall model structure
and form the basic building blocks of any IDAES process model. An introduction to each of the
IDAES modeling components can be found
:ref:`here<user_guide/components/index:Components>`.

Model Libraries
---------------
To provide a starting point for modelers in using the process modeling tools, the IDAES Process
Modeling Framework contains a library of models for common unit operations and thermophysical
properties. Modelers can use these models out-of-the-box to model their process applications or
can used these as building blocks for developing their own models. All models within IDAES are
designed to be fully open and extensible, allowing users to inspect and modify them to suit
their needs. Documentation of the available model libraries can be found
:ref:`here<user_guide/model_libraries/index:IDAES Model Libraries>`.

Modeling Extensions
--------------------
The IDAES tools set also provides users with access to a number of cutting edge tools not
directly related to process modeling. These tools are collected under the heading of Modeling
Extensions, and information on them can be found
:ref:`here<user_guide/modeling_extensions/index:Modeling Extensions>`.
