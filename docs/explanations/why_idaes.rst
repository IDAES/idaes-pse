Why IDAES
=========

The National Energy Technology Laboratory's Institute for the Design of Advanced Energy Systems 
(IDAES) is a powerful and versatile computational platform offering next-generation engineering 
capabilities for optimizing the design and operation of innovative chemical process 
and energy systems beyond current constraints on complexity, uncertainty, and scales ranging 
from materials to process to market.

The IDAES Integrated Platform was conceived in 2016 to specifically 
address the gaps between state-of-the-art simulation packages and algebraic modeling languages.

Major strengths of commercial simulation packages are their libraries of unit models and 
thermophysical properties. However, such simulation packages often have difficulty optimizing 
flowsheets and have limited support for incorporating models of non-standard, dynamic units, 
such as solids handling, and uncertainty quantification. On the other hand, AMLs are eminently 
flexible and readily support large-scale optimization, but considerable work is required to 
construct process models, which are often only useful for a one-time application.

The IDAES Integrated Platform represents an innovative approach for the design and 
optimization of chemical and energy processes by integrating an extensible, equation-oriented 
process model library with Pyomo (a Python-based AML). Built  specifically to enable rigorous 
large-scale mathematical optimization, the platform includes capabilities for conceptual 
design, steady-state and dynamic optimization, multi-scale modeling, uncertainty quantification, 
and the automated development of thermodynamic, physical property, and kinetic sub-models from 
experimental data.

Key Features
------------

Open Source
^^^^^^^^^^^

All IDAES Code is completely free and redistributable, the license is available
:ref:`here<explanations/license:License>`. Users are free to modify and redistribute code, and community 
development is encouraged.

Equation Oriented
^^^^^^^^^^^^^^^^^

By using an equation-oriented platform, users gain access to a wide range of highly efficient, 
derivative-based numerical solvers for a wide range of problem types, including support for 
both linear and non-linear problems, ordinary and partial differential equations, and problems 
involving binary and integer variables.

Fully-Featured Programming Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By building off of Python, a fully-featured programming environment, users gain access to a 
wide range of libraries for tools such as data visualization and management.

Extensible
^^^^^^^^^^

The source code for all models and tools is fully-open and visible to the user. This allows 
users to both see and understand what is happening in each model, but also modify and extend 
models to suit their needs.

Flexible Form
^^^^^^^^^^^^^

No single model form is best suited to all applications, thus the IDAES Integrated Platform
is built to provide users with access to a range of different model forms. This allows 
users to easily pick-and-choose from the available model forms to find the one best suited to 
their particular application.

Access to Advanced Capabilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

IDAES aims to provide an integrated platform for development of not just process models but also 
tools for solving and analyzing these problems. The platform supports conceptual 
design, parameter estimation, model predictive control, uncertainty quantification, and 
surrogate modeling.