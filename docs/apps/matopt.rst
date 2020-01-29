======================================================
MatOpt : Material Design via Mathematical Optimization
======================================================

The MatOpt (Material design via Mathematical Optimization) module provides simple tools for
creating Pyomo model objects for nanomaterial optimization-based design problems.
There are two main sub-modules contained in the package targeting different purposes.
The matopt.materials module contains objects and methods for efficiently representing and manipulating a nanomaterial and its design space.
The matopt.opt module contains objects and methods for speeding up the formulation of a Mixed-integer Linear Programming (MILP) model with simplified
modeling syntax and automatic model formulation. The main goals of this package are as follows.

* To simplify the representation of nanostructured materials, streamlining the creation of materials optimization problems.
* To provide a simple interface so that users do not need to understand the details of building mathematical optimization models or the syntax of the Pyomo package.
* To automate many of the necessary steps of materials optimization, speeding up the development of new models.

.. warning::
   MatOpt depends on Pyomo and Numpy. User's access to MILP solvers (*i.e.*, CPLEX, Gurobi, GLPK) through Pyomo is assumed

Basic Usage
-----------

**Define design canvas**

For formulating a material optimization problem, several pieces of information about the material and design space need to be specified.
To fulfill the need, matopt.materials module defines generic and simple objects for describing the type of material to be designed and design space (canvas) of
the intended material

Key objects

.. module:: apps.matopt.materials.lattices.lattice

.. autoclass:: Lattice

.. module:: apps.matopt.materials.canvas

.. autoclass:: Canvas

.. module:: apps.matopt.materials.atom

.. autoclass:: Atom

.. module:: apps.matopt.materials.design

.. autoclass:: Design


For a detailed description of the material representation system, please see **PAPER**

**Build model via descriptors**

Material type and design space specified provide indexes, sets, and parameters for an optimization model.
Users could in principle formulate a mathematical optimization model based on that information via and python-based
modeling language. MatOpt can accelerate this process by defining a simplified set of objects (matopt.opt.mat_modeling) for users to interact with.
With simplified syntax, users define a MatOptModel object, which will be translated into a pyomo model object automatically by MatOpt (matopt.opt.pyomo_modeling).

MatOpt uses MaterialDescriptor objects to represent variables, constraints, and objectives into
A MatOptModel object holds lists of MaterialDescriptor objects. By default, several universal site descriptors are pre-defined in the model.
From these, all other material descriptors can be defined.

+------------+-------------------------------------------------------------------------------------------+
| Descriptor | Explanation                                                                               |
+============+===========================================================================================+
| Yik        | Presence of a building block of type k at site i                                          |
+------------+-------------------------------------------------------------------------------------------+
| Yi         | Presence of any type of building block at site i                                          |
+------------+-------------------------------------------------------------------------------------------+
| Xijkl      | Presence of a building block of type k at site i and a building block of type l at site j |
+------------+-------------------------------------------------------------------------------------------+
| Xij        | Presence of any building block at site i and any building block at site j                 |
+------------+-------------------------------------------------------------------------------------------+
| Cikl       | Count of neighbors of type l next to a building block of type k at site i                 |
+------------+-------------------------------------------------------------------------------------------+
| Ci         | Count of any type of neighbors next to a building block at site i                         |
+------------+-------------------------------------------------------------------------------------------+

User-specified descriptors are defined by DescriptorRule objects in conjunction with Expr expression objects.
Available expressions include:

+--------------------+-----------------------------------------------------------------------------+
| Expression         | Explanation                                                                 |
+====================+=============================================================================+
| LinearExpr         | Multiplication and addition of coefficients to distinct MaterialDescriptors |
+--------------------+-----------------------------------------------------------------------------+
| SiteCombination    | Summation of site contributions from two sites                              |
+--------------------+-----------------------------------------------------------------------------+
| SumNeighborSites   | Summation of site contributions from all neighboring sites                  |
+--------------------+-----------------------------------------------------------------------------+
| SumNeighborBonds   | Summation of bond contributions to all neighboring sites                    |
+--------------------+-----------------------------------------------------------------------------+
| SumSites           | Summation across sites                                                      |
+--------------------+-----------------------------------------------------------------------------+
| SumBonds           | Summation across bonds                                                      |
+--------------------+-----------------------------------------------------------------------------+
| SumSiteTypes       | Summation across site types                                                 |
+--------------------+-----------------------------------------------------------------------------+
| SumBondTypes       | Summation across bond types                                                 |
+--------------------+-----------------------------------------------------------------------------+
| SumSitesAndTypes   | Summation across sites and site types                                       |
+--------------------+-----------------------------------------------------------------------------+
| SumBondsAndTypes   | Summation across bonds and bond types                                       |
+--------------------+-----------------------------------------------------------------------------+
| SumConfs           | Summation across conformation types                                         |
+--------------------+-----------------------------------------------------------------------------+
| SumSitesAndConfs   | Summation across sites and conformation types                               |
+--------------------+-----------------------------------------------------------------------------+

Several types of DescriptorRules are available.

+-------------------------+---------------------------------------------------------------------------------+
| Rule                    | Explanation                                                                     |
+=========================+=================================================================================+
| LessThan                | Descriptor less than or equal to an expression                                  |
+-------------------------+---------------------------------------------------------------------------------+
| EqualTo                 | Descriptor equal to an expression                                               |
+-------------------------+---------------------------------------------------------------------------------+
| GreaterThan             | Descriptor greater than or equal to an expression                               |
+-------------------------+---------------------------------------------------------------------------------+
| FixedTo                 | Descriptor fixed to a scalar value                                              |
+-------------------------+---------------------------------------------------------------------------------+
| PiecewiseLinear         | Descriptor equal to the evaluation of a piecewise linear function               |
+-------------------------+---------------------------------------------------------------------------------+
| Implies                 | Indicator descriptor that imposes other constraints if equal to 1               |
+-------------------------+---------------------------------------------------------------------------------+
| NegImplies              | Indicator descriptor that imposes other constraints if equal to 0               |
+-------------------------+---------------------------------------------------------------------------------+
| ImpliesSiteCombination  | Indicator bond-indexed descriptor that imposes constraints on the two sites     |
+-------------------------+---------------------------------------------------------------------------------+
| ImpliesNeighbors        | Indicator site-indexed descriptor that imposes constraints on neighboring sites |
+-------------------------+---------------------------------------------------------------------------------+

From the combination of pre-defined descriptors, expressions, and rules, we can specify a wide variety of other descriptors.

**Solve optimization model**

Once the model is fully specified, the user can optimize in light of a descriptor.
Several functions are provided for users to choose from.

.. module:: apps.matopt.opt.mat_modeling

.. autoclass:: MatOptModel
   :members: optimize, populate, maximize, minimize

MatOpt Output
-------------
Optimization results will be loaded into design objects automatically.
Users can then write material design(s) into files for further analysis and visualization.
MatOpt provides interfaces to several standard crystal structure file formats, including CFG, PDB, POSCAR, and XYZ (matopt.materials.parsers).

MatOpt Examples
---------------
Five case studies are provided in idaes.docs.app.matopt
In each case, a Jupyter notebook with explanation as well as an equivalent Python script is provided.
