=======================================================
MatOpt : Materials Design via Mathematical Optimization
=======================================================

The MatOpt (Materials Design via Mathematical Optimization) module provides simple tools for creating Pyomo model objects for optimization-based nanomaterials design. MatOpt can be used to design crystalline nanostructured materials, including but not limited to particles, surfaces, and periodic bulk structures.

The main goals of this package are as follows:

* To simplify the representation of nanostructured materials, streamlining the creation of materials optimization problems.
* To automate many of the necessary steps of materials optimization, speeding up the development of new models and accelerating new materials discovery.
* To provide a simple interface so that users do not need to handle the details of building efficient mathematical optimization models or the specific Pyomo syntax to do this.

Thank you for your interest in MatOpt. We would love to hear your feedback! Please report any thoughts, questions or bugs to: gounaris@cmu.edu

If you are using MatOpt, please consider citing:

* Hanselman, C.L., Yin, X., Miller, D.C. and Gounaris, C.E., 2020. MatOpt: A Python package for nanomaterials design using discrete optimization. *In preparation*.


Basic Usage
-----------

There are two main sub-modules contained in the package serving two disctinct purposes:

* The ``matopt.materials`` module contains objects and methods for efficiently representing and manipulating a nanomaterial and its design space.
* The ``matopt.opt`` module contains objects and methods for speeding up the casting of a Mixed-integer Linear Programming (MILP) model with simplified modeling syntax and automatic model formulation.

.. warning::
   MatOpt depends on Pyomo and Numpy. User access to the MILP solver CPLEX through Pyomo is assumed. For users who do not have access to CPLEX, the use of NEOS-CPLEX is suggested as an alternative.

**Define design canvas**

Several pieces of information about the material and design space need to be specified in order to formulate a materials optimization problem. To fulfill this need, the ``matopt.materials`` module defines generic and simple objects for describing the type of material to be designed and its design space, also referred to as a "canvas".

Some key objects are listed as follows:

.. module:: apps.matopt.materials.lattices.lattice

.. autoclass:: Lattice

.. module:: apps.matopt.materials.canvas

.. autoclass:: Canvas

.. module:: apps.matopt.materials.design

.. autoclass:: Design

**Build model via descriptors**

The material type and design space specified provide indices, sets, and parameters for the optimization model. Using simple syntax, inspired by materials-related terminology, MatOpt users define a ``MatOptModel`` object, which will be translated into a Pyomo ``ConcreteModel`` object automatically.

MatOpt uses ``MaterialDescriptor`` objects to represent variables, constraints, and objectives. A ``MatOptModel`` object holds lists of ``MaterialDescriptor`` objects. By default, several universal site descriptors are pre-defined in the model.

+----------------+-------------------------------------------------------------------------------------------+
| Descriptor     | Explanation                                                                               |
+================+===========================================================================================+
| ``Yik``        | Presence of a building block of type k at site i                                          |
+----------------+-------------------------------------------------------------------------------------------+
| ``Yi``         | Presence of any type of building block at site i                                          |
+----------------+-------------------------------------------------------------------------------------------+
| ``Xijkl``      | Presence of a building block of type k at site i and a building block of type l at site j |
+----------------+-------------------------------------------------------------------------------------------+
| ``Xij``        | Presence of any building block at site i and any building block at site j                 |
+----------------+-------------------------------------------------------------------------------------------+
| ``Cikl``       | Count of neighbors of type l next to a building block of type k at site i                 |
+----------------+-------------------------------------------------------------------------------------------+
| ``Ci``         | Count of any type of neighbors next to a building block at site i                         |
+----------------+-------------------------------------------------------------------------------------------+

User-specified descriptors are defined by ``DescriptorRule`` objects in conjunction with ``Expr`` expression objects. Available expressions include:

+------------------------+-----------------------------------------------------------------------------+
| Expression             | Explanation                                                                 |
+========================+=============================================================================+
| ``LinearExpr``         | Multiplication and addition of coefficients to distinct descriptors         |
+------------------------+-----------------------------------------------------------------------------+
| ``SiteCombination``    | Summation of site contributions from two sites                              |
+------------------------+-----------------------------------------------------------------------------+
| ``SumNeighborSites``   | Summation of site contributions from all neighboring sites                  |
+------------------------+-----------------------------------------------------------------------------+
| ``SumNeighborBonds``   | Summation of bond contributions to all neighboring sites                    |
+------------------------+-----------------------------------------------------------------------------+
| ``SumSites``           | Summation across sites                                                      |
+------------------------+-----------------------------------------------------------------------------+
| ``SumBonds``           | Summation across bonds                                                      |
+------------------------+-----------------------------------------------------------------------------+
| ``SumSiteTypes``       | Summation across site types                                                 |
+------------------------+-----------------------------------------------------------------------------+
| ``SumBondTypes``       | Summation across bond types                                                 |
+------------------------+-----------------------------------------------------------------------------+
| ``SumSitesAndTypes``   | Summation across sites and site types                                       |
+------------------------+-----------------------------------------------------------------------------+
| ``SumBondsAndTypes``   | Summation across bonds and bond types                                       |
+------------------------+-----------------------------------------------------------------------------+
| ``SumConfs``           | Summation across conformation types                                         |
+------------------------+-----------------------------------------------------------------------------+
| ``SumSitesAndConfs``   | Summation across sites and conformation types                               |
+------------------------+-----------------------------------------------------------------------------+

Several types of ``DescriptorRules`` are available.

+-----------------------------+---------------------------------------------------------------------------------+
| Rule                        | Explanation                                                                     |
+=============================+=================================================================================+
| ``LessThan``                | Descriptor less than or equal to an expression                                  |
+-----------------------------+---------------------------------------------------------------------------------+
| ``EqualTo``                 | Descriptor equal to an expression                                               |
+-----------------------------+---------------------------------------------------------------------------------+
| ``GreaterThan``             | Descriptor greater than or equal to an expression                               |
+-----------------------------+---------------------------------------------------------------------------------+
| ``FixedTo``                 | Descriptor fixed to a scalar value                                              |
+-----------------------------+---------------------------------------------------------------------------------+
| ``PiecewiseLinear``         | Descriptor equal to the evaluation of a piecewise linear function               |
+-----------------------------+---------------------------------------------------------------------------------+
| ``Implies``                 | Indicator descriptor that imposes other constraints if equal to 1               |
+-----------------------------+---------------------------------------------------------------------------------+
| ``NegImplies``              | Indicator descriptor that imposes other constraints if equal to 0               |
+-----------------------------+---------------------------------------------------------------------------------+
| ``ImpliesSiteCombination``  | Indicator bond-indexed descriptor that imposes constraints on the two sites     |
+-----------------------------+---------------------------------------------------------------------------------+
| ``ImpliesNeighbors``        | Indicator site-indexed descriptor that imposes constraints on neighboring sites |
+-----------------------------+---------------------------------------------------------------------------------+

From the combination of the above pre-defined descriptors, expressions, and rules, a user can specify a wide variety of other descriptors, as necessary.

.. currentmodule:: apps.matopt.opt.mat_modeling

.. autoclass:: MaterialDescriptor

**Solve optimization model**

Once the model is fully specified, the user can optimize it in light of a chosen descriptor to serve as the objective to be maximized or minimized, as appropriate. Several functions are provided for users to choose from.

.. module:: apps.matopt.opt.mat_modeling

.. autoclass:: MatOptModel
   :members: optimize, populate, maximize, minimize

MatOpt Output
-------------
The results of the optimization process will be loaded into ``Design`` objects automatically. Users can then save material design(s) into files for further analysis and visualization using suitable functions provided. MatOpt provides interfaces to several standard crystal structure file formats, including CFG, PDB, POSCAR, and XYZ.

MatOpt Examples
---------------
Five `case studies
<https://github.com/IDAES/examples-dev/tree/master/src/matopt>`_ are provided to illustrate the detailed usage of MatOpt. In each case, a Jupyter notebook with explanations as well as an equivalent Python script is provided.

References
----------
* Hanselman, C.L. and Gounaris, C.E., 2016. A mathematical optimization framework for the design of nanopatterned surfaces. *AIChE Journal*, 62(9), pp.3250-3263.
* Hanselman, C.L., Alfonso, D.R., Lekse, J.W., Matranga, C., Miller, D.C. and Gounaris, C.E., 2019. A framework for optimizing oxygen vacancy formation in doped perovskites. *Computers & Chemical Engineering*, 126, pp.168-177.
* Hanselman, C.L., Zhong, W., Tran, K., Ulissi, Z.W. and Gounaris, C.E., 2019. Optimization-based design of active and stable nanostructured surfaces. *The Journal of Physical Chemistry C*, 123(48), pp.29209-29218.
* Isenberg, N.M., Taylor, M.G., Yan, Z., Hanselman, C.L., Mpourmpakis, G. and Gounaris, C.E., 2020. Identification of optimally stable nanocluster geometries via mathematical optimization and density-functional theory. *Molecular Systems Design & Engineering*.
* Yin, X., Isenberg, N.M., Hanselman, C.L., Mpourmpakis, G. and Gounaris, C.E., 2020. A mathematical optimization-based design framework for identifying stable bimetallic nanoclusters. *In preparation*.
* Hanselman, C.L., Yin, X., Miller, D.C. and Gounaris, C.E., 2020. MatOpt: A Python package for nanomaterials design using discrete optimization. *In preparation*.
