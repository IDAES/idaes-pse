==================================
MatOpt: Nanomaterials Optimization
==================================

The MatOpt module provides tools for nanomaterials design using Mathematical Optimization. MatOpt can be used to design crystalline nanostructured materials, including but not limited to particles, wires, surfaces, and periodic bulk structures.

The main goals of this package are as follows:

* To automate many of the steps that are necessary for utilizing mathematical optimization to design materials, speeding up the development of new mathematical models and accelerating new materials discovery.
* To simplify the representation of nanostructured materials and their structure-function relationships as Pyomo objects, streamlining the creation of materials optimization problems in the Pyomo modeling language.
* To provide a simple interface so that users need not handle the details of casting efficient mathematical optimization models, invoking mathematical optimization solvers, or utilizing specialized Pyomo syntax to do this.

Thank you for your interest in MatOpt. We would love to hear your feedback! Please report any thoughts, questions or bugs to: gounaris@cmu.edu

If you are using MatOpt, please consider citing:

* Hanselman, C.L., Yin, X., Miller, D.C. and Gounaris, C.E., 2021. `MatOpt: A Python package for nanomaterials discrete optimization. <http://gounaris.cheme.cmu.edu/drafts/Draft_MATOPT.pdf>`_


MatOpt Examples
---------------
Several `case studies
<https://idaes.github.io/examples-pse/latest/Examples/MatOpt/index.html>`_ are provided to illustrate the detailed usage of MatOpt. In each case, a Jupyter notebook with explanations as well as an equivalent Python script is provided.

References
----------
* Hanselman, C.L. and Gounaris, C.E., 2016. `A mathematical optimization framework for the design of nanopatterned surfaces. <https://aiche.onlinelibrary.wiley.com/doi/full/10.1002/aic.15359>`_ *AIChE Journal*, 62(9), pp.3250-3263.
* Hanselman, C.L., Alfonso, D.R., Lekse, J.W., Matranga, C., Miller, D.C. and Gounaris, C.E., 2019. `A framework for optimizing oxygen vacancy formation in doped perovskites. <https://www.sciencedirect.com/science/article/pii/S0098135418310998>`_ *Computers & Chemical Engineering*, 126, pp.168-177.
* Hanselman, C.L., Zhong, W., Tran, K., Ulissi, Z.W. and Gounaris, C.E., 2019. `Optimization-based design of active and stable nanostructured surfaces. <https://pubs.acs.org/doi/abs/10.1021/acs.jpcc.9b08431>`_ *The Journal of Physical Chemistry C*, 123(48), pp.29209-29218.
* Isenberg, N.M., Taylor, M.G., Yan, Z., Hanselman, C.L., Mpourmpakis, G. and Gounaris, C.E., 2020. `Identification of optimally stable nanocluster geometries via mathematical optimization and density-functional theory. <https://pubs.rsc.org/en/content/articlelanding/2020/me/c9me00108e#!divAbstract>`_ *Molecular Systems Design & Engineering*, 5, pp.232-244.
* Yin, X., Isenberg, N.M., Hanselman, C.L., Mpourmpakis, G. and Gounaris, C.E., 2021. A mathematical optimization-based design framework for identifying stable bimetallic nanoclusters. *Molecular Systems Design & Engineering*, accepted for publication
* Hanselman, C.L., Yin, X., Miller, D.C. and Gounaris, C.E., 2021. `MatOpt: A Python package for nanomaterials discrete optimization. <http://gounaris.cheme.cmu.edu/drafts/Draft_MATOPT.pdf>`_


Basic Usage
-----------

There are two main sub-modules contained in the package serving two distinct purposes:

* The ``matopt.materials`` module contains objects and methods for efficiently representing and manipulating a nanomaterial and its design space.
* The ``matopt.opt`` module contains objects and methods for speeding up the casting of a Mixed-integer Linear Programming (MILP) model with simplified modeling syntax and automatic model formulation.

**Dependencies**

User access to the MILP solver CPLEX through Pyomo is assumed. For users who do not have access to CPLEX, the use of `NEOS-CPLEX <https://neos-guide.org/neos-interfaces#pyomo>`_ is suggested as an alternative.

**Define design canvas**

Several pieces of information about the material and design space need to be specified in order to formulate a materials optimization problem. To fulfill this need, the ``matopt.materials`` module defines generic and simple objects for describing the type of material to be designed and its design space, also referred to as a "canvas".

Some key objects are listed as follows:

.. module:: idaes.apps.matopt.materials.lattices.lattice

.. autoclass:: Lattice

.. module:: idaes.apps.matopt.materials.canvas

.. autoclass:: Canvas

.. module:: idaes.apps.matopt.materials.design

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

.. currentmodule:: idaes.apps.matopt.opt.mat_modeling

.. autoclass:: MaterialDescriptor

**Solve optimization model**

Once the model is fully specified, the user can optimize it in light of a chosen descriptor to serve as the objective to be maximized or minimized, as appropriate. Several functions are provided for users to choose from.

.. module:: idaes.apps.matopt.opt.mat_modeling

.. autoclass:: MatOptModel
   :members: optimize, populate, maximize, minimize

MatOpt Output
-------------
The results of the optimization process will be loaded into ``Design`` objects automatically. Users can then save material design(s) into files for further analysis and visualization using suitable functions provided. MatOpt provides interfaces to several standard crystal structure file formats, including CFG, PDB, POSCAR, and XYZ.