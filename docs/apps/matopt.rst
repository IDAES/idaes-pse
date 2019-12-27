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

* matopt.materials.lattices

 * matopt.materials.lattices.lattice
 * matopt.materials.lattices.cubic_lattice
 * matopt.materials.lattices.fcc_lattice
 * matopt.materials.lattices.perovskite_lattice
 * matopt.materials.lattices.unit_cell_lattice

* matopt.materials.canvas
* matopt.materials.atom
* matopt.materials.design

Advanced

* matopt.materials.geometry
* matopt.materials.tiling
* matopt.materials.transform_func
* matopt.materials.motifs

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
| m.Yik      | Presence of a building block of type k at site i                                          |
+------------+-------------------------------------------------------------------------------------------+
| m.Yi       | Presence of any type of building block at site i                                          |
+------------+-------------------------------------------------------------------------------------------+
| m.Xijkl    | Presence of a building block of type k at site i and a building block of type l at site j |
+------------+-------------------------------------------------------------------------------------------+
| m.Xij      | Presence of any building block at site i and any building block at site j                 |
+------------+-------------------------------------------------------------------------------------------+
| m.Cikl     | Count of neighbors of type l next to a building block of type k at site i                 |
+------------+-------------------------------------------------------------------------------------------+
| m.Ci       | Count of any type of neighbors next to a building block at site i                         |
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

* Method to create and optimize the materials design problem.

.. code-block:: python

 optimize(self,func,sense,nSolns=1,tee=True,disp=1,keepfiles=False, tilim=3600,trelim=None, solver='cplex')
        """Method to create and optimize the materials design problem.

        This method automatically creates a new optimization model every time it is called. Then, it solves the model via Pyomo with the
        CPLEX solver.

        If multiple solutions (called a 'solution pool') are desired, then
        the nSolns argument can be provided and the populate method will
        be called instead.

        Args:
        func (MaterialDescriptor/Expr): Material functionality to optimize.
        sense (int): flag to indicate the choice to minimize or maximize the
            functionality of interest.
            Choices: minimize/maximize (Pyomo constants 1,-1 respectively)
        nSolns (int): Optional, number of Design objects to return.
            Default: 1 (See MatOptModel.populate for more information)
        tee (bool): Optional, flag to turn on solver output.
            Default: True
        disp (int): Optional, flag to control level of MatOpt output.
            Choices: 0: No MatOpt output (other than solver tee)
                     1: MatOpt output for outer level method
                     2: MatOpt output for solution pool & individual solns.
            Default: 1
        keepfiles (bool): Optional, flag to save temporary pyomo files.
            Default: True
        tilim (float): Optional, solver time limit (in seconds).
            Default: 3600
        trelim (float): Optional, solver tree memeory limit (in MB).
            Default: None (i.e., Pyomo/CPLEX default)
        solver (str): Solver choice. Currently only cplex or
            neos-cplex are supported
            Default: cplex

        Returns:
        (Design/list<Design>) Optimal design or designs, depending on the
            number of solutions requested by argument nSolns.
        """

* If multiple solutions (called a 'solution pool') are desired, then the nSolns argument can be provided, and the populate method will be called instead.

.. code-block:: python

 populate(self,func,sense,nSolns,tee=True,disp=1,keepfiles=False, tilim=3600,trelim=None, solver='cplex')
        """Method to a pool of solutions that optimize the material model.

        This method automatically creates a new optimization model every time it is called. Then, it solves the model via Pyomo with the
        CPLEX solver.

        The populate method iteratively solves the model, interprets the
        solution as a Design object, creates a constraint to disallow that
        design, and resolves to find the next best design. We build a pool
        of Designs that are gauranteed to be the nSolns-best solutions in the
        material design space.

        Args:
        func (MaterialDescriptor/Expr): Material functionality to optimize.
        sense (int): flag to indicate the choice to minimize or maximize the
            functionality of interest.
            Choices: minimize/maximize (Pyomo constants 1,-1 respectively)
        nSolns (int): Optional, number of Design objects to return.
            Default: 1 (See MatOptModel.populate for more information)
        tee (bool): Optional, flag to turn on solver output.
            Default: True
        disp (int): Optional, flag to control level of MatOpt output.
            Choices: 0: No MatOpt output (other than solver tee)
                     1: MatOpt output for outer level method
                     2: MatOpt output for solution pool & individual solns.
            Default: 1
        keepfiles (bool): Optional, flag to save temporary pyomo files.
            Default: True
        tilim (float): Optional, solver time limit (in seconds).
            Default: 3600
        trelim (float): Optional, solver tree memeory limit (in MB).
            Default: None (i.e., Pyomo/CPLEX default)
        solver (str): Solver choice. Currently only cplex or
            neos-cplex are supported
            Default: cplex

        Returns:
        (list<Design>) A list of optimal Designs in order of decreasing
            optimality.
        """

* Wrapper functions for maximization/minimization.

.. code-block:: python

 maximize(self,func,**kwargs)
 minimize(self,func,**kwargs)

MatOpt Output
-------------
Optimization results will be loaded into design objects automatically.
Users can then write material design(s) into files for further analysis and visualization.
MatOpt provides interfaces to several standard crystal structure file formats, including CFG, PDB, POSCAR, and XYZ (matopt.materials.parsers).

MatOpt Examples
---------------
Five case studies are provided in idaes.docs.app.matopt
In each case, a Jupyter notebook with explanation as well as an equivalent Python script is provided.
