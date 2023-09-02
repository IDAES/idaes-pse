MatOpt : Materials Design via Mathematical Optimization
=======================================================

The MatOpt (Materials Design via Mathematical Optimization) module provides simple tools for creating Pyomo model objects for optimization-based nanomaterials design. MatOpt can be used to design crystalline nanostructured materials, including but not limited to particles, surfaces, and periodic bulk structures.

The main goals of this package are as follows:

-   To simplify the representation of nanostructured materials, streamlining the creation of materials optimization problems.
-   To automate many of the necessary steps of materials optimization, speeding up the development of new models and accelerating new materials discovery.
-   To provide a simple interface so that users do not need to handle the details of building efficient mathematical optimization models or the specific Pyomo syntax to do this.

Thank you for your interest in MatOpt. We would love to hear your feedback! Please report any thoughts, questions or bugs to: <gounaris@cmu.edu>

If you are using MatOpt, please consider citing:

-   Hanselman, C.L., Yin, X., Miller, D.C. and Gounaris, C.E., 2022. [MatOpt: A Python Package for Nanomaterials Design Using Discrete Optimization.](https://doi.org/10.1021/acs.jcim.1c00984) *Journal of Chemical Information and Modeling*, 62(2), pp.295-308.

Basic Usage
-----------

There are two main sub-modules contained in the package serving two distinct purposes:

-   The `matopt.materials` module contains objects and methods for efficiently representing and manipulating a nanomaterial and its design space.
-   The `matopt.opt` module contains objects and methods for speeding up the casting of a Mixed-integer Linear Programming (MILP) model with simplified modeling syntax and automatic model formulation.

**Dependencies**

User access to the MILP solver CPLEX through Pyomo is assumed. For users who do not have access to CPLEX, the use of [NEOS-CPLEX](https://neos-guide.org/neos-interfaces#pyomo) is suggested as an alternative.

**Define design canvas**

Several pieces of information about the material and design space need to be specified in order to formulate a materials optimization problem. To fulfill this need, the `matopt.materials` module defines generic and simple objects for describing the type of material to be designed and its design space, also referred to as a "canvas".

Some key objects are listed as follows:

**Build model via descriptors**

The material type and design space specified provide indices, sets, and parameters for the optimization model. Using simple syntax, inspired by materials-related terminology, MatOpt users define a `MatOptModel` object, which will be translated into a Pyomo `ConcreteModel` object automatically.

MatOpt uses `MaterialDescriptor` objects to represent variables, constraints, and objectives. A `MatOptModel` object holds lists of `MaterialDescriptor` objects. By default, several universal site descriptors are pre-defined in the model.

<table>
<col width="15%" />
<col width="84%" />
<thead>
<tr class="header">
<th align="left">Descriptor</th>
<th align="left">Explanation</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><code>Yik</code></td>
<td align="left">Presence of a building block of type k at site i</td>
</tr>
<tr class="even">
<td align="left"><code>Yi</code></td>
<td align="left">Presence of any type of building block at site i</td>
</tr>
<tr class="odd">
<td align="left"><code>Xijkl</code></td>
<td align="left">Presence of a building block of type k at site i and a building block of type l at site j</td>
</tr>
<tr class="even">
<td align="left"><code>Xij</code></td>
<td align="left">Presence of any building block at site i and any building block at site j</td>
</tr>
<tr class="odd">
<td align="left"><code>Cikl</code></td>
<td align="left">Count of neighbors of type l next to a building block of type k at site i</td>
</tr>
<tr class="even">
<td align="left"><code>Ci</code></td>
<td align="left">Count of any type of neighbors next to a building block at site i</td>
</tr>
</tbody>
</table>

User-specified descriptors are defined by `DescriptorRule` objects in conjunction with `Expr` expression objects. Available expressions include:

<table>
<col width="24%" />
<col width="75%" />
<thead>
<tr class="header">
<th align="left">Expression</th>
<th align="left">Explanation</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><code>LinearExpr</code></td>
<td align="left">Multiplication and addition of coefficients to distinct descriptors</td>
</tr>
<tr class="even">
<td align="left"><code>SiteCombination</code></td>
<td align="left">Summation of site contributions from two sites</td>
</tr>
<tr class="odd">
<td align="left"><code>SumNeighborSites</code></td>
<td align="left">Summation of site contributions from all neighboring sites</td>
</tr>
<tr class="even">
<td align="left"><code>SumNeighborBonds</code></td>
<td align="left">Summation of bond contributions to all neighboring sites</td>
</tr>
<tr class="odd">
<td align="left"><code>SumSites</code></td>
<td align="left">Summation across sites</td>
</tr>
<tr class="even">
<td align="left"><code>SumBonds</code></td>
<td align="left">Summation across bonds</td>
</tr>
<tr class="odd">
<td align="left"><code>SumSiteTypes</code></td>
<td align="left">Summation across site types</td>
</tr>
<tr class="even">
<td align="left"><code>SumBondTypes</code></td>
<td align="left">Summation across bond types</td>
</tr>
<tr class="odd">
<td align="left"><code>SumSitesAndTypes</code></td>
<td align="left">Summation across sites and site types</td>
</tr>
<tr class="even">
<td align="left"><code>SumBondsAndTypes</code></td>
<td align="left">Summation across bonds and bond types</td>
</tr>
<tr class="odd">
<td align="left"><code>SumConfs</code></td>
<td align="left">Summation across conformation types</td>
</tr>
<tr class="even">
<td align="left"><code>SumSitesAndConfs</code></td>
<td align="left">Summation across sites and conformation types</td>
</tr>
</tbody>
</table>

Several types of `DescriptorRules` are available.

<table>
<col width="26%" />
<col width="73%" />
<thead>
<tr class="header">
<th align="left">Rule</th>
<th align="left">Explanation</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><code>LessThan</code></td>
<td align="left">Descriptor less than or equal to an expression</td>
</tr>
<tr class="even">
<td align="left"><code>EqualTo</code></td>
<td align="left">Descriptor equal to an expression</td>
</tr>
<tr class="odd">
<td align="left"><code>GreaterThan</code></td>
<td align="left">Descriptor greater than or equal to an expression</td>
</tr>
<tr class="even">
<td align="left"><code>FixedTo</code></td>
<td align="left">Descriptor fixed to a scalar value</td>
</tr>
<tr class="odd">
<td align="left"><code>PiecewiseLinear</code></td>
<td align="left">Descriptor equal to the evaluation of a piecewise linear function</td>
</tr>
<tr class="even">
<td align="left"><code>Implies</code></td>
<td align="left">Indicator descriptor that imposes other constraints if equal to 1</td>
</tr>
<tr class="odd">
<td align="left"><code>NegImplies</code></td>
<td align="left">Indicator descriptor that imposes other constraints if equal to 0</td>
</tr>
<tr class="even">
<td align="left"><code>ImpliesSiteCombination</code></td>
<td align="left">Indicator bond-indexed descriptor that imposes constraints on the two sites</td>
</tr>
<tr class="odd">
<td align="left"><code>ImpliesNeighbors</code></td>
<td align="left">Indicator site-indexed descriptor that imposes constraints on neighboring sites</td>
</tr>
</tbody>
</table>

From the combination of the above pre-defined descriptors, expressions, and rules, a user can specify a wide variety of other descriptors, as necessary.

**Solve optimization model**

Once the model is fully specified, the user can optimize it in light of a chosen descriptor to serve as the objective to be maximized or minimized, as appropriate. Several functions are provided for users to choose from.

MatOpt Output
-------------

The results of the optimization process will be loaded into `Design` objects automatically. Users can then save material design(s) into files for further analysis and visualization using suitable functions provided. MatOpt provides interfaces to several standard crystal structure file formats, including CFG, PDB, POSCAR, and XYZ.

MatOpt Examples
---------------
Several [case studies](https://idaes.github.io/examples-pse/latest/Examples/MatOpt/index.html) are provided to illustrate the detailed usage of MatOpt. In each case, a Jupyter notebook with explanations as well as an equivalent Python script is provided.

References
----------

-   Hanselman, C.L. and Gounaris, C.E., 2016. [A mathematical optimization framework for the design of nanopatterned surfaces.](https://aiche.onlinelibrary.wiley.com/doi/full/10.1002/aic.15359) *AIChE Journal*, 62(9), pp.3250-3263.
-   Hanselman, C.L., Alfonso, D.R., Lekse, J.W., Matranga, C., Miller, D.C. and Gounaris, C.E., 2019. [A framework for optimizing oxygen vacancy formation in doped perovskites.](https://www.sciencedirect.com/science/article/pii/S0098135418310998) *Computers & Chemical Engineering*, 126, pp.168-177.
-   Hanselman, C.L., Zhong, W., Tran, K., Ulissi, Z.W. and Gounaris, C.E., 2019. [Optimization-based design of active and stable nanostructured surfaces.](https://pubs.acs.org/doi/abs/10.1021/acs.jpcc.9b08431) *The Journal of Physical Chemistry C*, 123(48), pp.29209-29218.
-   Isenberg, N.M., Taylor, M.G., Yan, Z., Hanselman, C.L., Mpourmpakis, G. and Gounaris, C.E., 2020. [Identification of optimally stable nanocluster geometries via mathematical optimization and density-functional theory.](https://pubs.rsc.org/en/content/articlelanding/2019/me/c9me00108e#!divAbstract) *Molecular Systems Design & Engineering*.
-   Yin, X., Isenberg, N.M., Hanselman, C.L., Dean, J.R., Mpourmpakis, G. and Gounaris, C.E., 2021. [Designing stable bimetallic nanoclusters via an iterative two-step optimization approach.](https://pubs.rsc.org/en/content/articlelanding/2021/ME/D1ME00027F) *Molecular Systems Design & Engineering*, 6(7), pp.545-557.
-   Yin, X. and Gounaris, C.E., 2022. [Search methods for inorganic materials crystal structure prediction.](https://doi.org/10.1016/j.coche.2021.100726) *Current Opinion in Chemical Engineering*, 35, p.100726.
-   Hanselman, C.L., Yin, X., Miller, D.C. and Gounaris, C.E., 2022. [MatOpt: A Python Package for Nanomaterials Design Using Discrete Optimization.](https://doi.org/10.1021/acs.jcim.1c00984) *Journal of Chemical Information and Modeling*, 62(2), pp.295-308.
