
ALAMOPY.ALAMO Options
=====================

This page lists in more detail the ALAMOPY options and the relation of ALAMO and ALAMOPY.

.. contents::
    :depth: 3


Basic ALAMOPY.ALAMO options
---------------------------

Data Arguments
^^^^^^^^^^^^^^

* **xmin**, **xmax**: minimum/maximum values of inputs, if not given they are calculated
* **zmin**, **zmax**: minimum/maximum values of outputs, if not given they are calculated
* **xlabels**: user-specified labels given to the inputs
* **zlabels**: user-specified labels given to the outputs

.. code-block:: python

  alamo(x_inputs, z_outputs, xlabels=['x1','x2'], zlabels=['z1','z2'])
  alamo(x_inputs, z_outputs, xmin=(-5,0),xmax=(10,15))


Available Basis Functions
^^^^^^^^^^^^^^^^^^^^^^^^^
* **linfcns, expfcns, logfcns, sinfcns, cosfcns**: 0-1 option to include linear, exponential, logarithmic, sine, and cosine transformations. For example

.. code-block:: python

  linfcns = 1, expfcns = 1, logfcns = 1, sinfcns = 1, cosfcns = 1


This results in basis functions =  x1, exp(x1), log(x1), sin(x1), cos(x1)
* **monomialpower, multi2power, multi3power**: list of monomial, binomial, and trinomial powers. For example

.. code-block:: python

  monomialpower = (2,3,4), multi2power = (1,2,3), multi3power = (1,2,3)


This results in the following basis functions: 

     * Monomial functions = x^2, x^3, x^4
     * Binomial functions = x1*x2, (x1*x2)^2, (x1*x2)^3
     * Trinomial functions = (x1*x2*x3), (x1*x2*x3)^2, (x1*x2*x3)^3

* **ratiopower**: list of ratio powers. For example

.. code-block:: python

  ratiopower = (1,2,3)

This results in basis functions = (x1/x2), (x1/x2)^2, (x1/x2)^3

.. code-block:: python

  alamo(x_inputs, z_outputs, linfcns=1, logfcns=1, expfcns=1)
  alamo(x_inputs, z_outputs, linfcns=1, multi2power=(2,3))


Note: Custom basis functions are discussed in the Advanced User Section.

ALAMO Regression Options
^^^^^^^^^^^^^^^^^^^^^^^^

* **showalm**: print ALAMO output to the screen
* **expandoutput**:  add a key to the output dictionary for multiple outputs
* **solvemip, builder, linearerror**:  A \01 indicator to solve with an optimizer (GAMSSOLVER), use a greedy heuristic, or use a linear objective instead of squared error.
* **modeler**:  Fitness metric to beused for model building (1-8)

  * 1. **BIC**: Bayesian infromation criterion
  * 2. **Cp**: Mallow's Cp
  * 3. **AICc**: the corrected Akaike's information criterion
  * 4. **HQC**: the Hannan-Quinn information criterion
  * 5. **MSE**: mean square error
  * 6. **SSEp**: sum of square error plus a penalty proportional to the model size (Note: convpen is the weight of the penalty)
  * 7. **RIC**: the risk information criterion
  * 8. **MADp**: the maximum absolute eviation plus a penalty proportional to  model size (Note: convpen is the weight of the penalty)

* **regularizer**: Regularization method used to reduce the number of potential basis functions before optimization of the selected fitness metric. Possible values are 0 and 1, corresponding to no regularization and regularization with the lasso, respectively.
* **maxterms**: Maximum number of terms to be fit in the model
* **convpen**: When MODELER is set to 6 or 8 the size of the model is weighted by CONVPEN.
* **almopt**: name of the alamo option file
* **simulator**: a python function to be used as a simulator for ALAMO, a variable that is a python function (not a string)
* **maxiter**: max iteration of runs


Validation Capabilities
^^^^^^^^^^^^^^^^^^^^^^^

* **xval, zval**: validation input/output variables
* **loo**: leave-one-out evaluation
* **lmo**: leave-many-out evaluation
* **cvfun**: cross-validation function (True/False)


File Options
^^^^^^^^^^^^
* **almname**: specify a name for the .alm file
* **savescratch**: saves .alm and .lst
* **savetrace**: saves tracefile
* **saveopt**: save .opt options file
* **savegams**: save the .gms gams file



ALAMOPY results dictionary
---------------------------


The results from alamopy.alamo are returned as a python dictionary.  The data can be accessed by using the dictionary keys listed below. For example

.. code-block:: python

  regression_results = doalamo(x_input, z_output, **kargs)
  model = regression_results['model']


Output models
^^^^^^^^^^^^^

* **f(model)**: A callable function
* **pymodel**: name of the python model written
* **model**: string of the regressed model

Note: A python script named after the output variables is written to the current directory. The model can be imported and used for further evaluation, for example to evaluate residuals:
        
.. code-block:: python

   import z1
   residuals = [y-z1.f(inputs[0],inputs[1]) for y,inputs in zip(z,x)]



Fitness metrics
^^^^^^^^^^^^^^^

* **size**: number of terms chosen in the regression
* **R2**: R2 value of the regression
* **Objective value metrics**: ssr, rmse, madp

Regression description
^^^^^^^^^^^^^^^^^^^^^^

* **version**: Version of ALAMO
* **xlabels, zlabels**: The labels used for the inputs/outputs
* **xdata, zdata**: array of xdata/zdata
* **ninputs, nbas**: number of inputs/basis functions


Performance specs
^^^^^^^^^^^^^^^^^
There are three types of regression problems that are used: ordinary linear regression (olr), classic linear regression (clr), and a mixed integer program (mip). Performance metrics include the number of each problems and the time spent on each type of problem. Additionally, the time spent on other operations and the total time are included.

* **numolr, olrtime, numclr, clrtime, nummip, miptime**: number of type of regression problems solved and time
* **othertime**: Time spent on other operations
* **totaltime**: Total time spent on the regression



Advanced user options in depth
------------------------------

Similar to ALAMO, there are advanced capabilities for customization and constrained regression facilitated by methods in ALAMOPY including custom basis functions, custom constraints on the response surface, and basis function groups. These methods interact with the regression using the alamo option file.

Custom Basis Functions
^^^^^^^^^^^^^^^^^^^^^^

Custom basis functions can be added to the built-in functions to expand the functional forms available. In ALAMO, this can be done with the following syntax

.. code-block:: python

  NCUSTOMBAS #
  BEGIN_CUSTOMBAS
  x1^2 * x2^2
  END_CUSTOMBAS

To use this advanced capability in ALAMOPY, the following function is called. Note it is necessary to use the xlabels assigned to the input parameters.

.. code-block:: python

  addCustomFunctions(fcn_list)
  addCustomFunctions(["x1^2 * x2^2", "...", "..." ...])



Custom Constraints
^^^^^^^^^^^^^^^^^^

Custom constraints can be placed on response surface or regressed function of the output variable. In ALAMO, this is controlled using custom constraints, CUSTOMCON. The constraints, a function **g(x_inputs, z_outputs)**  are applied to a specific output variable, which is the index of the output variable, and are less than or equal to 0 (g <= 0).

.. code-block:: python

  CRNCUSTOM #
  BEGIN_CUSTOMCON
  1 z1 - x1 + x2 + 1
  END_CUSTOMCON

To use this advanced capability in ALAMOPY, the following function is called. Note it is necessary to use the xlabels assigned to the input parameters.

.. code-block:: python

  addCustomConstraints(custom_constraint_list, **kargs)
  addCustomConstraints(["1 z1 - x1 + x2 +1", "...", "..." ...])



Basis Function Groups and Constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to imposing constraints on the response surface it produces, ALAMO has the ability to enforce constraints on groups of selected basis functions. This can be accomplished  using NGROUPS and identifying groups of basis functions. For ALAMO, this is achieved by first defining the groups with

.. code-block:: python

  NGROUPS 3
  BEGIN_GROUPS
  # Group-id Member-type Member-indices <Powers>
  1 LIN 1 2
  2 MONO 1 2
  3 GRP 1 2
  END_GROUPS


To add groups to ALAMOPY, you can use the following methods. Each Basis group has an index number that will be used as reference in the group constraints. The groups are defined by three or four parameters. Options for Member-type are LIN, LOG, EXP, SIN, COS, MONO, MULTI2, MULTI3, RATIO, GRP, RBF, and CUST.

.. code-block:: python

  addBasisGroup(type_of_function, input_indices, powers)
  addBasisGroups(groups)

  addBasisGroup("MONO", "1", "2")
  addBasisGroups([["LIN","1 2"],["MONO","1","2"],["GRP","1 2"]])

With the groups defined, constraints can be placed on the groups using the constraint-types NMT (no-more-than), ATL (at-least), REQ (requires), and XCL (exclude). For NMT and ATL the integer-parameter is the number of members in the group that should be selected based on the constraint. For REQ and XCL the integer-parameter is the group-id number of excluded or required basis functions.

.. code-block:: python

  BEGIN_GROUPCON
  # Group-id Output-id Constraint-type Integer-parameter
  3 1 NMT 1
  END_GROUPCON


To add the basis constraints to alamopy, you can use the following methods.

.. code-block:: python

  addBasisConstraint(group_id, output_id, constraint_type, intParam)
  addBasisConstraints(groups_constraint_list)

  addBasisConstraint(3,1,"NMT",1)
  addBasisConstraints([[3,1,"NMT",1]])

