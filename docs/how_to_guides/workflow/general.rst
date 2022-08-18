General Workflow
================

While IDAES offers significant freedom in how users write their models, they
are encouraged to follow this general workflow in order to make it easier for others to follow
their code.

This workflow is used throughout the tutorials and examples on the |examples-site|.

.. note::

    It is important to note that IDAES models are constructed upon execution of each line of
    code, and that most user defined options are only processed on model construction. This
    means that if the user wishes to make changes to any model construction option, it is
    necessary to rebuild the model from the beginning. Users should not be put off by this
    however, as model construction is generally very quick.

The general workflow for working with a model in IDAES is shown below:

.. contents:: :local:

1. Importing Modules
--------------------

IDAES is built upon a modular, object-oriented platform using Python, which requires users to
import the components from the appropriate model libraries. The necessary components and
libraries will vary from application to application, and were discussed earlier in this User
Guide, however some common components users will need include:

* Pyomo environment components (e.g. ConcreteModel, SolverFactory, TransformationFactory, Var, Constraint, objective) imported from `pyomo.environ`
* Pyomo network components (e.g. Arc, expand_arcs) from `pyomo.network`
* IDAES FlowsheetBlock, from `idaes.core`
* :ref:`Property packages<explanations/components/property_package/index:Property Package>` for materials of interest
* Unit models for process equipment, drawn from either the IDAES model libraries and/or user-defined models
* Data visualization and analysis tools. Common tools include degrees of freedom and scaling, a full list is provided :ref:`here<reference_guides/core/util/index:Utility Methods>`.
* External packages of interest to the user. Being built upon Python, users have access to the full range of Python libraries for working with and analyzing their models.

2. Building a Model
-------------------

The next step in the workflow is to create a model object which represents the problem to be
solved. The steps involved in this may vary depending on the problem being solved, but the
general procedure is as follows:

2.1 Create a Model Object
^^^^^^^^^^^^^^^^^^^^^^^^^

The foundation of any model in IDAES is a Pyomo `ConcreteModel` object, which is created as
follows:

.. code-block:: python

    m = ConcreteModel()

.. note::

    IDAES does not support the use of Pyomo `AbstractModels`

2.2 Add a Flowsheet to the Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The foundation of a process model within IDAES is the `FlowsheetBlock`, which forms the canvas
upon which the process will be constructed. A key aspect of the `FlowsheetBlock` is to define
whether the model will be steady-state or dynamic, and to define the time domain as appropriate.

.. code-block:: python

    m.fs = FlowsheetBlock(dynamic=False)

.. note::

    IDAES supports nested flowsheets to allow complex processes
    to be broken down into smaller sub-processes.

2.3 Add Property Packages to Flowsheet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All process models depend on calculations of thermophysical and chemical reaction properties,
which are represented in IDAES using property packages. Users need to add the property packages
they intend to use to the flowsheet.

.. code-block:: python

    m.fs.properties_1 = MyPropertyPackage.PhysicalParameterBlock()

.. note::

    Users can add as many property packages as they need to a flowsheet, and can determine which
    property package will be used for each unit operation as it is created.

2.4 Add Unit Models to Flowsheet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, the user can add Unit Models to their flowsheet to represent each unit operation in the
process.

.. code-block:: python

    m.fs.unit01 = UnitModel(property_package=m.fs.properties_1)

2.5 Define Unit Model Connectivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to describe the flow of material between unit operations, users must declare `Arcs`
(or streams) which connect the outlet of each unit operation to the inlet of the next.

.. code-block:: python

    m.fs.arc_1 = Arc(source=m.fs.unit01.outlet, destination=m.fs.unit02.inlet)

2.6 Expand Arcs
^^^^^^^^^^^^^^^

It is important to note that `Arcs` only define the connectivity between unit operations, but
do not create the actual model constraints needed to describe this. Once all `Arcs` in a
flowsheet have been defined, it is necessary to expand these `Arcs` using the Pyomo
`TransformationFactory`.

.. code-block:: python

    TransformationFactory("network.expand_arcs").apply_to(m)

.. note::

    Pyomo provides a number of other Transformations and tools that may be useful to the user
    depending on the application. Examples include the `gdp` and `dae` transformations.

2.7 Add Variables, Constraints and Objectives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, users can add any additional variables, constraints and objectives to their model.
These could include the objective function for which they wish to optimize, additional
constraints that provide limits on process performance, or simply additional quantities that
the user wishes to use in analyzing or visualizing the results.

3. Scaling the Model
--------------------

.. note::

    The IDAES scaling tools are currently under development.

Ensuring that a model is well scaled is important for increasing the efficiency and reliability
of solvers, and users should consider model scaling as an integral part of the modeling process.
IDAES provides a number of tool for assisting users with scaling their models, and details on
these can be found :ref:`here<reference_guides/core/util/scaling:Scaling Methods>`.

4. Specifying the Model
-----------------------

.. note::

    IDAES is in the process of developing a set of tools to assist users with working with units
    of measurement when fixing and displaying values.

The next step is to specify the model by fixing variables. which can be done using the form
`variable_name.fix(value)`. The variables that need to be fixed are application dependent,
but commonly include the feed state variables.

In order to prepare the model for initialization, it is necessary to fully specify the model,
such that there are no degrees of freedom. IDAES provides a tools for counting and reporting
the degrees of freedom in any model (or sub-model/block):

.. code-block:: python

    from idaes.core.util.model_statistics import degrees_of_freedom

    print(degrees_of_freedom(m))

.. note::

    Whilst it is not always necessary to fully define a model before initialization, it is much
    safer to do so as it ensures the model is well-defined. Most IDAES initialization tools
    check that the model is well-defined before proceeding, and will raise an Exception if it is
    not.

.. note::

    Depending on the solver to be used during initialization, it can be better to avoid putting
    bounds on variables and adding inequality constraints at this stage. For solving square
    problems (i.e. zero degrees of freedom), some solvers (e.g. IPOPT) perform better without
    bounds on the problem. These bounds and constraints can be added later when it comes time to
    optimize the problem.

5. Initializing the Model
-------------------------

The next step is to initialize the model. All IDAES models have established initialization
methods that can be called using `model.initialize()` which can be expected to take a model
from its initial state to a feasible solution for a set of initial guesses (within the models
expected operating range).

IDEAS workflows generally use a sequential-modular approach to
initialize flowsheets, where unit models are initialized sequentially, passing the outlet
state from one unit as the initial state for the next. An automated sequential-modular tool is
available through Pyomo and demonstrated in the tutorials.

6. Solving the Model
--------------------

.. important::

    The sequential-modular approach initializes each unit model individually, thus it is
    important to do a final solve of the overall flowsheet/model in order to complete the
    initialization process. In most cases, this final solve should only take a few iterations,
    as the state of each unit model should be at or near the final solution already.

In order to solve the model, it is necessary to create a solve object and set any desired solver
options (such as tolerances, iteration limits etc.).

.. code-block:: python

    solver = SolverFactory('solver_name')
    solver.options = {'tol': 1d-6}

    results = solver.solve(m)

Users should check the output from the solver to ensure a feasible solution was found using
the following:

.. code-block:: python

    print(results.solver.termination_condition)

Different problems will require different solvers, and users will need to experiment to find
those that work best for their problems. The default solver for most IDAES applications is
IPOPT, which can be downloaded using the ``idaes get-extensions`` command line.

7. Optimizing the Model
-----------------------

Once an initial solution has been found, users can proceed to solving the optimization problem
of interest. This procedure will vary by application but generally involves the following steps:

7.1) Unfix some degrees of freedom to provide the problem with decision variables, `variable_name.unfix()`.

7.2) Add bounds to variables and inequality constraints to constrain solution space, `variable_name.setlb(value)` and `var_name.setub(value)`

7.3) Call a solver and check the termination conditions, see step 6 Solving the Model.

.. note::

    Users may wish/need to use different solvers for initialization and optimization. IDAES and
    Pyomo support the use of multiple solvers as part of the same workflow for solving different
    types of problems.

8. Analyzing and Visualizing the Results
----------------------------------------

One of the benefits of the IDAES Integrated Platform is that it operates in a fully featured
programming language, which provides users a high degree of flexibility in analyzing their
models. For example, users can automate the simulation of the model across multiple objectives
or a range of parameters, store and save results from one or multiple solutions. Users also have
access to a wide range of tools for manipulating, plotting and visualizing the results.
