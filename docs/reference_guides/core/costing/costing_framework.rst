IDAES Process Costing Framework
===============================

.. contents:: Contents 
    :depth: 4

Introduction
------------

.. rubric:: Available Costing Libraries

This document outlines IDAES support for general process costing and integrating custom costing packages into IDAES models. The base costing framework supports attaching externally-defined blocks and methods for flowsheet and unit model cost calculations. Users must provide a costing package to use for cost calculations; IDAES provides a set of preexisting :ref:`Available Costing Libraries<reference_guides/core/costing/costing_libraries:Available Costing Libraries>` for use as examples and to guide users in developing their own costing packages.

When an instance of a unit model is passed to a `UnitModelCostingBlock`, a new sub-block is created 
on that unit named `costing` (i.e. `flowsheet.unit.costing`). All variables and constraints related to costing will be constructed within this new block.

Note that a `FlowsheetCostingBlock` must already exist, and users must also pass the correct costing method as an argument (see the example below for details on proper implementation).

Detailed Example
----------------

Below is a step-by-step example of the syntax required to integrate a custom costing package with an IDAES flowsheet.

Defining a Costing Package
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To start, users may pick a custom costing package and add it to a flowsheet or model via a costing parameter block.

* IDAES provides a default `FlowsheetCostingBlock` that will be overloaded with the selected costing package, and notify the user if any expected methods are missing.
* The `register_idaes_currency_units` method defines conversion rates for US Dollars based on the CE (Chemical Engineering) Cost Index. The CE Index assigns relative values to relate costs across years; for example, costs are converted from a `USD_2016` basis (`CE = 541.7`) to a `USD_2018` basis (`CE = 603.1`) by the relation `USD_2018 = 603.1/541.7 * USD_2016`. To standardize conversion, IDAES converts costs to a CE Index of 500 before converting to the final, user-selected cost units.
* Users may select cost units and set other base parameters via the `build_global_params` method, which scopes parameters for the entire costing package.
* `build_process_costs` allows integration of global cost methods such as ways to automatically aggregate (sum) total costs. It is not required to add anything to this method, but users may utilize existing aggregate cost variables in any relations they define.
* The `initialize_build` method allows for optional initialization steps for costs relations; note that unit model blocks and aggregate costs will already be initialized and don't need to be added here.

Users should define their costing package as below to properly overload the base methods, and may add their own methods to replace `cost_example_method`:

.. code:: python

    from idaes.core import FlowsheetCostingBlockData, register_idaes_currency_units

    @declare_process_block_class("CustomCosting")
    class CustomCostingData(FlowsheetCostingBlockData):

        # Register currency and conversion rates based on CE Index
        register_idaes_currency_units()

        def build_global_params(self):
            """
            This is where we can declare any global parameters we need, such as
            Lang factors, or coefficients for costing methods that should be
            shared across the process.
            You can do what you want here, so you could have e.g. sub-Blocks
            for each costing method to separate the parameters for each method.
            """
            # Set the base year for all costs
            self.base_currency = pyo.units.USD_2018
            # Set a base period for all operating costs
            self.base_period = pyo.units.year

        def build_process_costs(self):
            """
            This is where you do all your process wide costing.
            This is completely up to you, but you will have access to the
            following aggregate costs:
            1. self.aggregate_capital_cost
            2. self.aggregate_fixed_operating_cost
            3. self.aggregate_variable_operating_cost
            4. self.aggregate_flow_costs (indexed by flow type)
            """
            # process level costs go here
            pass

        def initialize_build(self):
            """
            Here we can add initialization steps for the things we built in
            build_process_costs.
            Note that the aggregate costs will be initialized by the framework.
            """
            # additional initialization steps go here
            pass

        @staticmethod
        def cost_example_method(
            blk,
            arg_1=Option1.option,
            arg_2=Option2.option,
            arg_3=Option3.option,
            arg_4=True,
        ):
            """
            An example costing method, use this to define arguments, variables and constraints for a unit model costing block.
            """
            # costing script goes here

Setting a Costing Package
^^^^^^^^^^^^^^^^^^^^^^^^^

Once the custom costing package (e.g. `idaes/models/costing/custom.py`) has been defined, it may be imported and set as a flowsheet-level costing block:

.. code:: python

    from pyomo.environ import ConcreteModel
    from idaes.core import FlowsheetBlock
    from my_costing_module import CustomCosting
    
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    
    m.fs.costing = CustomCosting()

Adding Capital Costing For a Specific Unit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From here, users may add capital costing for a specific unit in the process, in this case a `Heater`.

First, define the `Heater` unit:

.. code:: python

    from idaes.models.unit_models.heater import Heater
    from idaes.models.properties import thermodynamic_properties
    
    m.fs.properties = thermodynamic_properties.ThermoParameterBlock()  # property package for this unit model
    
    m.fs.unit = Heater(
        property_package=m.fs.properties,
        has_pressure_change=False,
        has_phase_equilibrium=True)

    # set inputs
    m.fs.unit.inlet.flow_mol[0].fix(100*pyunits.mol/pyunits.s)
    m.fs.unit.inlet.temperature[0].fix(300*pyunits.K)
    m.fs.unit.inlet.pressure[0].fix(101325*pyunits.Pa)
    m.fs.unit.inlet.mole_fraction[0, "A"] = 0.5  # fraction, dimensionless
    m.fs.unit.inlet.mole_fraction[0, "B"] = 0.5  # fraction, dimensionless
    
    m.fs.unit.outlet.temperature.fix(600*pyunits.K)

Using an existing IDAES costing library (e.g. SSLW), users may populate a costing block per a passed costing method. Suppose we are costing the heater using the SSLW `cost_fired_heater` method. IDAES provides a `UnitModelCostingBlock` to pass specific costing methods and argument to the base framework:

.. code:: python

    from idaes.models.costing.SSLW import SSLWCosting, SSLWCostingData
    from idaes.core import UnitModelCostingBlock
    from idaes.models.costing.SSLW import HeaterMaterial, HeaterSource

    # add a flowsheet costing block
    m.fs.costing = SSLWCosting()

    # add a heater costing block
    m.fs.H101.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=SSLWCostingData.cost_fired_heater,
        costing_method_arguments={
            "material_type": HeaterMaterial.CarbonSteel,
            "heat_source": HeaterSource.Fuel
        }
    )

The SSLW package contains a `unit_mapping` method to automatically lookup the proper method based on the unit model class (see the script `idaes.models.costing.sslw.py` for more details on implementing a unit mapping method). Omitting a `costing_method` argument or passing `None` will automatically trigger the SSLW module to lookup the `Heater` class and select the `cost_fired_heater` method:

.. code:: python

    from idaes.models.costing.SSLW import SSLWCosting, SSLWCostingData
    from idaes.core import UnitModelCostingBlock
    from idaes.models.costing.SSLW import HeaterMaterial, HeaterSource

    # add a flowsheet costing block
    m.fs.costing = SSLWCosting()

    # add a heater costing block
    m.fs.H101.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments={
            "material_type": HeaterMaterial.CarbonSteel,
            "heat_source": HeaterSource.Fuel
        }
    )

Alternatively, users may pass their own custom costing methods after writing their own costing packages:

.. code:: python

    from my_costing_module import CustomCosting, CustomCostingData
    from idaes.core import UnitModelCostingBlock
    from my_costing_module import Option1, Option2

    # add a flowsheet costing block
    m.fs.costing = CustomCosting()

    # add a heater costing block
    m.fs.H101.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=CustomCostingData.custom_cost_heater_method,
        costing_method_arguments={
            "option1": Option1.option,
            "option2": Option2.option
        }
    )

Adding Process Material and Utility Costs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once capital costing has been defined for desired units and the flowsheet as a whole, users may add process material and utility cost calculations to the flowsheet. Similar to capital cost methods, users may define their own process and utility cost methods.

Suppose our custom costing package defines the following costing method for fixed operating and maintenance costs. The contents are empty below; users should add their own variable and constraint definitions:

.. code:: python

    def get_fixed_OM_costs(
        b,
        nameplate_capacity=650,
        labor_rate=38.50,
        labor_burden=30,
        operators_per_shift=6,
        tech=1,
        fixed_TPC=None,
    ):
        """
        Creates constraints for the following fixed O&M costs in $MM/yr:
            1. Annual operating labor
            2. Maintenance labor
            3. Admin and support labor
            4. Property taxes and insurance
            5. Other fixed costs
            6. Total fixed O&M cost
            7. Maintenance materials (actually a variable cost, but scales off TPC)
        These costs apply to the project as a whole and are scaled based on the
        total TPC.
        Args:
            b: pyomo concrete model or flowsheet block
            nameplate_capacity: rated plant output in MW
            labor_rate: hourly rate of plant operators in project dollar year
            labor_burden: a percentage multiplier used to estimate non-salary
                labor expenses
            operators_per_shift: average number of operators per shift
            tech: int 1-7 representing the categories in get_PP_costing, used to
                determine maintenance costs
            TPC_value: The TPC in $MM that will be used to determine fixed O&M
            costs. If the value is None, the function will try to use the TPC
                calculated from the individual units.
        Returns:
            None
        """
        # code to define the method
        # builds and returns fixed OM costing variables and constraints

Calling the method above as `block.get_fixed_OM_costs(*args)` will add variables and constraints as attributes of the parent block, e.g. a `Var` named `block.total_fixed_OM_cost`, a `Constraint` named `block.total_fixed_OM_cost_rule` according to the method contents.

To build our variable set for fixed O&M costs, we call:

.. code:: python

    m.fs.costing.get_fixed_OM_costs(
                nameplate_capacity=nameplate_capacity,
                labor_rate=labor_rate,
                labor_burden=labor_burden,
                operators_per_shift=operators_per_shift,
                tech=tech
            )

As an example of the resulting block structure, the elements above would exist as `m.fs.costing.total_fixed_OM_cost` and `m.fs.costing.total_fixed_OM_cost_rule`. Users can leverage similar methods to add variable O&M costs, other generic plant-level costs, or custom expressions to calculate flowsheet-level attributes.

Aggregate Plant Costs
^^^^^^^^^^^^^^^^^^^^^

The IDAES Process Costing Framework will automatically compute the following costs for a flowsheet level costing block `m.fs.costing`:

* `m.fs.costing.aggregate_capital_cost`
* `m.fs.costing.aggregate_fixed_operating_cost`
* `m.fs.costing.aggregate_variable_operating_cost`
* `m.fs.costing.aggregate_flow_costs` (indexed by flow type)

These costs are the sums of their respective quantities, for example `m.fs.costing.aggregate_capital_cost` is the sum of all units for which `m.fs.unit.costing.capital_cost` exists. Therefore, to add these costs to the flowsheet users should call the following after defining all unit model costing blocks:

.. code:: python

    m.fs.costing.cost_process()

Units of Measurement
--------------------

It is important to highlight that the units of measurement of the model are part of the variables and expressions themselves and are converted as necessary. For example, if a costing package defines area in :math:`m^2`, while the cost correlations for heat exchangers require units to be in :math:`ft^2`, the costing method will convert the units to :math:`ft^2`. See the Pyomo Unit Containers documentation for further information on the subject: https://pyomo.readthedocs.io/en/stable/advanced_topics/units_container.html.

The base costing framework provides a method to automatically register currency units for the costing block based on Chemical Engineering (CE) Cost Index conversion rates for US Dollars. Users may define their own units for currency via methods discussed in the Pyomo units documentation linked above.
