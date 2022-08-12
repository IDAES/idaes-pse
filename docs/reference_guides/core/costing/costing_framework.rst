IDAES Process Costing Framework
===============================

.. contents:: Contents 
    :depth: 4

Introduction
------------

.. rubric:: Available Costing Libraries

This document outlines IDAES support for general process costing and integrating custom costing packages into IDAES models. The base costing framework supports attaching externally-defined blocks and methods for flowsheet and unit model cost calculations. Users must provide a costing package to use for cost calculations; IDAES provides a set of preexisting :ref:`Available Costing Libraries<reference_guides/core/costing/costing_libraries:Available Costing Libraries>`.

When an instance of a unit model is passed to a `UnitModelCostingBlock`, a new sub-block is created 
on that unit named `costing` (i.e. `flowsheet.unit.costing`). All variables and constraints related to costing will be constructed within this new block.

Note that a `FlowsheetCostingBlock` must already exist, and users must also pass the correct costing method as an argument (see the example below for details on proper implementation).

Example
-------
Below is a simple example of the syntax required to integrate a custom costing package with an IDAES flowsheet. After importing and building the relevant unit model and property package blocks, users may easily attach a flowsheet-level costing block and extract any required data class components:

* `UnitModel` stands for a unit model block - users should supply an IDAES or custom unit model
* `ParameterBlock` stands for a property package parameter block - users should supply an IDAES or custom property package
* `UnitModelCostingBlock` is the required class to define the costing sub-block
* `FlowsheetCostingBlock` is the required class to define the flowsheet-level costing level; when costing individual unit models with a custom costing package, inheritance allows the costing class to stand in for this block. Below, `CustomPackageCosting` inherits all flowsheet-level costing block properties and may be invoked instead.
* `FlowsheetCostingBlock` is the required costing data class, which likewise may be defined by invoking a `CustomPackageCostingData` via inheritance.
* `unit_costing_method` is a specific costing method that IDAES will pass to define costing variables and constraints. This is completely up to the user and grants freedom and flexibility to integrate custom costing for specific unit models.
* `CustomPackageOption1` and `CustomPackageOption2` are example costing method arguments, which pass an `option` to the costing block.


.. code:: python

    from pyomo.environ import ConcreteModel
    from idaes.core import FlowsheetBlock
    from custom_properties import ParameterBlock
    from idaes.models.unit_models import UnitModel

    from idaes.core import UnitModelCostingBlock
    from idaes.models.costing.CustomPackage import CustomPackageCosting, CustomPackageCostingData
    from idaes.models.costing.CustomPackage import CustomPackageOption1, CustomPackageOption2
    
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    
    m.fs.properties = ParameterBlock()  # this is a relevant thermodynamic property package
    
    m.fs.unit = UnitModel()  # this is a unit model supported by the custom costing package
    # set inputs
    # ...
    # adding flow, temperature, pressure, other state and unit specifications
    
    m.fs.costing = CustomPackageCosting()
    m.fs.unit.costing = UnitModelCostingBlock(
        default={
            "flowsheet_costing_block": m.fs.costing,
            "costing_method": CustomPackageCostingData.unit_costing_method,
            "costing_method_arguments": {
                    "option_1": CustomPackageOption1.option,
                    "option_2": CustomPackageOption2.option
            }
        }
    )

Units of Measurement
--------------------

It is important to highlight that the costing method interrogates the property package to determine the units of measurement of this model. For example, if a costing package defines area in :math:`m^2`, while the cost correlations for heat exchangers require units to be in :math:`ft^2`, the costing method will convert the units to :math:`ft^2`. The use of Pyomo-unit conversion tools is under development.

The base costing framework provides a method to automatically register currency units for the costing block based on Chemical Engineering (CE) Cost Index conversion rates for US Dollars. IDAES will check is Pyomo has standard currency units registered, and if not will load unit conversions between USD at CE indices of 500 and 394, as well as USD annually from 1990 to 2020.
