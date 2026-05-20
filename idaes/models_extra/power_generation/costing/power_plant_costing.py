#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Power Plant costing library
This method leverages NETL costing capabilities. Two main methods have been
developed to calculate the capital cost of power generation plants:

1. Fossil fueled power plants (from SCPC to IGCC) (get_PP_costing)
2. Supercritical CO2 power cycles (direct and indirect) (get_sCO2_unit_cost)

Other methods:

    - check_sCO2_costing_bounds() to display a warning if costing model have
      been used outside the range that where designed for
    - get_ASU_cost() to cost air separation units
"""

# TODO: Missing docstrings
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

__author__ = (
    "Costing Team (A. Noring, A. Deshpande, B. Paul, D. Caballero, and M. Zamarripa)"
)
__version__ = "1.0.0"

from pyomo.environ import (
    Param,
    Var,
    Constraint,
    value,
    Expr_if,
    units as pyunits,
)

import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetCostingBlockData
from idaes.models_extra.power_generation.costing.costing_dictionaries import (
    load_sCO2_costing_dictionary,
    register_power_plant_currency_units,
)

import idaes.logger as idaeslog
from idaes.core import declare_process_block_class

_log = idaeslog.getLogger(__name__)

# -----------------------------------------------------------------------------
# Power Plant Costing Library
# -----------------------------------------------------------------------------


@declare_process_block_class("PowerPlantCosting")
class PowerPlantCostingData(FlowsheetCostingBlockData):
    # Register currency and conversion rates based on CE Index
    # register_idaes_currency_units()
    if (
        not hasattr(pyunits, "USD_2008_Nov")
        and not hasattr(pyunits, "USD_2019_Sep")
        and not hasattr(pyunits, "USD_2018_Dec")
    ):
        register_power_plant_currency_units()

    def build_global_params(self):
        """
        This is where we can declare any global parameters we need, such as
        Lang factors, or coefficients for costing methods that should be
        shared across the process.

        You can do what you want here, so you could have e.g. sub-Blocks
        for each costing method to separate the parameters for each method.
        """
        # Set the base year for all costs
        self.base_currency = pyunits.USD_2018
        # Set a base period for all operating costs
        self.base_period = pyunits.year

    # -----------------------------------------------------------------------------
    # Supercritical CO2 Costing Library
    # -----------------------------------------------------------------------------
    def get_sCO2_unit_cost(
        self,
        equipment,
        scaled_param,
        temp_C=None,
        n_equip=1,
        CE_index_year="2018",
        custom_accounts=None,
    ):
        """
        Args:
            self: pyomo Block where constraints will be made
            equipment: the name of the SCO2 equipment to cost
            scaled_param: the scaling parameter (in appropriate units) for the
                elected equipment
            temp: the maximum temperature of the equipment. Not all types of
                equipment use a temperature correction factor, so it is
                optional
            n_equip: the number of pieces of equipment to cost
            CE_index_year: year for cost basis, e.g. "2018" to use 2018 dollars
            custom_accounts: user-defined dictionary to append to
                existing cost accounts dictionary

        Cost is in M$
        """
        # check to see if a costing block already exists
        if (
            self.parent_block().name
            in self.config.flowsheet_costing_block._registered_unit_costing  # pylint: disable=protected-access
        ):
            raise AttributeError(
                "{} already has an attribute costing. "
                "Check that you are not calling get_costing"
                " twice on the same model".format(self.name)
            )

        # define costing library
        self.library = "sCO2"
        self.equipment = equipment

        try:
            CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
            )

        # load sCO2 costing dictionary
        sCO2_costing_params = load_sCO2_costing_dictionary()

        param_dict = sCO2_costing_params[equipment]

        if custom_accounts is not None:
            # adding custom params to the dictionary
            for key, val in custom_accounts.items():
                param_dict[key] = val

        # define parameters
        self.ref_cost = Param(
            mutable=True,
            initialize=param_dict["a"],
            doc="Reference cost",
            units=pyunits.USD_2017,
        )

        self.exp = Param(
            mutable=True,
            initialize=param_dict["b"],
            doc="Scaling exponent",
            units=pyunits.dimensionless,
        )

        self.c = Param(
            mutable=True,
            initialize=param_dict["c"],
            doc="coefficient for temperature correction",
            units=pyunits.C**-1,
        )

        self.d = Param(
            mutable=True,
            initialize=param_dict["d"],
            doc="coefficient for temperature correction",
            units=pyunits.C**-2,
        )

        self.material_cost = Param(
            mutable=True,
            doc="material installation cost factor",
            initialize=param_dict["Material Cost"],
            units=pyunits.dimensionless,
        )

        self.labor_cost = Param(
            mutable=True,
            initialize=param_dict["Labor Cost"],
            doc="labor installation cost factor",
            units=pyunits.dimensionless,
        )

        # estimates for the percentages of TPC will be added later

        self.eng_fee = Param(
            mutable=True,
            initialize=0,
            doc="engineering fee percentage",
            units=pyunits.dimensionless,
        )

        self.process_conting = Param(
            mutable=True,
            initialize=0,
            doc="process contingency percentage",
            units=pyunits.dimensionless,
        )

        self.project_conting = Param(
            mutable=True,
            initialize=0,
            doc="project contingency percentage",
            units=pyunits.dimensionless,
        )

        # define variables
        # n_equip is left as a fixed variable to support MINLP optimization
        self.n_equip = Var(
            initialize=n_equip,
            doc="number of pieces of equipment",
            units=pyunits.dimensionless,
        )

        self.n_equip.fix(n_equip)

        self.scaled_param = Var(
            initialize=scaled_param,
            bounds=(0, 1e12),
            doc="scaled parameter",
            units=scaled_param.get_units(),
        )

        self.temp_factor = Var(
            initialize=1,
            bounds=(0.9, 100),
            doc="temperature correction factor",
            units=pyunits.dimensionless,
        )

        self.equipment_cost = Var(
            initialize=self.ref_cost,
            bounds=(0, 1e4),
            doc="equipment cost of sCO2 unit in $MM",
            units=CE_index_units,
        )

        self.bare_erected_cost = Var(
            initialize=self.ref_cost,
            bounds=(0, 1e4),
            doc="bare erected cost of sCO2 unit" "in $MM",
            units=CE_index_units,
        )

        self.total_plant_cost = Var(
            initialize=self.ref_cost,
            bounds=(0, 1e4),
            doc="total plant cost of sCO2 unit" "in $MM",
            units=CE_index_units,
        )

        # divides the scaled parameter by the number of pieces of equipment
        def scaled_param_rule(costing):
            return costing.scaled_param * costing.n_equip == scaled_param

        self.scaled_param_eq = Constraint(rule=scaled_param_rule)

        # check if equipment requires a temperature correction factor
        if equipment in [
            "Axial turbine",
            "Radial turbine",
            "Coal-fired heater",
            "Natural gas-fired heater",
            "Recuperator",
        ]:
            if temp_C is None:
                raise ValueError(
                    "Temperature argument is "
                    "required to cost %s equipment" % equipment
                )

            else:
                self.temperature = Var(
                    initialize=500,
                    bounds=(0, 1e6),
                    doc="dummy var for temperature",
                    units=pyunits.C,
                )

                self.temp_eq = Constraint(expr=self.temperature == temp_C)

                def temp_correction_rule(costing):  # rule for temp correction
                    return (
                        Expr_if(
                            costing.temperature < 550 * pyunits.C,
                            1e-6 * costing.temperature / pyunits.C + 1,
                            1
                            + costing.c * (costing.temperature - 550 * pyunits.C)
                            + costing.d * (costing.temperature - 550 * pyunits.C) ** 2,
                        )
                        == costing.temp_factor
                    )

                self.temp_correction_eq = Constraint(rule=temp_correction_rule)
        else:
            self.temp_factor.fix(1)

        # rule for equipment cost
        def equipment_cost_rule(costing):
            return (
                costing.equipment_cost
                == costing.n_equip
                * pyunits.convert(costing.ref_cost, CE_index_units)
                * (
                    (costing.scaled_param / costing.scaled_param.get_units())
                    ** costing.exp
                )
                * costing.temp_factor
            )

        self.equipment_cost_eq = Constraint(rule=equipment_cost_rule)

        # rule for bare erected cost
        def bare_erected_cost_rule(costing):
            return costing.bare_erected_cost == costing.equipment_cost * (
                1 + costing.material_cost + costing.labor_cost
            )

        self.bare_erected_cost_eq = Constraint(rule=bare_erected_cost_rule)

        # rule for calculating total plant cost
        def total_plant_cost_rule(costing):
            return costing.total_plant_cost == costing.bare_erected_cost * (
                1 + costing.eng_fee + costing.process_conting + costing.project_conting
            )

        self.total_plant_cost_eq = Constraint(rule=total_plant_cost_rule)

        # add variable and constraint scaling
        if equipment in ["Recuperator", "Direct air cooler"]:
            iscale.set_scaling_factor(self.scaled_param, 1e-5)
        else:
            iscale.set_scaling_factor(self.scaled_param, 1)

        iscale.set_scaling_factor(self.equipment_cost, 1e3)
        iscale.set_scaling_factor(self.bare_erected_cost, 1e3)
        iscale.set_scaling_factor(self.total_plant_cost, 1e3)
        iscale.constraint_scaling_transform(
            self.equipment_cost_eq, 1e-6, overwrite=False
        )
        iscale.constraint_scaling_transform(
            self.bare_erected_cost_eq, 1e3, overwrite=False
        )
        iscale.constraint_scaling_transform(
            self.bare_erected_cost_eq, 1e3, overwrite=False
        )

    def check_sCO2_costing_bounds(b):
        # load sCO2 costing dictionary
        sCO2_costing_params = load_sCO2_costing_dictionary()
        # This method accepts a flowsheet-level costing block
        for o in b.parent_block().component_objects(descend_into=False):
            # look for costing blocks
            if o.name in [block.name for block in b._registered_unit_costing]:
                if o.library == "sCO2":
                    equipment = o.equipment
                    lower_bound = sCO2_costing_params[equipment]["Lower Bound"]
                    upper_bound = sCO2_costing_params[equipment]["Upper Bound"]
                    if value(o.scaled_param) < lower_bound:
                        print(
                            """%s: The scaled parameter (%f) is below the lower
                            bound (%f)."""
                            % (value(o.name), value(o.scaled_param), lower_bound)
                        )
                    elif value(o.scaled_param) > upper_bound:
                        print(
                            """%s: The scaled parameter (%f) is above the upper
                            bound (%f)."""
                            % (value(o.name), value(o.scaled_param), upper_bound)
                        )
                    else:
                        print(
                            """%s: The scaled parameter is within the
                            bounds."""
                            % value(o.name)
                        )

    # -----------------------------------------------------------------------------
    # Air Separation Unit Costing Library
    # -----------------------------------------------------------------------------
    def get_ASU_cost(self, scaled_param, CE_index_year="2018"):
        # scaled parameter is O2 flowrate in TPD
        # only one set of parameters used, ref is hard coded for TPD and 2017

        params = {
            "Reference Cost": 3.26e6,
            "Reference Parameter": 13078,
            "Exponent": 0.7,
            "Eng Fee": 0.097,
            "Process": 0,
            "Project": 0.110,
        }

        # check to see if a costing block already exists
        if (
            self.parent_block().name
            in self.config.flowsheet_costing_block._registered_unit_costing  # pylint: disable=protected-access
        ):
            raise AttributeError(
                "{} already has an attribute costing. "
                "Check that you are not calling get_costing"
                " twice on the same model".format(self.name)
            )

        # define costing library
        self.library = "ASU"

        try:
            CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
            )

        # define parameters
        self.ref_cost = Param(
            initialize=params["Reference Cost"],
            mutable=True,
            doc="ASU reference cost",
            units=pyunits.USD_2008_Nov,
        )

        self.ref_param = Param(
            initialize=params["Reference Parameter"],
            mutable=True,
            doc="ASU reference parameter value",
            units=pyunits.ton / pyunits.d,
        )

        self.exp = Param(
            initialize=params["Exponent"],
            mutable=True,
            doc="ASU scaling exponent",
            units=pyunits.dimensionless,
        )

        self.eng_fee = Param(
            mutable=True,
            initialize=params["Eng Fee"],
            doc="engineering fee percentage",
            units=pyunits.dimensionless,
        )

        self.process_conting = Param(
            mutable=True,
            initialize=params["Process"],
            doc="process contingency percentage",
            units=pyunits.dimensionless,
        )

        self.project_conting = Param(
            mutable=True,
            initialize=params["Project"],
            doc="project contingency percentage",
            units=pyunits.dimensionless,
        )

        # define variables
        self.bare_erected_cost = Var(
            initialize=params["Reference Cost"],
            bounds=(0, 1e4),
            doc="scaled bare erected cost in $MM",
            units=CE_index_units,
        )

        self.total_plant_cost = Var(
            initialize=params["Reference Cost"],
            bounds=(0, 1e4),
            doc="total plant cost in $MM",
            units=CE_index_units,
        )

        # rule for scaling BEC

        def bare_erected_cost_rule(costing):
            return (
                costing.bare_erected_cost
                == pyunits.convert(
                    costing.ref_cost,
                    CE_index_units,
                )
                * (scaled_param / costing.ref_param) ** costing.exp
            )

        self.bare_erected_cost_eq = Constraint(rule=bare_erected_cost_rule)

        # rule for calculating TPC
        def total_plant_cost_rule(costing):
            return costing.total_plant_cost == costing.bare_erected_cost * (
                1 + costing.eng_fee + costing.process_conting + costing.project_conting
            )

        self.total_plant_cost_eq = Constraint(rule=total_plant_cost_rule)

        # add variable and constraint scaling
        iscale.set_scaling_factor(self.bare_erected_cost, 1)
        iscale.set_scaling_factor(self.total_plant_cost, 1)

        iscale.constraint_scaling_transform(
            self.bare_erected_cost_eq, 1e-3, overwrite=False
        )
        iscale.constraint_scaling_transform(
            self.total_plant_cost_eq, 1, overwrite=False
        )
