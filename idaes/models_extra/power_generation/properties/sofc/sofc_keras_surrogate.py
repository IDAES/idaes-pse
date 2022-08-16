#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################

import os

from pyomo.environ import SolverFactory, Var, units as pyunits

from idaes.core import declare_process_block_class, UnitModelBlockData
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError

from idaes.core.surrogate.keras_surrogate import KerasSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock

import idaes.logger as idaeslog


path = os.path.dirname(os.path.abspath(__file__))


@declare_process_block_class("SofcSurrogate")
class SofcSurrogateData(UnitModelBlockData):
    def build(self):
        # create flowsheet input variables
        self.current_density = Var(
            self.flowsheet().time,
            initialize=4000,
            units=pyunits.A / pyunits.m**2,
            bounds=[1500, 6000],
        )

        self.fuel_temperature = Var(
            self.flowsheet().time,
            initialize=773.15,
            units=pyunits.kelvin,
            bounds=[288.15, 873.15],
        )

        self.internal_reforming = Var(
            self.flowsheet().time,
            initialize=0.4,
            units=pyunits.dimensionless,
            bounds=[0, 1],
            doc="fraction of natural gas internally reformed",
        )

        self.air_temperature = Var(
            self.flowsheet().time,
            initialize=973.15,
            units=pyunits.kelvin,
            bounds=[823.15, 1073.15],
        )

        self.air_recirculation = Var(
            self.flowsheet().time,
            initialize=0.5,
            units=pyunits.kelvin,
            bounds=[0, 0.8],
            doc="fraction of cathode air recirculated",
        )

        self.otc_ratio = Var(
            self.flowsheet().time,
            initialize=2.1,
            units=pyunits.dimensionless,
            bounds=[1.5, 3],
            doc="oxygen to carbon molar ratio at anode inlet",
        )

        self.fuel_utilization = Var(
            self.flowsheet().time,
            initialize=0.85,
            units=pyunits.dimensionless,
            bounds=[0.4, 0.95],
        )

        self.air_utilization = Var(
            self.flowsheet().time,
            initialize=0.5,
            units=pyunits.dimensionless,
            bounds=[0.125, 0.833],
        )

        self.pressure = Var(
            self.flowsheet().time,
            initialize=1,
            units=pyunits.atm,
            bounds=[1, 5],
            doc="SOFC operating pressure",
        )

        # create flowsheet output variables
        self.stack_voltage = Var(self.flowsheet().time, initialize=0.8, units=pyunits.V)

        self.max_cell_temperature = Var(
            self.flowsheet().time, initialize=1023.15, units=pyunits.kelvin
        )

        self.delta_cell_temperature = Var(
            self.flowsheet().time, initialize=100, units=pyunits.kelvin
        )

        self.anode_outlet_temperature = Var(
            self.flowsheet().time, initialize=800, units=pyunits.kelvin
        )

        self.cathode_outlet_temperature = Var(
            self.flowsheet().time, initialize=800, units=pyunits.kelvin
        )

        # create input and output variable object lists for flowsheet
        inputs = [
            self.current_density,
            self.fuel_temperature,
            self.internal_reforming,
            self.air_temperature,
            self.air_recirculation,
            self.otc_ratio,
            self.fuel_utilization,
            self.air_utilization,
            self.pressure,
        ]

        outputs = [
            self.stack_voltage,
            self.max_cell_temperature,
            self.delta_cell_temperature,
            self.anode_outlet_temperature,
            self.cathode_outlet_temperature,
        ]

        # create the Pyomo/IDAES block that corresponds to the surrogate
        keras_surrogate = KerasSurrogate.load_from_folder(
            os.sep.join([path, "sofc_surrogate_data"])
        )
        self.surrogate = SurrogateBlock()
        self.surrogate.build_model(
            keras_surrogate,
            formulation=KerasSurrogate.Formulation.FULL_SPACE,
            input_vars=inputs,
            output_vars=outputs,
        )

    def initialize(self, outlvl=idaeslog.NOTSET, solver="ipopt", optarg={"tol": 1e-6}):

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        opt = SolverFactory(solver)
        opt.options = optarg

        init_log.info_low("Starting initialization...")

        if degrees_of_freedom(self) != 0:
            raise ConfigurationError("User needs to check degrees of freedom")

        # solve model
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 1 {}.".format(idaeslog.condition(res)))
        init_log.info_high("Initialization Step 1 Complete.")
        init_log.info("Initialization Complete.")
