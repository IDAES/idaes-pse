#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

"""
Initializer for the fixed bed TSA unit model.
"""

# Import Pyomo libraries
from pyomo.environ import (
    value,
    check_optimal_termination,
    Block,
    units,
    Var,
    Constraint,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# import IDAES core libraries
from idaes.core.initialization import ModularInitializerBase
from idaes.core.initialization.initializer_base import StoreState
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.model_serializer import to_json, from_json
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver

from idaes.models_extra.temperature_swing_adsorption import SteamCalculationType


import idaes.logger as idaeslog

__author__ = "Daison Yancy Caballero, Alex Noring"

# Set up logger
_log = idaeslog.getLogger(__name__)


class FixedBedTSA0DInitializer(ModularInitializerBase):
    """
    Initializer for 0D Fixed Bed TSA units.

    """

    def initialize(
        self,
        model: Block,
        initial_guesses: dict = None,
        json_file: str = None,
        output_level=None,
        exclude_unused_vars: bool = True,
        heating_time_guess=1000,
        cooling_time_guess=500,
    ):

        if not exclude_unused_vars:
            exclude_unused_vars = True
            _log.warning(
                "The FixedBedTSA0DInitializer was run with "
                "exclude_unused_vars=False, but the FixedBedTSA0D model "
                "contains unused Vars. exclude_unused_vars has been set "
                "to True to avoid raising an InitializationError."
            )

        self.heating_time_guess = heating_time_guess
        self.cooling_time_guess = cooling_time_guess

        super().initialize(
            model=model,
            initial_guesses=initial_guesses,
            json_file=json_file,
            output_level=output_level,
            exclude_unused_vars=exclude_unused_vars,
        )

    def initialization_routine(
        self,
        blk: Block,
    ):
        """
        Initialization routine for fixed bed TSA unit.

        Args:
            blk: model to be initialized

        Returns:
            None

        Raises:
            InitializationError: If degrees of freedom is not zero at the start
            of each initialization step.

        """
        # set up logger for initialization and solve
        init_log = idaeslog.getInitLogger(blk.name, self.get_output_level(), tag="unit")
        solve_log = idaeslog.getSolveLogger(
            blk.name, self.get_output_level(), tag="unit"
        )

        # create solver
        if self.config.solver is None:
            self.opt = get_solver(self.config.solver, self.config.solver_options)
        else:
            self.opt = self.config.solver

        # initialization of fixed bed TSA model unit
        init_log.info("Starting fixed bed TSA initialization")

        tsa_state = to_json(blk, wts=StoreState, return_dict=True)

        # 1 - solve heating step

        # 1.1) fix states in the fixed bed TSA inlet  ("flow_mol_in_total",
        # "mol  e_frac_in", and "pressure_adsorption"). These are equal
        # to those states coming from the exhaust gas stream in the CCS system

        vars_lst_heating = ["flow_mol_in_total", "pressure_adsorption", "mole_frac_in"]

        cons_lst_heating = ["flow_mol_in_total_eq", "pressure_in_eq", "mole_frac_in_eq"]

        self._calculate_and_fix_variable_from_constraint(
            blk, variable_list=vars_lst_heating, constraint_list=cons_lst_heating
        )

        # 1.2) initial solution using false position method

        # deactivate final condition constraint and fix time
        blk.heating.fc_temperature_eq.deactivate()
        blk.heating.time.fix(self.heating_time_guess * units.s)

        # check degrees of freedom and solve
        if degrees_of_freedom(blk.heating) == 0:
            self._false_position_method(
                blk,
                cycle_step=blk.heating,
                t_guess=self.heating_time_guess,
            )
        else:
            raise InitializationError(
                "Degrees of freedom is not zero during initialization of "
                "heating step. Fix/unfix appropriate number of variables "
                "to result in zero degrees of freedom for initialization."
            )

        # 1.3) activate final condition constraint and solve entire step

        blk.heating.fc_temperature_eq.activate()
        blk.heating.time.unfix()

        # check degrees of freedom and solve
        if degrees_of_freedom(blk.heating) == 0:
            self._step_initialize(cycle_step=blk.heating)
        else:
            raise InitializationError(
                "Degrees of freedom is not zero during initialization of "
                "heating step. Fix/unfix appropriate number of variables "
                "to result in zero degrees of freedom for initialization."
            )

        # 2 - solve cooling step

        # 2.1) fix mole fraction at end of heating step

        vars_lst_cooling = ["mole_frac_heating_end"]
        cons_lst_cooling = ["mole_frac_heating_end_eq"]

        self._calculate_and_fix_variable_from_constraint(
            blk, variable_list=vars_lst_cooling, constraint_list=cons_lst_cooling
        )

        # 2.2) initial solution using false position method

        # deactivate final condition constraint and fix time
        blk.cooling.fc_temperature_eq.deactivate()
        blk.cooling.time.fix(self.cooling_time_guess * units.s)

        # check degrees of freedom and solve
        if degrees_of_freedom(blk.cooling) == 0:
            self._false_position_method(
                blk,
                cycle_step=blk.cooling,
                t_guess=self.cooling_time_guess,
            )
        else:
            raise InitializationError(
                "Degrees of freedom is not zero during initialization of "
                "cooling step. Fix/unfix appropriate number of variables "
                "to result in zero degrees of freedom for initialization."
            )

        # 2.3) activate final condition constraint and solve entire model

        blk.cooling.fc_temperature_eq.activate()
        blk.cooling.time.unfix()

        # check degrees of freedom and solve
        if degrees_of_freedom(blk.cooling) == 0:
            self._step_initialize(cycle_step=blk.cooling)
        else:
            raise InitializationError(
                "Degrees of freedom is not zero during initialization of "
                "cooling step. Fix/unfix appropriate number of variables "
                "to result in zero degrees of freedom for initialization."
            )

        # 3 - solve pressurization step

        # 3.1) fix mole fraction, temperature, pressure and loadings at
        # end of cooling step

        vars_lst_pressurization = [
            "mole_frac_cooling_end",
            "pressure_cooling_end",
            "loading_cooling_end",
        ]

        cons_lst_pressurization = [
            "mole_frac_cooling_end_eq",
            "pressure_cooling_end_eq",
            "loading_cooling_end_eq",
        ]

        self._calculate_and_fix_variable_from_constraint(
            blk,
            variable_list=vars_lst_pressurization,
            constraint_list=cons_lst_pressurization,
        )

        if blk.calculate_beds:
            # if "calculate_beds" is True, there is an extra variable for
            # "velocity_in" and it needs to be fixed to initialize the
            # pressurization and adsorption step
            velocity_fixed = blk.velocity_in.fixed
            if not velocity_fixed:
                blk.velocity_in.fix()
                blk.pressure_drop.unfix()

        # 3.2) check degrees of freedom and solve

        if degrees_of_freedom(blk.pressurization) == 0:
            self._step_initialize(
                cycle_step=blk.pressurization,
            )
        else:
            raise InitializationError(
                "Degrees of freedom is not zero during initialization of "
                "pressurization step. Fix/unfix appropriate number of "
                "variables to result in zero degrees of freedom for "
                "initialization."
            )

        # 4 - solve adsorption step

        # 4.1) fix mole fraction and loadings at end of pressurization step

        vars_lst_adsorption = [
            "mole_frac_pressurization_end",
            "loading_pressurization_end",
        ]

        cons_lst_adsorption = [
            "mole_frac_pressurization_end_eq",
            "loading_pressurization_end_eq",
        ]

        self._calculate_and_fix_variable_from_constraint(
            blk, variable_list=vars_lst_adsorption, constraint_list=cons_lst_adsorption
        )

        # 4.2) check degrees of freedom and solve

        if degrees_of_freedom(blk.adsorption) == 0:
            self._step_initialize(cycle_step=blk.adsorption)
        else:
            raise InitializationError(
                "Degrees of freedom is not zero during initialization of "
                "adsorption step. Fix/unfix appropriate number of variables "
                "to result in zero degrees of freedom for initialization."
            )

        # 5 - solve entire fixed bed TSA model

        # 5.1) unfix variables and activate constraints that were fixed and
        # deactivated during individual steps
        from_json(blk, sd=tsa_state, wts=StoreState)

        if blk.calculate_beds:
            calculate_variable_from_constraint(blk.velocity_in, blk.pressure_drop_eq)
        else:
            calculate_variable_from_constraint(blk.pressure_drop, blk.pressure_drop_eq)

        # 5.2) deactivate compressor
        if blk.config.compressor:
            blk.compressor.deactivate()

        # 5.3) deactivate steam calculation constraints
        if blk.config.steam_calculation != SteamCalculationType.none:
            if blk.config.steam_calculation == SteamCalculationType.rigorous:
                blk.steam_heater.deactivate()
            blk.flow_mass_steam_eq.deactivate()
            blk.flow_mass_steam.fix()

        # 5.4) check degrees of freedom and solve
        if degrees_of_freedom(blk) == 0:

            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = self.opt.solve(blk, tee=slc.tee)

            if check_optimal_termination(res):
                init_log.info(
                    "Initialization of fixed bed TSA model "
                    "completed {}.".format(idaeslog.condition(res))
                )
            else:
                _log.warning(
                    "Initialization of fixed bed TSA model "
                    "Failed {}.".format(blk.name)
                )
        else:
            raise InitializationError(
                "Degrees of freedom is not zero during initialization of "
                "fixed bed TSA model. Fix/unfix appropriate number of "
                "variables to result in zero degrees of freedom for "
                "initialization."
            )

        # 6 - solve compressor unit
        if blk.config.compressor:

            # initialization of compressor
            init_log.info("Starting initialization of compressor.")

            # activate compressor model
            blk.compressor.activate()

            # 6.1) fix state at inlet of compressor to match states in
            # exhaust gas stream in the CCS system. Fix pressure drop in
            # compressor unit.
            for t in blk.flowsheet().time:
                for i in blk.config.compressor_properties.component_list:
                    calculate_variable_from_constraint(
                        blk.compressor.unit.inlet.flow_mol_comp[t, i],
                        blk.compressor.flow_mol_in_compressor_eq[t, i],
                    )
                calculate_variable_from_constraint(
                    blk.compressor.unit.inlet.temperature[t],
                    blk.compressor.temperature_in_compressor_eq[t],
                )
                calculate_variable_from_constraint(
                    blk.compressor.unit.inlet.pressure[t],
                    blk.compressor.pressure_in_compressor_eq[t],
                )
                calculate_variable_from_constraint(
                    blk.compressor.unit.deltaP[t],
                    blk.compressor.pressure_drop_tsa_compressor_eqn[t],
                )

            blk.compressor.unit.inlet.flow_mol_comp[:, :].fix()
            blk.compressor.unit.inlet.temperature[:].fix()
            blk.compressor.unit.inlet.pressure[:].fix()
            blk.compressor.unit.deltaP[:].fix()

            # deactivate related constraints for fixed variables.
            blk.compressor.flow_mol_in_compressor_eq.deactivate()
            blk.compressor.temperature_in_compressor_eq.deactivate()
            blk.compressor.pressure_in_compressor_eq.deactivate()
            blk.compressor.pressure_drop_tsa_compressor_eqn.deactivate()

            # 6.2) check degrees of freedom and solve
            if degrees_of_freedom(blk.compressor) == 0:

                # TODO: switch to new initialization method when implemented
                # for FlueGasStateBlock
                blk.compressor.unit.initialize()

                # re-solve compressor model
                with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                    res = self.opt.solve(blk.compressor, tee=slc.tee)

                if check_optimal_termination(res):
                    init_log.info(
                        "Initialization of compressor completed {}.".format(
                            idaeslog.condition(res)
                        )
                    )
                else:
                    _log.warning(
                        "Initialization of compressor Failed {}.".format(
                            blk.compressor.unit.name
                        )
                    )
            else:
                raise InitializationError(
                    "Degrees of freedom is not zero during initialization of "
                    "compressor model. Fix/unfix appropriate number of "
                    "variables to result in zero degrees of freedom for "
                    "initialization."
                )

        # 7 - solve steam calculation
        if blk.config.steam_calculation != SteamCalculationType.none:

            if blk.config.steam_calculation == SteamCalculationType.rigorous:

                # initialization of steam heater
                init_log.info(
                    "Starting initialization of heater model for steam " "calculation."
                )

                # activate steam heater model
                blk.steam_heater.activate()

                # 7.1) deactivate constraints for total saturation condition
                #      and heat_duty_heater_eq
                blk.steam_heater.unit.vapor_frac_out_eq.deactivate()
                blk.steam_heater.unit.heat_duty_heater_eq.deactivate()

                # 7.2) assume a dumy inlet flow rate and heat duty and fix them
                Fin = 100 * units.mol / units.s
                Q = -500000 * units.W
                blk.steam_heater.unit.inlet.flow_mol[0].fix(Fin)
                blk.steam_heater.unit.heat_duty.fix(Q)

                # 7.3) check degrees of freedom and solve
                if degrees_of_freedom(blk.steam_heater) == 0:

                    # initialize steam heater
                    heater_initializer = self.get_submodel_initializer(
                        blk.steam_heater.unit
                    )
                    heater_initializer.initialize(blk.steam_heater.unit)

                else:
                    raise InitializationError(
                        "Degrees of freedom is not zero during initialization "
                        "of heater model. Fix/unfix appropriate number of "
                        "variables to result in zero degrees of freedom for "
                        "initialization."
                    )

                # 7.4) activate constraint for total saturation condition and
                #      unfix heat duty
                blk.steam_heater.unit.vapor_frac_out_eq.activate()
                blk.steam_heater.unit.heat_duty.unfix()

                # 7.5) solve model for total saturation at outlet
                if degrees_of_freedom(blk.steam_heater) == 0:

                    init_log.info_high(
                        "Starting initialization of heater model "
                        "for total saturation at outlet."
                    )

                    with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                        res = self.opt.solve(blk.steam_heater, tee=slc.tee)

                    if check_optimal_termination(res):
                        init_log.info_high(
                            "Initialization of heater model "
                            "for total saturation at outlet "
                            "completed {}.".format(idaeslog.condition(res))
                        )
                    else:
                        _log.warning(
                            "Initialization of heater model for "
                            "total saturation at outlet Failed {}.".format(
                                blk.steam_heater.unit.name
                            )
                        )
                else:
                    raise InitializationError(
                        "Degrees of freedom is not zero during initialization "
                        "of heater model for total saturation at outlet. "
                        "Fix/unfix appropriate number of variables to result "
                        "in zero degrees of freedom for initialization."
                    )

                # 7.6) unfix inlet flow rate and fix heat duty
                blk.steam_heater.unit.inlet.flow_mol[0].unfix()
                calculate_variable_from_constraint(
                    blk.steam_heater.unit.heat_duty[0],
                    blk.steam_heater.unit.heat_duty_heater_eq,
                )
                blk.steam_heater.unit.heat_duty.fix()

                # 7.7) solve model for steam flow rate
                if degrees_of_freedom(blk.steam_heater) == 0:

                    init_log.info_high(
                        "Starting initialization of heater model "
                        "for total steam flow rate."
                    )

                    with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                        res = self.opt.solve(blk.steam_heater, tee=slc.tee)

                    if check_optimal_termination(res):
                        init_log.info_high(
                            "Initialization of heater model "
                            "for total steam flow rate "
                            "completed: {}.".format(idaeslog.condition(res))
                        )
                        init_log.info(
                            "Initialization of heater model "
                            "for steam calculation completed {}.".format(
                                idaeslog.condition(res)
                            )
                        )
                    else:
                        _log.warning(
                            "Initialization of heater model for "
                            "total steam flow rate Failed {}.".format(
                                blk.steam_heater.unit.name
                            )
                        )
                else:
                    raise InitializationError(
                        "Degrees of freedom is not zero during initialization "
                        "of heater model for total steam flow rate. "
                        "Fix/unfix appropriate number of variables to result "
                        "in zero degrees of freedom for initialization."
                    )

            calculate_variable_from_constraint(
                blk.flow_mass_steam, blk.flow_mass_steam_eq
            )

        # 8 - solve fixed bed TSA model, steam calculation constraints and
        #     compressor unit simultaneously

        if blk.config.compressor:

            # 8.1) unfix state that were fixed in 6.1
            blk.compressor.unit.inlet.flow_mol_comp[:, :].unfix()
            blk.compressor.unit.inlet.temperature[:].unfix()
            blk.compressor.unit.inlet.pressure[:].unfix()
            blk.compressor.unit.deltaP[:].unfix()

            # activate constraints that were deactivated in 6.1
            blk.compressor.flow_mol_in_compressor_eq.activate()
            blk.compressor.temperature_in_compressor_eq.activate()
            blk.compressor.pressure_in_compressor_eq.activate()
            blk.compressor.pressure_drop_tsa_compressor_eqn.activate()

        if blk.config.steam_calculation != SteamCalculationType.none:

            # 8.2 unfix variables and activate constraints that were fixed
            #     and deactivated in 7
            if blk.config.steam_calculation == SteamCalculationType.rigorous:
                blk.steam_heater.unit.heat_duty_heater_eq.activate()
                blk.steam_heater.unit.heat_duty.unfix()

            blk.flow_mass_steam_eq.activate()
            blk.flow_mass_steam.unfix()

        # 8.3) check degrees of freedom and solve
        if (
            blk.config.compressor
            or blk.config.steam_calculation != SteamCalculationType.none
        ):
            if degrees_of_freedom(blk) == 0:

                with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                    res = self.opt.solve(blk, tee=slc.tee)
                if (
                    blk.config.compressor
                    and blk.config.steam_calculation != SteamCalculationType.none
                ):
                    if check_optimal_termination(res):
                        init_log.info(
                            "Initialization of fixed bed TSA, steam "
                            "calculation and compressor models completed {}"
                            ".".format(idaeslog.condition(res))
                        )
                    else:
                        _log.warning(
                            "Initialization of fixed bed TSA, "
                            "steam calculation and compressor "
                            "models Failed {}.".format(blk.name)
                        )
                if (
                    blk.config.compressor
                    and blk.config.steam_calculation != SteamCalculationType.none
                ):
                    if check_optimal_termination(res):
                        init_log.info(
                            "Initialization of fixed bed TSA "
                            "and compressor models completed {}.".format(
                                idaeslog.condition(res)
                            )
                        )
                    else:
                        _log.warning(
                            "Initialization of fixed bed TSA "
                            "and compressor models Failed {}.".format(blk.name)
                        )
                if (
                    not blk.config.compressor
                    and blk.config.steam_calculation != SteamCalculationType.none
                ):
                    if check_optimal_termination(res):
                        init_log.info(
                            "Initialization of fixed bed TSA and steam "
                            "calculation models completed {}"
                            ".".format(idaeslog.condition(res))
                        )
                    else:
                        _log.warning(
                            "Initialization of fixed bed TSA "
                            "and steam calculation "
                            "models Failed {}.".format(blk.name)
                        )
            else:
                if (
                    blk.config.compressor
                    and blk.config.steam_calculation != SteamCalculationType.none
                ):
                    raise InitializationError(
                        "Degrees of freedom is not zero during initialization "
                        "of fixed bed TSA, steam calculation and compressor "
                        "models. Fix/unfix appropriate number of variables to "
                        "result in zero degrees of freedom for "
                        "initialization."
                    )
                if (
                    blk.config.compressor
                    and blk.config.steam_calculation == SteamCalculationType.none
                ):
                    raise InitializationError(
                        "Degrees of freedom is not zero during initialization "
                        "of fixed bed TSA and compressor "
                        "models. Fix/unfix appropriate number of variables to "
                        "result in zero degrees of freedom for "
                        "initialization."
                    )
                if (
                    not blk.config.compressor
                    and blk.config.steam_calculation != SteamCalculationType.none
                ):
                    raise InitializationError(
                        "Degrees of freedom is not zero during initialization "
                        "of fixed bed TSA and steam calculation "
                        "models. Fix/unfix appropriate number of variables to "
                        "result in zero degrees of freedom for "
                        "initialization."
                    )

    def _step_initialize(self, cycle_step=None):
        """
        Initialization routine for TSA cycle steps.

        Keyword Arguments:
            outlvl    : output level of initialisation routine
            solver    : str indicating which solver to use during
                        initialization
            optarg    : dictionary with solver options
            cycle_step: block model for cycle step

        """
        # set up logger
        init_log = idaeslog.getInitLogger(
            cycle_step.name, self.get_output_level(), tag="unit"
        )
        solve_log = idaeslog.getSolveLogger(
            cycle_step.name, self.get_output_level(), tag="unit"
        )

        # initialization of cycle steps
        init_log.info(
            "Starting initialization of "
            + str(cycle_step).rsplit(".", maxsplit=1)[-1]
            + " step."
        )

        # solve cycle step
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = self.opt.solve(cycle_step, tee=slc.tee)

        if check_optimal_termination(res):
            init_log.info(
                "Initialization of "
                + str(cycle_step).rsplit(".", maxsplit=1)[-1]
                + " step completed {}.".format(idaeslog.condition(res))
            )
        else:
            _log.warning(
                "Initialization of "
                + str(cycle_step).rsplit(".", maxsplit=1)[-1]
                + " step Failed {}.".format(cycle_step.name)
            )

    def _false_position_method(self, blk, cycle_step=None, t_guess=None):
        """
        False position method to provide initial solution for TSA cycle steps.

        Keyword Arguments:
            outlvl    : output level of initialisation routine
            solver    : str indicating which solver to use during
                        initialization
            optarg    : dictionary with solver options
            cycle_step: block model for cycle step
            x0        : initial guess for time

        """
        # set up logger
        init_log = idaeslog.getInitLogger(
            cycle_step.name, self.get_output_level(), tag="unit"
        )
        solve_log = idaeslog.getSolveLogger(
            cycle_step.name, self.get_output_level(), tag="unit"
        )

        # initial interval containing a root to apply false position method
        init_log.info_high(
            "Initialization of "
            + str(cycle_step).rsplit(".", maxsplit=1)[-1]
            + " step: step 1a: finding initial interval containing a root"
            " to apply false position method"
        )

        # fix time to initial guess
        x0 = t_guess
        cycle_step.time.fix(x0)

        # counter
        count = 1

        # solve model
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = self.opt.solve(cycle_step, tee=slc.tee)

        if check_optimal_termination(res):
            init_log.info_high(
                "Initialization of "
                + str(cycle_step).rsplit(".", maxsplit=1)[-1]
                + " step: step 1a - iteration {0}, completed {1}.".format(
                    count, idaeslog.condition(res)
                )
            )
        else:
            _log.warning(
                "Initialization of "
                + str(cycle_step).rsplit(".", maxsplit=1)[-1]
                + " step: step 1a - iteration {0}, Failed {1}.".format(
                    count, cycle_step.name
                )
            )

        # save solution in initial guess, f(x0)
        if str(cycle_step).rsplit(".", maxsplit=1)[-1] == "heating":
            f_x0 = value(blk.temperature_desorption - cycle_step.temperature[1])
        else:
            f_x0 = value(cycle_step.temperature[1] - blk.temperature_adsorption)
        f_x0_start = f_x0

        # iterate to find an interval with the root
        condition = True
        while condition:
            # update counter and x0
            count += 1

            if f_x0_start >= 0:
                x_new = 1.2 * x0
            else:
                x_new = x0 / 2.0

            # fix time to x_new guess
            cycle_step.time.fix(x_new)

            # solve model
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = self.opt.solve(cycle_step, tee=slc.tee)

            if check_optimal_termination(res):
                init_log.info_high(
                    "Initialization of "
                    + str(cycle_step).rsplit(".", maxsplit=1)[-1]
                    + " step: step 1a - iteration {0}, completed {1}.".format(
                        count, idaeslog.condition(res)
                    )
                )
            else:
                _log.warning(
                    "Initialization of "
                    + str(cycle_step).rsplit(".", maxsplit=1)[-1]
                    + " step: step 1a - iteration {0}, Failed {1}.".format(
                        count, cycle_step.name
                    )
                )

            # save solution in new initial guess, f(x_new)
            if str(cycle_step).rsplit(".", maxsplit=1)[-1] == "heating":
                f_x_new = value(blk.temperature_desorption - cycle_step.temperature[1])
            else:
                f_x_new = value(cycle_step.temperature[1] - blk.temperature_adsorption)

            # set up new interval until find one containing the root
            if f_x0_start >= 0:
                if f_x_new >= 0:
                    x0 = x_new
                    f_x0 = f_x_new
                else:
                    x1 = x_new
                    f_x1 = f_x_new
                    condition = False
            else:
                if f_x_new >= 0:
                    x1 = x_new
                    f_x1 = f_x_new
                    condition = False
                else:
                    x0 = x_new
                    f_x0 = f_x_new

        # implementing false position method
        init_log.info_high(
            "Initialization of "
            + str(cycle_step).rsplit(".", maxsplit=1)[-1]
            + " step: step 1b: implementing false position method"
        )

        # counter
        count = 1
        condition = True

        # check condition to stop
        while condition:

            # compute new approximated root as x2
            x2 = x0 - (x1 - x0) * f_x0 / (f_x1 - f_x0)

            # fix time to x_2 guess
            cycle_step.time.fix(x2)

            # solve model
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = self.opt.solve(cycle_step, tee=slc.tee)

            if check_optimal_termination(res):
                init_log.info_high(
                    "Initialization of "
                    + str(cycle_step).rsplit(".", maxsplit=1)[-1]
                    + " step: step 1b - iteration {0}, completed {1}.".format(
                        count, idaeslog.condition(res)
                    )
                )
            else:
                _log.warning(
                    "Initialization of "
                    + str(cycle_step).rsplit(".", maxsplit=1)[-1]
                    + " step: step 1b - iteration {0}, Failed {1}.".format(
                        count, cycle_step.name
                    )
                )

            # save solution in new x_2, f(x_2)
            if str(cycle_step).rsplit(".", maxsplit=1)[-1] == "heating":
                f_x2 = value(blk.temperature_desorption - cycle_step.temperature[1])
            else:
                f_x2 = value(cycle_step.temperature[1] - blk.temperature_adsorption)

            # check if f(x_0)*f(x_2) is negative
            if f_x0 * f_x2 < 0:
                x1 = x2
            else:
                x0 = x2

            # update counter and set up new condition |f(x_2)| > error
            count += 1
            condition = (f_x2**2) ** 0.5 > 1

    def _calculate_and_fix_variable_from_constraint(
        self, obj, variable_list=None, constraint_list=None
    ):
        """
        Method to calculate variables from a constraints, then fix the
        variables and deactivate the constraints.

        Keyword Arguments:
            variable_list: a list of str with names of variables to be
                calculated.
            constraint_list: a list of str with names of constraints
                used to calculate the variables in variable_list.

        """
        if variable_list is None:
            variable_list = []
        if constraint_list is None:
            constraint_list = []

        v_list = []
        c_list = []

        for v in obj.component_objects(Var, descend_into=True):
            if v.local_name in variable_list:
                v_list.append(v)
        for c in obj.component_objects(Constraint, descend_into=True):
            if c.local_name in constraint_list:
                c_list.append(c)

        v_c_tuple = list(zip(v_list, c_list))

        for k in v_c_tuple:
            if k[1].dim() == 0:
                calculate_variable_from_constraint(k[0], k[1])
            elif k[1].dim() == 1:
                for var_index in k[1].index_set():
                    calculate_variable_from_constraint(k[0][var_index], k[1][var_index])

        for v in v_list:
            v.fix()

        for c in c_list:
            c.deactivate()
