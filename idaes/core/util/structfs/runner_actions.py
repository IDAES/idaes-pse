###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
###############################################################################
"""
Predefined Actions for the generic Runner.
"""
# stdlib
from collections import defaultdict
from collections.abc import Callable
import time
from typing import Union, Optional

# third-party
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.base.unit_model import ProcessBlockData
import pandas as pd
from pyomo.network.port import ScalarPort

# package
from .runner import Action
from .fsrunner import FlowsheetRunner


class Timer(Action):
    """Simple step/run timer action."""

    def __init__(self, runner, **kwargs):
        """Constructor.

        Args:
            runner: Associated Runner object
            kwargs: Additional optional arguments for Action constructor

        Attributes:
            step_times: Dict with key step name and value a list of
                        timings for that step
            run_times: List of timings for a run (sequence of steps)
        """
        super().__init__(runner, **kwargs)
        self.step_times = defaultdict(list)
        self.run_times = []
        self._step_begin = {}
        self._run_begin = None
        self._step_order = []

    def before_step(self, step_name):
        self._step_begin[step_name] = time.time()
        if len(self.run_times) == 0:
            self._step_order.append(step_name)

    def after_step(self, step_name):
        t1 = time.time()
        t0 = self._step_begin.get(step_name, None)
        if t0 is None:
            self.log.warning(f"Timer: step {step_name} end without begin")
        else:
            dt = t1 - t0
            self.step_times[step_name].append(dt)
            self._step_begin[step_name] = None

    def before_run(self):
        self._run_begin = time.time()

    def after_run(self):
        t1 = time.time()
        if self._run_begin is None:
            self.log.warning("Timer: run end without begin")
        else:
            dt = t1 - self._run_begin
            self.run_times.append(dt)
            self._run_begin = None

    def summary(self):
        data = []
        for i, run_time in enumerate(self.run_times):
            step_times = [(k, self.step_times[k][i]) for k in self._step_order]
            step_total = sum((item[1] for item in step_times))
            data.append(
                {
                    "run": run_time,
                    "steps": step_times,
                    "inclusive": step_total,
                    "exclusive": run_time - step_total,
                }
            )
        return data


# Hold degrees of freedom for one FlowsheetRunner 'step'
# {key=component: value=dof}
UnitDofType = dict[str, int]


class UnitDofChecker(Action):
    """Check degrees of freedom on unit models.

    After a (caller-named) step or steps, check the degrees
    of freedom on each unit model by the method of
    fixing the inlet, applying the `degrees_of_freedom()` function,
    and unfixing the inlet again. The calculated values are
    saved and passed to an optional caller-provided function.

    At the end of a run, the degrees of freedom for the entire
    model are checked, saved, and passed to an optional function.
    """

    def __init__(
        self,
        runner: FlowsheetRunner,
        flowsheet: str,
        steps: Union[str, list[str]],
        step_func: Optional[Callable[[str, UnitDofType], None]] = None,
        run_func: Optional[Callable[[dict[str, UnitDofType], int], None]] = None,
        **kwargs,
    ):
        """Constructor.

        Args:
            runner: Associated Runner object (provided by `add_action`)
            flowsheet: Variable name for flowsheet, e.g. "fs"
            steps: Step or steps at which to run the checking action
            step_func: Function to call with calculated DoF values for one step.
                  Takes name of step and dictionary with per-unit degrees of freedom
                  (see `UnitDofType` alias).
            run_func: Function to call with calculated DoF values for each step, as well
                  as overall model DoF.
            kwargs: Additional optional arguments for Action constructor

        Raises:
            ValueError: if `steps` list is empty, or no callback functions provided
        """
        super().__init__(runner, **kwargs)
        if hasattr(steps, "lower"):  # string-like
            self._steps = {steps}
        else:  # assume it is list-like
            if len(steps) == 0:
                raise ValueError("At least one step name must be provided")
            self._steps = set(steps)
        self._steps_dof: dict[str, UnitDofType] = {}
        self._model_dof = None
        self._step_func, self._run_func = step_func, run_func
        self._fs = flowsheet

    def after_step(self, step_name: str):
        step_name = self._runner._norm_name(step_name)
        if step_name not in self._steps:
            return

        fs = self._get_flowsheet()

        units_dof = {}
        for unit in fs.component_objects(descend_into=True):
            if self._is_unit_model(unit):
                units_dof[unit.name] = self._get_dof(unit)
        self._steps_dof[step_name] = units_dof  # save
        if self._step_func:
            self._step_func(step_name, units_dof)

    def after_run(self):
        fs = self._get_flowsheet()
        model_dof = degrees_of_freedom(fs)
        self._model_dof = model_dof
        if self._run_func:
            self._run_func(self._steps_dof, model_dof)

    def _get_flowsheet(self):
        m = self._runner.model
        if self._fs:
            return getattr(m, self._fs)
        return m

    @staticmethod
    def _is_unit_model(block):
        return isinstance(block, ProcessBlockData)

    def as_dataframe(self) -> pd.DataFrame:
        """Format per-step DoF as a Pandas `DataFrame`.

        Returns:
            DataFrame: Step (str), Unit (str), DoF (int)
        """
        step_names, unit_names, dofs = [], [], []

        # add DoF for each step
        for sn, data in self._steps_dof.items():
            for un, dof in data.items():
                step_names.append(sn)
                unit_names.append(un)
                dofs.append(dof)

        # add model DoF
        step_names.append("RUN")
        unit_names.append(self._fs)
        dofs.append(self._model_dof)

        return pd.DataFrame(
            {"after_step": step_names, "unit_name": unit_names, "dof": dofs}
        )

    def get_unit_dof(self, step_name: str) -> UnitDofType:
        """Get DoF for each unit, as measured after the given step.

        Args:
            step_name: Step for which to get the per-unit degrees of freedom.

        Returns:
            UnitDofType

        Raises:
            KeyError: If `step_name` is unknown, or has no data
            ValueError: There is no degrees_of_freedom data at all
        """
        if not self._steps_dof:
            raise ValueError("No degrees of freedom have been calculated")
        if step_name not in self._steps:
            raise KeyError(
                f"Unknown step. name={step_name} known={','.join(self._steps)}"
            )
        return self._steps_dof[step_name]

    def steps(self, only_with_data: bool = False) -> list[str]:
        """Get list of steps for which unit degrees of freedom are calculated.

        Args:
            only_with_data: If True, do not return steps with no data

        Returns:
            list of step names
        """
        if only_with_data:
            return [s for s in self._steps if s in self._steps_dof]
        return list(self._steps)

    @staticmethod
    def _get_dof(block, fix_inlets: bool = True):
        name = block.name
        if fix_inlets:
            inlets = [
                c
                for c in block.component_objects(descend_into=False)
                if isinstance(c, ScalarPort)
                and (c.name.endswith("inlet") or c.name.endswith("recycle"))
            ]
            free_me = []
            for inlet in inlets:
                if not inlet.is_fixed():
                    inlet.fix()
                    free_me.append(inlet)

        dof = degrees_of_freedom(block)

        if fix_inlets:
            for inlet in free_me:
                inlet.free()

        return dof
