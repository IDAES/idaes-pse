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
Predefined Actions for the generic Runner.
"""
# stdlib
from collections import defaultdict
from collections.abc import Callable
from io import StringIO
import sys
import time
from typing import Union, Optional

# third-party
import pandas as pd
from pyomo.network.port import ScalarPort

# package
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.base.unit_model import ProcessBlockData
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
        self.step_times: list[dict[str, float]] = []
        self.run_times: list[float] = []
        self._run_begin, self._step_begin = None, {}
        self._step_order = runner.list_steps()

    def before_step(self, step_name):
        self._step_begin[step_name] = time.time()

    def after_step(self, step_name):
        t1 = time.time()
        t0 = self._step_begin.get(step_name, None)
        if t0 is None:
            self.log.warning(f"Timer: step {step_name} end without begin")
        else:
            self._cur_step_times[step_name] = t1 - t0
            self._step_begin[step_name] = None

    def before_run(self):
        self._run_begin = time.time()
        self._cur_step_times = {}
        self._step_begin = {}

    def after_run(self):
        t1 = time.time()
        if self._run_begin is None:
            self.log.warning("Timer: run end without begin")
        else:
            self.run_times.append(t1 - self._run_begin)
            self._run_begin = None
            filled_times = {}
            for step in self._runner.list_steps():
                filled_times[step] = self._cur_step_times.get(step, -1)
            self.step_times.append(filled_times)

    def __len__(self):
        return len(self.run_times)

    def get_history(self) -> list[dict]:
        """Summarize timings

        Returns:
            Summary of timings (in seconds) for each run in `run_times`:
              - 'run': Time for the run
              - 'steps': dict of `{<step_name>: <time(sec)>}`
              - 'inclusive': total time spent in the steps
              - 'exclusive': difference between run time and inclusive time
        """
        return [self._get_summary(i) for i in range(0, len(self.run_times))]

    def _get_summary(self, i):
        rt, st = self.run_times[i], self.step_times[i]
        step_total = sum((max(t, 0) for t in st.values()))
        return {
            "run": rt,
            "steps": st,
            "inclusive": step_total,
            "exclusive": rt - step_total,
        }

    def summary(self, stream=None, run_idx=-1) -> str:
        if stream is None:
            stream = StringIO()

        if len(self.run_times) == 0:
            return ""  # nothing to summarize

        d = self._get_summary(run_idx)

        stream.write("Time per step:\n\n")
        slen, ttot = -1, 0
        for s, t in d["steps"].items():
            if t >= 0:
                slen = max(slen, len(s))
                ttot += t
        sfmt = "  {{s:{slen}s}} : {{t:8.3f}}  {{p:4.1f}}%\n"
        for s, t in d["steps"].items():
            if t >= 0:
                fmt = sfmt.format(slen=slen)
                stream.write(fmt.format(s=s, t=t, p=(t / ttot * 100)))

        stream.write(f"\nTotal time: {d['run']:.3f} s\n")

        if isinstance(stream, StringIO):
            return stream.getvalue()

    def _ipython_display_(self):
        print(self.summary())

    def report(self) -> dict:
        return self.step_times[-1].copy()  # most recent run


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
        step_name = self._runner.normalize_name(step_name)
        if step_name not in self._steps:
            self.log.debug(f"Do not check DoF for step: {step_name}")
            return

        fs = self._get_flowsheet()

        model_dof = degrees_of_freedom(self._get_flowsheet())
        units_dof = {self._fs: model_dof}
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

    def summary(self, stream=sys.stdout, step=None):
        if stream is None:
            stream = StringIO()

        def write_step(sdof, indent=4):
            sdof = self._steps_dof[step]
            istr = " " * indent
            unit_names = list(sdof.keys())
            ulen = max((len(u) for u in unit_names))
            dfmt = f"{istr}{{u:{ulen}s}} : {{d}}\n"
            unit_names.sort()
            for unit in unit_names:
                dof = sdof[unit]
                stream.write(dfmt.format(u=unit, d=dof))

        stream.write(f"Degrees of freedom: {self._model_dof}\n\n")
        if step is None:
            stream.write("Degrees of freedom after steps:\n")
            for step in self._runner._steps:
                if step in self._steps_dof:
                    stream.write(f"  {step}:\n")
                    write_step(self._steps_dof[step])
        else:
            write_step(self._steps_dof[step], indent=0)

        if isinstance(stream, StringIO):
            return stream.getvalue()

    def _ipython_display_(self):
        self.summary()

    def get_dof(self) -> dict[str, UnitDofType]:
        """Get degrees of freedom

        Returns:
            dict[str, UnitDofType]: Mapping of step name to per-unit DoF when
               the step completed.
        """
        return self._steps_dof.copy()

    def get_dof_model(self) -> int:
        """Get degrees of freedom for the model.

        Returns:
            int: Last calculated DoF for the model.
        """
        return self._model_dof

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

    def report(self) -> dict:
        return {"steps": self.get_dof(), "model": self.get_dof_model()}

    @staticmethod
    def _get_dof(block, fix_inlets: bool = True):
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
