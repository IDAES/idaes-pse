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

from collections import defaultdict
import time

from pyomo.network.port import ScalarPort
from idaes.core.util.model_statistics import degrees_of_freedom

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
            self.log.warning(f"Timer: run end without begin")
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


class DOFChecker(Action):
    def __init__(self, runner: FlowsheetRunner, **kwargs):
        """Constructor.

        Args:
            runner: Associated Runner object
            kwargs: Additional optional arguments for Action constructor
        """
        super().__init__(runner, **kwargs)

    def after_step(self, step_name: str):
        m = self._runner.model
        if step_name == "special step":
            unit_dof = {}
            for unit in m.component_objects(descend_into=False):
                # XXX: Check whether really IDAES unit model
                unit_dof[unit.name] = self._get_dof(unit)
        # XXX: check that DOF is what was wanted?


def _get_dof(self, block, fix_inlets: bool = True):
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
    self._component_dof[name] = dof
    if fix_inlets:
        for inlet in free_me:
            inlet.free()
    # if expected is not None:
    #     assert dof == expected, (
    #         f"Degrees of freedom for component "
    #         f"{name} ({dof}) does not match "
    #         f"expected ({expected})"
    #     )
