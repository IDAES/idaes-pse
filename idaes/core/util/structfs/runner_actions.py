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
from collections.abc import Callable
from io import StringIO
import logging
import re
import sys
import time
from typing import Union, Optional

# third-party
from pyomo.network.port import ScalarPort
from pyomo.core.base.var import IndexedVar
from pyomo.core.base.param import IndexedParam
import pyomo.environ as pyo
from pydantic import BaseModel, Field

try:
    from idaes_connectivity.base import Connectivity, Mermaid
except ImportError:
    Connectivity = None

# package
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.base.unit_model import ProcessBlockData
from .runner import Action
from .fsrunner import FlowsheetRunner


class Timer(Action):
    """Simple step/run timer action."""

    class Report(BaseModel):
        """Report returned by report() method."""

        # {"step_name": <float time>, ..} for each step
        timings: dict[str, float] = Field(default={})

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

    def summary(self, stream=None, run_idx=-1) -> str | None:
        """Summary of the timings.

        Args:
            stream: Output stream, with `write()` method. Return a string if None.
            run_idx: Index of run, -1 meaning "last one"

        Returns:
            str: If output stream was None, the text summary; otherwise None
        """
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
                stream.write(fmt.format(s=s, t=t, p=t / ttot * 100))

        stream.write(f"\nTotal time: {d['run']:.3f} s\n")

        if isinstance(stream, StringIO):
            return stream.getvalue()

        return None

    def _ipython_display_(self):
        print(self.summary())

    def report(self) -> Report:
        """Report the timings.

        Returns:
            The report object
        """
        rpt = self.Report(timings=self.step_times[-1].copy())
        return rpt


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

    class Report(BaseModel):
        """Report on degrees of freedom in a model."""

        steps: dict[str, UnitDofType] = Field(default={})
        model: int = Field(default=0)

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
        """Actions performed after a run."""
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
        """Readable summary of the degrees of freedom.

        Args:
            stream: Output stream, with `write()` method. Return a string if None.
            step: Specific step to summarize, otherwise all steps.

        Returns:
            The summary as a string if `stream` was None, otherwise None
        """
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
            for step in self._runner.list_steps():
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

    def report(self) -> Report:
        """Machine-readable report of degrees of freedom.

        Returns:
            Report object
        """
        return self.Report(steps=self.get_dof(), model=self.get_dof_model())

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


class CaptureSolverOutput(Action):
    """Capture the solver output."""

    def __init__(self, runner, **kwargs):
        super().__init__(runner, **kwargs)
        self._logs = {}
        self._solver_out = None

    def before_step(self, step_name: str):
        """Action performed before the step."""
        if self._is_solve_step(step_name):
            self._solver_out = StringIO()
            self._save_stdout, sys.stdout = sys.stdout, self._solver_out

    def after_step(self, step_name: str):
        """Action performed after the step."""
        if self._solver_out is not None:
            self._logs[step_name] = self._solver_out.getvalue()
            self._solver_out = None
            sys.stdout = self._save_stdout

    def _is_solve_step(self, name: str):
        return name.startswith("solve")

    def report(self) -> dict:
        """Machine-readable report with solver output.

        Returns:
            Report dict, {'solver_logs': "<text-log>"}
        """
        return {"solver_logs": self._logs}


class ModelVariables(Action):
    """Extract and format model variables."""

    VAR_TYPE, PARAM_TYPE = "V", "P"

    class Report(BaseModel):
        """Report for ModelVariables.

        The value of `tree` is a tree represented as a nested dict,
        where each sub-component has a class and sub-component key, and
        values (variables or parameters), which are leaves of the tree
        (i.e., do not have sub-components), have a type and value key.

        A given component is represented a dict like:
        `{'name': {'t': 'class.name', 'sub': { ...sub-components... }}}`

        For the parameter/variable value node, the 'sub' key will be gone and
        the value for 't' will be a code 'P' or 'V' for parameter or variable.
        The values are given as a list of items under the 'v' key:
        - For a parameter, each is `[index, value]`.
        - For a variable, each is `[index, value, fixed, stale, lb, ub]`,
          where 'lb' means lower bound and 'ub' means upper bound.
          The lb and ub can be None (or, in JSON, null) if no bound exists.
          Fixed and stale are boolean values corresponding to the Pyomo variable
          attributes of the same name.
        For scalar parameters/variables, there is only 1 item
        and its index is None/null. Otherwise, it is an indexed parameter or variable.
        Thus, the value dict will look like one of:

        Parameter: `{'name': {'t': 'P', 'v': [[<index>, <value], ...]}}`
        Variable: `{'name': {'t': 'V', 'v': [[<index>, <value>, True/False,
                                              True/False, <lb>, <ub>], ...]}}`

        This structure is designed so the value of 't' can determine the function
        to call to process the contents of a node.
        """

        #: Tree of variables
        tree: dict = Field(default={})

    def __init__(self, runner, **kwargs):
        assert isinstance(runner, FlowsheetRunner)  # makes no sense otherwise
        super().__init__(runner, **kwargs)

    def after_run(self):
        """Actions performed after the run."""
        self._saved_paths = {}  # fast lookup used in _add_block()
        self.log = logging.getLogger(self.log.name)
        self._extract_vars(self._runner.model)

    def _extract_vars(self, m):
        var_tree = {}

        for c in m.component_objects():
            # get component type
            if self._is_var(c):
                subtype = self.VAR_TYPE
            elif self._is_param(c):
                subtype = self.PARAM_TYPE
            else:
                continue  # ignore other components
            if self._dbg:
                self.log.debug(f"_extract_vars: component {c} ({subtype})")
            # start new block
            b = []
            # add its variables:
            #   - add each value from an indexed var/param,
            #   - this also works for scalars
            for index in c:
                v = c[index]
                value = pyo.value(v, exception=False)
                if subtype == self.VAR_TYPE:
                    # index, value, is-fixed, is-stale, lower-bound, upper-bound
                    item = (index, value, v.fixed, v.stale, v.lb, v.ub)
                else:
                    # index, value
                    item = (index, value)
                b.append(item)
            # add block to tree
            self._add_block(m, var_tree, c.name, b, subtype)

        self._vars = var_tree

    @staticmethod
    def _is_var(c):
        return c.is_variable_type() or isinstance(c, IndexedVar)

    @staticmethod
    def _is_param(c):
        return c.is_parameter_type() or isinstance(c, IndexedParam)

    def _add_block(self, m, tree: dict, name: str, block, subtype):
        assert name
        # get parts of the name, accounting for indexes
        tok_pat = r"[^.\[\]]+(?:\[[^\]]*\])?"
        parts = re.findall(tok_pat, name)
        assert parts
        # insert in tree
        t, prev = tree, None
        cur_cmp = m
        type_key, child_key = "t", "sub"
        # walk down tree to where this belongs
        if self._dbg:
            self.log.debug(f"_add_block: parts={parts}, name={name}")
        cur_path, p, prev = "m", None, {}
        for p in parts:
            cur_path += "." + p
            if p in t:
                cur_cmp = self._saved_paths[cur_path]
            else:
                # get the child component
                if p.endswith("]"):  # indexed
                    cur_cmp, cur_cmp_s = self._indexed_sub(cur_cmp, p)
                else:  # scalar
                    cur_cmp_s = cur_cmp = getattr(cur_cmp, p)
                self._saved_paths[cur_path] = cur_cmp
                try:
                    # prefer IDAES process block class name
                    clazz = str(cur_cmp_s.process_block_class())
                except AttributeError:
                    # fall back to component class name
                    clazz = cur_cmp_s.__class__.__name__
                t[p] = {type_key: self._bare_class(clazz), child_key: {}}
            prev, t = t, t[p][child_key]
        # put value in leaf node
        type_code = "P" if subtype == self.PARAM_TYPE else "V"
        prev[p] = {"v": block, type_key: type_code}

    @classmethod
    def _indexed_sub(cls, block, part):
        """Get both indexed and scalar parts, e.g.,
        if 'part' is "foo[0]" then return (block.foo[0], block.foo).
        """
        pos = part.find("[")
        # scalar part of this subcomponent
        part_s = part[:pos]
        sub_s = getattr(block, part_s)
        # extract index as tuple
        str_idx = part[pos + 1 : -1]
        idx_parts = str_idx.split(",")
        idx_tuple = []
        for ip in idx_parts:
            try:
                idx = int(ip)  # try int
            except ValueError:
                try:
                    idx = float(ip)  # try float
                except ValueError:
                    idx = ip  # use string
            idx_tuple.append(idx)
        # use the index tuple to get the correct component
        sub_i = sub_s[idx_tuple]

        return sub_i, sub_s

    @staticmethod
    def _bare_class(s):
        m = re.match(r"<class '(.*)'>", str(s))
        return m.group(1) if m else s

    def report(self) -> Report:
        """Report containing model variable values."""
        return self.Report(tree=self._vars)


class MermaidDiagram(Action):
    """Action to generate a Mermaid diagram after the run."""

    class Report(BaseModel):
        """Report containing a Mermaid diagram."""

        diagram: list[str]  #: each item is one line

    def __init__(self, runner, **kwargs):
        super().__init__(runner, **kwargs)
        self._images = True
        self._model_root_split = []

    def show_unit_images(self, value: bool):
        """Whether Mermaid displays images for units.

        Args:
            value: If true, display images. Otherwise, don't.
        """
        self._images = bool(value)

    def set_model_root(self, path: str):
        """Set path to root of model to display (default is model itself).

        Args:
            path: Dotted path like "fs" or "fs.component"
        """
        self._model_root_split = path.split(".")

    def after_run(self):
        """Build Mermaid diagram after the run."""
        if Connectivity is None:
            self.diagram = None
        else:
            root = self._runner.model
            for p in self._model_root_split:
                root = getattr(root, p)
            conn = Connectivity(input_model=root)
            self.diagram = Mermaid(conn, component_images=self._images)

    def report(self) -> Report | dict:
        """Report containing the Mermaid diagram.

        Returns:
            Report object if idaes_connectivity is active, otherwise
            an empty dictionary
        """
        if self.diagram is None:
            return {}
        mermaid_lines = self.diagram.write(None).split("\n")
        return self.Report(diagram=mermaid_lines)
