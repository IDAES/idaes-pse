# {py:mod}`structfs.runner_actions`

```{py:module} structfs.runner_actions
:noindex:
```

```{autodoc2-docstring} structfs.runner_actions
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Timer <structfs.runner_actions.Timer>`
  - ```{autodoc2-docstring} structfs.runner_actions.Timer
    :summary:
    ```
* - {py:obj}`UnitDofChecker <structfs.runner_actions.UnitDofChecker>`
  - ```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker
    :summary:
    ```
* - {py:obj}`CaptureSolverOutput <structfs.runner_actions.CaptureSolverOutput>`
  - ```{autodoc2-docstring} structfs.runner_actions.CaptureSolverOutput
    :summary:
    ```
* - {py:obj}`ModelVariables <structfs.runner_actions.ModelVariables>`
  - ```{autodoc2-docstring} structfs.runner_actions.ModelVariables
    :summary:
    ```
* - {py:obj}`MermaidDiagram <structfs.runner_actions.MermaidDiagram>`
  - ```{autodoc2-docstring} structfs.runner_actions.MermaidDiagram
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`UnitDofType <structfs.runner_actions.UnitDofType>`
  - ```{autodoc2-docstring} structfs.runner_actions.UnitDofType
    :summary:
    ```
````

### API

``````{py:class} Timer(runner, **kwargs)
:canonical: structfs.runner_actions.Timer
:noindex:

Bases: {py:obj}`structfs.runner.Action`

```{autodoc2-docstring} structfs.runner_actions.Timer
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.runner_actions.Timer.__init__
```

`````{py:class} Report(/, **data: typing.Any)
:canonical: structfs.runner_actions.Timer.Report
:noindex:

Bases: {py:obj}`pydantic.BaseModel`

```{autodoc2-docstring} structfs.runner_actions.Timer.Report
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.runner_actions.Timer.Report.__init__
```

````{py:attribute} timings
:canonical: structfs.runner_actions.Timer.Report.timings
:noindex:
:type: dict[str, float]
:value: >
   'Field(...)'

```{autodoc2-docstring} structfs.runner_actions.Timer.Report.timings
```

````

`````

````{py:method} before_step(step_name)
:canonical: structfs.runner_actions.Timer.before_step
:noindex:

```{autodoc2-docstring} structfs.runner_actions.Timer.before_step
```

````

````{py:method} after_step(step_name)
:canonical: structfs.runner_actions.Timer.after_step
:noindex:

```{autodoc2-docstring} structfs.runner_actions.Timer.after_step
```

````

````{py:method} before_run()
:canonical: structfs.runner_actions.Timer.before_run
:noindex:

```{autodoc2-docstring} structfs.runner_actions.Timer.before_run
```

````

````{py:method} after_run()
:canonical: structfs.runner_actions.Timer.after_run
:noindex:

```{autodoc2-docstring} structfs.runner_actions.Timer.after_run
```

````

````{py:method} __len__()
:canonical: structfs.runner_actions.Timer.__len__
:noindex:

```{autodoc2-docstring} structfs.runner_actions.Timer.__len__
```

````

````{py:method} get_history() -> list[dict]
:canonical: structfs.runner_actions.Timer.get_history
:noindex:

```{autodoc2-docstring} structfs.runner_actions.Timer.get_history
```

````

````{py:method} _get_summary(i)
:canonical: structfs.runner_actions.Timer._get_summary
:noindex:

```{autodoc2-docstring} structfs.runner_actions.Timer._get_summary
```

````

````{py:method} summary(stream=None, run_idx=-1) -> str | None
:canonical: structfs.runner_actions.Timer.summary
:noindex:

```{autodoc2-docstring} structfs.runner_actions.Timer.summary
```

````

````{py:method} _ipython_display_()
:canonical: structfs.runner_actions.Timer._ipython_display_
:noindex:

```{autodoc2-docstring} structfs.runner_actions.Timer._ipython_display_
```

````

````{py:method} report() -> Report
:canonical: structfs.runner_actions.Timer.report
:noindex:

```{autodoc2-docstring} structfs.runner_actions.Timer.report
```

````

``````

````{py:data} UnitDofType
:canonical: structfs.runner_actions.UnitDofType
:noindex:
:value: >
   None

```{autodoc2-docstring} structfs.runner_actions.UnitDofType
```

````

``````{py:class} UnitDofChecker(runner: structfs.fsrunner.FlowsheetRunner, flowsheet: str, steps: typing.Union[str, list[str]], step_func: typing.Optional[collections.abc.Callable[[str, structfs.runner_actions.UnitDofType], None]] = None, run_func: typing.Optional[collections.abc.Callable[[dict[str, structfs.runner_actions.UnitDofType], int], None]] = None, **kwargs)
:canonical: structfs.runner_actions.UnitDofChecker
:noindex:

Bases: {py:obj}`structfs.runner.Action`

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker.__init__
```

`````{py:class} Report(/, **data: typing.Any)
:canonical: structfs.runner_actions.UnitDofChecker.Report
:noindex:

Bases: {py:obj}`pydantic.BaseModel`

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker.Report
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker.Report.__init__
```

````{py:attribute} steps
:canonical: structfs.runner_actions.UnitDofChecker.Report.steps
:noindex:
:type: dict[str, structfs.runner_actions.UnitDofType]
:value: >
   'Field(...)'

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker.Report.steps
```

````

````{py:attribute} model
:canonical: structfs.runner_actions.UnitDofChecker.Report.model
:noindex:
:type: int
:value: >
   'Field(...)'

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker.Report.model
```

````

`````

````{py:method} after_step(step_name: str)
:canonical: structfs.runner_actions.UnitDofChecker.after_step
:noindex:

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker.after_step
```

````

````{py:method} after_run()
:canonical: structfs.runner_actions.UnitDofChecker.after_run
:noindex:

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker.after_run
```

````

````{py:method} _get_flowsheet()
:canonical: structfs.runner_actions.UnitDofChecker._get_flowsheet
:noindex:

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker._get_flowsheet
```

````

````{py:method} _is_unit_model(block)
:canonical: structfs.runner_actions.UnitDofChecker._is_unit_model
:noindex:
:staticmethod:

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker._is_unit_model
```

````

````{py:method} summary(stream=sys.stdout, step=None)
:canonical: structfs.runner_actions.UnitDofChecker.summary
:noindex:

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker.summary
```

````

````{py:method} _ipython_display_()
:canonical: structfs.runner_actions.UnitDofChecker._ipython_display_
:noindex:

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker._ipython_display_
```

````

````{py:method} get_dof() -> dict[str, structfs.runner_actions.UnitDofType]
:canonical: structfs.runner_actions.UnitDofChecker.get_dof
:noindex:

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker.get_dof
```

````

````{py:method} get_dof_model() -> int
:canonical: structfs.runner_actions.UnitDofChecker.get_dof_model
:noindex:

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker.get_dof_model
```

````

````{py:method} steps(only_with_data: bool = False) -> list[str]
:canonical: structfs.runner_actions.UnitDofChecker.steps
:noindex:

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker.steps
```

````

````{py:method} report() -> Report
:canonical: structfs.runner_actions.UnitDofChecker.report
:noindex:

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker.report
```

````

````{py:method} _get_dof(block, fix_inlets: bool = True)
:canonical: structfs.runner_actions.UnitDofChecker._get_dof
:noindex:
:staticmethod:

```{autodoc2-docstring} structfs.runner_actions.UnitDofChecker._get_dof
```

````

``````

`````{py:class} CaptureSolverOutput(runner, **kwargs)
:canonical: structfs.runner_actions.CaptureSolverOutput
:noindex:

Bases: {py:obj}`structfs.runner.Action`

```{autodoc2-docstring} structfs.runner_actions.CaptureSolverOutput
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.runner_actions.CaptureSolverOutput.__init__
```

````{py:method} before_step(step_name: str)
:canonical: structfs.runner_actions.CaptureSolverOutput.before_step
:noindex:

```{autodoc2-docstring} structfs.runner_actions.CaptureSolverOutput.before_step
```

````

````{py:method} after_step(step_name: str)
:canonical: structfs.runner_actions.CaptureSolverOutput.after_step
:noindex:

```{autodoc2-docstring} structfs.runner_actions.CaptureSolverOutput.after_step
```

````

````{py:method} _is_solve_step(name: str)
:canonical: structfs.runner_actions.CaptureSolverOutput._is_solve_step
:noindex:

```{autodoc2-docstring} structfs.runner_actions.CaptureSolverOutput._is_solve_step
```

````

````{py:method} report() -> dict
:canonical: structfs.runner_actions.CaptureSolverOutput.report
:noindex:

```{autodoc2-docstring} structfs.runner_actions.CaptureSolverOutput.report
```

````

`````

``````{py:class} ModelVariables(runner, **kwargs)
:canonical: structfs.runner_actions.ModelVariables
:noindex:

Bases: {py:obj}`structfs.runner.Action`

```{autodoc2-docstring} structfs.runner_actions.ModelVariables
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.runner_actions.ModelVariables.__init__
```

`````{py:class} Report(/, **data: typing.Any)
:canonical: structfs.runner_actions.ModelVariables.Report
:noindex:

Bases: {py:obj}`pydantic.BaseModel`

```{autodoc2-docstring} structfs.runner_actions.ModelVariables.Report
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.runner_actions.ModelVariables.Report.__init__
```

````{py:attribute} variables
:canonical: structfs.runner_actions.ModelVariables.Report.variables
:noindex:
:type: dict
:value: >
   'Field(...)'

```{autodoc2-docstring} structfs.runner_actions.ModelVariables.Report.variables
```

````

`````

````{py:method} after_run()
:canonical: structfs.runner_actions.ModelVariables.after_run
:noindex:

```{autodoc2-docstring} structfs.runner_actions.ModelVariables.after_run
```

````

````{py:method} _extract_vars(m)
:canonical: structfs.runner_actions.ModelVariables._extract_vars
:noindex:

```{autodoc2-docstring} structfs.runner_actions.ModelVariables._extract_vars
```

````

````{py:method} _is_var(c)
:canonical: structfs.runner_actions.ModelVariables._is_var
:noindex:
:staticmethod:

```{autodoc2-docstring} structfs.runner_actions.ModelVariables._is_var
```

````

````{py:method} _is_param(c)
:canonical: structfs.runner_actions.ModelVariables._is_param
:noindex:
:staticmethod:

```{autodoc2-docstring} structfs.runner_actions.ModelVariables._is_param
```

````

````{py:method} _add_block(tree: dict, name: str, block)
:canonical: structfs.runner_actions.ModelVariables._add_block
:noindex:
:staticmethod:

```{autodoc2-docstring} structfs.runner_actions.ModelVariables._add_block
```

````

````{py:method} report() -> Report
:canonical: structfs.runner_actions.ModelVariables.report
:noindex:

```{autodoc2-docstring} structfs.runner_actions.ModelVariables.report
```

````

``````

``````{py:class} MermaidDiagram
:canonical: structfs.runner_actions.MermaidDiagram
:noindex:

Bases: {py:obj}`structfs.runner.Action`

```{autodoc2-docstring} structfs.runner_actions.MermaidDiagram
```

`````{py:class} Report(/, **data: typing.Any)
:canonical: structfs.runner_actions.MermaidDiagram.Report
:noindex:

Bases: {py:obj}`pydantic.BaseModel`

```{autodoc2-docstring} structfs.runner_actions.MermaidDiagram.Report
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.runner_actions.MermaidDiagram.Report.__init__
```

````{py:attribute} diagram
:canonical: structfs.runner_actions.MermaidDiagram.Report.diagram
:noindex:
:type: list[str]
:value: >
   None

```{autodoc2-docstring} structfs.runner_actions.MermaidDiagram.Report.diagram
```

````

`````

````{py:method} after_run()
:canonical: structfs.runner_actions.MermaidDiagram.after_run
:noindex:

```{autodoc2-docstring} structfs.runner_actions.MermaidDiagram.after_run
```

````

````{py:method} report() -> Report | dict
:canonical: structfs.runner_actions.MermaidDiagram.report
:noindex:

```{autodoc2-docstring} structfs.runner_actions.MermaidDiagram.report
```

````

``````
