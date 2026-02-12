# {py:mod}`structfs.fsrunner`

```{py:module} structfs.fsrunner
:noindex:
```

```{autodoc2-docstring} structfs.fsrunner
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Context <structfs.fsrunner.Context>`
  - ```{autodoc2-docstring} structfs.fsrunner.Context
    :summary:
    ```
* - {py:obj}`BaseFlowsheetRunner <structfs.fsrunner.BaseFlowsheetRunner>`
  - ```{autodoc2-docstring} structfs.fsrunner.BaseFlowsheetRunner
    :summary:
    ```
* - {py:obj}`FlowsheetRunner <structfs.fsrunner.FlowsheetRunner>`
  - ```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner
    :summary:
    ```
````

### API

`````{py:class} Context()
:canonical: structfs.fsrunner.Context
:noindex:

Bases: {py:obj}`dict`

```{autodoc2-docstring} structfs.fsrunner.Context
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.fsrunner.Context.__init__
```

````{py:property} model
:canonical: structfs.fsrunner.Context.model
:noindex:

```{autodoc2-docstring} structfs.fsrunner.Context.model
```

````

````{py:property} solver
:canonical: structfs.fsrunner.Context.solver
:noindex:

```{autodoc2-docstring} structfs.fsrunner.Context.solver
```

````

`````

`````{py:class} BaseFlowsheetRunner(solver=None, tee=True)
:canonical: structfs.fsrunner.BaseFlowsheetRunner
:noindex:

Bases: {py:obj}`structfs.runner.Runner`

```{autodoc2-docstring} structfs.fsrunner.BaseFlowsheetRunner
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.fsrunner.BaseFlowsheetRunner.__init__
```

````{py:attribute} STEPS
:canonical: structfs.fsrunner.BaseFlowsheetRunner.STEPS
:noindex:
:value: >
   ('build', 'set_operating_conditions', 'set_scaling', 'initialize', 'set_solver', 'solve_initial', 'a...

```{autodoc2-docstring} structfs.fsrunner.BaseFlowsheetRunner.STEPS
```

````

````{py:method} run_steps(first: str = Runner.STEP_ANY, last: str = Runner.STEP_ANY, before=None, after=None, **kwargs)
:canonical: structfs.fsrunner.BaseFlowsheetRunner.run_steps
:noindex:

```{autodoc2-docstring} structfs.fsrunner.BaseFlowsheetRunner.run_steps
```

````

````{py:method} reset()
:canonical: structfs.fsrunner.BaseFlowsheetRunner.reset
:noindex:

````

````{py:method} _create_model()
:canonical: structfs.fsrunner.BaseFlowsheetRunner._create_model
:noindex:

```{autodoc2-docstring} structfs.fsrunner.BaseFlowsheetRunner._create_model
```

````

````{py:property} model
:canonical: structfs.fsrunner.BaseFlowsheetRunner.model
:noindex:

```{autodoc2-docstring} structfs.fsrunner.BaseFlowsheetRunner.model
```

````

````{py:property} results
:canonical: structfs.fsrunner.BaseFlowsheetRunner.results
:noindex:

```{autodoc2-docstring} structfs.fsrunner.BaseFlowsheetRunner.results
```

````

````{py:method} annotate_var(variable: object, key: str = None, title: str = None, desc: str = None, units: str = None, rounding: int = 3, is_input: bool = True, is_output: bool = True, input_category: str = 'main', output_category: str = 'main') -> object
:canonical: structfs.fsrunner.BaseFlowsheetRunner.annotate_var
:noindex:

```{autodoc2-docstring} structfs.fsrunner.BaseFlowsheetRunner.annotate_var
```

````

````{py:property} annotated_vars
:canonical: structfs.fsrunner.BaseFlowsheetRunner.annotated_vars
:noindex:
:type: dict[str]

```{autodoc2-docstring} structfs.fsrunner.BaseFlowsheetRunner.annotated_vars
```

````

`````

``````{py:class} FlowsheetRunner(**kwargs)
:canonical: structfs.fsrunner.FlowsheetRunner
:noindex:

Bases: {py:obj}`structfs.fsrunner.BaseFlowsheetRunner`

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.__init__
```

`````{py:class} DegreesOfFreedom(runner)
:canonical: structfs.fsrunner.FlowsheetRunner.DegreesOfFreedom
:noindex:

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.DegreesOfFreedom
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.DegreesOfFreedom.__init__
```

````{py:method} model()
:canonical: structfs.fsrunner.FlowsheetRunner.DegreesOfFreedom.model
:noindex:

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.DegreesOfFreedom.model
```

````

````{py:method} __getattr__(name)
:canonical: structfs.fsrunner.FlowsheetRunner.DegreesOfFreedom.__getattr__
:noindex:

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.DegreesOfFreedom.__getattr__
```

````

````{py:method} __str__()
:canonical: structfs.fsrunner.FlowsheetRunner.DegreesOfFreedom.__str__
:noindex:

````

````{py:method} _ipython_display_()
:canonical: structfs.fsrunner.FlowsheetRunner.DegreesOfFreedom._ipython_display_
:noindex:

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.DegreesOfFreedom._ipython_display_
```

````

`````

`````{py:class} Timings(runner)
:canonical: structfs.fsrunner.FlowsheetRunner.Timings
:noindex:

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.Timings
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.Timings.__init__
```

````{py:property} values
:canonical: structfs.fsrunner.FlowsheetRunner.Timings.values
:noindex:
:type: list[dict]

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.Timings.values
```

````

````{py:property} history
:canonical: structfs.fsrunner.FlowsheetRunner.Timings.history
:noindex:
:type: str

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.Timings.history
```

````

````{py:method} __str__()
:canonical: structfs.fsrunner.FlowsheetRunner.Timings.__str__
:noindex:

````

````{py:method} _ipython_display_()
:canonical: structfs.fsrunner.FlowsheetRunner.Timings._ipython_display_
:noindex:

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.Timings._ipython_display_
```

````

`````

````{py:method} build()
:canonical: structfs.fsrunner.FlowsheetRunner.build
:noindex:

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.build
```

````

````{py:method} solve_initial()
:canonical: structfs.fsrunner.FlowsheetRunner.solve_initial
:noindex:

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.solve_initial
```

````

````{py:method} show_diagram()
:canonical: structfs.fsrunner.FlowsheetRunner.show_diagram
:noindex:

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.show_diagram
```

````

``````
