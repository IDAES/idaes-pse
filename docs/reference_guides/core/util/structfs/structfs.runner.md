# {py:mod}`structfs.runner`

```{py:module} structfs.runner
:noindex:
```

```{autodoc2-docstring} structfs.runner
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Step <structfs.runner.Step>`
  - ```{autodoc2-docstring} structfs.runner.Step
    :summary:
    ```
* - {py:obj}`Runner <structfs.runner.Runner>`
  - ```{autodoc2-docstring} structfs.runner.Runner
    :summary:
    ```
* - {py:obj}`Action <structfs.runner.Action>`
  - ```{autodoc2-docstring} structfs.runner.Action
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`__author__ <structfs.runner.__author__>`
  - ```{autodoc2-docstring} structfs.runner.__author__
    :summary:
    ```
* - {py:obj}`_log <structfs.runner._log>`
  - ```{autodoc2-docstring} structfs.runner._log
    :summary:
    ```
* - {py:obj}`ActionType <structfs.runner.ActionType>`
  - ```{autodoc2-docstring} structfs.runner.ActionType
    :summary:
    ```
````

### API

````{py:data} __author__
:canonical: structfs.runner.__author__
:noindex:
:value: >
   'Dan Gunter (LBNL)'

```{autodoc2-docstring} structfs.runner.__author__
```

````

````{py:data} _log
:canonical: structfs.runner._log
:noindex:
:value: >
   'Logger(...)'

```{autodoc2-docstring} structfs.runner._log
```

````

`````{py:class} Step(name: str, func: typing.Callable)
:canonical: structfs.runner.Step
:noindex:

```{autodoc2-docstring} structfs.runner.Step
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.runner.Step.__init__
```

````{py:attribute} SEP
:canonical: structfs.runner.Step.SEP
:noindex:
:value: >
   '::'

```{autodoc2-docstring} structfs.runner.Step.SEP
```

````

````{py:method} add_substep(name: str, func: typing.Callable)
:canonical: structfs.runner.Step.add_substep
:noindex:

```{autodoc2-docstring} structfs.runner.Step.add_substep
```

````

`````

````{py:data} ActionType
:canonical: structfs.runner.ActionType
:noindex:
:value: >
   'TypeVar(...)'

```{autodoc2-docstring} structfs.runner.ActionType
```

````

`````{py:class} Runner(steps: typing.Sequence[str])
:canonical: structfs.runner.Runner
:noindex:

```{autodoc2-docstring} structfs.runner.Runner
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.runner.Runner.__init__
```

````{py:attribute} STEP_ANY
:canonical: structfs.runner.Runner.STEP_ANY
:noindex:
:value: >
   '-'

```{autodoc2-docstring} structfs.runner.Runner.STEP_ANY
```

````

````{py:method} __getitem__(key)
:canonical: structfs.runner.Runner.__getitem__
:noindex:

```{autodoc2-docstring} structfs.runner.Runner.__getitem__
```

````

````{py:method} add_step(name: str, func: typing.Callable)
:canonical: structfs.runner.Runner.add_step
:noindex:

```{autodoc2-docstring} structfs.runner.Runner.add_step
```

````

````{py:method} add_substep(base_name, name, func)
:canonical: structfs.runner.Runner.add_substep
:noindex:

```{autodoc2-docstring} structfs.runner.Runner.add_substep
```

````

````{py:method} run_step(name)
:canonical: structfs.runner.Runner.run_step
:noindex:

```{autodoc2-docstring} structfs.runner.Runner.run_step
```

````

````{py:method} run_steps(first: str = '', last: str = '', after: str = '', before: str = '')
:canonical: structfs.runner.Runner.run_steps
:noindex:

```{autodoc2-docstring} structfs.runner.Runner.run_steps
```

````

````{py:method} _run_steps(first: str, last: str, endpoints: tuple[bool, bool])
:canonical: structfs.runner.Runner._run_steps
:noindex:

```{autodoc2-docstring} structfs.runner.Runner._run_steps
```

````

````{py:method} reset()
:canonical: structfs.runner.Runner.reset
:noindex:

```{autodoc2-docstring} structfs.runner.Runner.reset
```

````

````{py:method} list_steps(all_steps=False) -> list[str]
:canonical: structfs.runner.Runner.list_steps
:noindex:

```{autodoc2-docstring} structfs.runner.Runner.list_steps
```

````

````{py:method} add_action(name: str, action_class: type, *args, **kwargs) -> object
:canonical: structfs.runner.Runner.add_action
:noindex:

```{autodoc2-docstring} structfs.runner.Runner.add_action
```

````

````{py:method} get_action(name: str) -> structfs.runner.ActionType
:canonical: structfs.runner.Runner.get_action
:noindex:

```{autodoc2-docstring} structfs.runner.Runner.get_action
```

````

````{py:method} remove_action(name: str)
:canonical: structfs.runner.Runner.remove_action
:noindex:

```{autodoc2-docstring} structfs.runner.Runner.remove_action
```

````

````{py:method} _find_step(reverse=False)
:canonical: structfs.runner.Runner._find_step
:noindex:

```{autodoc2-docstring} structfs.runner.Runner._find_step
```

````

````{py:method} normalize_name(s: typing.Optional[str]) -> str
:canonical: structfs.runner.Runner.normalize_name
:noindex:
:classmethod:

```{autodoc2-docstring} structfs.runner.Runner.normalize_name
```

````

````{py:method} _step_begin(name: str)
:canonical: structfs.runner.Runner._step_begin
:noindex:

```{autodoc2-docstring} structfs.runner.Runner._step_begin
```

````

````{py:method} _substep_begin(base: str, name: str)
:canonical: structfs.runner.Runner._substep_begin
:noindex:

```{autodoc2-docstring} structfs.runner.Runner._substep_begin
```

````

````{py:method} _step_end(name: str)
:canonical: structfs.runner.Runner._step_end
:noindex:

```{autodoc2-docstring} structfs.runner.Runner._step_end
```

````

````{py:method} _substep_end(base: str, name: str)
:canonical: structfs.runner.Runner._substep_end
:noindex:

```{autodoc2-docstring} structfs.runner.Runner._substep_end
```

````

````{py:method} step(name: str)
:canonical: structfs.runner.Runner.step
:noindex:

```{autodoc2-docstring} structfs.runner.Runner.step
```

````

````{py:method} substep(base: str, name: str)
:canonical: structfs.runner.Runner.substep
:noindex:

```{autodoc2-docstring} structfs.runner.Runner.substep
```

````

````{py:method} report() -> dict[str, dict]
:canonical: structfs.runner.Runner.report
:noindex:

```{autodoc2-docstring} structfs.runner.Runner.report
```

````

`````

`````{py:class} Action(runner: structfs.runner.Runner, log: typing.Optional[logging.Logger] = None)
:canonical: structfs.runner.Action
:noindex:

Bases: {py:obj}`abc.ABC`

```{autodoc2-docstring} structfs.runner.Action
```

```{rubric} Initialization
```

```{autodoc2-docstring} structfs.runner.Action.__init__
```

````{py:method} before_step(step_name: str)
:canonical: structfs.runner.Action.before_step
:noindex:

```{autodoc2-docstring} structfs.runner.Action.before_step
```

````

````{py:method} before_substep(step_name: str, substep_name: str)
:canonical: structfs.runner.Action.before_substep
:noindex:

```{autodoc2-docstring} structfs.runner.Action.before_substep
```

````

````{py:method} after_step(step_name: str)
:canonical: structfs.runner.Action.after_step
:noindex:

```{autodoc2-docstring} structfs.runner.Action.after_step
```

````

````{py:method} after_substep(step_name: str, substep_name: str)
:canonical: structfs.runner.Action.after_substep
:noindex:

```{autodoc2-docstring} structfs.runner.Action.after_substep
```

````

````{py:method} before_run()
:canonical: structfs.runner.Action.before_run
:noindex:

```{autodoc2-docstring} structfs.runner.Action.before_run
```

````

````{py:method} after_run()
:canonical: structfs.runner.Action.after_run
:noindex:

```{autodoc2-docstring} structfs.runner.Action.after_run
```

````

````{py:method} report() -> pydantic.BaseModel | dict
:canonical: structfs.runner.Action.report
:noindex:
:abstractmethod:

```{autodoc2-docstring} structfs.runner.Action.report
```

````

`````
