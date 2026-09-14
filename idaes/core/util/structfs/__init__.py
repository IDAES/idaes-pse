
__generated_with = "0.23.16"

# %%
from unittest.mock import Mock

# %%
class MockFlowsheetRunner:
    def __init__(self, *args, **kwargs):
        self._args, self._kwargs = args, kwargs
        self.ctx = Mock()
        
    def __getattr__(self, name):
        return Mock()
        
    def step(self, *args, **kwargs):
        def decorator(func):
            def wrapper(ctx=self.ctx):
                func(ctx)
            return wrapper
        return decorator

# %%
try:
    from idaes_fi.structfs import FlowsheetRunner
except ImportError:
    FlowsheetRunner = MockFlowsheetRunner

FS = FlowsheetRunner(foo=1)
@FS.step("name")
def foo(ctx):
    ctx.model = "hello"

@FS.step("next")
def bar(ctx):
    m = ctx.model
    print(f"model = {m}")

# %%
FS.run_steps()

# %%
