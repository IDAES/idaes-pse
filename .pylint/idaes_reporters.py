from pathlib import Path
import time

from pylint.message import Message
from pylint.reporters import BaseReporter, text, JSONReporter


class DisplayProgress(BaseReporter):
    """
    Display analyzed modules with total elapsed time since beginning of run.
    """

    name = "progress"
    time_format = "07.3f"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._start = None

    @property
    def elapsed(self):
        f = time.perf_counter
        if self._start is None:
            self._start = f()
            return 0.0
        return f() - self._start

    def on_set_current_module(self, module, filepath):
        if filepath is None: return
        self.writeln(string=f"{self.elapsed:{self.time_format}} {filepath}")

    def _display(self, layout):
        pass


class GHACheckAnnotations(text.TextReporter):
    name = "gha"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._map_pylint_to_gha = {
            "critical": "error",
            "warning": "warning",
            "info": "notice",
        }

    def write_message(self, m: Message) -> None:
        gha_type = self._map_pylint_to_gha.get(m.category, "notice")
        gha_title = f"{m.msg_id} ({m.symbol})"
        line = f"::{gha_type} file={m.path},line={m.line},endLine={m.end_line},title={gha_title}::{m.msg}"
        self.writeln(line)

    def display_reports(self, *args, **kwargs):
        # avoid displaying summary sections for this reporter
        pass


def register(linter):
    "This function needs to be defined for the plugin to be picked up by pylint"
    linter.register_reporter(DisplayProgress)
    linter.register_reporter(GHACheckAnnotations)
