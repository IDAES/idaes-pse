from pathlib import Path
import time

from pylint.message import Message
from pylint.reporters import BaseReporter, text, JSONReporter


class MultiReporter(BaseReporter):
    """
    Combines multiple reporters so that the pylint message stream can be displayed in real time
    and simultaneously saved in a machine-readable format for further "offline" processing.
    """

    name = "multi"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.console = text.ColorizedTextReporter()
        self.json = JSONReporter(output=open("pylint.json", "w"))
        self._start = None

    def __iter__(self):
        yield self.console
        yield self.json

    def handle_message(self, msg):
        for rep in self:
            rep.handle_message(msg)

    def writeln(self, **kwargs):
        self.console.writeln(**kwargs)

    @property
    def elapsed(self):
        f = time.perf_counter
        if self._start is None:
            self._start = f()
            return 0.0
        return f() - self._start

    def on_set_current_module(self, module, filepath):
        self.writeln(string=f"{self.elapsed:07.3f} {filepath}")
        self.console._template = self.console.line_format

    def display_messages(self, layout):
        for rep in [self.json]:
            rep.display_messages(layout)

    def _display(self, layout):
        ...


class GHAReporter(text.TextReporter):
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


def register(linter):
    "This function needs to be defined for the plugin to be picked up by pylint"
    linter.register_reporter(MultiReporter)
