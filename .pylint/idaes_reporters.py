from pathlib import Path
import time

from pylint.reporters import BaseReporter, text, JSONReporter


class MultiReporter(BaseReporter):
    """
    Combines multiple reporters so that the pylint message stream can be displayed in real time
    and simultaneously saved in a machine-readable format for further "offline" processing.
    """
    name = 'multi'
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.console = text.ColorizedTextReporter()
        self.json = JSONReporter(output=open('pylint.json', 'w'))
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
        f = time.clock
        if self._start is None:
            self._start = f()
            return 0.
        return f() - self._start

    def on_set_current_module(self, module, filepath):
        self.writeln(string=f'{self.elapsed:07.3f} {filepath}')
        self.console._template = self.console.line_format

    def display_messages(self, layout):
        for rep in [self.json]:
            rep.display_messages(layout)

    def _display(self, layout):
        ...


def register(linter):
    "This function needs to be defined for the plugin to be picked up by pylint"
    linter.register_reporter(MultiReporter)
