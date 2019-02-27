"""
This module implements pytest plugin for Sphinx doc tests.

In a nutshell, it uses the pytest `pytest_collect_file()` plugin hook
to recognize the Sphinx Makefile. Then it does a quick and dirty parse
of that Makefile to extract the command Sphinx is using to run the
doctests, which it recognizes by being the first command in the
Makefile target named by `SPHINX_DOCTEST_TARGET`. The parser is
able to handle simple Makefile variable expansion, though not currently
nested variables so don't do that.

The mechanics of the pytest plugin mechanism are such that the Makefile
is wrapped with a subclass of :class:`pytest.File`, :class:`SphinxMakefile`,
which implements the `collect` method to yield a subclass of :class:`pytest.Item`
called :class:`SphinxItem`, that in turn implements a few methods to run the
test and report the result. The bulk of the code in running the test is parsing
the output to look for errors, and thus decide whether all the doctests passed,
or not.

The drawback of this whole setup is of course some extra complexity.
The advantage is that (a) whatever the Makefile does is what this plugin
should do, for running the command, as long as the command is the first
(and only significant) thing that occurs in the target, and (b) if there ends up
being more than one Makefile, it should all continue to work.
"""
import os
import re
import subprocess

#
import pytest

__author__ = 'Dan Gunter'

# Modify this to match the target in the
# Sphinx Makefile for running the doctests
SPHINX_DOCTEST_TARGET = "doctest"


def pytest_collect_file(parent, path):
    if path.ext == "" and path.basename == "Makefile":
        if file_contains("sphinx", str(path)):
            return SphinxMakefile(path, parent)


def file_contains(s, path):
    result = False
    with open(path) as f:
        for line in f:
            if s in line:
                result = True
                break
    return result


class SphinxMakefile(pytest.File):

    # simple way to find variables in a Makefile
    makefile_var_defn = re.compile(r"\s*([A-Za-z_]+)\s*=(.*)")

    def collect(self):
        cmd = self._get_doctest_command()
        wd = os.path.dirname(self.fspath)
        if cmd is not None:
            yield SphinxItem('doctests', self, wd, cmd)

    def _get_doctest_command(self):
        result = None
        mkvars = {}  # variables defined in Makefile
        with self.fspath.open() as makefile:
            cmd = ""
            state = "pre"
            for line in makefile:
                s = line.strip()
                # look for 'doctest' target
                if state == "pre":
                    if s.startswith(SPHINX_DOCTEST_TARGET):
                        state = "next"
                    else:
                        # look for a NAME = value definition, save it in 'mkvars'
                        m = self.makefile_var_defn.match(s)
                        if m is not None:
                            value = m.group(2).strip()
                            expanded_value = self._expand(value, mkvars)
                            mkvars[m.group(1)] = expanded_value
                # primed for this to be the doctest command
                elif state == "next":
                    # ignore print statements
                    if "echo " in s or "printf " in s:
                        continue
                    # use 'mkvars' to expand out command
                    if s.endswith("\\"):
                        cmd += s[:-1] + " "
                    else:
                        cmd += s
                        result = self._expand(cmd, mkvars)
                        break
            if result is None and state == "next":
                raise ValueError(f"EOF while parsing doctest command '{cmd}'")
        return result

    @staticmethod
    def _expand(s, mkvars):
        """Expand out the variable refs by making them look like format vars, then
        simply using string formatting.
        """
        fstr, n = re.subn(r"\$\(([A-Za-z_]+)\)", "{\\1}", s)
        return fstr.format(**mkvars)


class SphinxItem(pytest.Item):
    def __init__(self, name, parent, wd, cmd):
        """New Sphinx doctest item.

        Args:
            name (str): Item name
            parent (pytest.File): Parent item
            wd (str): working directory
            cmd (str): Command to run
        """
        super().__init__(name, parent)
        # working directory and command
        self.wd, self.cmd = wd, cmd
        self.successes, self.failures = None, None
        self.total_successes, self.total_failures = -1, -1
        self.failure_list = []

    def runtest(self):
        """Run the Sphinx doctest.
        """
        old_d = os.getcwd()
        print(f"doctest command: {self.cmd}")
        try:
            os.chdir(self.wd)
            args = self.cmd.split()
            try:
                # print(f"Running [{self.cmd}] from dir {self.wd}")
                proc = subprocess.Popen(
                    args, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                self.successes, self.failures = self._parse_output(proc.stdout)
                rc = proc.wait()
                if rc != 0:
                    raise RuntimeError(f"non-zero exit code: {rc}")
                self.total_successes = sum(self.successes.values())
                self.total_failures = sum((x[0] for x in self.failures.values()))
            except Exception as exc:
                print(f"Unexpected failure in Sphinx command: {exc}")
                raise SphinxCommandFailed(self.cmd, str(exc))
            # report a test failure
            if self.total_failures > 0:
                self.failure_list = [
                    self.failures[k][1]
                    for k in self.failures
                    if self.failures[k][0] > 0
                ]
                raise SphinxHadErrors()
        finally:
            os.chdir(old_d)

    def repr_failure(self, excinfo):
        """This is called when self.runtest() raises an exception.
        """
        if not isinstance(excinfo.value, SphinxHadErrors):
            return None
        header = f"doctest execution failed: {self.total_failures} failures"
        body = "\n".join(self.failure_list)
        return header + "\n" + body

    def reportinfo(self):
        summary = f"Sphinx doctests in {self.wd}"
        return f"doctests in {self.wd}", 0, summary

    @staticmethod
    def _parse_output(input_stream):
        """Parse output of Sphinx doctest.
        """

        def _msgblock(messages):
            lines, indent, marker = [], None, "  - "
            for m in messages:
                first_line = marker
                indent = " " * len(first_line)
                for s in m.split("\n"):
                    if first_line:
                        first_line += s
                        lines.append(first_line)
                        first_line = None
                    else:
                        lines.append(indent + s)
            return "\n".join(lines)

        state = "ok"
        failures, successes = {}, {}
        cur_doc, n_passed, n_failed = "", -1, -1
        failed_msg, failed_msgs = None, []
        for line_bytes in input_stream:
            line = line_bytes.decode('utf-8').rstrip()
            try:
                if state == "ok":
                    if line.startswith('****'):
                        # start of a failure message
                        state = "fail"
                        failed_msg = ""
                    elif line.startswith("Document:"):
                        # start of a new doc
                        _, cur_doc = line.split(None, 1)
                        n_passed, n_failed, failed_msgs = 0, 0, []
                    elif re.match(r"\d+ passed and \d+ failed", line):
                        # summary of pass/fail for a doc
                        m = re.match(r"(\d+) passed and (\d+) failed", line)
                        n_passed, n_failed = int(m.group(1)), int(m.group(2))
                    elif ("Test passed" in line) or ("Test Failed" in line):
                        # final line of report for a doc
                        successes[cur_doc] = n_passed
                        if n_failed > 0:
                            failures[cur_doc] = (n_failed, _msgblock(failed_msgs))
                        else:
                            failures[cur_doc] = (0, [])
                        cur_doc = None
                elif state == "fail":
                    if line.startswith("****"):
                        state = "ok"
                        failed_msgs.append(failed_msg)
                    else:
                        if failed_msg:
                            failed_msg += "\n"
                        failed_msg += line
            except Exception as exc:
                print(f"Parse error on line '{line}' :: {exc}")
                raise
        # print(f"DONE: successes={successes}, failures={failures}")
        return successes, failures


class SphinxCommandFailed(Exception):
    def __init_(self, cmd, msg):
        super().__init__(f"Sphinx doctest command ({cmd}) failed: {msg}")


class SphinxHadErrors(Exception):
    pass
