##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
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
import pathlib
import subprocess

#
import pytest

__author__ = 'Dan Gunter'

# Modify this to match the target in the
# Sphinx Makefile for running the doctests
SPHINX_DOCTEST_TARGET = "doctest"

# Modify this for non-standard sphinx-build commands.
# Knowing the Sphinx build command is needed to replicate it
SPHINX_BUILD = "sphinx-build"

g_sphinx_warn_file = None


def pytest_collect_file(parent, path):
    global g_sphinx_warn_file
    if path.ext == "" and path.basename == "Makefile":
        if file_contains("sphinx", str(path)):
            makefile = SphinxMakefile(path, parent)
            g_sphinx_warn_file = makefile.warnings_file
            return makefile
    elif g_sphinx_warn_file is not None and path.basename == g_sphinx_warn_file:
        return SphinxWarnings(path, parent)


def file_contains(s, path):
    result = False
    with open(path) as f:
        for line in f:
            if s in line:
                result = True
                break
    return result


class SphinxWarnings(pytest.File):
    def collect(self):
        path = pathlib.Path(self.fspath)  # convert to std
        yield SphinxWarningsItem('Sphinx warnings', self, path)


class SphinxWarningsItem(pytest.Item):
    def __init__(self, name, parent, path: pathlib.Path):
        super().__init__(name, parent)
        self._path = path
        self._warnings = None

    def runtest(self):
        if self._path.stat().st_size > 0:
            with self._path.open() as f:
                self._warnings = f.read()
            raise SphinxHadErrors

    def repr_failure(self, excinfo):
        """This is called when self.runtest() raises an exception.
        """
        return f"Build had warnings and/or errors:\n{self._warnings}"

    def reportinfo(self):
        summary = f"Sphinx doc warnings in {self._path}"
        return f"doc build warnings in {self._path}", 0, summary


class SphinxMakefile(pytest.File):

    # simple way to find variables in a Makefile
    makefile_var_defn = re.compile(r"\s*([A-Za-z_]+)\s*=(.*)")

    def collect(self):
        cmd = self._get_doctest_command()
        wd = os.path.dirname(self.fspath)
        if cmd is not None:
            yield SphinxDoctestItem('Sphinx doctests', self, wd, cmd)

    @property
    def warnings_file(self):
        """Get warnings and errors output file, if any, from the Sphinx Makefile.
        """
        mkvars = {}  # variables defined in Makefile
        with self.fspath.open() as makefile:
            for line in makefile:
                s = line.strip()
                # look for definition of Sphinx options
                if s.startswith("SPHINXOPTS"):
                    opts = s[s.find('=') + 1 :].strip()
                    # return value of '-w' option, if found
                    m = re.match("-w\s*['\"]?([^\"' \t]+)", opts)
                    return m.group(1) if m else None

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
                    # use 'mkvars' to expand out command
                    if s.endswith("\\"):
                        cmd += s[:-1] + " "
                    else:
                        cmd += s
                        expanded = self._expand(cmd, mkvars)
                        if SPHINX_BUILD not in expanded:
                            cmd = ""
                        else:
                            result = expanded
                            state = "done"
                elif state == "done":
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


class SphinxDoctestItem(pytest.Item):

    __executed = False

    def __init__(self, name, parent, wd, cmd):
        """New Sphinx doctest item.

        Args:
            name (str): Item name
            parent (pytest.File): Parent item
            wd (str): working directory
            cmd (str): Command to run
        """
        super().__init__(name, parent)
        self._sess = parent.session
        # working directory and command
        self.wd, self.cmd = wd, cmd
        self.successes, self.failures = None, None
        self.total_successes, self.total_failures = -1, -1
        self.failure_list = []

    def runtest(self):
        """Run the Sphinx doctest.
        """
        self._insert_items_at = self._sess.items.index(self) + 1
        #  print(f"\n@@ SphinxDoctestItem.run(): session items: {self._sess.items}.\n"
        #      f"Insert at: {self._insert_items_at}\n")
        if self.__executed:
            return
        self.__executed = True
        old_d = os.getcwd()
        print(f"doctest command: {self.cmd}")
        try:
            os.chdir(self.wd)
            args = self.cmd.split()
            try:
                proc = subprocess.run(
                    args,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    timeout=300,
                    encoding='utf-8',
                )
                self._parse_output(proc.stdout)
            except Exception as exc:
                print(f"Unexpected failure in Sphinx doctest command: {exc}")
                raise SphinxCommandFailed(self.cmd, str(exc))
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

    def _parse_output(self, output):
        """Parse output of Sphinx doctest.
        """

        def _msgblock(messages):
            lines, indent, marker = [], None, " ** "
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

        state, failed_msg = "ok", ""
        cur_doc, n_passed, n_failed = "", -1, -1
        for line in output.split('\n'):
            line = line.rstrip()
            # print(f"@@ doctest: {line}")
            try:
                if state == "ok":
                    if line.startswith('****'):
                        # start of a failure message
                        state = "fail"
                        failed_msg = ""
                    elif line.startswith("Document:"):
                        # start of a new doc
                        _, cur_doc = line.split(None, 1)
                        n_passed, n_failed = 0, 0
                    elif re.match(r"\d+ passed and \d+ failed", line):
                        # summary of pass/fail for a doc
                        m = re.match(r"(\d+) passed and (\d+) failed", line)
                        n_passed, n_failed = int(m.group(1)), int(m.group(2))
                    elif ("Test passed" in line) or ("Test Failed" in line):
                        # final line of report for a doc
                        if n_passed > 0:
                            self.tests_ok(cur_doc, n_passed)
                        cur_doc = None
                elif state == "fail":
                    if line.startswith("****"):
                        # end of failure record
                        self.test_failed(cur_doc, failed_msg)
                        state = "ok"
                        # do not reset 'cur_doc', there may be more..
                    else:
                        # middle of failure record
                        if failed_msg:
                            failed_msg += "\n"
                        failed_msg += line
            except Exception as exc:
                print(f"Parse error on line '{line}' :: {exc}")
                raise
        # print(f"@@ done doctest output")

    def tests_ok(self, context, num):
        # print(f"Tests in {context}: {num} OK")
        test_name = "Sphinx doctest in: {context}"
        if num > 1:
            test_name += " ({i})"
        for i in range(num, 0, -1):
            success_item = SphinxDoctestSuccess(
                test_name.format(**locals()), self.parent
            )
            self._sess.items.insert(self._insert_items_at, success_item)
        self._sess.testscollected += num

    def test_failed(self, context, failure):
        test_name = "Sphinx doctest in: {context}.rst"
        item = SphinxDoctestFailure(test_name.format(**locals()), self.parent, failure)
        self._sess.items.insert(self._insert_items_at, item)
        self._sess.testscollected += 1


class SphinxDoctestSuccess(pytest.Item):
    def runtest(self):
        return


class SphinxDoctestFailure(pytest.Item):
    def __init__(self, name, parent, details):
        super().__init__(name, parent)
        self.details = details

    def runtest(self):
        raise SphinxHadErrors()

    def repr_failure(self, excinfo):
        return "doctest failure:\n" + self.details

    def reportinfo(self):
        return f"doctest in {self.name}", 0, f"FAILED {self.name}"


class SphinxCommandFailed(Exception):
    def __init_(self, cmd, msg):
        super().__init__(f"Sphinx doctest command ({cmd}) failed: {msg}")


class SphinxHadErrors(Exception):
    pass
