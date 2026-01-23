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
Utility code for documentation-related tests.
"""

__author__ = "Dan Gunter (LBNL)"

import re


class Docstring:
    """Pick the marked code sections out of a docstring.

    This is useful for (at least) writing tests that use the code
    snippets in a markdown (e.g. myst and with autodoc2) docstring,
    without duplicating those code snippets across the original module
    and test file.
    """

    LABEL_OPTION = ":name:"

    def __init__(self, text: str, style: str = "markdown"):
        self._code = {}
        if style == "markdown":
            self._init_markdown(text)
        # TODO: also handle RST
        else:
            raise ValueError(
                f"Unknown docstring style: {style}. " f"Must be one of: markdown"
            )

    def code(self, section: str, func_prefix: str = None) -> str:
        """Get code section.

        Args:
            section: Label for the section
            func_prefix: Add this prefix to all top-level
               function definitions. If not provided,
               the prefix will be the section name and an
               underscore. Give an empty string to
               have no prefixes.

        Returns:
            Text of the code section, for exec()
        """
        if func_prefix is None:
            func_prefix = f"{section}_"
        section_lines = self._prefix_def(self._code[section], func_prefix)
        return "\n".join(section_lines)

    def _prefix_def(self, lines: list[str], prefix: str):
        # Note that this will only match function definitions at
        # module scope (not methods in a class, which are already,
        # presumably, in a unique-enough namespace).
        expr = re.compile(r"def\s+(?P<name>\w+)")
        result = []
        for line in lines:
            m = expr.match(line)
            if m:
                # prepend prefix to function name
                pos = m.start(1)
                line = line[:pos] + prefix + line[pos:]
            result.append(line)
        return result

    def _init_markdown(self, text: str):
        lines = text.split("\n")
        state = 0
        for line in lines:
            ls_line = line.lstrip()
            if state == 0 and ls_line.startswith("```{code}"):
                indent = len(line) - len(ls_line)  # spaces to left
                state = 1
                section_name = None
                section_lines = []
            elif state > 0:
                if ls_line.startswith("```"):
                    if section_name is None:
                        pass  # silently ignore
                    elif len(section_lines) == 0:
                        pass  # silently ignore
                    else:
                        self._code[section_name] = section_lines
                    state = 0
                elif state == 1:
                    if ls_line.startswith(self.LABEL_OPTION):
                        section_name = ls_line[7:].strip()
                    elif ls_line == "" or ls_line.startswith(":"):
                        pass
                    else:
                        state = 2  # got past metadata
                        section_lines.append(line[indent:].rstrip())
                else:
                    section_lines.append(line[indent:].rstrip())
