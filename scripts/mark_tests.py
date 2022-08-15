#!/usr/bin/env python
"""
Fill in missing pytest mark decorators to test functions.

Adds 'import pytest' to all test scripts that lack it, then
adds '@pytest.mark.unit' to all tests that are not already marked
with either 'component', 'integration', or 'unit'.

The user must provide a list of test files to process,
in a text file, listed by by relative or absolute path.
"""
TEST_FILE_LIST = "testfiles.txt"

import logging
import re

_log = logging.getLogger("mark_tests")
_h = logging.StreamHandler()
_h.setFormatter(logging.Formatter(fmt="%(asctime)s [%(levelname)s] %(message)s"))
_log.addHandler(_h)

_log.setLevel(logging.DEBUG)


def process_testfile(filename):
    _log.debug(f"Processing {filename}")
    with open(filename, "r") as testfile:
        lines = testfile.readlines()

    insert_indices = []
    pytest_imported = False
    insert_pytest_import_at = -1
    in_docstring = False  # hacky hacky hacky
    in_multiline_import = False  # more hacky hacky ugh

    for line_ind in range(len(lines)):
        # skip any comments
        if re.match(r" *#", lines[line_ind]):
            continue

        if re.match(r' *"""', lines[line_ind]):
            if not in_docstring:
                in_docstring = True
                continue
            else:
                in_docstring = False

        if in_multiline_import:
            # TODO
            _log.debug("")

        if not pytest_imported:  # we're at/near the top of the file still
            # look for import statments; when we hit something that isn't
            # an import, whitepsace, comment, or docstring
            if re.match(r"import pytest", lines[line_ind]):
                pytest_imported = True
            elif re.match(r"(import|from).*,\n$", lines[line_ind]):
                in_multiline_import = True
            elif re.match(r"import|from", lines[line_ind]):
                insert_pytest_import_at = line_ind + 1
            elif re.match(r" *($|#)", lines[line_ind]) and insert_pytest_import_at < 0:
                insert_pytest_import_at = line_ind + 1

        def_ind = lines[line_ind].find("def test_")
        if def_ind > -1:
            # We're at a test function definition.
            #
            # Search backwards for one of (ignoring indent spaces!):
            #   0 - head of file
            #   1 - @pytest.mark.integration
            #   2 - @pytest.mark.component
            #   3 - @pytest.mark.unit
            #   4 - an empty line, docstring, or class decl
            #   5 - a comment or another unspecified decorator
            #   6 - something else
            #
            # 0: add @pytest.mark.unit, then continue, but warn about the unusual situation
            # 1,2,3: we're done, continue
            # 4: add @pytest.mark.unit, then continue
            # 5: keep looking backwards
            # 6: error; abort! There should always be whitespace before a def
            backtrack_index = line_ind - 1
            while backtrack_index > 0:
                if re.match(r' *($|"""|class)', lines[backtrack_index]):  # case 4
                    insert_indices.append((line_ind, def_ind))
                    _log.debug(
                        f"\tAdding `@pytest.mark.unit` to {lines[line_ind][def_ind+4:-2]}"
                    )
                    break
                elif re.match(
                    r" *@pytest.mark.(integration|component|unit)",
                    lines[backtrack_index],
                ):  # case 1,2,3
                    break
                elif re.match(r" *(#|@)", lines[backtrack_index]):  # case 5
                    backtrack_index -= 1
                else:  # case 5: ABORT!
                    # HACKY HACKY HACKY EWW
                    if re.match(r" *reason", lines[backtrack_index]) and re.match(
                        r" *@pytest.mark.skipif", lines[backtrack_index - 1]
                    ):
                        backtrack_index -= 2
                        continue
                    _log.warning(
                        f"Error at line {line_ind} (function {lines[line_ind][def_ind+4:-2]})"
                        f" in file {filename}! Skipping file."
                    )
                    return -1
            else:  # case 0, head of file: unusual situation!
                insert_indices.append((0, def_ind))
                _log.debug(
                    f"\tAdding `@pytest.mark.unit` to {lines[line_ind][def_ind+4:-2]}"
                )
                _log.debug("\tNOTE: this is the first line of the file!")

    if len(insert_indices) == 0 and pytest_imported:
        return 0

    while len(insert_indices) > 0:
        mark_index, indent_level = insert_indices.pop()
        lines.insert(mark_index, " " * indent_level + "@pytest.mark.unit\n")

    if not pytest_imported:
        lines.insert(insert_pytest_import_at, "import pytest\n")
        _log.debug(f"Inserted `import pytest` at line {insert_pytest_import_at}.")

    with open(filename, "w") as testfile:
        for line in lines:
            testfile.write(line)

    return 1


if __name__ == "__main__":
    with open(TEST_FILE_LIST, "r") as filelist:
        summary = ""
        for filename in filelist:
            success = process_testfile(filename.strip())
            if success == -1:
                summary += f"Manual intervention required in file {filename}!\n"

        summary += "mark_tests completed.\n"
        _log.info(summary)
