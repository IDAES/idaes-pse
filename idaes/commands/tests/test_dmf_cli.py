"""
Tests for DMF CLI
"""
# stdlib
import json
import os
from pathlib import Path
# third-party
from click.testing import CliRunner
import pytest
# package
from idaes.dmf.cli import init, register, info
from idaes.dmf.dmfbase import DMFConfig
from idaes.dmf.workspace import Workspace
from . import random_tempdir

@pytest.fixture(scope="module")
def runner():
    return CliRunner()

DATAFILE = "foo.txt"

@pytest.fixture()
def dmfctx(random_tempdir):
    """Switch DMF context to a random subdir, then switch back when done.
    """
    path = (random_tempdir / '.dmf').absolute()
    DMFConfig._filename = str(path)
    print("@@ dmfconfig.filename = " + str(path))
    origdir = os.getcwd()
    os.chdir(random_tempdir)
    with open(DATAFILE, 'w') as fp:
        fp.write("This is some sample data")
    yield path
    os.unlink(DATAFILE)
    DMFConfig._filename = str(Path('~/.dmf').expanduser())
    os.chdir(origdir)

def test_dmf_find(dmfctx, runner):
    filename = DATAFILE
    # register some objects in a workspace
    runner.invoke(init, ['ws', '--create', '--name', 'foo', '--desc', 'foo desc'])
    result = runner.invoke(register, [filename])
    id_all = result.output.strip()
    id_4 = id_all[:4]
    #
    result = runner.invoke(info, ["--no-color", id_4], catch_exceptions=False)
    assert result.exit_code == 0
    assert filename in result.output
    #
    result = runner.invoke(info, ["--no-color", "--format", "json", id_4],
                           catch_exceptions=False)
    assert result.exit_code == 0
    assert filename in result.output
    out = result.output.strip()
    assert out.startswith("{") and out.endswith("}")
    assert '"relations"' in out
    #
    result = runner.invoke(info, ["--no-color", "--format", "jsonc", id_4],
                           catch_exceptions=False)
    assert result.exit_code == 0
    assert filename in result.output
    out = result.output.strip()
    j = json.loads(out)
    assert len(j['datafiles']) == 1

def test_dmf_init(dmfctx, runner):
    result = runner.invoke(init, ['ws', '--create', '--name', 'foo',
        '--desc', 'foo workspace description'])
    assert result.exit_code == 0
    assert (Path('ws') / Workspace.WORKSPACE_CONFIG).exists()
    #
    result = runner.invoke(init, ['ws2', '--create'], input='foo\nfoo desc\n')
    assert result.exit_code == 0
    assert (Path('ws2') / Workspace.WORKSPACE_CONFIG).exists()
    #
    result = runner.invoke(init, ['doesnotexist'])
    assert result.exit_code != 0
    assert 'not found' in result.output
    #
    os.mkdir('some_random_directory')
    result = runner.invoke(init, ['some_random_directory'])
    assert result.exit_code != 0
    #
    result = runner.invoke(init, ['ws', '--create'])
    assert result.exit_code != 0
    assert 'exists' in result.output
 

# # setup:: dmf-ls

#     from pathlib import Path
#     from click.testing import CliRunner
#     from idaes.dmf.cli import init, ls, register
#     from idaes.dmf.dmfbase import DMFConfig
#     runner = CliRunner()

#     fsctx = runner.isolated_filesystem()
#     fsctx.__enter__()
#     DMFConfig._filename = str(Path('.dmf').absolute())
#     runner.invoke(init, ['ws', '--create', '--name', 'foo', '--desc', 'foo desc'])
#     files = [f"foo1{n}" for n in range(5)]
#     files.append("bar1")
#     for f in files:
#         with open(f, 'w') as fp:
#             fp.write("This is some sample data")
#         runner.invoke(register, [f])  # add file to DMF


#     fsctx.__exit__(None, None, None)
#     DMFConfig._filename = str(Path('~/.dmf').expanduser())


# # code:: dmf-ls
#     :hide:

#     result = runner.invoke(ls, ['--no-color'])
#     assert result.exit_code == 0
#     output1 = result.output
#     result = runner.invoke(ls, ['--no-color'])
#     assert result.output == output1

# # code:: dmf-ls
#     :hide:

#     result = runner.invoke(ls, ['--no-color', '-S', 'modified'])
#     assert result.exit_code == 0
#     output1 = result.output
#     result = runner.invoke(ls, ['--no-color', '--sort', 'modified'])
#     assert result.output == output1


# # setup:: dmf-register

#     from pathlib import Path
#     import re, json
#     from click.testing import CliRunner
#     from idaes.dmf.cli import init, register, info
#     from idaes.dmf.dmfbase import DMFConfig
#     runner = CliRunner()

#     fsctx = runner.isolated_filesystem()
#     fsctx.__enter__()
#     DMFConfig._filename = str(Path('.dmf').absolute())
#     runner.invoke(init, ['ws', '--create', '--name', 'foo', '--desc', 'foo desc'])
#     filename = "file.csv"
#     with open(filename, 'w') as fp:
#         fp.write("index,time,value\n1,0.1,1.0\n2,0.2,1.3\n")



#     fsctx.__exit__(None, None, None)
#     DMFConfig._filename = str(Path('~/.dmf').expanduser())


# # code:: dmf-register
#     :hide:

#     result = runner.invoke(register, ["file.csv", "--info"], catch_exceptions=False)
#     assert result.exit_code == 0
#     assert filename in result.output
#     assert "version" in result.output

# # code:: dmf-register
#     :hide:

#     result = runner.invoke(register, ["file.csv",], catch_exceptions=False)
#     assert result.exit_code != 0
#     result = runner.invoke(register, ["file.csv", "--no-unique"], catch_exceptions=False)
#     assert result.exit_code == 0

# # code:: dmf-register
#     :hide:

#     not_json = "notreally.json"
#     with open(not_json, "w") as fp:
#         fp.write("totally bogus\n")
#     result = runner.invoke(register, [not_json], catch_exceptions=False)
#     assert result.exit_code == 0
#     result = runner.invoke(register, [not_json, "--strict", "--no-unique"], catch_exceptions=False)
#     assert result.exit_code != 0

# # code:: dmf-register
#     :hide:

#     not_nb = "my.ipynb"
#     with open(not_nb, "w") as fp:
#         fp.write("foo\n")
#     result = runner.invoke(register, [not_nb, '-t', 'notebook'])
#     assert result.exit_code != 0
#     result = runner.invoke(register, [not_nb, '-t', 'data'])
#     assert result.exit_code == 0


# # code:: dmf-register
#     :hide:

#     for text_file in "shoebox", "shoes", "closet":
#         open(f"{text_file}.txt", "w")
#     result = runner.invoke(register, ["shoebox.txt"], catch_exceptions=False)
#     assert result.exit_code == 0
#     shoebox_id = result.output.strip()
#     result = runner.invoke(register, ["shoes.txt", "--contained", shoebox_id], catch_exceptions=False)
#     assert result.exit_code == 0
#     shoe_id = result.output.strip()
#     result = runner.invoke(info, [shoe_id, "--format", "jsonc"])
#     assert result.exit_code == 0

# # code:: dmf-register
#     :hide:

#     for text_file in "shoebox", "shoes", "closet":
#         open(f"{text_file}.txt", "w")
#     result = runner.invoke(register, ["closet.txt", "--is-subject",
#                            "--contained", shoebox_id], catch_exceptions=False)
#     assert result.exit_code == 0
#     closet_id = result.output.strip()
#     result = runner.invoke(info, [shoebox_id, "--format", "jsonc"])
#     assert result.exit_code == 0
#     data = json.loads(result.output)
#     assert len(data["relations"]) == 2
#     for rel in data["relations"]:
#         if rel["role"] == "subject":
#             assert rel["identifier"] == shoe_id
#         else:
#             assert rel["identifier"] == closet_id

# # setup:: dmf-related

#     from pathlib import Path
#     import re, json
#     from click.testing import CliRunner
#     from idaes.dmf import DMF, resource
#     from idaes.dmf.cli import init, related
#     from idaes.dmf.dmfbase import DMFConfig
#     runner = CliRunner()

#     fsctx = runner.isolated_filesystem()
#     fsctx.__enter__()
#     DMFConfig._filename = str(Path('.dmf').absolute())
#     runner.invoke(init, ['ws', '--create', '--name', 'foo', '--desc', 'foo desc'])
#     # add the fully-connected 4 resources
#     dmf = DMF()
#     rlist = [resource.Resource(value={"desc": ltr, "aliases": [ltr],
#                                "tags": ["graph"]})
#              for ltr in "ABCD"]
#     A_id = rlist[0].id  # root resource id, used in testcode
#     relation = resource.PR_USES
#     for r in rlist:
#         for r2 in rlist:
#             if r is r2:
#                 continue
#             resource.create_relation_args(r, relation, r2)
#     for r in rlist:
#         dmf.add(r)


#     fsctx.__exit__(None, None, None)
#     DMFConfig._filename = str(Path('~/.dmf').expanduser())


# # code:: dmf-related
#     :hide:

#     result = runner.invoke(related, [A_id, '--no-unicode', '--no-color'],
#                            catch_exceptions=False)
#     assert result.exit_code == 0
#     rlines = result.output.split('\n')
#     nrelations = sum(1 for _ in filter(lambda s: resource.PR_USES in s, rlines))
#     assert nrelations == 12  # 3 blocks of (1 + 3)


# # setup:: dmf-rm

#     from pathlib import Path
#     import re, json
#     from click.testing import CliRunner
#     from idaes.dmf.cli import init, register, ls, rm
#     from idaes.dmf.dmfbase import DMFConfig
#     runner = CliRunner()
#     # logging
#     import logging
#     log = logging.getLogger("cli")
#     _h = logging.FileHandler("/tmp/sphinx-dmf-cli.log")
#     log.addHandler(_h)
#     log.setLevel(logging.INFO)
#     # setup workspace
#     fsctx = runner.isolated_filesystem()
#     fsctx.__enter__()
#     DMFConfig._filename = str(Path('.dmf').absolute())
#     runner.invoke(init, ['ws', '--create', '--name', 'foo', '--desc', 'foo desc'])
#     # add some files named `file{1-5}`
#     for i in range(1,6):
#         filename = f"file{i}.txt"
#         with open(filename, 'w') as fp:
#             fp.write(f"file #{i}")
#         runner.invoke(register, [filename])


#     fsctx.__exit__(None, None, None)
#     DMFConfig._filename = str(Path('~/.dmf').expanduser())
#     # Comment to save log for debugging:
#     try:
#         Path("/tmp/sphinx-dmf-cli.log").unlink()
#     except FileNotFoundError:
#         pass

# # code:: dmf-rm
#     :hide:

#     result = runner.invoke(ls, ['--no-color', '--no-prefix'])
#     assert result.exit_code == 0
#     output1 = result.output
#     output1_lines = output1.split('\n')
#     rsrc_id = output1_lines[1].split()[0] # first token
#     log.info(f"resource id=`{rsrc_id}`")
#     result = runner.invoke(rm, [rsrc_id, "--yes", "--no-list"], catch_exceptions=False)
#     assert result.exit_code == 0
#     result = runner.invoke(ls, ['--no-color', '--no-prefix'])
#     assert result.exit_code == 0
#     output2 = result.output
#     output2_lines = output2.split('\n')
#     assert len(output2_lines) == len(output1_lines) - 1

# # code:: dmf-rm
#     :hide:

#     result = runner.invoke(ls, ['--no-color'])
#     assert result.exit_code == 0
#     output1 = result.output
#     output1_lines = output1.split('\n')
#     rsrc_id = output1_lines[1].split()[0] # first token
#     log.info(f"resource id=`{rsrc_id}`")
#     result = runner.invoke(rm, [rsrc_id, "--yes", "--no-list"], catch_exceptions=False)
#     assert result.exit_code == 0
#     result = runner.invoke(ls, ['--no-color', '--no-prefix'])
#     assert result.exit_code == 0
#     output2 = result.output
#     output2_lines = output2.split('\n')
#     assert len(output2_lines) == len(output1_lines) - 1

# # setup:: dmf-status

#     from pathlib import Path
#     from click.testing import CliRunner
#     from idaes.dmf.cli import init, status
#     from idaes.dmf.dmfbase import DMFConfig
#     runner = CliRunner()

#     fsctx = runner.isolated_filesystem()
#     fsctx.__enter__()
#     DMFConfig._filename = str(Path('.dmf').absolute())
#     runner.invoke(init, ['ws', '--create', '--name', 'foo', '--desc', 'foo desc'])


#     fsctx.__exit__(None, None, None)
#     DMFConfig._filename = str(Path('~/.dmf').expanduser())

# # code:: dmf-status
#     :hide:

#     result = runner.invoke(status, ['--no-color'])
#     assert result.exit_code == 0
#     assert "settings" in result.output
#     assert "name: foo" in result.output

# # code:: dmf-status
#     :hide:

#     result = runner.invoke(status, ['--no-color', '--show', 'files'])
#     assert result.exit_code == 0
#     assert "settings" in result.output
#     assert "name: foo" in result.output
#     assert "files:" in result.output

# # code:: dmf-status
#     :hide:

#     result = runner.invoke(status, ['--no-color', '--show', 'files',
#         '--show', 'htmldocs'])
#     assert result.exit_code == 0
#     assert "settings" in result.output
#     assert "name: foo" in result.output
#     assert "html" in result.output

# # code:: dmf-status
#     :hide:

#     result = runner.invoke(status, ['--no-color', '-a'])
#     assert result.exit_code == 0
#     assert "settings" in result.output
#     assert "name: foo" in result.output
#     assert "html" in result.output
#     assert "logging:" in result.output
