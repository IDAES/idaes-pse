#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Jupyter magics for the DMF.
"""
# stdlib
import inspect
import logging
import os
from urllib.parse import urlparse
import webbrowser

# third-party
from IPython.core.magic import Magics, magics_class, line_magic
from IPython.display import display_markdown

# local
from . import dmfbase, errors, help, workspace

__author__ = "Dan Gunter"

# Logging

_log = logging.getLogger(__name__)


# Custom exception classes


class DMFMagicError(errors.DMFError):
    def __init__(self, errmsg, usermsg=None):
        super(DMFMagicError, self).__init__(errmsg)
        if usermsg is True:
            usermsg = errmsg
        self.message = usermsg


@magics_class
class DmfMagics(Magics):
    """Implement "magic" commands in Jupyter/IPython for
    interacting with the DMF and IDAES more generally.

    In order to allow easier testing, the functionality is
    broken into two classes. This class has the decorated method(s)
    for invoking the 'magics', and :class:`DmfMagicsImpl` has the
    state and functionality.
    """

    def __init__(self, shell):
        super(DmfMagics, self).__init__(shell)
        self._impl = DmfMagicsImpl(shell)
        self._last_ok = None

    @property
    def last_ok(self):
        return self._last_ok is not False  # None or True

    @line_magic
    def dmf(self, line):
        """DMF outer command.

        Example::

            %dmf <subcommand> [subcommand args..]
        """
        try:
            result = self._impl.dmf(line)
            self._last_ok = True
        except DMFMagicError as err:
            _log.error("Error with magic command: {}".format(err))
            result = None  # err.message
        return result


class DmfMagicsImpl(object):
    """State and implementation called by DmfMagics.

    On failure of any method, a `DMFMagicError` is raised, that
    should be handled by the line or cell magic that invoked it.
    """

    # Encoding which commands need 'init'.
    # Also encode args required to that command.
    # Key is command, value is #args: '*'=any, '+'=1 or more.
    _NEED_INIT = {"info": "*", "help": "+"}

    def __init__(self, shell):
        self._dmf = None
        self._shell = shell

    @property
    def initialized(self):
        return self._dmf is not None

    def dmf(self, line):
        """DMF outer command"""
        # Parse input into a subcommand and other tokens
        line = line.strip()
        if line == "":
            # No command will invoke "help"
            tokens = ["help"]
        else:
            tokens = line.split()
        # Find sub-method matching subcommand
        subcmd = tokens[0]
        submeth = getattr(self, "dmf_" + subcmd, None)
        if submeth is None:
            txt = "Unrecognized command `{}`".format(subcmd)
            # self._dmf_markdown(txt)
            raise DMFMagicError(txt)
        # Invoke sub-method
        params = tokens[1:]
        try:
            return submeth(*params)
        except Exception as err:
            msg = "Error in `%dmf {}`: {}".format(subcmd, err)
            # self._dmf_markdown(msg)
            raise DMFMagicError(msg, usermsg=msg)

    def dmf_init(self, path, *extra):
        """Initialize DMF (do this before most other commands).
        *Arguments*: `path ["create"]`

        Args:
            path (str): Full path to DMF home
            extra (str): Extra tokens. If 'create', then try to create the
                         path if it is not found.
        """
        kwargs, create = {}, False
        path = os.path.expanduser(path)
        if len(extra) > 0:
            if extra[0].lower() == "create":
                if os.path.exists(os.path.join(path, dmfbase.DMF.WORKSPACE_CONFIG)):
                    self._dmf_markdown(
                        'Cannot create new configuration at "{}": '
                        'file "{}" exists. Will try without '
                        '"create" option'.format(path, dmfbase.DMF.WORKSPACE_CONFIG)
                    )
                else:
                    kwargs = {"create": True, "add_defaults": True}
                    create = True
            else:
                _log.warning('Ignoring extra argument to "init"')
        try:
            self._dmf = dmfbase.DMF(path, **kwargs)
        except errors.WorkspaceNotFoundError:
            assert not create
            msg = (
                'Workspace not found at path "{}". '
                "If you want to create a new workspace, add the word "
                '"create", after the path, to the command.'.format(path)
            )
            # self._dmf_markdown(msg)
            raise DMFMagicError(msg)
        except errors.WorkspaceCannotCreateError as err:
            txt = 'Workspace could not be created at path "{}": {}'.format(path, err)
            raise DMFMagicError(txt, usermsg=txt)
        except errors.DMFError as err:
            raise DMFMagicError("Initializing workspace: {}".format(err), usermsg=True)

        self._dmf.set_meta({"name": os.path.basename(path)})
        if create:
            self._dmf_markdown('*Success!* Created new workspace at "{}"'.format(path))
        else:
            self._dmf_markdown('*Success!* Using workspace at "{}"'.format(path))

    def dmf_workspaces(self, *paths):
        """List DMF workspaces. Optionally takes one or more paths
        to use as a starting point. By default, start from current
        directory. *Arguments*: `[paths..]`

        Args:
            paths (List[str]): Paths to search, use "." by default
        """
        wslist = []
        if len(paths) == 0:
            paths = ["."]
        for root in paths:
            wslist.extend(workspace.find_workspaces(root))
        any_good_workspaces = False
        output_table = []
        for w in sorted(wslist):
            try:
                ws = workspace.Workspace(w)
                output_table.append((w, ws.name, ws.description))
                any_good_workspaces = True
            except errors.WorkspaceError:
                pass  # XXX: Should we print a warning?
        if not any_good_workspaces:
            # either no paths, or all paths raised an error
            return "No valid workspaces found\n"
        else:
            lines = ["| Path | Name | Description |", "| ---- | ---- | ----------- |"]
            for row in output_table:
                rowstr = "| {} | {} | {} |".format(row[0], row[1], row[2])
                lines.append(rowstr)
            listing = "\n".join(lines)
            self._dmf_markdown(listing)
            return "{:d} workspace(s) found".format(len(output_table))

    def dmf_list(self):
        """List resources in the current workspace. *Arguments*: none."""
        self._init_required("list")
        lines = [
            "| ID | Name(s) | Type | Modified | Description | ",
            "| -- | ------- | ---- | -------- | ----------- |",
        ]
        for rsrc in self._dmf.find():
            msince = rsrc.v["modified"]
            rowstr = "| {id} | {names} | {type} | {mdate} | {desc} |".format(
                id=rsrc.id,
                names=",".join(rsrc.v["aliases"]),
                type=rsrc.type,
                mdate=msince,
                desc=rsrc.v["desc"],
            )
            lines.append(rowstr)
        listing = "\n".join(lines)
        self._dmf_markdown(listing)
        return True

    def dmf_status(self, *topics):
        """Provide information about DMF current state. *Arguments*: none

        Args:
            (List[str]) topics: List of topics

        Returns:
            None
        """
        self._init_required("info")
        if topics:
            self._dmf_markdown("Sorry, no topic-specific info yet available")
            raise DMFMagicError("Topics not supported")
        # configuration info
        text_lines = ["## Configuration"]
        for key, value in self._dmf.meta.items():
            hdr = "  * {}".format(key)
            if isinstance(value, list):
                text_lines.append("{}:".format(hdr))
                for v in value:
                    text_lines.append("    - {}".format(v))
            elif isinstance(value, dict):
                text_lines.append("{}:".format(hdr))
                for k2, v2 in value.items():
                    text_lines.append("    - {}: {}".format(k2, v2))
            else:
                text_lines.append("{}: {}".format(hdr, value))
        conf_info = "\n".join(text_lines)
        # all info
        all_info = "\n".join((conf_info,))
        self._dmf_markdown(all_info)

    def dmf_help(self, *names):
        """Provide help on IDAES objects and classes.
        Invoking with no arguments gives general help.
        Invoking with one argument looks for help in the docs
        on the given object or class.
        *Arguments*: `[name]`.
        """
        if len(names) == 0:
            # give some general help for magics
            return self._magics_help()
        self._init_required("help")
        if len(names) > 1:
            _log.warning(
                "DMF Help is restricted to only one object or class "
                "at a time. Ignoring trailing arguments ({})".format(
                    " ".join(names[1:])
                )
            )
        name = names[0]
        # Check some special names first.
        # To re-use the <module>.<class> mechanism, translate them into
        # some pseudo-classes in the "dmf.help" pseudo-module.
        p_module, p_class = None, None
        if name.lower() == "dmf":
            p_module = "dmf"
            p_class = "DMF"
        elif name.lower() in ("help",):
            p_module = "dmf"
            p_class = "Help"
        elif name.lower() in ("idaes", "*", ""):
            p_module = "idaes"
            p_class = "Home"
        # Get the help
        help_locations = None
        if p_module is not None:
            try:
                help_locations = help.get_html_docs(self._dmf, p_module, p_class)
            except DMFMagicError as err:
                _log.debug(
                    "Getting help for pseudo-class {}::{}, error: {}".format(
                        p_module, p_class, err
                    )
                )
                # for Result
                name = "pseudo-module {}.{}".format(p_module, p_class)
        else:
            try:
                help_locations = self._find_help_for_object(name)
            except DMFMagicError as err:
                _log.debug("Getting help for object {}, error: {}".format(name, err))
        # Result
        if help_locations:
            self._show_help_in_browser(help_locations)
        else:
            self._dmf_markdown('No Sphinx docs found for "{}"'.format(name))
        return None

    def _init_required(self, subcmd):
        """If no DMF, display init-required message and
        raise an exception. Otherwise, do nothing and return False.
        """
        if self._dmf is not None:
            return False
        msg = "Must call `init` before command `{}`".format(subcmd)
        self._dmf_markdown(msg)
        raise DMFMagicError(msg)

    @staticmethod
    def _dmf_markdown(text):
        display_markdown(text, raw=True)

    def _find_help_for_object(self, name):
        """Get object by evaluating name as an expression.
        If that doesn't work, try just using it as a string.
        """
        obj, oname = None, None
        try:
            obj = self._shell.ev(name)
            _log.debug("Looking for HTML docs for object: {}".format(obj))
        except Exception:
            oname = name
            _log.debug("Looking for HTML docs for object: {}".format(oname))
        try:
            result = help.find_html_docs(self._dmf, obj=obj, obj_name=oname)
        except (ValueError, AttributeError):
            if obj is not None:
                raise DMFMagicError("Cannot find help for {}".format(obj))
            else:
                raise DMFMagicError("Cannot find help for {}".format(oname))
        return result

    @staticmethod
    def _show_help_in_browser(help_locations):
        """Open help docs in the browser."""
        for url in help_locations:
            if urlparse(url).scheme == "":
                url = "file://" + url
            _log.debug('Opening URL "{}"'.format(url))
            webbrowser.open_new(url)

    def _magics_help(self):
        """Introspect to give a list of commands."""
        # Build a dictionary of magic methods and their
        # descriptions (from their docstrings).
        help_dict = {"dmf": {}, "idaes": {}}
        for name, meth in inspect.getmembers(self, inspect.ismethod):
            if name.startswith("dmf_"):
                help_dict["dmf"][name] = self._extract_help_text(meth)
            elif name.startswith("idaes_"):
                help_dict["idaes"][name] = self._extract_help_text(meth)
        # build help text from dict
        txt = ""
        for mname in "idaes", "dmf":
            if len(help_dict[mname]) == 0:
                continue
            txt_list = []
            for name in sorted(help_dict[mname].keys()):
                cmd = name[len(mname) + 1 :]
                cmd_help = "-  `%{} {}` - {}".format(mname, cmd, help_dict[mname][name])
                txt_list.append(cmd_help)
            section = "{}## {} magic commands".format(
                "\n\n" if txt else "", mname.upper()
            )
            txt += section + "\n\n" + "\n".join(txt_list)
        return self._dmf_markdown(txt)

    @staticmethod
    def _extract_help_text(meth):
        doc = inspect.getdoc(meth)
        # extract first sentence (strip off period)
        sentence = doc[: doc.find("\n\n")].strip(".")
        # take out line breaks from sentence
        hdr = sentence.replace("\n", " ")
        return hdr


###############################################################################
# Register at import time
###############################################################################

_registered = False


def register():
    """Register with IPython on import (once)."""
    global _registered
    if _registered:
        return
    try:
        ip = get_ipython()  # noqa: F821  pylint: disable=undefined-variable
        _log.debug("Registering DMF magics")
        ip.register_magics(DmfMagics)
    except:  # noqa: E722
        pass
    _registered = True


register()
