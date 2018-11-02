##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
# 
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Jupyter magics for the DMF.
"""
# stdlib
import inspect
import os
import webbrowser
# third-party
import pendulum
from IPython.core.magic import Magics, magics_class, line_magic
from IPython.display import display_markdown
import six
# local
from . import dmfbase, util, errors, help, workspace

__author__ = 'Dan Gunter <dkgunter@lbl.gov>'

_log = util.get_logger('magics')


@magics_class
class DmfMagics(Magics):

    # Encoding which commands need 'init'.
    # Key is command, value is #args: '*'=any, '+'=1 or more.
    NEED_INIT_CMD = {'info': '*', 'help': '+'}

    def __init__(self, shell):
        super(DmfMagics, self).__init__(shell)
        self._dmf = None

    @line_magic
    def dmf(self, line):
        """DMF outer command
        """
        line = line.strip()
        if line == '':
            tokens = ['help']
        else:
            tokens = line.split()
        subcmd = tokens[0]
        if subcmd == 'workspaces':
            pass
        elif self._dmf is None:
                required, why = self._init_required(subcmd, tokens)
                if required:
                    return self._dmf_markdown(why)
        submeth = getattr(self, 'dmf_' + subcmd, None)
        if submeth is None:
            txt = "Unrecognized command `{}`".format(subcmd)
            return self._dmf_markdown(txt)
        params = tokens[1:]
        try:
            return submeth(*params)
        except Exception as err:
            msg = 'Error in `%dmf {}`: {}'.format(subcmd, err)
            return self._dmf_markdown(msg)

    @line_magic
    def idaes(self, line):
        """%idaes magic
        """
        line = line.strip()
        if line == '':
            tokens = ['help']
        else:
            tokens = line.split()
        subcmd = tokens[0]
        if subcmd == 'help':
            return self.idaes_help(*tokens[1:])
        else:
            txt = "Unrecognized command `{}`".format(subcmd)
            return self._dmf_markdown(txt)

    def _init_required(self, subcmd, tokens):
        """Return whether init is required before this particular
        subcommand invocation and, if so, why.
        """
        req, why = False, ''
        if subcmd in self.NEED_INIT_CMD:
            nargs_code = self.NEED_INIT_CMD[subcmd]
            if nargs_code == '*' or (nargs_code == '+' and len(tokens) > 1):
                req = True
                nwhy = '' if nargs_code == '*' else ' with 1 or more args'
                why = 'Must call `init` before command `{}`{}'.format(
                    subcmd, nwhy)
        return req, why

    def dmf_init(self, path, *extra):
        """Initialize DMF (do this before most other commands).

        Args:
            path (str): Full path to DMF home
        """
        kwargs, create = {}, True
        if len(extra) > 0:
            if extra[0].lower() == 'create':
                kwargs = {'create': True, 'add_defaults': True}
                create = True
            else:
                _log.warn('Ignoring extra argument to "init"')

        try:
            self._dmf = dmfbase.DMF(path, **kwargs)
        except errors.DMFWorkspaceNotFoundError:
            if not create:
                msg = 'Workspace not found at path "{}". ' \
                      'If you want to create a new workspace, add the word ' \
                      '"create", after the path, to the command.'.format(path)
                return msg
            else:
                return 'Workspace could not be created at path "{}"'\
                       .format(path)
        except errors.DMFBadWorkspaceError as err:
            return 'Error initializing workspace: {}'.format(err)

        self._dmf.set_meta({'name': os.path.basename(path)})
        return None

    def dmf_workspaces(self, *paths):
        """List DMF workspaces.

        Args:
            paths (List[str]): Paths to search, use "." by default
        """
        wslist = []
        if len(paths) == 0:
            paths = ['.']
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
            return('ERROR: No valid workspaces found\n')
        else:
            lines = ['| Path | Name | Description |',
                     '| ---- | ---- | ----------- |']
            for row in output_table:
                rowstr = '| {} | {} | {} |'.format(row[0], row[1], row[2])
                lines.append(rowstr)
            listing = '\n'.join(lines)
            self._dmf_markdown(listing)

    def dmf_list(self):
        """List resources in the current workspace.
        """
        lines = ['| ID | Name(s) | Type | Modified | Description | ',
                 '| -- | ------- | ---- | -------- | ----------- |']
        for rsrc in self._dmf.find():
            msince = pendulum.from_timestamp(rsrc.modified).diff_for_humans()
            rowstr = '| {id} | {names} | {type} | {mdate} | {desc} |'.format(
                id=rsrc.id_, names=','.join(rsrc.aliases), type=rsrc.type,
                mdate=msince, desc=rsrc.desc)
            lines.append(rowstr)
        listing = '\n'.join(lines)
        return self._dmf_markdown(listing)

    def dmf_info(self, *topics):
        """Provide information about DMF current state for whatever
        'topics' are provided. With no topic, provide general information
        about the configuration.

        Args:
            (List[str]) topics: List of topics

        Returns:
            None
        """
        if topics:
            self._dmf_markdown('Sorry, no topic-specific info yet available')
            return
        # configuration info
        text_lines = ['## Configuration']
        for key, value in six.iteritems(self._dmf.meta):
            hdr = '  * {}'.format(key)
            if isinstance(value, list):
                text_lines.append('{}:'.format(hdr))
                for v in value:
                    text_lines.append('    - {}'.format(v))
            elif isinstance(value, dict):
                text_lines.append('{}:'.format(hdr))
                for k2, v2 in six.iteritems(value):
                    text_lines.append('    - {}: {}'.format(k2, v2))
            else:
                text_lines.append('{}: {}'.format(hdr, value))
        conf_info = '\n'.join(text_lines)
        # all info
        all_info = '\n'.join((conf_info, ))
        self._dmf_markdown(all_info)

    def _dmf_markdown(self, text):
        display_markdown(text, raw=True)

    def dmf_help(self, *names):
        """Provide help on IDAES objects and classes.

        Invoking with no arguments gives general help.
        Invoking with one or more arguments looks for help in the docs
        on the given objects or classes.
        """
        if len(names) == 0:
            # give some general help for magics
            return self._magics_help()
        if len(names) > 1:
            raise ValueError('Only one object or class at a time')
        name = names[0]
        # Check some special names first.
        # To re-use the <module>.<class> mechanism, translate them into
        # some pseudo-classes in the "dmf.help" pseudo-module.
        if name.lower() in ('help', 'dmf'):
            p_module = 'dmf.help'
            p_class = name.title()
            helpfiles = help.get_html_docs(self._dmf, p_module, p_class)
        else:
            helpfiles = self._find_help_for_object(name)
        if helpfiles:
            self._show_help_in_browser(helpfiles)
        else:
            return 'No Sphinx docs found for "{}"'.format(name)

    def idaes_help(self, *names):
        """Provide help on IDAES objects and classes.

        Invoking with no arguments gives general help.
        Invoking with one or more arguments looks for help in the docs
        on the given objects or classes.
        """
        if len(names) == 0:
            # give some general help for magics
            return self._magics_help()
        if len(names) > 1:
            raise ValueError('Only one object or class at a time')
        name = names[0]
        # Check some special names first.
        # To re-use the <module>.<class> mechanism, translate them into
        # some pseudo-classes in the "idaes.help" pseudo-module.
        if name.lower() in ('help', 'idaes'):
            p_module = 'idaes.help'
            p_class = name.title()
            helpfiles = help.get_html_docs(self._dmf, p_module, p_class)
        else:
            helpfiles = self._find_help_for_object(name)
        if helpfiles:
            self._show_help_in_browser(helpfiles)
        else:
            return 'No Sphinx docs found for "{}"'.format(name)

    def _find_help_for_object(self, name):
        """Get object by evaluating name as an expression."""
        try:
            obj = self.shell.ev(name)
        except Exception as err:
            raise errors.DmfError('Cannot evaluate object/class for help: '
                                  '{}'.format(err))
        _log.debug('Looking for HTML docs for object: {}'.format(obj))
        return help.find_html_docs(self._dmf, obj)

    @staticmethod
    def _show_help_in_browser(helpfiles):
        """Open help docs in the browser."""
        first_option = webbrowser._tryorder[0]
        if first_option in ('xdg-open', 'chromium'):
            # for these browsers, prefixing with file:// allows
            # the anchors in the URL to work
            url = 'file://' + helpfiles[0]
        else:
            url = helpfiles[0]
        _log.debug('Opening URL "{}"'.format(url))
        webbrowser.open_new(url)

    def _magics_help(self):
        """Introspect to give a list of commands."""
        # Build a dictionary of magic methods and their
        # descriptions (from their docstrings).
        help_dict = {'dmf': {}, 'idaes': {}}
        for name, meth in inspect.getmembers(self, inspect.ismethod):
            if name.startswith('dmf_'):
                help_dict['dmf'][name] = self._extract_help_text(meth)
            elif name.startswith('idaes_'):
                help_dict['idaes'][name] = self._extract_help_text(meth)
        # build help text from dict
        txt = ''
        for mname in 'idaes', 'dmf':
            txt_list = []
            for name in sorted(help_dict[mname].keys()):
                cmd = name[len(mname) + 1:]
                cmd_help = '* `%{} {}` - {}'.format(
                    mname, cmd, help_dict[mname][name])
                txt_list.append(cmd_help)
            section = '{}{} magic commands:'.format(
                '\n\n' if txt else '', mname.upper())
            txt += section + '\n' + '\n'.join(txt_list)
        return self._dmf_markdown(txt)

    @staticmethod
    def _extract_help_text(meth):
        doc = inspect.getdoc(meth)
        # extract first sentence (strip off period)
        sentence = doc[:doc.find('.')]
        # take out line breaks from sentence
        hdr = sentence.replace('\n', ' ')
        return hdr

# def dmf_add(self, name, *params):
#     """Add one resource
#     """
#     try:
#         obj = self.shell.ev(name)
#     except Exception as err:
#         raise errors.DmfError('Cannot evaluate object/class for "add": {}'
#                               .format(err))
#     if isinstance(obj, idaes_models.core.FlowsheetBlockData):
#         rsrc = self._rfactory.create_flowsheet(obj)
#         self._dmf.add(rsrc)
#     else:
#         return 'Unknown object type "{}". Cannot add as a DMF resource.'\
#             .format(type(obj))


_registered = False


def register():
    global _registered
    if _registered:
        return
    try:
        ip = get_ipython()  # noqa: F821
        _log.info('Registering DMF magics')
        ip.register_magics(DmfMagics)
    except:                 # noqa: E722
        pass
    _registered = True


register()
