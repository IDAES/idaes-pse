##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
# stdlib
import webbrowser
# package
from idaes import logger
from .model_server import FlowsheetServer

_log = logger.getLogger(__name__)

web_server = None


def visualize(
    flowsheet, name: str = "flowsheet", save_as=None, browser: bool = True, port=None,
        log_level=logger.WARNING
):
    """Visualizes the flowsheet in a web application.
    
    Opens a browser window to display the visualization app, as well as
    directly showing the URL in case the browser fails to open.

    Args:
        flowsheet: IDAES flowsheet to visualize
        name: Name of flowsheet to display as the title of the visualization
        save_as: If a string or path then save to a file.
        browser: If true, open a browser
        log_level: An IDAES logging level, see :mod:`idaes.logger`, to set for all the visualiztion

    Returns:
        None.

    Raises:
        ValueError if the data storage at 'save_as' can't be opened
    """
    global web_server

    _init_logging(log_level)

    if web_server is None:
        web_server = FlowsheetServer(port=port)
        web_server.start()
    else:
        _log.info(f"Using HTTP server on localhost, port {web_server.port}")

    web_server.add_flowsheet(name, flowsheet, save_as)

    # Open a browser window for the UI
    url = f"http://localhost:{web_server.port}/app"
    if browser:
        success = webbrowser.open(url + f"?id={name}")
        _log.debug(f"Opened in browser window: {success}")


def _init_logging(lvl):
    ui_logger = logger.getIdaesLogger("ui", level=lvl, tag="ui")
    ui_logger.setLevel(lvl)
