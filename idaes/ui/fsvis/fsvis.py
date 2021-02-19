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
    flowsheet,
    name: str = "flowsheet",
    save_as=None,
    browser: bool = True,
    port: int = None,
    log_level: int = logger.WARNING,
    quiet: bool = False,
) -> int:
    """Visualize the flowsheet in a web application.

    The web application is started in a separate thread and this function returns immediately.

    Also open a browser window to display the visualization app.
    The URL is also printed (unless ``quiet`` is True).

    **Note**: The visualization server runs in its own thread. If the program that it is running in stops,
    the visualization UI will not be able to save or refresh its view. This is not an issue in a
    `REPL <https://en.wikipedia.org/wiki/Read%E2%80%93eval%E2%80%93print_loop>`_
    like the Python console, IPython, or Jupyter Notebook,
    since these all run until the user explicitly closes them. But if you are running from a script, you need to do
    something to avoid having the program exit after the `visualize()` method returns (which happens very quickly).
    For example, loop forever in a try/catch clause that will handle KeyboardInterrupt exceptions::

        # Example code for a script, to keep program running after starting visualize() thread
        my_model.fs.visualize()  # this returns immediately
        try:
            print("Type ^C to stop the program")
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            print("Program stopped")

    Args:
        flowsheet: IDAES flowsheet to visualize
        name: Name of flowsheet to display as the title of the visualization
        save_as: If a string or path then save to a file.
        browser: If true, open a browser
        port: Start listening on this port. If not given, find an open port.
        log_level: An IDAES logging level, which is a superset of the built-in :mod:`logging` module levels.
                  See :mod:`idaes.logger` for details
        quiet: If True, suppress printing any messages to standard output (console)

    Returns:
        Port number where server is listening

    Raises:
        ValueError: if the data storage at 'save_as' can't be opened
    """
    global web_server

    _init_logging(log_level)

    if web_server is None:
        web_server = FlowsheetServer(port=port)
        web_server.start()
        if not quiet:
            print("Started visualization server")
    else:
        _log.info(f"Using HTTP server on localhost, port {web_server.port}")

    new_name = web_server.add_flowsheet(name, flowsheet, save_as)
    if new_name != name:
        _log.warning(f"Flowsheet name changed: old='{name}' new='{new_name}'")
        if not quiet:
            print(f"Flowsheet name changed to '{new_name}'")
        name = new_name

    # Open a browser window for the UI
    url = f"http://localhost:{web_server.port}/app?id={name}"
    if browser:
        success = webbrowser.open(url)
        if success:
            _log.debug(f"Flowsheet opened in browser window")
        else:
            _log.warning(f"Could not open flowsheet URL '{url}' in browser")
            if not quiet:
                print("Error: Unable to open flowsheet in web browser.")

    if not quiet:
        print(f"Flowsheet visualization at: {url}")

    return web_server.port


def _init_logging(lvl):
    ui_logger = logger.getIdaesLogger("ui", level=lvl, tag="ui")
    ui_logger.setLevel(lvl)
