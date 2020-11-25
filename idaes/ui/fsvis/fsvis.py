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
# server backend code for fsvis
import json
import logging
import os
import requests
from requests.exceptions import ConnectionError
from slugify import slugify
import time
import webbrowser

from .app import App, find_free_port
from .persist import DataStore, MemoryDataStore
from .model_server import ModelServer
from ..flowsheet_serializer import FlowsheetSerializer

_log = logging.getLogger(__name__)


def visualize(
    flowsheet,
    name: str = "flowsheet",
    save_as=None,
    browser: bool = True
):
    """Visualizes the flowsheet in a web application.
    
    Opens a browser window to display the visualization app, as well as
    directly showing the URL in case the browser fails to open.

    Args:
        flowsheet: IDAES flowsheet to visualize
        name: Name of flowsheet to display as the title of the visualization
        save_as: If a string or path then save to a file.
        browser: If true, open a browser

    Returns:
        None.

    Raises:
        ValueError if the data storage at 'save_as' can't be opened
    """
    # Get singleton for the web app
    web_app = App()

    # Add data store
    if save_as is None:
        store = MemoryDataStore()
    else:
        try:
            store = DataStore.create(dest=save_as)
        except ValueError as err:
            raise ValueError(f"Cannot create new data store: {err}")
    web_app.data_storage_manager.add(name, store)

    # Start model server, and add its address to flask
    model_server = ModelServer(flowsheet, name)
    model_server.start()
    web_app.model_server_manager.add(name, model_server.addr)

    # Open a browser window for the UI
    url = f"http://{web_app.host}:{web_app.port}/app"
    if browser:
        success = webbrowser.open(url + f"?id={name}")
        _log.debug(f"Opened in browser window: {success}")



    # # Check if the {name}.viz file exists and overwrite is not true. If it was True
    # # then we want to serialize the flowsheet and reset to the original
    # file_path = os.path.expandvars(
    #     os.path.join(os.path.expanduser("~"), ".idaes", "viz", f"{name}.viz")
    # )
    # if os.path.isfile(file_path) and not overwrite:
    #     print(f"Model {name} visualization exists. Reloading the visualization.")
    #     print(
    #         "If you don't want to load the existing visualization specify overwrite=True "
    #         "when calling visualize"
    #     )
    #     with open(file_path, "r") as viz_file:
    #         serialized_flowsheet = json.load(viz_file)
    # else:
    #     serialized_flowsheet = FlowsheetSerializer().serialize(flowsheet, name)
    # serialized_flowsheet = FlowsheetSerializer().serialize(flowsheet, name)
    #
    # # Set up the server URL
    # url = f"http://{server.host}:{server.port}/app"
    # model_server = ModelServer.getInstance(
    #     flowsheet,
    #     name,
    #     f"http://{server.host}:{server.port}/fs?id={name}",
    #     model_server_host,
    # )
    # model_server_url = f"http://{model_server_host}:{model_server.port}"
    # try_to_connect(
    #     requests.post,
    #     url,
    #     json=serialized_flowsheet,
    #     params={"id": slugify(name), "modelurl": model_server_url},
    # )


def try_to_connect(f, *args, retries=127, **kwargs):
    for i in range(retries):
        try:
            _log.debug(f"attempt {i} of {retries}")
            return f(*args, **kwargs)
        except ConnectionError as e:
            time.sleep(0.1)
            _log.info(f"connection error: attempt {i}; {e}")
            continue
