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

from idaes.ui.fsvis.app import App as fsvis_server
from idaes.ui.fsvis.app import find_free_port
from idaes.ui.fsvis.model_server import ModelServer
from idaes.ui.flowsheet_serializer import FlowsheetSerializer


# See model_server.py for more information about this module's role in the visualizer


# serialize flowsheet and launch the app
def visualize(flowsheet, name, browser=True, overwrite=False, model_server_host='127.0.0.1'):
    """Visualizes the flowsheet, assigning it to the given name. 
    
    Attempts to
    open a browser window to display the visualization app, as well as
    directly showing the URL in case the browser fails to open.

    Usage example:
    
    m = ConcreteModel()
    m.fs = FlowsheetBlock(...)
    ...
    visualize(m.fs, "draftview")
        
    Args: 
        flowsheet: An IDAES flowsheetBlock to be visualized.
        name: A name string to assign to the visualization. This name will be
            attached to the visualization instance as long as the visualization app
            server stays running (e.g. until the parent python kernel is shut down)

    Returns:   
        None.

    Raises:
        None.#TODO
    """
    # Start the model server that contains a reference to the model so that the flask 
    # server can ping it when refresh is called in order to get the updated model
    server = fsvis_server()

    url = f"http://{server.host}:{server.port}/app"

    model_server = ModelServer.getInstance(flowsheet, name, f"http://{server.host}:{server.port}/fs?id={name}", model_server_host)

    model_server_url = f"http://{model_server_host}:{model_server.port}"

    # Check if the {name}.viz file exists and overwrite is not true. If it was True
    # then we want to serialize the flowsheet and reset to the original
    file_path = os.path.expandvars(os.path.join(os.path.expanduser("~"), ".idaes", "viz", f"{name}.viz"))
    if os.path.isfile(file_path) and not overwrite:
        print(f"Model {name} visualization exists. Reloading the visualization.")
        print("If you don't want to load the existing visualization specify overwrite=True "
              "when calling visualize")
        with open(file_path, "r") as viz_file:
            serialized_flowsheet = json.load(viz_file)
    else:
        serialized_flowsheet = FlowsheetSerializer().serialize(flowsheet, name)

    repeat_until_connection_available(requests.post, url, json=serialized_flowsheet, 
                                      params={'id': slugify(name), "modelurl": model_server_url})
    if browser:
        success = webbrowser.open(url + f"?id={name}&modelurl={model_server_url}")
        print(f'Opened in browser window: {success}')
        print(f'{url}?id={name}&modelurl={model_server_url}')

    return server


# possibly should be changed to _repeat_until_connection_available()
def repeat_until_connection_available(f, *args, retries=127, **kwargs):
    for i in range(retries):
        try:
            print(f'attempt {i} of {retries}')
            return f(*args, **kwargs)
        except ConnectionError as e:
            time.sleep(0.1)
            print(f'connection error: attempt {i}; {e}')
            continue # consider logging
            
    # raise ConnectionRefusedError?? how
    # or maybe just print to debug and stop?
        