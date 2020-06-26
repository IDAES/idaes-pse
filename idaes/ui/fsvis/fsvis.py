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
import os
import requests
from requests.exceptions import ConnectionError
from slugify import slugify
import time
import webbrowser

from idaes.ui.fsvis.app import App as fsvis_server
from idaes.ui.flowsheet_serializer import FlowsheetSerializer

# serialize flowsheet and launch the app
def visualize(flowsheet, name, browser=True, overwrite=False):
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
    server = fsvis_server()
    url = f"http://{server.host}:{server.port}/app"
        
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
                                      params={'id': slugify(name)})
    if browser:
        success = webbrowser.open(url + f"?id={name}")
        print(f'Opened in browser window: {success}')
        print(f'{url}?id={name}')
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

#
def get_stored_fs(fs_id):
    return

# compare serialized model with stored to see if changes need to be made at all
def _model_has_changed(fs_id, new):
    return True

# record structural changes from the flowsheet
def update_stored_model(fs_id, new):
    return

# record manual layout changes from jointjs
def update_layout(fs_id, new):
    return

# sort the elements of the json representation of the flowsheet
# to facilitate comparison between versions
def _canonicalize_jointjs_output(json):
    return
