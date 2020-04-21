# server backend code for fsvis
import requests
from requests.exceptions import ConnectionError
import time
import webbrowser

from idaes.ui.fsvis.app import App as fsvis_server
from idaes.ui.flowsheet_serializer import FlowsheetSerializer

# serialize flowsheet and launch the app
def visualize(flowsheet, name):
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
    
    # TODO: make sure `name` is valid for use as a URL query string!! E.g. no spaces...

    serialized_flowsheet = FlowsheetSerializer().serialize(flowsheet, name)
    
    repeat_until_connection_available(requests.post, url, json=serialized_flowsheet, 
                        params={'id': name})
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
