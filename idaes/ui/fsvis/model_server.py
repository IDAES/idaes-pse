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

# Communication between the difference pieces of the visualization
# app.py - Flask server
# model_server.py - Server that has a reference to the model and allows the javascript to make
# requests to get the updated model
# fsvis.py - Starts the flask server and the model server and opens the web browser to the 
# visualization
# 
# Example communication:
# 1. User calls m.fs.visualize("test")
# 2. flowsheet_model.py calls idaes.ui.fsvis.visualize()
# 3. fsvis.py starts the flask server (app.py) and the model server (model_server.py)
# 4. fsvis.py opens a web browser to localhost:port/app?id=test&modelurl=localhost:model_port
# 5. Web browser requests the serialized model from the flask server and displays it using rappid
# 6. User makes a change to the model. model_server.py's model gets updated because the model reference
#    was updated
# 7. User clicks Refresh Graph on the javascript
# 8. Javscript makes an ajax GET request to the model server.
# 9. model_server.py's SimpleModelServerHandler do_GET gets called and serializes the updated model 
#    and returns it
# 10. Javascript receives the serialized model and displays it


import http.server
import json
import requests
import threading

from idaes.ui.flowsheet_comparer import compare_models, model_jointjs_conversion
from idaes.ui.flowsheet_serializer import FlowsheetSerializer
from idaes.ui.fsvis.app import find_free_port
from idaes.ui.fsvis.server import DataStorage


class ServerVariablesNotSetError(Exception):
    pass


class SimpleModelServerHandler(http.server.BaseHTTPRequestHandler):
    """
    This is the server handler for the model server. This takes care of the communcation of the model server
    """
    flowsheet = None
    name = None
    flask_url = None

    # Overrides the BaseHTTPRequestHandler's do_OPTIONS to allow cross origin requests
    def do_OPTIONS(self):           
        self.send_response(200, "ok")       
        self.send_header("Access-Control-Allow-Origin", "*")                
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "X-Requested-With") 

    # Overrides the BaseHTTPRequestHandler's do_GET in order to serve the new serialized model
    def do_GET(self):
        if None in [self.flowsheet, self.name, self.flask_url]:
            raise ServerVariablesNotSetError

        serialized_flowsheet = FlowsheetSerializer().serialize(self.flowsheet, self.name)
        r = requests.get(self.flask_url)
        
        diff_model, model_json = compare_models(r.json(), serialized_flowsheet, keep_old_model=True)
        new_flowsheet = model_jointjs_conversion(diff_model, model_json)

        # Need to set the response and the headers or else you get CORS errors
        self.send_response(200)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Content-type", "application/json")
        http.server.SimpleHTTPRequestHandler.end_headers(self)

        self.wfile.write(json.dumps(new_flowsheet).encode(encoding="utf_8"))


class ModelServer():
    """Singleton to run a server with a reference to the IDAES model in a thread.
    """
    __instance = None

    @staticmethod
    def getInstance(flowsheet, name, flask_url, host, port=0):
       """ Static access method. """
       if ModelServer.__instance == None:
          ModelServer(flowsheet, name, flask_url, host, port)
       return ModelServer.__instance

    def __init__(self, flowsheet, name, flask_url, host, port=0):
        """ Virtually private constructor. """
        if ModelServer.__instance != None:
            return
        else:
            ModelServer.__instance = self
            self.port = port
            self._start_model_server(flowsheet, name, flask_url, host)

    def _start_model_server_thread(self, flowsheet, name, flask_url, host):
        daemon = threading.Thread(name='daemon_model_server',
                                  target=self._setup_model_server,
                                  args=(flowsheet, name, flask_url, host))
        daemon.setDaemon(True) # Set as a daemon so it will be killed once the main thread is dead.
        daemon.start()

    # Start a daemon thread for the model server

    def _start_model_server(self, flowsheet, name, flask_url, host):
        if self.port == 0:
            # Find a free port on the host if there isn't an existing port
            self.port = find_free_port(host)

        # Create a thread and start the server in the thread
        self._start_model_server_thread(flowsheet, name, flask_url, host)
        

    # Set up a http server for serving the model to the flask server refresh button
    def _setup_model_server(self, flowsheet, name, flask_url, host):
        SimpleModelServerHandler.flowsheet = flowsheet
        SimpleModelServerHandler.name = name
        SimpleModelServerHandler.flask_url = flask_url

        with http.server.HTTPServer((host, self.port), SimpleModelServerHandler) as httpd:
            try:
                httpd.serve_forever()
            except Exception:
                httpd.shutdown()
