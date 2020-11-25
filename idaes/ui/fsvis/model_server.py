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
from idaes.ui.fsvis import persist


class ServerVariablesNotSetError(Exception):
    pass


class ModelServerHandler(http.server.BaseHTTPRequestHandler):
    """This is the server handler for the model server. This takes care of the communication of the model server
    """
    flowsheet, name = None, None

    # def do_OPTIONS(self):
    #     self.send_response(200, "ok")
    #     self.send_header("Access-Control-Allow-Origin", "*")
    #     self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
    #     self.send_header("Access-Control-Allow-Headers", "X-Requested-With")
    #     self.end_headers()

    def do_GET(self):
        """Get a model from the app.
        """
        value = json.dumps(self.flowsheet).encode(encoding="utf-8")
        self.send_response(200)
        #self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        #self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Content-type", "application/json")
        self.end_headers()
        self.wfile.write(value)
        #
        # if None in [self.flowsheet, self.name, self.flask_url]:
        #     raise ServerVariablesNotSetError
        #
        # serialized_flowsheet = FlowsheetSerializer().serialize(self.flowsheet, self.name)
        # r = requests.get(self.flask_url)
        #
        # diff_model, model_json = compare_models(r.json(), serialized_flowsheet, keep_old_model=True)
        # new_flowsheet = model_jointjs_conversion(diff_model, model_json)
        #
        # # Need to set the response and the headers or else you get CORS errors
        # self.send_response(200)
        # self.send_header("Access-Control-Allow-Origin", "*")
        # self.send_header("Access-Control-Allow-Headers", "Content-Type")
        # self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        # self.send_header("Content-type", "application/json")
        # http.server.SimpleHTTPRequestHandler.end_headers(self)
        #
        # self.wfile.write(json.dumps(new_flowsheet).encode(encoding="utf_8"))


class ModelServer(threading.Thread):
    def __init__(self, flowsheet, name, port=0):
        self._host = "127.0.0.1"
        if port == 0:
            # Find a free port on the host if there isn't an existing port
            self._port = find_free_port(host)
        else:
            self._port = port
        self._fs, self._name = flowsheet, name
        super().__init__()

    @property
    def addr(self):
        return self._host, self._port

    def run(self):
        ModelServerHandler.flowsheet, ModelServerHandler.name = self._fs, self._name
        with http.server.HTTPServer((self._host, self._port), ModelServerHandler) as httpd:
            try:
                httpd.serve_forever()
            except Exception:
                httpd.shutdown()
