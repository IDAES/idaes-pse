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
from pathlib import Path
import socket
import threading
from urllib.parse import urlparse

from idaes.logger import getLogger
# from idaes.ui.flowsheet_comparer import compare_models, model_jointjs_conversion
from . import persist
from ..flowsheet_serializer import FlowsheetSerializer

_log = getLogger(__name__)

# Directories
_this_dir = Path(__file__).parent.absolute()
_static_dir = _this_dir / "static"
_template_dir = _this_dir / "templates"


class FlowsheetServer(http.server.HTTPServer):
    """A simple HTTP server that runs in its own thread.

    This server is used for *all* models for a given process, so every request needs to contain
    the ID of the model that should be used in that transaction.
    """

    def __init__(self, port=None):
        """Create HTTP server
        """
        self._port = port or find_free_port()
        print(f"@@ found free port {self._port}")
        super().__init__(("127.0.0.1", self._port), FlowsheetServerHandler)
        self._dsm = persist.DataStoreManager()
        self._flowsheets = {}
        self._fss = FlowsheetSerializer()

    @property
    def port(self):
        return self._port

    def start(self):
        """Start the server, which will spawn a thread.
        """
        self._thr = threading.Thread(target=self._run)
        self._thr.setDaemon(True)
        self._thr.start()

    def _run(self):
        """Run in a separate thread.
        """
        _log.info(f"Serve forever on localhost:{self._port}")
        try:
            self.serve_forever()
        except Exception:
            _log.info("Shutting down server")
            self.shutdown()

    def add_flowsheet(self, id_, flowsheet, save_as):
        """Add a flowsheet, and also the method of saving it.
        """
        self._flowsheets[id_] = flowsheet
        store = persist.DataStore.create(save_as)
        self._dsm.add(id_, store)
        store.save(flowsheet)

    def save_flowsheet_data(self, id_, data):
        """Save the flowsheet data to the appropriate store.
        """
        self._dsm.save(id_, data)

    def get_flowsheet_obj(self, id_):
        """Get a flowsheet with the given ID.
        """
        return self._flowsheets[id_]

    def serialize_flowsheet(self, flowsheet, id_):
        try:
            result = self._fss.serialize(flowsheet, id_)
        except (AttributeError, KeyError) as err:
            raise ValueError(f"Error serializing flowsheet: {err}")
        return result


class FlowsheetServerHandler(http.server.SimpleHTTPRequestHandler):
    """Handle requests from the IDAES flowsheet visualization (IFV) web page.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def do_GET(self):
        """Get a model from the app.
        """
        u, id_ = self._parse_flowsheet_url(self.path)
        _log.debug(f"do_GET: path={self.path} id=={id_}")
        if u.path == "/app":
            if id_ is None:
                self.send_error(500, message="Query parameter 'id' is required for '/app'")
            else:
                self._get_app(u, id_)
        elif u.path == "/fs":
            if id_ is None:
                self.send_error(500, message="Query parameter 'id' is required for '/fs'")
            else:
                self._get_fs(u, id_)
        else:
            # Try to serve a file
            self.directory = _static_dir  # keep here: would be overwritten if set earlier
            super().do_GET()

    def _parse_flowsheet_url(self, path):
        u, id_ = urlparse(self.path), None
        if u.query:
            queries = dict([q.split("=") for q in u.query.split("&")])
            id_ = queries.get("id", None)
        return u, id_

    def do_PUT(self):
        """Receive an updated flowsheet from the web app.
        """
        u, id_ = self._parse_flowsheet_url(self.path)
        _log.debug(f"do_PUT: path={self.path} query={u.query}")
        if u.path == "/fs":
            if id_ is None:
                self.send_error(500, message="Query parameter 'id' is required for '/fs'")
            else:
                self._put_fs(u, id_)

    def _get_app(self, url, id_):
        p = Path(_template_dir / "index.html")
        with open(p, "r") as fp:
            s = fp.read()
            page = s.format(flowsheet_id=id_)
        self._write_html(200, page)

    def _get_fs(self, url, id_):
        try:
            flowsheet = self.server.get_flowsheet_obj(id_)
        except KeyError:
            self._bad_id_error(id_)
            return
        try:
            flowsheet_json = self.server.serialize_flowsheet(flowsheet, id_)
        except ValueError as err:
            self.send_error(500, message="Serialization error", explain=str(err))
        else:
            self._write_json(200, flowsheet_json)

    def _put_fs(self, url, id_):
        # read and parse flowsheet sent from application
        self.server.rfile.read()
        app_flowsheet = json.loads(utf8_decode(bytes))
        # save application flowsheet
        self.server.save_flowsheet_data(app_flowsheet)
        self._write_json(200, {"message": "Success"})

    def _no_id_error(self):
        self._write_json(404, {"message": "Identifier missing from request"})

    def _bad_id_error(self, id_):
        self._write_json(404, {"message": "No flowsheet found for identifier", "id": id_})

    def _write_json(self, code, data):
        str_json = json.dumps(data)
        value = utf8_encode(str_json)
        self.send_response(code)
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Content-type", "application/json")
        self.send_header("Content-length", str(len(value)))
        self.end_headers()
        self.wfile.write(value)

    def _write_html(self, code, page):
        value = utf8_encode(page)
        self.send_response(code)
        self.send_header("Content-type", "text/html")
        self.send_header("Content-length", str(len(value)))
        self.end_headers()
        self.wfile.write(value)


def utf8_encode(s: str):
    return s.encode(encoding="utf-8")

def utf8_decode(b: bytes):
    return b.decode(encoding="utf-8")

def find_free_port():
    import time
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    time.sleep(1)  # wait for socket cleanup!!!
    return port
