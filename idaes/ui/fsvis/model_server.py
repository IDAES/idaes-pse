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

# stdlib
import http.server
import json
from pathlib import Path
import socket
import threading
from typing import Dict, Union
from urllib.parse import urlparse

# package
from idaes import logger
from idaes.ui.flowsheet_comparer import compare_models, model_jointjs_conversion
from . import persist
from ..flowsheet_serializer import FlowsheetSerializer

_log = logger.getLogger(__name__)

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
        _log.info(f"Starting HTTP server on localhost, port {self._port}")
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
        _log.debug(f"Serve forever on localhost:{self._port}")
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
        _log.debug(f"Flowsheet '{id_}' storage is {store}")
        self._dsm.add(id_, store)
        # Serialize as JSON string
        fs_dict = FlowsheetSerializer().serialize(flowsheet, id_)
        store.save(fs_dict)

    def save_flowsheet(self, id_, flowsheet: Union[Dict, str]):
        """Save the flowsheet to the appropriate store.
        """
        self._dsm.save(id_, flowsheet)

    def load_flowsheet(self, id_) -> Union[Dict, str]:
        return self._dsm.load(id_)

    def get_flowsheet_obj(self, id_):
        """Get a flowsheet with the given ID.
        """
        return self._flowsheets[id_]

    def serialize_flowsheet(self, flowsheet, id_):
        try:
            result = FlowsheetSerializer().serialize(flowsheet, id_)
        except (AttributeError, KeyError) as err:
            raise ValueError(f"Error serializing flowsheet: {err}")
        return result


class FlowsheetServerHandler(http.server.SimpleHTTPRequestHandler):
    """Handle requests from the IDAES flowsheet visualization (IFV) web page.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    # === GET ===

    def do_GET(self):
        """Process a request to receive data.

        Routes:
          * `/app`: Return the web page
          * `/fs`: Retrieve an updated flowsheet.
          * `/path/to/file`: Retrieve file stored static directory
        """
        u, id_ = self._parse_flowsheet_url(self.path)
        _log.debug(f"do_GET: path={self.path} id=={id_}")
        if u.path in ("/app", "/fs") and id_ is None:
            self.send_error(
                400, message=f"Query parameter 'id' is required for '{u.path}'"
            )
        if u.path == "/app":
            self._get_app(id_)
        elif u.path == "/fs":
            self._get_fs(id_)
        else:
            # Try to serve a file
            self.directory = _static_dir  # keep here: overwritten if set earlier
            super().do_GET()

    def _get_app(self, id_):
        """Read index file, process to insert flowsheet identifier, and return it.
        """
        p = Path(_template_dir / "index.html")
        with open(p, "r") as fp:
            s = fp.read()
            page = s.format(flowsheet_id=id_)
        self._write_html(200, page)

    def _get_fs(self, id_: str):
        """Get updated flowsheet.

        Algorithm:
            1. get saved flowsheet (if any) from datastore
            2. get flowsheet object from memory
            3. if both flowsheets exist, merge them (otherwise use the memory flowsheet)
            4. save merged flowsheet to the datastore (this way datastore mirrors what app shows)
            5. return merged flowsheet to app

        Args:
            id_: Flowsheet identifier

        Returns:
            None
        """
        # Get flowsheet from datastore
        fs_store_data, fs_store_dict = None, None
        try:
            fs_store_data = self.server.load_flowsheet(id_)
            _log.debug("Found stored flowsheet")
        except KeyError:
            _log.debug("No stored flowsheet found, using flowsheet from memory")
            pass  # this is ok; could happen initially
        if fs_store_data is not None:
            if isinstance(fs_store_data, dict):
                fs_store_dict = fs_store_data
            else:
                try:
                    fs_store_dict = json.loads(fs_store_data)
                except Exception as err:
                    self.send_error(
                        500, message="Cannot parse stored JSON data", explain=str(err)
                    )
                    return
        # Get flowsheet from memory
        try:
            fs_obj = self.server.get_flowsheet_obj(id_)
        except KeyError:
            self.send_error(404, message=f"No flowsheet found in memory for id={id_}")
            return
        fs_obj_dict = FlowsheetSerializer().serialize(fs_obj, id_)
        # Merge datastore and memory flowsheets
        if fs_store_dict is None:
            # no stored value, so use value from memory
            fs_merged = fs_obj_dict
        else:
            # update the 'cells' part to reflect the new 'model' part
            model_diff, _ = compare_models(fs_store_dict, fs_obj_dict)
            if not model_diff:  # model parts are the same
                _log.debug("Stored flowsheet and model in memory are the same")
                # Still need to update the model with the values from the object, since changes
                # in the *values* inside the streams and units may have occurred
                fs_merged = {
                    "model": fs_obj_dict["model"],
                    "cells": fs_store_dict["cells"],
                }
            else:
                if _log.isEnabledFor(logger.DEBUG):
                    num_diffs, plural = (
                        len(model_diff),
                        "s" if len(model_diff) > 1 else "",
                    )
                    _log.debug(
                        f"Stored flowsheet and model in memory differ by {num_diffs} change{plural}; "
                        f"updating JointJS display data"
                    )
                # modify the JointJS display info in 'cells' based on the diff
                merged_cells = model_jointjs_conversion(model_diff, fs_store_dict)
                # the merged value has the new model info plus the updated 'cells'
                fs_merged = {
                    "model": fs_obj_dict["model"],
                    "cells": merged_cells["cells"],
                }
        assert fs_merged
        # Save merged flowsheet
        if fs_merged is fs_store_dict:
            _log.debug("Flowsheet has not changed so skipping save of merged value")
        else:
            _log.debug("Storing merged flowsheet")
            self.server.save_flowsheet(id_, fs_merged)
        # Return merged flowsheet to app
        self._write_json(200, fs_merged)

    # === PUT ===

    def do_PUT(self):
        """Process a request to store data.
        """
        u, id_ = self._parse_flowsheet_url(self.path)
        _log.debug(f"do_PUT: route={u} id={id_}")
        if u.path in ("/fs",) and id_ is None:
            self.send_error(
                400, message=f"Query parameter 'id' is required for '{u.path}'"
            )
        if u.path == "/fs":
            self._put_fs(id_)

    def _put_fs(self, id_):
        # read  flowsheet from request (read(LENGTH) is required to avoid hanging)
        read_len = int(self.headers.get("Content-Length", "-1"))
        data = utf8_decode(self.rfile.read(read_len))
        # save flowsheet
        self.server.save_flowsheet(id_, data)
        self.send_response(200, message="success")

    # === Internal methods ===

    def _write_json(self, code, data):
        str_json = json.dumps(data)
        value = utf8_encode(str_json)
        self.send_response(code)
        # self.send_header("Access-Control-Allow-Headers", "Content-Type")
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

    def _parse_flowsheet_url(self, path):
        u, id_ = urlparse(self.path), None
        if u.query:
            queries = dict([q.split("=") for q in u.query.split("&")])
            id_ = queries.get("id", None)
        return u, id_

    # === Logging ===

    def log_message(self, fmt, *args):
        """Override to send messages to our module logger instead of stderr
        """
        msg = "%s - - [%s] %s" % (
            self.address_string(),
            self.log_date_time_string(),
            fmt % args,
        )
        _log.debug(msg)


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
