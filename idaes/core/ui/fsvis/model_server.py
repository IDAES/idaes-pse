#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Visualization server back-end.

The main class is `FlowsheetServer`, which is instantiated from the `visualize()` function.
"""

# stdlib
import http.server
import json
from pathlib import Path
import re
import socket
import threading
from typing import Dict, Union
from urllib.parse import urlparse

# package
from idaes import logger
from ..flowsheet import FlowsheetDiff, FlowsheetSerializer
from . import persist, errors

_log = logger.getLogger(__name__)

# Directories
_this_dir = Path(__file__).parent.absolute()
_static_dir = _this_dir / "static"
_template_dir = _this_dir / "templates"


class FlowsheetServer(http.server.HTTPServer):
    """A simple HTTP server that runs in its own thread.

    This server is used for *all* models for a given process, so every request needs to contain
    the ID of the model that should be used in that transaction.

    The only methods that the visualization function needs to call are the constructor, `start()` to
     start running the server, and `add_flowsheet()`, to a add a new flowsheet.
    """

    def __init__(self, port=None):
        """Create HTTP server"""
        self._port = port or find_free_port()
        _log.info(f"Starting HTTP server on localhost, port {self._port}")
        super().__init__(("127.0.0.1", self._port), FlowsheetServerHandler)
        self._dsm = persist.DataStoreManager()
        self._flowsheets = {}
        self._thr = None
        self._settings_block = {}

    @property
    def port(self):
        return self._port

    def start(self):
        """Start the server, which will spawn a thread."""
        self._thr = threading.Thread(target=self._run)
        self._thr.setDaemon(True)
        self._thr.start()

    def add_setting(self, key: str, value):
        """Add a setting to the flowsheet's settings block. Settings block is
        a dict that has general setting values related to the UI server. Such
        values could be retrieved to set some settings in the UI.

        An example setting value is the `save_model_time_interval` which sets
        the time interval at which the model checks if the graph has changed
        or not for the model to be saved.

        Args:
            key: Setting name
            value: Setting value
        """
        self._settings_block[key] = value

    def get_setting(self, key: str):
        """Get a setting value from the flowsheet's settings block.

        Args:
            key: Setting name

        Returns:
            Setting value. None if Setting name (key) doesn't exist
        """
        if key not in self._settings_block:
            _log.warning(f"key '{key}' is not set in the flowsheet settings block")
            return None
        return self._settings_block[key]

    def add_flowsheet(self, id_, flowsheet, store: persist.DataStore) -> str:
        """Add a flowsheet, and also the method of saving it.

        Args:
            id_: Name of flowsheet
            flowsheet: Flowsheet object
            store: A DataStore for saving the flowsheet

        Returns:
            Name of flowsheet, modified as necessary to be URL friendly

        Raises:
            ProcessingError: if the flowsheet can't be serialized
            DatastoreError: If the flowsheet can't be saved
        """
        # replace all but 'unreserved' (RFC 3896) chars with a dash; remove duplicate dashes
        id_ = self.canonical_flowsheet_name(id_)
        self._flowsheets[id_] = flowsheet
        _log.debug(f"Flowsheet '{id_}' storage is {store}")
        self._dsm.add(id_, store)
        # First try to update, so as not to overwrite saved value
        try:
            self.update_flowsheet(id_)
        except errors.FlowsheetNotFoundInDatastore:
            _log.debug(f"No existing flowsheet found in {store}: saving new value")
            # If not found in datastore, save new value
            fs_dict = FlowsheetSerializer(flowsheet, id_).as_dict()
            store.save(fs_dict)
        else:
            _log.debug(f"Existing flowsheet found in {store}: saving merged value")
        return id_

    @staticmethod
    def canonical_flowsheet_name(name: str) -> str:
        """Create a canonical flowsheet name from the name provided by the user.

        Replace all but 'unreserved' (RFC 3896) chars plus '~' with a dash and remove duplicate dashes.
        The result will not have whitespace, slashes, punctuation, or any special characters.

        Args:
            name: User-provided name

        Returns:
            New name
        """
        return re.sub(r"-+", "-", re.sub(r"[^a-zA-Z0-9-._]", "-", name))

    # === Public methods called only by HTTP handler ===

    def save_flowsheet(self, id_, flowsheet: Union[Dict, str]):
        """Save the flowsheet to the appropriate store.

        Raises:
            ProcessingError, if parsing of JSON failed (see :meth:`DataStoreManager.save()`)
        """
        try:
            self._dsm.save(id_, flowsheet)
        except errors.DatastoreError as err:
            raise errors.ProcessingError(f"While saving flowsheet: {err}")
        except KeyError as err:
            raise errors.ProcessingError(f"While saving flowsheet: {err}")

    def update_flowsheet(self, id_: str) -> Dict:
        """Update flowsheet.

        The returned flowsheet is also saved to the datastore.

        Args:
            id_: Identifier of flowsheet to update.

        Returns:
            Merged value of flowsheets in datastore and current value in memory

        Raises:
            FlowsheetUnknown if the flowsheet id is not known
            FlowsheetNotFound (subclass) if the flowsheet id is known, but it can't be retrieved
            ProcessingError for internal errors
        """
        # Get saved flowsheet from datastore
        try:
            saved = self._load_flowsheet(id_)
        except KeyError:
            raise errors.FlowsheetUnknown(id_)
        except ValueError:
            raise errors.FlowsheetNotFoundInDatastore(id_)
        # Get current value from memory
        try:
            obj = self._get_flowsheet_obj(id_)
        except KeyError:
            raise errors.FlowsheetNotFoundInMemory(id_)
        try:
            obj_dict = self._serialize_flowsheet(id_, obj)
        except ValueError as err:
            raise errors.ProcessingError(f"Cannot serialize flowsheet: {err}")
        # Compare saved and current value
        diff = FlowsheetDiff(saved, obj_dict)
        _log.debug(f"diff: {diff}")
        if not diff:
            # If no difference do nothing
            _log.debug("Stored flowsheet is the same as the flowsheet in memory")
            merged = saved
        else:
            # Otherwise, save this merged value before returning it
            num, pl = len(diff), "s" if len(diff) > 1 else ""
            _log.debug(f"Stored flowsheet and model in memory differ by {num} item{pl}")
            self.save_flowsheet(id_, diff.merged())
        # Return [a copy of the] merged value
        return diff.merged(do_copy=True)

    # === Internal methods ===

    def _load_flowsheet(self, id_) -> Union[Dict, str]:
        return self._dsm.load(id_)

    def _get_flowsheet_obj(self, id_):
        """Get a flowsheet with the given ID."""
        return self._flowsheets[id_]

    @staticmethod
    def _serialize_flowsheet(id_, flowsheet):
        try:
            result = FlowsheetSerializer(flowsheet, id_).as_dict()
        except (AttributeError, KeyError) as err:
            raise ValueError(f"Error serializing flowsheet: {err}")
        return result

    def _run(self):
        """Run in a separate thread."""
        _log.debug(f"Serve forever on localhost:{self._port}")
        try:
            self.serve_forever()
        except Exception as err:
            _log.info(f"Shutting down server due to error: {err}")
            self.shutdown()


class FlowsheetServerHandler(http.server.SimpleHTTPRequestHandler):
    """Handle requests from the IDAES flowsheet visualization (IFV) web page."""

    def __init__(self, *args, **kwargs):
        self.directory = (
            None  # silence warning about initialization outside constructor
        )
        super().__init__(*args, **kwargs)
        # Server should return text/javascript MIME type for served JS files (issue 259)
        self.extensions_map[".js"] = "text/javascript"

    # === GET ===

    def do_GET(self):
        """Process a request to receive data.

        Routes:
          * `/app`: Return the web page
          * `/fs`: Retrieve an updated flowsheet.
          * `/setting`: Retrieve a setting value.
          * `/path/to/file`: Retrieve file stored static directory
        """
        u, queries = self._parse_flowsheet_url(self.path)
        id_ = queries.get("id", None) if queries else None

        _log.debug(f"do_GET: path={self.path} id=={id_}")
        if u.path in ("/app", "/fs") and id_ is None:
            self.send_error(
                400, message=f"Query parameter 'id' is required for '{u.path}'"
            )
            return

        if u.path == "/app":
            self._get_app(id_)
        elif u.path == "/fs":
            self._get_fs(id_)
        elif u.path == "/setting":
            setting_key_ = queries.get("setting_key", None) if queries else None
            if setting_key_ is None:
                self.send_error(
                    400,
                    message=f"Query parameter 'setting_key' is required for '{u.path}'",
                )
                return
            self._get_setting(setting_key_)
        else:
            # Try to serve a file
            self.directory = _static_dir  # keep here: overwritten if set earlier
            super().do_GET()

    def _get_app(self, id_):
        """Read index file, process to insert flowsheet identifier, and return it."""
        p = Path(_template_dir / "index.html")
        with open(p, "r", encoding="utf-8") as fp:
            s = fp.read()
            page = s.format(flowsheet_id=id_)
        self._write_html(200, page)

    def _get_fs(self, id_: str):
        """Get updated flowsheet.

        Args:
            id_: Flowsheet identifier

        Returns:
            None
        """
        try:
            merged = self.server.update_flowsheet(id_)
        except errors.FlowsheetUnknown as err:
            # User error: user asked for a flowsheet by an unknown ID
            self.send_error(404, message=str(err))
            return
        except (errors.FlowsheetNotFound, errors.ProcessingError) as err:
            # Internal error: flowsheet ID is found, but other things are missing
            self.send_error(500, message=str(err))
            return
        # Return merged flowsheet
        self._write_json(200, merged)

    def _get_setting(self, setting_key_: str):
        """Get setting value.

        Args:
            id_: Flowsheet identifier
            setting_key_: Setting name (key)

        Returns:
            Setting value
        """
        self._write_json(200, {"setting_value": self.server.get_setting(setting_key_)})

    # === PUT ===

    def do_PUT(self):
        """Process a request to store data."""
        u, queries = self._parse_flowsheet_url(self.path)
        id_ = queries.get("id", None) if queries else None
        _log.info(f"do_PUT: route={u} id={id_}")
        if u.path in ("/fs",) and id_ is None:
            self._write_text(
                400, message=f"Query parameter 'id' is required for '{u.path}'"
            )
            return
        if u.path == "/fs":
            self._put_fs(id_)

    def _put_fs(self, id_):
        # read  flowsheet from request (read(LENGTH) is required to avoid hanging)
        read_len = int(self.headers.get("Content-Length", "-1"))
        data = utf8_decode(self.rfile.read(read_len))
        # save flowsheet
        try:
            self.server.save_flowsheet(id_, data)
        except errors.ProcessingError as err:
            self._write_text(400, message=str(err))
            return
        except Exception as err:
            self._write_text(500, message=str(err))
            return
        self._write_text(200, message="success")

    # === Internal methods ===

    def _write_text(self, code, message: str):
        value = utf8_encode(message)
        self.send_response(code)
        self.send_header("Content-type", "application/text")
        self.send_header("Content-length", str(len(value)))
        self.end_headers()
        self.wfile.write(value)

    def _write_json(self, code, data):
        str_json = json.dumps(data)
        value = utf8_encode(str_json)
        self.send_response(code)
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
        u, queries = urlparse(path), None
        if u.query:
            queries = dict([q.split("=") for q in u.query.split("&")])
        return u, queries

    # === Logging ===

    def log_message(self, fmt, *args):
        """Override to send messages to our module logger instead of stderr"""
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
