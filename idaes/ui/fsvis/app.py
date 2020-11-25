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
"""
Flask app for flowsheet viewer server.

This should not be invoked directly by the user. Inside the code to
visualize a flowsheet, the app will be started like:

    try:
        self._app = App(port=5678)
    except RuntimeError as err:
        # could not start server..

The server defines one endpoint, "/fs", to which you can post or from which you
can get a JSON blob. For either operation, supply a unique identifier for the blob
with the "?id=<value>" query parameter. So, for example, getting JSON blob "abc123"
would look like: "http://127.0.0.1:5555/fs?id=abc123", and POST would use the same
URL but provide JSON content. See the tests for programmatic examples using the
`requests` library.
"""
# standard library
from multiprocessing import Process
import json
import os
import requests
import threading
import socket

# third party
from flask import Flask, jsonify, request, render_template  # , send_static_file
from flask_cors import CORS, cross_origin
from werkzeug.serving import make_server

# local
from idaes.ui.fsvis.persist import DataStoreManager

# globals
app = Flask(
    __name__, static_url_path="", static_folder="static", template_folder="templates"
)
cors = CORS(app)
app.config["CORS_HEADERS"] = "Content-Type"

# -----------------
# Exception classes
# -----------------


class IdLookupError(Exception):
    def __init__(self, id_):
        self.id_ = id_
        super().__init__()


class NoDataError(Exception):
    pass


class NoIdError(Exception):
    pass


class NoModelServerUrlError(Exception):
    pass


# -----------------


class ServerThread(threading.Thread):
    def __init__(self, app, host, pport):
        threading.Thread.__init__(self)
        self.srv = make_server(host, pport, app)
        self.ctx = app.app_context()
        self.ctx.push()

    def run(self):
        self.srv.serve_forever()

    def shutdown(self):
        self.srv.shutdown()


def find_free_port(host):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind((host, 0))
    return s.getsockname()[1]


class ModelServerManager:

    def __init__(self):
        self.data = {}

    def add(self, name, info):
        self.data[name] = info

    def get(self, name):
        return self.data[name]


class App:
    """Singleton to run flask server as a forked process.
    """
    instance = None

    class __App:
        def __init__(self, **kwargs):
            self._server, self.host, self.port = None, None, None
            self.data_storage_manager = DataStoreManager()
            self.model_server_manager = ModelServerManager()
            self.start(**kwargs)

        def is_running(self) -> bool:
            return self._server and self._server.is_alive()

        def start(self, host: str = "127.0.0.1"):
            """Start server.

            Does nothing if server is already running.
            """
            if self._server:
                return
            port = find_free_port(host)
            self._server = ServerThread(app, host, port)
            self._server.start()
            self.host, self.port = host, port

        def stop(self) -> None:
            """Stop server.
            """
            if self.is_running():
                self._server.shutdown()

    def __init__(self, **kwargs):
        if App.instance is None:
            try:
                App.instance = App.__App(**kwargs)
            except RuntimeError as err:
                raise RuntimeError(f"Could not start Flask server: {err}")

    def __getattr__(self, name):
        """Delegate all access to instance.
        """
        return getattr(self.instance, name)


# ------------
# Flask routes
# ------------


@app.route("/fs", methods=["GET", "POST"])
@cross_origin()
def flowsheet():
    """Store/retrieve flowsheet JSON
    """
    db = App().data_storage_manager
    try:
        id_ = request.args["id"]
    except KeyError:
        raise NoIdError()
    if request.method == "POST":
        # TODO: 1) Get model from model_server
        # TODO: 2) Merge model with layout in 'data_in'
        # TODO: 3) Save merged data to the db
        # TODO: 4) Return merged data to caller
        data_in = request.get_json()
        data_out = db.update(id_, data_in)
        response = jsonify(data_out)
        response.headers.add("Access-Control-Allow-Origin", "*")
        return response

    elif request.method == "GET":
        try:
            id_ = request.args["id"]
        except KeyError:
            raise NoIdError()
        try:
            host, port = App().model_server_manager.get(id_)
        except KeyError:
            raise IdLookupError(id_)
        model = get_model(host, port)
        response = jsonify(model)
        #response.headers.add("Access-Control-Allow-Origin", "*")
        return response


def get_model(host, port):
    r = requests.get(f"http://{host}:{port}")
    return r.text


@app.route("/app", methods=["GET", "POST"])
@cross_origin()
def display():
    """Display the web app.
    """
    db = App().data_storage_manager
    try:
        id_ = request.args["id"]
    except KeyError:
        raise NoIdError()
    # try:
    #     model_server_url_ = request.args["modelurl"]
    # except KeyError:
    #     raise NoModelServerUrlError()
    # if request.method == "POST":
    #     data_in = request.get_json()
    #     #data = db.update(id_, data_in)
    if request.method == "GET":
        try:
            host, port = App().model_server_manager.get(id_)
        except KeyError:
            raise IdLookupError(id_)
        model = get_model(host, port)
#        data = db.load(id_)
    else:
        raise RuntimeError("Unexpected HTTP method")
    return render_template("index.html", model=model)


# --------------------
# Flask error handlers
# --------------------


@app.errorhandler(IdLookupError)
def lookup_exception_handler(error):
    return f"Flowsheet id='{error.id_}' not found", 404


@app.errorhandler(NoDataError)
def nodata_exception_handler(error):
    return "No flowsheet data", 500


@app.errorhandler(NoIdError)
def noid_exception_handler(error):
    return "Flowsheet id is missing from request", 500
