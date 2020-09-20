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


# See model_server.py for more information about this module's role in the visualizer


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
import threading
import socket

# third party
from flask import Flask, jsonify, request, render_template#, send_static_file
from flask_cors import CORS, cross_origin
from werkzeug.serving import make_server

# local
from idaes.ui.fsvis.server import DataStorage

# globals
app = Flask(__name__, static_url_path='', 
        static_folder='static', 
        template_folder='templates')
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'
db = DataStorage()


# Custom exceptions (handled by Flask error handlers)
class LookupError(Exception):
    def __init__(self, id_):
        self.id_ = id_
        super().__init__()


class NoDataError(Exception):
    pass


class NoIdError(Exception):
    pass


class NoModelServerUrlError(Exception):
    pass


# classes/functions
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


class App:
    """Singleton to run flask server as a forked process.
    """
    instance = None

    class __App:
        def __init__(self, **kwargs):
            self._server, self.host, self.port = None, None, None
            self.start(**kwargs)

        def is_running(self) -> bool:
            return self._server and self._server.is_alive()

        def start(self, host: str = '127.0.0.1'):
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


def find_free_port(host):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind((host, 0))
    return s.getsockname()[1]


def save(data, id_):
    path = os.path.expandvars(os.path.join(os.path.expanduser("~"), ".idaes", "viz"))
    if not os.path.exists(path):
        os.makedirs(path)
    file_path = os.path.join(path, f"{id_}.viz")
    with open(file_path, "w") as viz_file:
        json.dump(data, viz_file)


def update(request, id_):
    data = request.get_json()
    if data is None:
        raise NoDataError()
    db.update(id_, data)
    save_header = request.headers.get("Source", None)
    if save_header == "save_button":
        save(data, id_)
    return db.fetch(id_)


def fetch(request, id_):
    data = db.fetch(id_)
    if data is None:
        raise LookupError(id_)
    return data


# Flask routes
@app.route("/fs", methods=["GET", "POST"])
@cross_origin()
def flowsheet():
    """Store/retrieve flowsheet JSON
    """
    id_ = request.args.get("id", None)
    if id_ is None:
        raise NoIdError()
    if request.method == "POST":
        data = update(request, id_)
        response = jsonify(data)
        response.headers.add("Access-Control-Allow-Origin", "*")
        return response

    elif request.method == "GET":
        data = fetch(request, id_)
        response = jsonify(data)
        response.headers.add("Access-Control-Allow-Origin", "*")
        return response


@app.route("/app", methods=["GET", "POST"])
@cross_origin()
def display():
    """Display the web app."""
    id_ = request.args.get("id", None)
    model_server_url_ = request.args.get("modelurl", None)
    if id_ is None:
        raise NoIdError()
    if model_server_url_ is None:
        raise NoModelServerUrlError()
    if request.method == "POST":
        data = update(request, id_)
    elif request.method == "GET":
        data = fetch(request, id_)
    return render_template('index.html', model=data, modelurl=model_server_url_)


# Flask error handlers
@app.errorhandler(LookupError)
def lookup_exception_handler(error):
    return f"Flowsheet id='{error.id_}' not found", 404


@app.errorhandler(NoDataError)
def nodata_exception_handler(error):
    return "No flowsheet data", 500


@app.errorhandler(NoIdError)
def noid_exception_handler(error):
    return "Flowsheet id is missing from request", 500
