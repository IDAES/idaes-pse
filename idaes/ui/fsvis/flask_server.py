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
import socket
# third party
from flask import Flask, request, render_template#, send_static_file
# local
from idaes.ui.fsvis.server import DataStorage

# globals

app = Flask(__name__, static_url_path='', 
        static_folder='draft/static', 
        template_folder='draft/templates')
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


# classes/functions

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

        def start(self, host: str = '127.0.0.1', port: int = 5555, probe: int = 10):
            """Start server.

            Does nothing if server is already running.
            """
            if self._server:
                return
            self._server = "starting"
            pport = self._probe_ports(host, port, probe)
            if pport == 0:
                p1, p2 = port, port + probe - 1
                raise RuntimeError(f"Could not find an open port in range {p1}..{p2}")
            self._server = Process(target=app.run, args=(host, pport))
            self._server.start()
            self.host, self.port = host, pport

        @staticmethod
        def _probe_ports(host, start, num):
            found_port = 0  # invalid
            for port in range(start, start + num):
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                try:
                    s.bind((host, port))
                    found_port = port  # it worked
                except socket.error as e:
                    pass
                finally:
                    s.close()
                if found_port > 0:
                    break
            return found_port

        def stop(self) -> None:
            """Stop server.
            """
            if self.is_running():
                self._server.terminate()
                self._server.join()

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


# Flask routes

@app.route("/fs", methods=["GET", "POST"])
def flowsheet():
    """Store/retrieve flowsheet JSON
    """
    id_ = request.args.get("id", None)
    if id_ is None:
        raise NoIdError()
    if request.method == "POST":
        data = request.get_json()
        if data is None:
            raise NoDataError()
        db.save(id_, data)
        return id_
    elif request.method == "GET":
        data = db.fetch(id_)
        if data is None:
            raise LookupError(id_)
        return data

#@app.route("/test", methods=["GET"])
#def testpage():
#    """Quick hello world"""
#    import testpages
#    page = testpages.helloworld
#    return page

@app.route("/app", methods=["GET"])
def display():
    """Display the web app."""
    return render_template('app.html')

#@app.route("/draft")
#def send_static_draft():
#    return app.send_static_file('/draft.html')


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
