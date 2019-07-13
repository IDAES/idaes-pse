from notebook.utils import url_path_join
from notebook.base.handlers import IPythonHandler
import asyncio

class IDAESHandler(IPythonHandler):
    async def get(self):
        with open('icons/mixer.svg', 'r') as image_file:
            self.write(image_file.read())

def _jupyter_server_extension_paths():
    """
    Set up the server extension for collecting metrics
    """
    return [{
        'module': 'IDAES',
    }]

def load_jupyter_server_extension(nbapp):
    """
    Called during notebook start
    """
    route_pattern = url_path_join(nbapp.web_app.settings['base_url'], '/icons/mixer.svg')
    nbapp.web_app.add_handlers('.*', [(route_pattern, IDAESHandler)])
