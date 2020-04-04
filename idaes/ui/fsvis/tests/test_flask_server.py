# #!/usr/bin/env python
# # coding: utf-8
# """

# """
# import time
# import multiprocessing
# import pytest
# import requests

# from idaes.ui.fsvis.flask_server import App


# @pytest.fixture()
# def app():
#     a = App()
#     time.sleep(2)
#     yield a
#     App().stop()


# def test_app(app):
#     App()  # should not run twice
#     j = {"test": 123}  # payload
#     jid = b"abcdef"  # id: use un-encoded bytes
#     url = f"http://{app.host}:{app.port}/fs"  # server url
#     # put the JSON payload into the server
#     resp = requests.post(url, json=j, params={'id': jid})
#     assert resp.status_code == 200
#     assert resp.content == jid
#     # get back the payload by its ID
#     resp = requests.get(url, params={'id': jid})
#     assert resp.json() == j
#     # try an invalid ID
#     resp = requests.get(url, params={'id': 'foo'})
#     assert resp.status_code == 404


# def test_app_stop(app):
#     # stop the app
#     app.stop()
#     # setup
#     j = {"test": 123}  # payload
#     jid = b"abcdef"  # id: use un-encoded bytes
#     url = f"http://{app.host}:{app.port}/fs"  # server url
#     # should fail (server was stopped)
#     pytest.raises((requests.ConnectionError, requests.exceptions.ReadTimeout),
#                   requests.post, url, json=j, params={'id': jid}, timeout=1)
