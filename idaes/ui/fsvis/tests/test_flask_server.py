#!/usr/bin/env python
# coding: utf-8
"""

"""
import pytest
import requests

from idaes.ui.fsvis.flask_server import App


@pytest.fixture()
def app():
    yield App()
    App().stop()


def test_app(app):
    App()  # should not run twice
    j = {"test": 123}
    jid = b"abcdef"
    h, p = App().host, App().port
    url = f"http://{h}:{p}/fs"
    resp = requests.post(url, json=j, params={'id': jid})
    assert resp.status_code == 200
    assert resp.content == jid
    resp = requests.get(url, params={'id': jid})
    assert resp.json() == j
    resp = requests.get(url, params={'id': 'foo'})
    assert resp.status_code == 404


def test_app_stop(app):
    app.stop()
    j = {"test": 123}
    jid = b"abcdef"
    h, p = app.host, app.port
    url = f"http://{h}:{p}/fs"
    # should fail (server stopped)
    pytest.raises((requests.ConnectionError, requests.exceptions.ReadTimeout),
                  requests.post, url, json=j, params={'id': jid}, timeout=1)



