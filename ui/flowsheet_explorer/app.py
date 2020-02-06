from glob import glob
from flask import Flask, render_template, request
from werkzeug.contrib.cache import SimpleCache

import json
import os

cache = SimpleCache()

app = Flask(__name__)
app.debug = True
app.config['TEMPLATES_AUTO_RELOAD'] = True


@app.route("/", methods=['GET', 'POST'])
def index():
    cached_files = cache.get("modeled_files")
    if not cached_files:
        path = os.path.join(os.getcwd(), "static/modeled/")
        cached_files = [os.path.basename(name) for name in sorted(glob(path + "*-*-*_*.svg"))]
        cache.set("modeled_files", cached_files, timeout=300)

    file_count = len(cached_files)
    cached_files_dict = dict(enumerate(cached_files))
    slider_val = 0
    if request.method == 'POST':
        # POST means that there is input from the slider or the forward/back buttons and we need to retrieve the value from date_slider
        # GET means that someone is visiting the page but they haven't interacted with the slider or buttons
        slider_val = request.form.get('date_slider')
        if not slider_val:
            slider_val = 0
        else:
            slider_val = int(slider_val)

    if file_count != 0:
        date = os.path.splitext(cached_files[slider_val])[0]
        day, hour = date.split("_")
        date = day + " " + hour + ":00"
        return render_template("index.html", date=date, slider_value=slider_val, range=file_count, raw_img="/static/raw/" + cached_files[slider_val], modeled_img="/static/modeled/" + cached_files[slider_val], slider_labels=json.dumps(cached_files_dict))


if __name__ == '__main__':
    app.run()
