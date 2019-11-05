from flask import Flask, render_template, request

import ast
import json
import os

app = Flask(__name__)
app.debug = True
app.config['TEMPLATES_AUTO_RELOAD'] = True


@app.route("/", methods=['GET', 'POST'])
def index():
    path = os.path.join(os.getcwd(), "static/modeled")
    files = [name for name in os.listdir(path) if os.path.isfile(os.path.join(path, name))]
    file_count = len(files)
    files_dict = dict(enumerate(files))
    slider_val = int(file_count / 2)
    if request.method == 'POST':
        print("Got here")
        slider_val = request.form.get('date_slider')
        if not slider_val:
            slider_val = int(file_count / 2)
        else:
            slider_val = int(slider_val)

    print(slider_val)
    date = os.path.splitext(files[slider_val])[0]
    day, hour = date.split("_")
    date = day + " " + hour + ":00"
    if file_count != 0:
        return render_template("index.html", date=date, slider_value=slider_val, range=file_count, raw_img="/static/raw/" + files[slider_val], modeled_img="/static/modeled/" + files[slider_val], slider_labels=json.dumps(files_dict))


if __name__ == '__main__':
    app.run()



    # return render_template("index.html", start_time=start_time, end_time=end_time)
