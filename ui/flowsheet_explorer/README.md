# Flowsheet Explorer

## Initial Setup

1. `conda activate <idaes env name>`


###### Skip this step if you've already installed flask

2. `pip install flask`

3. `cd` into the directory that has `app.py`

4. `mkdir static/raw`

5. Copy all of your raw images to the `raw` directory

6. `mkdir static/modeled`

7. Copy all of your modeled images to the `modeled` directory


## Running the flowsheet explorer

1. `cd` into the directory that has `app.py`

2. `conda activate <idaes env name>`

3. `flask run`

## Navigating

On your browser, go to the URL that is output by the flask app. (Generally `http://127.0.0.1:5000/`)

Use the slider and arrow buttons to navigate through the images. 

The right image is the raw image and the left is the modeled image.