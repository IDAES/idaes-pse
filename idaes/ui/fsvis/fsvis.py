# server backend code for fsvis

# from notebook end, call m.visualize, which serializes m and launches the app


#
def get_stored_fs(fs_id):
    return

# compare serialized model with stored to see if changes need to be made at all
def _model_has_changed(fs_id, new):
    return True

# record structural changes from the flowsheet
def update_stored_model(fs_id, new):
    return

# record manual layout changes from jointjs
def update_layout(fs_id, new):
    return

# sort the elements of the json representation of the flowsheet
# to facilitate comparison between versions
def _canonicalize_jointjs_output(json):
    return
