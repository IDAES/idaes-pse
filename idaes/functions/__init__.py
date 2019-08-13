from pyomo.common.fileutils import this_file_dir as _this_file_dir

def functions_lib():
    import os
    plib = os.path.join(_this_file_dir(), "functions.so")
    return plib

def functions_available():
    import os
    return os.path.isfile(functions_lib())
