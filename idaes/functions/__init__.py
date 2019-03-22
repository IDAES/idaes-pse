def functions_lib():
    import os
    plib = os.path.join(os.path.dirname(__file__), "functions.so")
    return plib

def functions_available():
    import os
    return os.path.isfile(functions_lib())
