def functions_lib():
    import os
    import idaes
    plib = os.path.join(idaes.lib_directory, "functions.so")
    return plib

def functions_available():
    import os
    return os.path.isfile(functions_lib())
