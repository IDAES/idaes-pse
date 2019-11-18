import os
import idaes
from shutil import copyfile
import urllib.request as ur

def download_binaries(local=None, url=None):
    """
    Download IDAES solvers and libraries and put them in the right location. Need
    to supply either local or url argument.

    Args:
        local (str): a local directory containting binary files to install
        url (str): a url to download binary files to install files

    Returns:
        None
    """
    idaes._create_lib_dir()
    idaes._create_bin_dir()
    solvers_zip = os.path.join(idaes.bin_directory, "idaes-solvers.zip")
    libs_zip = os.path.join(idaes.lib_directory, "idaes-lib.zip")
    if url is not None:
        if not url.endswith("/"):
            c = "/"
        else:
            c = ""
        solvers_from = c.join([url, "idaes-solvers.zip"])
        libs_from = c.join([url, "idaes-lib.zip"])
        ur.urlretrieve(solvers_from, solvers_zip)
        ur.urlretrieve(libs_from, libs_zip)
    elif local is not None:
        solvers_from = os.path.join(local, "idaes-solvers.zip")
        libs_from =  os.path.join(local, "idaes-lib.zip")
        copyfile(solvers_from, solvers_zip)
        copyfile(libs_from, libs_zip)
    else:
        raise Exception("Must provide a location to download binaries")

    with zipfile.ZipFile(solvers_zip, 'r') as f:
        f.extractall(idaes.bin_directory)
    with zipfile.ZipFile(libs_zip, 'r') as f:
        f.extractall(idaes.lib_directory)
