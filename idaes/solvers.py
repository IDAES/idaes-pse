import os
import zipfile
import idaes
from shutil import copyfile
import urllib.request as ur
from pyomo.common.download import FileDownloader

def download_binaries(url=None):
    """
    Download IDAES solvers and libraries and put them in the right location. Need
    to supply either local or url argument.

    Args:
        url (str): a url to download binary files to install files

    Returns:
        None
    """
    idaes._create_lib_dir()
    idaes._create_bin_dir()
    solvers_zip = os.path.join(idaes.bin_directory, "idaes-solvers.zip")
    libs_zip = os.path.join(idaes.lib_directory, "idaes-lib.zip")
    fd = FileDownloader()
    arch = fd.get_sysinfo()
    if url is not None:
        if not url.endswith("/"):
            c = "/"
        else:
            c = ""
        solvers_from = c.join([url, "idaes-solvers-{}-{}.zip".format(arch[0], arch[1])])
        libs_from = c.join([url, "idaes-lib-{}-{}.zip".format(arch[0], arch[1])])
        fd.set_destination_filename(solvers_zip)
        fd.get_binary_file(solvers_from)
        fd.set_destination_filename(libs_zip)
        fd.get_binary_file(libs_from)
    else:
        raise Exception("Must provide a location to download binaries")

    with zipfile.ZipFile(solvers_zip, 'r') as f:
        f.extractall(idaes.bin_directory)
    with zipfile.ZipFile(libs_zip, 'r') as f:
        f.extractall(idaes.lib_directory)
