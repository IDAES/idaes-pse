#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
import re


def almpywriter(data, debug):
    """
    This function writes the file
    - <almname>alm.py
    y=<fname>.f(X)
    handle the multiple output
    """
    if data["opts"]["noutputs"] > 1 or debug["expandoutput"]:
        for output_name, mod_res in data["results"]["model"].items():
            almpywriter_help(data, mod_res, output_name)
    else:
        almpywriter_help(data, data["results"]["model"], data["labs"]["savezlabels"][0])


def almpywriter_help(data, mod_res, output_name):
    """
    This function writes the file
    - <almname>alm.py
    y=<fname>.f(X)
    preliminary formatting to get the model ready to write
    """
    model = mod_res.split("=")[1]
    model = model + " "
    tlist = ("sin", "cos", "log", "exp", "ln")
    for tok in tlist:
        if tok == "ln":
            model = model.replace(tok, "np.log")
        else:
            model = model.replace(tok, "np." + tok)
    model = model.replace("^", "**")

    with open(output_name + ".py", "w") as r:
        r.write("import numpy as np\n")
        r.write("def f(*X):\n")
        i = 0
        for lab in data["labs"]["savexlabels"]:
            r.write("    " + lab + "= X[" + str(i) + "]\n")
            i = i + 1
        r.write("    return " + model + "\n")


def almcvwriter(data):
    """
    This function writes the file
    - <almname>alm.py
    y=<fname>.f(X)
    handle the multiple output
    """
    if data["opts"]["noutputs"] > 1:
        for output_name, mod_res in data["results"]["model"].items():
            almcvwriter_help(data, mod_res, output_name)
    else:
        almcvwriter_help(data, data["results"]["model"], data["labs"]["savezlabels"][0])


def almcvwriter_help(data, mod_res, output_name):
    """
    This function writes the file
    - <almname>cv.py
    y=<fname>.f(X,params)
    """
    model = mod_res.split("=")[1]
    model = model + " "
    tlist = ("sin", "cos", "log", "exp")
    for tok in tlist:
        model = model.replace(tok, "np." + tok)
    model = model.replace("^", "**")

    # remove numerical values of coefficients
    model = re.sub(r"(?:.\w)?([0-9]+\.[0-9]+)", "b", model)
    model = re.sub(r"[E]{1}.\d{3}", "", model)
    model = model.replace("-", "+")

    tind = 0
    tstr = ""
    while "b" in model:
        model = model.replace("b", "B[" + str(tind) + "]", 1)
        tstr = tstr + ",B" + str(tind)
        tind = tind + 1
    with open(data["stropts"]["almname"].split(".")[0] + "cv.py", "w") as r:
        line = "def f(X, B):\n"
        r.write(line)
        r.write("    import numpy as np\n")
        i = 0
        for lab in data["labs"]["savexlabels"]:
            r.write("    " + lab + "= X[" + str(i) + "]\n")
            i = i + 1
        r.write("    return " + model + "\n")


def wrapwriter(sim):
    """
    This subroutine writes a temporary python file that is used
    to wrap python-based function that lack ALAMO's input/output.txt I/O
    """
    import os
    import inspect

    name = "simwrapper"
    name = name + ".py"
    with open(name, "w") as r:
        r.write("#!/usr/bin/python\n")
        r.write("def main():\n")
        r.write("    import " + inspect.getmodule(sim).__name__ + "\n")
        r.write("    infile = 'input.txt'\n")
        r.write("    outfile = 'output.txt'\n")
        r.write("    fin = open(infile, 'r')\n")
        r.write("    fout = open(outfile, 'w')\n")
        r.write("    newline = fin.readline()\n")
        r.write("    newlist = newline.split()\n")
        r.write("    n = int(newlist[0])\n")
        r.write("    for p in range(0, n):\n")
        r.write("        newline = fin.readline()\n")
        r.write("        newlist = newline.split()\n")
        r.write("        ninputs = len(newlist)\n")
        r.write("        x = [0] * (ninputs + 1)\n")
        r.write("        for k in range(0, ninputs):\n")
        r.write("            x[k] = float(newlist[k])\n")
        # r.write("        x[ninputs] = "+sim.__name__+"(x[:-1])\n")
        r.write(
            "        x[ninputs] = "
            + inspect.getmodule(sim).__name__
            + "."
            + sim.__name__
            + "(*x[:-1])\n"
        )
        # r.write("        x[ninputs] = "+str(sim)+"(x[:-1])\n")
        r.write("        for k in range(0, len(x)):\n")
        r.write("            fout.write(str(x[k]) + ' ')\n")
        r.write("        fout.write(' \\n')\n")
        r.write("\n\n")
        r.write("if __name__ == '__main__':\n")
        r.write("    main()\n")

    os.chmod(name, 509)
    return name
