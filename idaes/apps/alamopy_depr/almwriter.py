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

from .multos import deletefile, movefile, catfile


def almwriter(data, debug, vargs, kwargs):
    """
    This function writes a .alm file for the given data
    """

    xdata = vargs[0]
    zdata = vargs[1]
    if debug["validation"]:
        xvaldata = vargs[2]
        zvaldata = vargs[3]

    with open(data["stropts"]["almname"], "w") as a:
        for arg in data["opts"].keys():
            if arg == "sigma" and data["opts"][arg] < 0:
                continue
            else:
                a.write(arg + " " + str(data["opts"][arg]) + "\n")

        for arg in list(["xlabels", "zlabels"]):
            a.write(
                arg
                + " "
                + str(data["labs"]["save" + arg])[1:-1]
                .replace(",", "")
                .replace("'", "")
                + "\n"
            )

        for arg in data["lstopts"].keys():
            if arg in kwargs.keys() and arg not in list(["xlabels", "zlabels"]):
                temp = arg + " "
                for value in data["lstopts"][arg]:
                    value = "".join(c for c in str(value) if c not in "(),")
                    temp = temp + str(value) + " "
                    temp = temp.replace("[", "")
                    temp = temp.replace("]", "")
                a.write(temp + "\n")

        a.write("trace 1\n")
        for arg in data["stropts"].keys():
            if arg not in list(["almopt", "almname"]):
                a.write(arg + " " + str(data["stropts"][arg]) + "\n")

        # xmax and xmin writing in set4
        if kwargs.get("xmin", None) is not None:
            for arg in data["set4"].keys():
                temp = arg + " "
                for value in data["set4"][arg]:
                    value = "".join(c for c in str(value) if c not in "(),[] ")
                    temp = temp + str(value) + " "
                a.write(temp + "\n")
        else:
            for arg in data["set4"].keys():
                a.write(
                    arg + " " + str(data["set4"][arg])[0:-1].replace(",", "") + "\n"
                )

        # Write data and valdata
        a.write("begin_data\n")
        for i in range(data["opts"]["ndata"]):
            temp = ""
            if data["opts"]["ninputs"] > 1:
                for j in range(data["opts"]["ninputs"]):
                    temp = temp + str(xdata[i][j]) + " "
            else:
                temp = temp + str(xdata[i]) + " "
            for j in range(data["opts"]["noutputs"]):
                temp = temp + str(zdata[i][j]) + " "
            a.write(temp + "\n")
        a.write("end_data\n")
        if debug["validation"]:
            a.write("begin_valdata\n")
            for i in range(data["opts"]["nvaldata"]):
                temp = ""
                for j in range(data["opts"]["ninputs"]):
                    temp = temp + str(xvaldata[i][j]) + " "
                for j in range(data["opts"]["noutputs"]):
                    temp = temp + str(zvaldata[i][j]) + " "
                a.write(temp + "\n")
            a.write("end_valdata\n")

        if "solvemip" in data["opts"].keys():
            if data["opts"]["solvemip"] == 1:
                a.write("gams " + debug["gamsloc"] + "\n")
        if debug["savegams"]:
            a.write("enum1 0\n")
            a.write("enumall 0\n")
            a.write("removescratch 0\n")

        # This option can be used to test the flow of alamo
        if debug["hardset"]:
            temp = "maxterms "
            for i in range(data["opts"]["noutputs"]):
                temp = temp + "1 "
            a.write(temp)

    # Append text file if specified
    if "almopt" in data["stropts"].keys():
        catfile(
            "almtemp", str(data["stropts"]["almname"]), str(data["stropts"]["almopt"])
        )
        deletefile(str(data["stropts"]["almname"]))
        movefile("almtemp " + str(data["stropts"]["almname"]))
