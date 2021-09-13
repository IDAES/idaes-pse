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
# This subroutine is no longer used in the pyomo-ripe implementation


def read(infile):
    # Arguemnts
    # infile : ripe input file to be parsed
    import ripe  # pylint: disable=import-error
    import sys
    import numpy as np

    data, rspace = ripe.data, ripe.rspace
    # debug = ripe.debug

    # Hard coded delimiter
    dlm = " "

    # thigns that needs checks - mechansisms, clc, pclabels,
    # define options to be parsed automatically
    strlist = ["simulator"]
    intlist = ["ccon", "eobj", "npc", "npcval", "ntimedata", "timeseries"]
    llist = [
        "mechanisms",
        "clc",
        "slabels",
        "pclabels",
        "slabels",
        "pclabels",
        "noise",
        "adsspec",
        "adslabel",
        "revlabel",
        "revspec",
    ]
    nlist = ["savescratch", "c0vary", "volume"]
    mlist = ["begin_tdata", "begin_tvaldata", "begin_vdata"]
    endlist = ["end_tdata", "end_tvaldata", "end_vdata"]
    templist = ["T", "Tval", "v"]
    temp2list = ["npct", "npcval", "npcv"]
    labs = {}
    labs2 = {}
    for item, lab in zip(mlist, templist):
        labs[item] = lab
    for item, lab in zip(mlist, temp2list):
        labs2[item] = lab

    # open and read .ripe
    f = open(infile)
    f1 = f.read()
    f.close()
    rd = f1.split("\n")
    ln = 0
    while ln < len(rd):
        if rd[ln] == "":
            ln = ln + 1
            continue
        if rd[ln][0] in list(["%", "#"]):
            ln = ln + 1
            continue
        line = rd[ln].split(dlm)
        while "" in line:
            line.remove("")
        while " " in line:
            line.remove(" ")
        inopt = line[0]
        if inopt not in rspace["inopt"]:
            sys.stdout.write(
                "Error : I don't understand the input file at line : "
                + str(ln + 1)
                + "\n"
            )
            exit
        # check for string inputs
        elif inopt in strlist:
            data["otherdata"][inopt] = str(line[-1])
            sys.stdout.write(data["otherdata"][inopt])
        #            ln+=1
        #            continue
        elif inopt in llist:
            data["otherdata"][inopt] = line[1:]
        elif inopt in intlist:
            data["otherdata"][inopt] = int(line[-1])
        elif inopt in nlist:
            data["otherdata"][inopt] = line[-1]
        elif inopt in mlist:
            n = 1
            nline = rd[ln + n].split(dlm)
            while nline[0] not in endlist:
                try:
                    nline = rd[ln + n].split(dlm)
                except Exception:
                    nline = rd[ln + n]
                n += 1
            if n > 3:
                dlist = rd[ln + 1 : ln + n].split("\n")
                dlist = [s.split(dlm) for s in dlist]
            else:
                dlist = rd[ln + 1]
            data["procdata"][labs[inopt]] = np.asarray(dlist)
            data["procdata"][labs2[inopt]] = n
            ln = ln + n + 1
            continue

        # concdata is read outside of the loop
        elif inopt == "begin_concdata":
            # timeseries data is read differently from data indexed over process conditions
            if "timeseries" in data["otherdata"]:
                ln = readtimeseries(rd, ln, dlm, "procdata", "concdata")
                continue
            else:
                for s in data["otherdata"]["slabels"]:
                    data["concdata"][s] = np.zeros([data["otherdata"]["npc"]])
                t2 = 0
                for i in range(ln + 1, ln + 1 + data["otherdata"]["npc"]):
                    rd2 = rd[i].split(dlm)
                    while "" in rd2:
                        rd2.remove("")
                    while " " in rd2:
                        rd2.remove(" ")
                    t3 = 0
                    for s in data["otherdata"]["slabels"]:
                        data["concdata"][s][t2] = float(rd2[t3])
                        t3 = t3 + 1
                    t2 = t2 + 1
                ln = ln + data["otherdata"]["npc"] + 1
                continue

        # read validation conc data
        elif inopt == "begin_valdata":
            if "timeseries" in data["otherdata"]:
                ln = readtimeseries(rd, ln, dlm, "procvaldata", "concvaldata")
                continue
            else:
                for s in data["otherdata"]["slabels"]:
                    data["concvaldata"][s] = np.zeros([data["otherdata"]["npcval"]])
                t2 = 0
                for i in range(ln + 1, ln + 1 + data["otherdata"]["npcval"]):
                    rd2 = rd[i].split(dlm)
                    while "" in rd2:
                        rd2.remove("")
                    while " " in rd2:
                        rd2.remove(" ")
                    t3 = 0
                    for s in data["otherdata"]["slabels"]:
                        data["concvaldata"][s][t2] = float(rd2[t3])
                        t3 = t3 + 1
                    t2 = t2 + 1
                ln = ln + data["otherdata"]["npcval"] + 1
                continue
        # Read initial concentrations
        elif inopt == "begin_c0data":
            data["procdata"]["c0data"] = {}
            if "pclabels" in data["otherdata"]:
                if "x0" in data["otherdata"]["pclabels"]:
                    if "npc" in data["otherdata"]:
                        temp = int(data["otherdata"]["npc"])
                    else:
                        sys.stdout.write(
                            "Must specify number of process conditions before c0data"
                        )
                elif "T" in data["otherdata"]["pclabels"]:
                    if "npc" in data["otherdata"]:
                        temp = int(data["otherdata"]["npc"])
                    else:
                        sys.stdout.write(
                            "Must specify number of process conditions before c0data"
                        )
                else:
                    temp = 1
            else:
                temp = 1
            for s in data["otherdata"]["slabels"]:
                data["procdata"]["c0data"][s] = np.zeros(temp)
            for i in range(temp):
                rd2 = rd[ln + i + 1].split(dlm)
                while "" in rd2:
                    line.remove("")
                while " " in rd2:
                    line.remove(" ")
                t1 = 0
                for s in data["otherdata"]["slabels"]:
                    data["procdata"]["c0data"][s][i] = float(rd2[t1])
                    t1 = t1 + 1
            ln = ln + 1 + temp
            continue

        # Read custom adsorption ~constants(T)
        # Current ads schemes must be applied to all adsspec
        elif inopt == "begin_ads":
            if "adsspec" not in data["otherdata"]:
                sys.stdout.write("Define species involved in adsorption\n")
                exit
            elif "adslabel" not in data["otherdata"]:
                sys.stdout.write("Define labels for adsorption terms\n")
                exit
            # Initilize the data dictionary
            data["procdata"]["adata"] = {}
            for i in data["otherdata"]["adslabel"]:
                data["procdata"]["adata"][i] = {}
                for j in data["otherdata"]["slabels"]:
                    data["procdata"]["adata"][i][j] = ""
            t1 = 0
            for i in range(
                ln + 1,
                ln
                + len(data["otherdata"]["adslabel"]) * len(data["otherdata"]["adsspec"])
                + 1,
            ):
                data["procdata"]["adata"][data["otherdata"]["adslabel"][t1]][
                    data["otherdata"]["slabels"][i - ln - 1]
                ] = str(rd[i])
                if i - ln - i == len(data["otherdata"]["adsspec"]):
                    t1 = t1 + 1
            ln = (
                ln
                + 1
                + len(data["otherdata"]["adslabel"]) * len(data["otherdata"]["adsspec"])
            )
            # Take care of species not included in adsorption
            for s in [
                val
                for val in data["otherdata"]["slabels"]
                if val not in data["otherdata"]["adsspec"]
            ]:
                for r in data["otherdata"]["adslabel"]:
                    data["procdata"]["adata"][r][s] = "0.0"
            continue

        # Read eq constant(T)
        # 1 label per reaction > one line in data section
        elif inopt == "begin_rev":
            if "revspec" not in data["otherdata"]:
                sys.stdout.write("Define species involved in reversible reaction\n")
                exit
            elif "revlabel" not in data["otherdata"]:
                sys.stdout.write("Define labels for eq constants\n")
                exit
            # Initilize the data dictionary
            data["procdata"]["revdata"] = {}
            for i in data["otherdata"]["revlabel"]:
                data["procdata"]["revdata"][i] = {}
                for j in data["otherdata"]["revspec"]:
                    data["procdata"]["revdata"][i][j] = ""
            t1 = 0
            for i in range(ln + 1, ln + len(data["otherdata"]["revspec"]) + 1):
                data["procdata"]["revdata"][data["otherdata"]["revlabel"][0]][
                    data["otherdata"]["revspec"][i - ln - 1]
                ] = str(rd[i])
            ln = ln + 1 + len(data["otherdata"]["revspec"])
            continue

        # End of the loop
        ln = ln + 1
    # perform some checks
    if "ntimedata" not in data["otherdata"]:
        data["otherdata"]["ntimedata"] = 0


def readtimeseries(rd, ln, dlm, st1, st2):
    import ripe  # pylint: disable=import-error

    # Read time series data
    import numpy as np

    data, debug = ripe.data, ripe.debug

    data[st1]["time"] = np.array(
        np.zeros([data["otherdata"]["timeseries"], data["otherdata"]["ntimedata"]])
    )
    for s in data["otherdata"]["slabels"]:
        data[st2][s] = np.zeros(
            [data["otherdata"]["timeseries"], data["otherdata"]["ntimedata"]]
        )
        t1 = 0
        t2 = 0
        for i in range(
            ln + 1, ln + 1 + data["otherdata"]["ntimedata"] * data["otherdata"]["npc"]
        ):
            rd2 = rd[i].split(dlm)
            if len(rd2) == 1:
                rd2 = rd[i].split("\t")
                while "" in rd2:
                    rd2.remove("")
                while " " in rd2:
                    rd2.remove(" ")
                data["procdata"]["time"][t1][t2] = float(rd2[0])
                t3 = 0
                for s in data["otherdata"]["slabels"]:
                    # Check for true model specified to debug
                    if not debug["true"]["use"]:
                        data[st2][s][t1][t2] = float(rd2[t3 + 1])
                    else:
                        data[st2][s][t1][t2] = 1 - np.exp(
                            -debug["true"]["k"]
                            * np.exp(-debug["true"]["e"] / (8.314 * data[st1]["T"][t1]))
                            * data[st1]["time"][t1][t2]
                        )
                    t3 = t3 + 1
                t2 = t2 + 1
                if t2 == data["otherdata"]["ntimedata"]:
                    t1 = t1 + 1
                    t2 = 0
        ln = ln + data["otherdata"]["ntimedata"] * data["otherdata"]["timeseries"] + 3
        return ln
