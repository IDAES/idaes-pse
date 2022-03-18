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
# This subroutine is no longer used in the ripe-pyomo implementation

# PYLINT-TODO-CHECK "ripe" and "alamopy" cause pylint to report import errors
# this should be checked to see if these are real errors and/or this module is still in use
# pylint: disable=import-error


def writeripe(aterm, ccon):
    import os
    import numpy as np
    import ripe

    data = ripe.data

    # wrapper for solving RIPE problems
    ndata = data["otherdata"]["npc"] * max(data["otherdata"]["ntimedata"], 1)
    if ccon < 0:
        clist = range(1, min(ndata, len(data["otherdata"]["slabels"])))
        [aterm, kils] = writemodel(aterm, ccon)
        if kils == 1:
            os.system("gams ripe.gms")
            f = open("rates.put")
            lf = f.read()
            f.close()
            lf2 = lf.split("\n")
            lf2 = lf2[0].split(" ")
            lf2 = [float(x) for x in lf2 if x != ""]
            sigmaest = [0] * len(lf2)
            bic = [0] * int(ndata)
            for i in range(len(lf2)):
                sigmaest[i] = lf2[i] / (ndata - 1)
        else:
            exit
    else:
        clist = list([ccon])
    for ccon in clist:
        [aterm, kils] = writemodel(aterm, ccon)
        if kils == 1:
            os.system("gams ripe.gms")
            f = open("rates.put")
            lf = f.read()
            f.close()
            lf2 = lf.split("\n")
            lf2 = lf2[0].split(" ")
            lf2 = [float(x) for x in lf2 if x != ""]
            sigmaest = [0] * len(lf2)
            bic = [0] * int(ndata)
            for i in range(len(lf2)):
                sigmaest[i] = lf2[i] / (ndata - 1)
            bic[ccon] = (
                sum(np.divide(lf2, sigmaest))
                + np.log(data["otherdata"]["npc"] * len(data["otherdata"]["slabels"]))
                * ccon
            )
            # ripeci(aterm, ccon)
        else:
            exit


def writemodel(aterm, ccon):
    import numpy as np
    import os
    import sys
    from idaes.apps import alamopy_depr as alamopy
    import ripe
    from alamopy import deletefile

    data, rspace = ripe.data, ripe.rspace

    # First an activity term matrix is defined to make assignment easier later
    if "clc" in data["otherdata"]:
        aterm = np.zeros(
            [
                data["otherdata"]["npc"],
                data["otherdata"]["ntimedata"],
                data["rxnet"]["mechanisms"],
            ]
        )
        for i in range(data["otherdata"]["npc"]):
            for j in range(data["otherdata"]["ntimedata"]):
                # t is an index on the current mechanism
                t = 0
                if "clc" in data["otherdata"]:
                    for m in range(len(data["otherdata"]["clc"])):
                        mech = data["otherdata"]["clc"][m]
                        # Check for CLC mechanisms
                        if mech == "c1":
                            aterm[i, j, t] = (
                                1.0 - data["concdata"]["X"][i][j]
                            )  # random nucleation
                        if mech == "c2":
                            aterm[i, j, t] = (
                                2.0
                                / 3.0
                                * data["concdata"]["X"][i][j] ** (2.0 / 3.0 - 3.0 / 2.0)
                            )  # Power law of oer 2/3
                        if mech == "c3":
                            aterm[i, j, t] = 2.0 * data["concdata"]["X"][i][j] ** (
                                2.0 - 1.0 / 2.0
                            )  # Power law model of order 2
                        if mech == "c4":
                            aterm[i, j, t] = 4.0 * data["concdata"]["X"][i][j] ** (
                                4.0 - 1.0 / 4.0
                            )  # power law model of order 4
                        if mech == "c5":
                            aterm[i, j, t] = 3.0 * data["concdata"]["X"][i][j] ** (
                                3.0 - 1.0 / 3.0
                            )  # power law model of order 4
                        if mech == "c6":
                            aterm[i, j, t] = 1.0  # power law model of order 1
                        if mech == "c7":
                            aterm[i, j, t] = (
                                2.0
                                * (1.0 - data["concdata"]["X"][i][j])
                                * np.power(
                                    -np.log(1.0 - data["concdata"]["X"][i][j]), 1.5
                                )
                            )  # nucleation mode - Avrami-Erofeev eq of order 2
                        if mech == "c8":
                            aterm[i, j, t] = (
                                3.0
                                * (1.0 - data["concdata"]["X"][i][j])
                                * np.power(
                                    -np.log(1.0 - data["concdata"]["X"][i][j]),
                                    3.0 - (1.0 / 3.0),
                                )
                            )  # '' order 3
                        if mech == "c9":
                            aterm[i, j, t] = (
                                0.5
                                * (1.0 - data["concdata"]["X"][i][j])
                                * (-np.log(1.0 - data["concdata"]["X"][i][j]))
                                ** (0.5 - 1.0 / 0.5)
                            )  # '' order 0.5
                        if mech == "c10":
                            aterm[i, j, t] = (
                                1.5
                                * (1.0 - data["concdata"]["X"][i][j])
                                * np.power(
                                    -np.log(1.0 - data["concdata"]["X"][i][j]),
                                    1.5 - (1.0 / 1.5),
                                )
                            )  # '' order 1.5
                        if mech == "c11":
                            aterm[i, j, t] = (
                                4.0
                                * (1.0 - data["concdata"]["X"][i][j])
                                * (-np.log(1.0 - data["concdata"]["X"][i][j]))
                                ** (4.0 - 1.0 / 4.0)
                            )  # '' order 4
                        if mech == "c12":
                            aterm[i, j, t] = data["concdata"]["X"][i][j] * (
                                1.0 - data["concdata"]["X"][i][j]
                            )  # Prout tompkins
                        # diffusion models
                        if mech == "c13":
                            aterm[i, j, t] = (
                                3.0
                                * (1.0 - data["concdata"]["X"][i][j]) ** (1.0 / 3.0)
                                * (
                                    2.0
                                    * (
                                        (1.0 - data["concdata"]["X"][i][j])
                                        ** (-1.0 / 3.0)
                                        - 1.0
                                    )
                                )
                                ** (-1)
                            )  # Jander equation (3d difsusion spherical sym)
                        if mech == "c14":
                            aterm[i, j, t] = (
                                3.0
                                / 2.0
                                * (1.0 + data["concdata"]["X"][i][j]) ** (2.0 / 3.0)
                                * (
                                    (
                                        (1.0 + data["concdata"]["X"][i][j])
                                        ** (-1.0 / 3.0)
                                        - 1.0
                                    )
                                )
                                ** (-1)
                            )  # antijander
                        if mech == "c15":
                            aterm[i, j, t] = 1.0 / (
                                -np.log(1.0 - data["concdata"]["X"][i][j])
                            )  # valensi
                        if mech == "c16":
                            aterm[i, j, t] = 1.0 / (
                                2.0 * data["concdata"]["X"][i][j]
                            )  # parabolic
                        if mech == "c17":
                            aterm[i, j, t] = (3.0 / 2.0) * (
                                (
                                    (1.0 - data["concdata"]["X"][i][j]) ** (-1.0 / 3.0)
                                    - 1.0
                                )
                            ) ** (
                                -1
                            )  # Ginstling-Brounstein eq 3d diff cylindrical
                        if mech == "c18":
                            aterm[i, j, t] = (
                                3.0
                                / 2.0
                                * (1.0 - data["concdata"]["X"][i][j]) ** (4.0 / 3.0)
                                * (
                                    (
                                        (1.0 - data["concdata"]["X"][i][j])
                                        ** (-1.0 / 3.0)
                                        - 1.0
                                    )
                                )
                                ** (-1)
                            )  # zhuralev, lesohkhin, tempelman
                        if mech == "c19":
                            aterm[i, j, t] = (1.0 - data["concdata"]["X"][i][j]) ** (
                                2.0 / 3.0
                            )
                        # if(mech=='c19'): aterm[i,j,t]=  (0.2)**(0.535) * 1.5 * (0.112)**(1.0/1.5)*((1.0)-data['concdata']['X'][i][j])*(-np.log(1.0-data['concdata']['X'][i][j]))**(1.0-1.0/1.5)                                                                                                                                                                      #                    if(mech=='c20'): aterm[i,j,t]=  (0.3)**(0.312) * 2.76 * (11.24)**(1.0/2.76)*((0.18 *1.479)-data['concdata']['X'][i][j])*(-np.log(1.0-data['concdata']['X'][i][j])/(0.18 *1.479))**(1.0-1.0/2.76)                                                                                                                                              #                    if(mech=='c21'): aterm[i,j,t]=  (0.3)**(2.13) * 1.63 * (28.4)**(1.0/1.63)*((0.135*1.479)-data['concdata']['X'][i][j])*(-np.log(1.0-data['concdata']['X'][i][j])/(0.135*1.479))**(1.0-1.0/1.63)
                        t = t + 1
    else:
        for m in range(data["rxnet"]["mechanisms"]):
            tempa = 1
            for s in data["otherdata"]["slabels"]:
                if data["rxnet"][i, s] < 0:
                    tempa = tempa * float(data["concdata"][s][i][j]) ** float(
                        abs(data["rxnet"][i, s])
                    )
                    aterm[i, j, m] = tempa

        for q in range(data["rxnet"]["mechanisms"]):
            if np.isinf(aterm[i, j, q]):
                aterm[i, j, q] = 0
                if np.isnan(aterm[i, j, q]):
                    aterm[i, j, q] = 0

        if "time" not in data["procdata"]:
            for q in range(data["rxnet"]["nrx"]):
                temp = 1
                for s in data["otherdata"]["slabels"]:
                    if data["rxnet"][q, s] < 0:
                        temp = temp * data["concdata"][s][i] ** abs(data["rxnet"][q, s])
                aterm[i, q] = temp

    # check for isothermal, single species observations
    if len(data["otherdata"]["slabels"]) == 1:
        if "T" not in data["otherdata"]["pclabels"]:
            sys.stdout = open(os.devnull, "w")
            r = alamopy.doalamo(
                np.asarray(aterm[0, :, :]),
                np.asarray(
                    data["ratedata"][str(data["otherdata"]["slabels"][0])][0, :]
                ),
                linfcns=1,
                constant=0,
                maxterms=1,
                expandoutput=1,
            )
            r = alamopy.almconfidence(r)
            sys.stdout = sys.__stdout__
            selected = r["model"].split(" * ")[-1]
            selected = selected[1:]
            coeff = r["model"].split(" * ")[0]
            coeff = coeff.split(" = ")[-1]
            cfi = r["conf_inv"][0].split(" : ")[-1]
            slist = [
                "Selected model is : " + rspace["clcnames"][int(selected) - 1],
                "Selected model form : " + rspace["clcforms"][int(selected) - 1],
                "Estimated rate constant : " + cfi,
            ]
            sys.stdout.write(slist[0] + "\n")
            sys.stdout.write(slist[1] + "\n")
            sys.stdout.write(slist[2] + "\n")
            with open("output.txt", "w") as o:
                o.write(slist[0] + "\n")
                o.write(slist[1] + "\n")
                o.write(slist[2] + "\n")
            deletefile("*logscratch*", "tempalm.py", "sv.alm", "tmpscratch")
            return [aterm, 0]
    else:
        writeminlp(aterm, ccon)
        return [aterm, 0]


def writeminlp(aterm, ccon):
    import ripe
    import numpy as np

    # Main ripe subroutine
    # initialize some stuff

    data, debug = ripe.data, ripe.debug

    ndatadic = {}
    for what in data["procdata"]:
        try:
            ndatadic[what] = len(data["procdata"][what])
        except Exception:
            ndatadic[what] = 1
        if "time" == what:
            ndatadic[what] = np.size(data["procdata"][what], 1)

    with open("ripe.gms", "w") as r:
        r.write("$offdigit\n")
        r.write("$offsymxref offsymlist\n")

        # write out species no matter what
        if "clc" in data["otherdata"]:
            r.write("set s /X/;\n")
        else:
            r.write("set s /a*" + data["otherdata"]["slabels"][-1] + "/;\n")
        r.write("alias(s,s1,s2,s3,s4);\n")

        # write npc set
        if data["otherdata"]["npc"] > 0:
            r.write("set pc /1*" + str(data["otherdata"]["npc"]) + "/;\n")

        # Is it isothermal?
        if "T" in data["procdata"]:
            r.write("parameters temp(pc);\n")
            if ndatadic["T"] > 1:
                #            r.write('set pc /1*'+str(ndatadic['T'])+'/;\n')
                for i in range(ndatadic["T"]):
                    r.write(
                        "temp('"
                        + str(i + 1)
                        + "')="
                        + str(data["procdata"]["T"][i])
                        + ";\n"
                    )
            else:
                r.write("temp('1')=" + str(data["procdata"]["T"]) + ";\n")
        if "time" in data["procdata"]:
            r.write("set time /1*" + str(ndatadic["time"]) + "/;\n")
            r.write("parameters t(pc,time);\n")
            for i in range(ndatadic["T"]):
                for j in range(ndatadic["time"]):
                    r.write(
                        "t('"
                        + str(i + 1)
                        + "','"
                        + str(j + 1)
                        + "')="
                        + str(data["procdata"]["time"][i][j])
                        + ";\n"
                    )

        # Some more equations and sets
        r.write("equations obj,cardcon,m1l,m1u,mc;\n")
        if "T" in data["otherdata"]["pclabels"]:
            r.write("equations m2l, m2u;\n")
        r.write("set m /1*" + str(data["rxnet"]["mechanisms"]) + "/;\n")
        r.write("set rx /1*" + str(data["rxnet"]["nrx"]) + "/;\n")

        if "clc" in data["otherdata"]:
            # iscomm = 0
            tstr = "set mclc(m) / "
            for m in data["otherdata"]["clc"]:
                tstr = tstr + str(data["gprxindex"][0][m])
                if m != data["otherdata"]["clc"][-1]:
                    tstr = tstr + ","
                else:
                    tstr = tstr + "/;\n"
            r.write(tstr)

        # Now we gotta write all those concentration values
        if data["otherdata"]["npc"] > 0:
            if "time" in data["procdata"]:
                r.write(
                    "parameters conc(pc,time,s),rate(pc,time,s),aterm(pc,time,m);\n"
                )
                #                r.write('parameters conc(pc,time,s),rate(pc,time,s),aterm(pc,time,rc);\n')
                for i in range(ndatadic["T"]):
                    for j in range(ndatadic["time"]):
                        for s in data["otherdata"]["slabels"]:
                            r.write(
                                "conc('"
                                + str(i + 1)
                                + "','"
                                + str(j + 1)
                                + "','"
                                + s
                                + "')="
                                + str(data["concdata"][s][i][j])
                                + ";\n"
                            )
                            r.write(
                                "rate('"
                                + str(i + 1)
                                + "','"
                                + str(j + 1)
                                + "','"
                                + s
                                + "')="
                                + str(data["ratedata"][s][i][j])
                                + ";\n"
                            )
                            for q in range(data["rxnet"]["mechanisms"]):
                                r.write(
                                    "aterm('"
                                    + str(i + 1)
                                    + "','"
                                    + str(j + 1)
                                    + "','"
                                    + str(q + 1)
                                    + "')="
                                    + str(aterm[i, j, q])
                                    + ";\n"
                                )

            else:
                r.write("parameters conc(pc,s),rate(pc,s);\n")
                for i in range(data["otherdata"]["npc"]):
                    for s in data["otherdata"]["slabels"]:
                        r.write(
                            "conc('"
                            + str(i + 1)
                            + "','"
                            + s
                            + "')="
                            + str(data["concdata"][s][i])
                            + ";\n"
                        )
                        r.write(
                            "rate('"
                            + str(i + 1)
                            + "','"
                            + s
                            + "')="
                            + str(data["ratedata"][s][i])
                            + ";\n"
                        )

        elif "time" in data["procdata"]:
            r.write("parameters conc(time,s),rate(time,s),aterm(time,m);\n")
            for j in range(ndatadic["time"]):
                for s in data["otherdata"]["slabels"]:
                    r.write(
                        "conc('"
                        + str(j + 1)
                        + "','"
                        + s
                        + "')="
                        + str(data["concdata"][s][0][j])
                        + ";\n"
                    )
                    r.write(
                        "rate('"
                        + str(j + 1)
                        + "','"
                        + s
                        + "')="
                        + str(data["ratedata"][s][0][j])
                        + ";\n"
                    )
                for m in range(len(data["otherdata"]["clc"])):
                    r.write(
                        "aterm('"
                        + str(j + 1)
                        + "','"
                        + str(m + 1)
                        + "')="
                        + str(aterm[0, j, m])
                        + ";\n"
                    )
        if "ads" in data["otherdata"]["mechanisms"]:
            for i in data["otherdata"]["adslabel"]:
                r.write("parameter " + i + "(pc,s);\n")
                for j in data["otherdata"]["slabels"]:
                    r.write(
                        i
                        + "(pc,'"
                        + j
                        + "')="
                        + data["procdata"]["adata"][i][j]
                        + ";\n"
                    )
        if 1 == 0:
            #        if ( 'adata' in data['procdata'] ):
            r.write("set adss(s) / ")
            for s in data["otherdata"]["adsspec"]:
                r.write(s)
                if s != data["otherdata"]["adsspec"][-1]:
                    r.write(",")
            r.write(" /;\n")
            r.write("set noads(s) / ")
            for s in [
                val
                for val in data["otherdata"]["slabels"]
                if val not in data["otherdata"]["adsspec"]
            ]:
                r.write(s)
                if (
                    s
                    != [
                        val
                        for val in data["otherdata"]["slabels"]
                        if val not in data["otherdata"]["adsspec"]
                    ][-1]
                ):
                    r.write(",")
            r.write(" /;\n")
            r.write("equations adeq, adyeq;\n")
            # This section write kads(pc,rx)=
            for i in data["otherdata"]["adslabel"]:
                r.write("parameter " + i + "(pc,s);\n")
                for j in data["otherdata"]["slabels"]:
                    r.write(
                        i
                        + "(pc,'"
                        + j
                        + "')="
                        + data["procdata"]["adata"][i][j]
                        + ";\n"
                    )

            if "time" in data["procdata"]:
                r.write("parameter aax(pc,rx,time);\n")
                for n in range(data["otherdata"]["npc"]):
                    for x in range(data["rxnet"]["nrx"]):
                        for t in range(data["otherdata"]["ntimedata"]):
                            for s in data["otherdata"]["slabels"]:
                                if (
                                    data["rxnet"][x, s] < 0
                                    and s in data["otherdata"]["adsspec"]
                                ):
                                    r.write(
                                        "aax('"
                                        + str(n + 1)
                                        + "','"
                                        + str(x + 1)
                                        + "','"
                                        + str(t + 1)
                                        + "')=kads('"
                                    )
                                    r.write(
                                        str(n + 1)
                                        + "','"
                                        + s
                                        + "')/(1+sum(s1,conc('"
                                        + str(n + 1)
                                        + "','"
                                        + str(t + 1)
                                        + "',s1)"
                                    )
                                    r.write("*kads('" + str(n + 1) + "',s1)));\n")
                                elif (
                                    data["rxnet"][x, s] < 0
                                    and s not in data["otherdata"]["adsspec"]
                                ):
                                    #  else:
                                    r.write(
                                        "aax('"
                                        + str(n + 1)
                                        + "','"
                                        + str(x + 1)
                                        + "','"
                                        + str(t + 1)
                                        + "')=1.0;\n"
                                    )
                r.write("variable adsaux(pc,rx,time);\n")
                r.write("adsaux.l(pc,rx,time)=1;\n")
            else:
                r.write("parameter aax(pc,rx);\n")
                for n in range(data["otherdata"]["npc"]):
                    for x in range(data["rxnet"]["nrx"]):
                        swtch = 0
                        for s in data["otherdata"]["slabels"]:
                            if (
                                (data["rxnet"][x, s] < 0)
                                and s in data["otherdata"]["adsspec"]
                                and swtch == 0
                            ):
                                r.write(
                                    "aax('"
                                    + str(n + 1)
                                    + "','"
                                    + str(x + 1)
                                    + "')=kads('"
                                )
                                r.write(
                                    str(n + 1)
                                    + "','"
                                    + s
                                    + "')/(1+sum(s1,conc('"
                                    + str(n + 1)
                                    + "',s1)"
                                )
                                r.write("*kads('" + str(n + 1) + "',s1)));\n")
                                swtch = 1
                            elif (
                                (data["rxnet"][x, s] < 0)
                                and s not in data["otherdata"]["adsspec"]
                                and swtch == 0
                            ):
                                r.write(
                                    "aax('"
                                    + str(n + 1)
                                    + "','"
                                    + str(x + 1)
                                    + "')=1.0;\n"
                                )
                                swtch = 1
                r.write("variable adsaux(pc,rx);\n")
                r.write("adsaux.l(pc,rx)=1;\n")

        # Get that equilibrium info in there too
        temp1 = 0
        if "revdata" in data["procdata"]:
            for i in data["otherdata"]["revlabel"]:
                jj = 0
                r.write("parameter " + i + "(pc,rx);\n")
                for rx in range(data["rxnet"]["nrx"]):
                    if rx in data["rxnet"]["rev"]:
                        j = data["otherdata"]["revspec"][jj]
                        r.write(
                            i
                            + "(pc,'"
                            + str(int(rx) + 1)
                            + "')="
                            + data["procdata"]["revdata"][i][j]
                            + ";\n"
                        )
                        jj = jj + 1
                    else:
                        r.write(i + "(pc,'" + str(int(rx) + 1) + "')=" + str(1) + ";\n")

        if 1 == 0:
            #        if ('revdata' in data['procdata']):
            temp1 = "set rvx1(rx) / 1"
            temp2 = "set rvx2(rx) / 2"
            for i in range(3, data["rxnet"]["nrx"]):
                if i % 2 == 0:
                    temp2 = temp2 + ", " + str(i)
                else:
                    temp1 = temp1 + ", " + str(i)
            temp1 = temp1 + " /;\n"
            temp2 = temp2 + " /;\n"
            r.write(temp1)
            r.write(temp2)

        r.write("parameter sigest(s);\n")
        sc = 0
        for s in data["otherdata"]["slabels"]:
            try:
                tempn = str(data["otherdata"]["noise"][sc])
            except Exception:
                tempn = "0"
            r.write("sigest('" + s + "')=" + tempn + ";\n")
            sc = sc + 1

        r.write("parameter rmat(m,rx,s),pmat(m,rx,s),revp(m,rx,s);\n")
        r.write("revp(m,rx,s) = 0;\n")
        for i in range(data["rxnet"]["nrx"]):
            for j in range(data["rxnet"]["mechanisms"]):
                if "clc" in data["otherdata"]:
                    mech = data["otherdata"]["clc"][j]
                else:
                    mech = data["otherdata"]["mechanisms"][j]
                for s in data["otherdata"]["slabels"]:
                    if 1 == 1:
                        # first handle reversible reactions - we are going to use the pmat for ma the reactants
                        if mech == "rev":
                            if data["rxnet"][i, s] < 0:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(data["rxnet"][i, s])
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                            elif data["rxnet"][i, s] > 0:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(data["rxnet"][i, s])
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(data["rxnet"][i, s])
                                    + ";\n"
                                )
                            else:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                        elif mech == "ma":
                            if data["rxnet"][i, s] < 0:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(data["rxnet"][i, s])
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(abs(data["rxnet"][i, s]))
                                    + ";\n"
                                )
                            elif data["rxnet"][i, s] > 0:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(data["rxnet"][i, s])
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(data["rxnet"][i, s])
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                            else:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                        elif mech == "homcat":
                            r.write(
                                "rmat('"
                                + str(j + 1)
                                + "','"
                                + str(i + 1)
                                + "','"
                                + s
                                + "')="
                                + str(data["rxnet"][i, s])
                                + ";\n"
                            )
                            r.write(
                                "pmat('"
                                + str(j + 1)
                                + "','"
                                + str(i + 1)
                                + "','"
                                + s
                                + "')="
                                + str(data["rxnet"][i, s])
                                + ";\n"
                            )
                        elif mech == "ads":
                            if data["rxnet"][i, s] < 0:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(data["rxnet"][i, s])
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(abs(data["rxnet"][i, s]))
                                    + ";\n"
                                )
                            elif data["rxnet"][i, s] > 0:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(data["rxnet"][i, s])
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(data["rxnet"][i, s])
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                            else:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                        elif mech == "iT":
                            if data["rxnet"][i, s] < 0:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(data["rxnet"][i, s])
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(abs(data["rxnet"][i, s]))
                                    + ";\n"
                                )
                            elif data["rxnet"][i, s] > 0:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(data["rxnet"][i, s])
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(data["rxnet"][i, s])
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                            else:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                        elif "clc" in data["otherdata"]:
                            #                            for q in data['otherdata']['clc']:
                            if data["rxnet"][i, s] < 0:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(data["rxnet"][i, s])
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(abs(data["rxnet"][i, s]))
                                    + ";\n"
                                )
                            elif data["rxnet"][i, s] > 0:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(data["rxnet"][i, s])
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                            else:
                                r.write(
                                    "rmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                                r.write(
                                    "pmat('"
                                    + str(j + 1)
                                    + "','"
                                    + str(i + 1)
                                    + "','"
                                    + s
                                    + "')="
                                    + str(0)
                                    + ";\n"
                                )
                                # j=j+1
                    else:
                        r.write(
                            "rmat('"
                            + str(j + 1)
                            + "','"
                            + str(i + 1)
                            + "','"
                            + s
                            + "')=0;\n"
                        )
                        r.write(
                            "pmat('"
                            + str(j + 1)
                            + "','"
                            + str(i + 1)
                            + "','"
                            + s
                            + "')=0;\n"
                        )

        r.write("binary variables y(rx,m);\n")
        if 1 == 0:
            #        if ( 'adata' in data['procdata']):
            r.write("binary variable yads;\n")
        r.write("positive variables k(rx,m);\n")
        if data["otherdata"]["npc"] > 0:
            if "T" in data["otherdata"]["pclabels"]:
                r.write("positive variables E(rx,m);\n")
        r.write("variable modelcomp;\n")
        if "time" in data["procdata"]:
            r.write("variables estrate(pc,time,s);\n")
            r.write("equations er(pc,time,s);\n")
        else:
            r.write("variables estrate(pc,s);\n")
            r.write("equations er(pc,s);\n")

        if data["otherdata"]["eobj"] == 1:
            if "time" in data["procdata"]:
                r.write("positive variables err(pc,time,s);\n")
            else:
                r.write("positive variables err(pc,s);\n")
            r.write("equations res1,res2;\n")

        if "time" in data["procdata"]:
            if 1 == data["otherdata"]["eobj"]:
                r.write(
                    "obj.. modelcomp =e= sum(pc,sum(s,sum(time,err(pc,time,s))));\n"
                )
            else:
                r.write(
                    "obj.. modelcomp =e= sum(pc,sum(s,sum(time,power(rate(pc,time,s)-estrate(pc,time,s),2))));\n"
                )
            # first write irrev model
            if "clc" in data["otherdata"]:
                if 1 == data["otherdata"]["eobj"]:
                    r.write(
                        "  res1(pc,time,s).. rate(pc,time,s) - estrate(pc,time,s) =l= err(pc,time,s);\n"
                    )
                    r.write(
                        "  res2(pc,time,s).. estrate(pc,time,s) - rate(pc,time,s) =l= err(pc,time,s);\n"
                    )
                r.write("  er(pc,time,s).. estrate(pc,time,s) =e=\n")
                if data["otherdata"]["npc"] > 1:
                    r.write(
                        "  sum(mclc,(sum(rx,(rmat(mclc,rx,s)*k(rx,mclc)*exp((-1*E(rx,mclc))/(8.314*(temp(pc))))*\n"
                    )
                else:
                    r.write("  sum(mclc,(sum(rx,(rmat(mclc,rx,s)*k(rx,mclc)*\n")
                r.write("  aterm(pc,time,mclc)))));\n")
            else:
                # Need to write independent T
                r.write("  er(pc,time,s).. estrate(pc,time,s) =e=\n")
                r.write(
                    "  (sum(rx,(rmat('1',rx,s)*k(rx,'1')*exp(-E(rx,'1')/(8.314*temp(pc)))*\n"
                )
                r.write("  aterm(pc,time,'1'))));\n")
            # Write the adsorption parameter is needed
            if 1 == 0:
                #            if ( 'adata' in data['procdata']):
                r.write("adsaux(pc,rx,time)*\n")
            if "rev" in data["otherdata"]["mechanisms"]:
                r.write("rmat('2',rx,s)*k(rx,'2')*exp(-E(rx,'2')/(8.314*temp(pc)))*\n")
                if "adata" in data["procdata"]:
                    r.write("adsaux(pc,rx,time)*\n")
                r.write("  (prod(s1,power(conc(pc,time,s1),pmat('2',rx,s1)))\n")
                r.write(
                    "-(prod(s3,power(conc(pc,time,s3),revp('2',rx,s3)))/keq(pc))))))),2))));\n"
                )
        else:
            mcount = 1
            # Input reactor is not dynamic
            r.write(
                "obj.. modelcomp =e= sum(pc,sum(s,(1/sigest(s))*power((rate(pc,s)-estrate(pc,s)),2)));\n"
            )
            # Estiamted rate
            r.write("  er(pc,s).. estrate(pc,s) =e=\n")
            if "iT" in data["otherdata"]["mechanisms"]:
                mcount = 1
                r.write(
                    "  sum(rx,(rmat('"
                    + str(mcount)
                    + "',rx,s)*k(rx,'"
                    + str(mcount)
                    + "')*\n"
                )
                if len(data["otherdata"]["mechanisms"]) > 1:
                    r.write(
                        "  (prod(s1,power(conc(pc,s1),pmat('"
                        + str(mcount)
                        + "',rx,s1)))))+(\n"
                    )
                else:
                    r.write(
                        "  (prod(s1,power(conc(pc,s1),pmat('"
                        + str(mcount)
                        + "',rx,s1))))));\n"
                    )
            if "ma" in data["otherdata"]["mechanisms"]:  # first write irrev models
                r.write(
                    "  sum(rx,(rmat('1',rx,s)*k(rx,'1')*exp(-E(rx,'1')/(8.314*temp(pc)))*\n"
                )
                r.write("  (prod(s1,power(conc(pc,s1),pmat('1',rx,s1)))))+(\n")
            #            if ( 1==0):
            #            if ( 'adata' in data['procdata']):
            #                r.write("adsaux(pc,rx)*\n")

            # now rev rxns
            if "rev" in data["otherdata"]["mechanisms"]:
                r.write("rmat('2',rx,s)*k(rx,'2')*exp(-E(rx,'2')/(8.314*temp(pc)))*\n")
                if 1 == 0:
                    #                if ( 'adata' in data['procdata']):
                    r.write("adsaux(pc,rx)*\n")
                r.write("  (prod(s1,power(conc(pc,s1),pmat('1',rx,s1)))\n")
                r.write(
                    "-(prod(s3,power(conc(pc,s3),pmat('2',rx,s3)))/keq(pc,rx))))+(\n"
                )
            # r.write("-(prod(s3,power(power(conc(pc,s3),revp('2',rx,s3)),pmat('2',rx,s3))/keq(pc,rx)))))+(\n")
            # Now the homcat models
            if "homcat" in data["otherdata"]["mechanisms"]:
                r.write("rmat('3',rx,s)*k(rx,'3')*exp(-E(rx,'3')/(8.314*temp(pc)))*\n")
                if 1 == 0:
                    # if ( 'adata' in data['procdata']):
                    r.write("adsaux(pc,rx)*\n")
                r.write("  (prod(s1,power(conc(pc,s1),pmat('3',rx,s1)))))+(\n")
            #            r.write("-(prod(s3,power(power(conc(pc,s3),revp('3',rx,s3)),pmat('2',rx,s3))/keq(pc,rx)))))))),2)));\n")
            # Now surface reactions
            if "ads" in data["otherdata"]["mechanisms"]:
                r.write("rmat('4',rx,s)*k(rx,'4')*exp(-E(rx,'4')/(8.314*temp(pc)))*\n")
                r.write(
                    "sum(s1,kads(pc,s1)*conc(pc,s1))/power(1+sum(s1,kads(pc,s1)*conc(pc,s1)),2)*\n"
                )
                r.write("  (prod(s1,power(conc(pc,s1),pmat('1',rx,s1))))));\n")
        #            r.write("))),2)));\n")

        r.write("cardcon.. sum(rx,sum(m,y(rx,m)))=e=" + str(ccon) + ";\n")
        # Set variable bounds here !!
        # Big M constraints big-M
        r.write("m1l(rx,m).. k(rx,m)=g=0*y(rx,m);\n")
        r.write("m1u(rx,m).. k(rx,m)=l=20*y(rx,m);\n")
        if data["otherdata"]["npc"] > 1:
            if "T" in data["otherdata"]["pclabels"]:
                r.write("m2l(rx,m).. E(rx,m)=g=500*y(rx,m);\n")
                r.write("m2u(rx,m).. E(rx,m)=l=10000*y(rx,m);\n")
        r.write("mc(rx).. sum(m,y(rx,m))=l=1;\n")
        # Disqualify irrev back rx if rev is chosen
        if "rev" in data["otherdata"]["mechanisms"]:
            for i in range(data["rxnet"]["nrx"]):
                if i in data["rxnet"]["rev"]:
                    r.write("equation req" + str(i + 1) + ";\n")
                    r.write(
                        "req"
                        + str(i + 1)
                        + ".. y('"
                        + str(i + 1)
                        + "','2')+y('"
                        + str(i + 2)
                        + "','1')=l=1;\n"
                    )
        # handle adsorption binary
        if 1 == 0:
            #        if ('adata' in data['procdata']):
            if "time" not in data["procdata"]:
                r.write("adeq(pc,rx).. adsaux(pc,rx) =e=")
                r.write(" yads*aax(pc,rx)+(1-yads)*1;\n")
                r.write("adyeq.. yads =l= 1;\n")
            else:
                r.write("adeq(pc,rx,time).. adsaux(pc,rx,time) =e=")
                r.write(" yads*aax(pc,rx,time)+(1-yads)*1;\n")
                r.write("adyeq.. yads =l= 1;\n")

        r.write("model ccmodel /all/;\n")
        # fix bad binaries
        #        if('clc' not in data['otherdata']):
        #            for i in range(data['rxnet']['nrx']):
        #                for j in range(len(data['rxnet']['mechanisms'])):
        #                    if ( i not in data['rxnet'][data['otherdata']['mechanisms'][j]]):
        #                        r.write("y.fx('"+str(i+1)+"','"+str(j+1)+"')=0;\n")

        #        r.write("option minlp=dicopt, reslim=1000,optca=1e-7,optcr=1e-7;\n")
        if 0:
            r.write("y.l('1','19')=1;\n")
            r.write("y.fx('1','20')=0;\n")
            r.write("y.fx('1','21')=0;\n")
            r.write("y.fx('1','18')=0;\n")
            r.write("y.fx('1','17')=0;\n")
            r.write("y.fx('1','16')=0;\n")
            r.write("y.fx('1','15')=0;\n")
        #            r.write("k.l('1','19')="+str(10926)+";\n")
        #            r.write("E.l('1','19')="+str(55319)+";\n")
        elif debug["true"]["fix"]:
            for i in range(1, len(debug["true"]["k"]) + 1):
                r.write("y.fx('" + str(i) + "','1')=1;\n")
                r.write(
                    "k.fx('"
                    + str(i)
                    + "','1')="
                    + str(debug["true"]["k"][i - 1])
                    + ";\n"
                )
                if data["otherdata"]["npc"] > 1:
                    r.write(
                        "E.fx('"
                        + str(i)
                        + "','1')="
                        + str(debug["true"]["e"][i - 1])
                        + ";\n"
                    )
        #    r.write("solve ccmodel minimizing modelcomp using minlp;\n")
        r.write("option minlp=baron, reslim=10000,optca=1e-7,optcr=1e-7;\n")
        r.write("solve ccmodel minimizing modelcomp using minlp;\n")
        if "time" in data["procdata"]:
            r.write("display estrate.l, rate;\n")
        else:
            r.write("display rate;\n")
        if "clc" in data["otherdata"]:
            # rates.put is estrate then rate_obs
            r.write("file rates;\n")
            r.write("put rates;\n")
            r.write("  loop(pc,loop(time,\n")
            r.write("  put ord(pc),ord(time),estrate.l(pc,time,'X'):20:10\n")
            r.write("  , rate(pc,time,'X'):20:10\n")
            r.write("  , rate(pc,time,'X'):20:10\n")
            r.write("/););\n")
        elif "T" not in data["otherdata"]["pclabels"]:
            r.write("parameter resid(s);\n")
            r.write(
                "resid(s) = sum(pc, power( rate(pc,s) - (sum(rx,(rmat('1',rx,s)*k.l(rx,'1')*(prod(s1,power(conc(pc,s1),pmat('1',rx,s1))))))),2));\n"
            )
            r.write("file rates;\n")
            r.write("put rates;\n")
            r.write("  loop(s,\n")
            r.write("  put resid(s):20:10);\n")
            r.write(" putclose rates;\n")
        else:
            r.write("parameter resid(s);\n")
            r.write("resid(s) = sum(pc, power( rate(pc,s) - estrate.l(pc,s),2));\n")
            r.write("file rates;\n")
            r.write("put rates;\n")
            r.write("  loop(s,\n")
            r.write("  put resid(s):20:10);\n")
            r.write(" putclose rates;\n")

        r.write("file res;\n")
        r.write("put res;\n")
        r.write("res.pc=4;\n")
        r.write("put modelcomp.l:20:10 / ;\n")
        r.write("  loop(rx,loop(m,\n")
        r.write("  put ord(rx),ord(m),y.l(rx,m):1:0\n")
        r.write("  , k.l(rx,m):20:10")
        if "T" in data["otherdata"]["pclabels"]:
            r.write(",E.l(rx,m):20:10")
        r.write("\n")
        r.write("/););\n")
        if 0:
            #        if( 'adata' in data['procdata']):
            if "time" in data["procdata"]:
                r.write("  loop(pc,loop(rx,loop(time,\n")
                r.write("  put ord(rx), adsaux.l(pc,rx,time):20:10/);););\n")
            else:
                r.write("  loop(pc,loop(rx,\n")
                r.write("  put ord(rx), adsaux.l(pc,rx):20:10/););\n")
        r.write("putclose res;\n")

    return aterm
