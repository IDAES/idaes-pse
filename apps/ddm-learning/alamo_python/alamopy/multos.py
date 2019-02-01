##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#                                                                          
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
# These subroutines ensure compatibility across linux and windows operating systems
import os
import platform

def deletefile(*fname):
    tos = platform.platform()
    if 'Windows' in tos:
        for name in fname:
            os.system("del "+name)
    else:
        for name in fname:
            os.system("rm "+name)

def movefile(*fname):
    tos = platform.platform()
    if 'Windows' in tos:
        for name in fname:
            os.system("move "+name)
    else:
        for name in fname:
            os.system("mv "+name)

def copyfile(outf,inf):
    tos = platform.platform()
    if 'Windows' in tos:
        os.system("copy "+inf+' '+outf)
    else:
        os.system("cp "+inf+' '+outf)


def catfile(outf,*fname):
    tos = platform.platform()
    if 'Windows' in tos:
        ostr = ''
        for name in fname:
            ostr += "echo "+name+' >> '+outf+' & '
        ostr = ostr[:-3]
        os.system(ostr)
    else:
        ostr = 'cat '
        for name in fname:
            ostr += name+' '
        ostr+='> '+outf
        os.system(ostr)

