##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
from idaes.surrogate import ripe
import numpy as np
from idaes.surrogate.ripe import mechs as mechs

def main():
    spec = ['X']
    # Import data from csv
    data = np.genfromtxt('clc.csv', delimiter=',')
    t = data[:,0]
    xdata = data[:,1]
    stoich = [1]

    # User pre-defined clc rate forms found in RIPE
    # mechs = ripe.clcforms
    clc_mechs = [mechs.powerlawp5, mechs.powerlaw2, mechs.powerlaw3, mechs.powerlaw4, mechs.avrami2, mechs.avrami3, mechs.avrami4, mechs.avrami5, mechs.randomnuc, mechs.ptompkins, mechs.jander, mechs.antijander, mechs.valensi, mechs.parabolic, mechs.gb3d, mechs.zlt, mechs.grain]

    # Identify optimal kinetic mechanism
    results = ripe.ripemodel(xdata,stoichiometry=stoich,mechanisms=clc_mechs,time=t)


if __name__ == "__main__":
    main()
