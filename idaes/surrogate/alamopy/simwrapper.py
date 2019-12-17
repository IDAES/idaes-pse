#!/usr/bin/python
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
def main():
    from idaes.surrogate.alamopy import examples
    infile = 'input.txt'
    outfile = 'output.txt'
    fin = open(infile, 'r')
    fout = open(outfile, 'w')
    newline = fin.readline()
    newlist = newline.split()
    n = int(newlist[0])
    for p in range(0, n):
        newline = fin.readline()
        newlist = newline.split()
        ninputs = len(newlist)
        x = [0] * (ninputs + 1)
        for k in range(0, ninputs):
            x[k] = float(newlist[k])
        x[ninputs] = examples.sixcamel(*x[:-1])
        for k in range(0, len(x)):
            fout.write(str(x[k]) + ' ')
        fout.write(' \n')


if __name__ == '__main__':
    main()
