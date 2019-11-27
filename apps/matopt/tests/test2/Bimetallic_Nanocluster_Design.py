
import os
from math import sqrt

from matopt import *

# === MONOMETALIC OPT === 
Lat = FCCLattice(IAD=1.0)
Canv = Canvas()
Canv.addLocation(np.array([0,0,0],dtype=float))
Canv.addShells(1,Lat.getNeighbors)
Atoms = [Atom('Cu')]
N = 9
m = MatOptModel(Canv,Atoms)
Vals = [sqrt(CN) for CN in range(0,13)]
BPs = [CN for CN in range(0,13)]
m.addSitesDescriptor('CNRi',bounds=(0,sqrt(12)),integer=False,
                     rules=PiecewiseLinear(values=Vals,
                                           breakpoints=BPs,
                                           input_desc=m.Ci))
m.addGlobalDescriptor('Ecoh',rules=EqualTo(SumSites(desc=m.CNRi,
                                                    coefs=(1/(N*sqrt(12))))))
m.addGlobalDescriptor('Size',bounds=(N,N),
                      rules=EqualTo(SumSites(desc=m.Yi)))
D = m.maximize(m.Ecoh,tilim=10,solver='neos-cplex')

# === PROCESSING SOLUTION === 
Canv = Canvas()
for i in range(len(D)):
    if(D.Contents[i] is not None):
        Canv.addLocation(D.Canvas.Points[i])
Canv.setNeighborsFromFunc(Lat.getNeighbors)

Atoms = [Atom('Cu'),Atom('Ag')]
CompBounds = {Atom('Cu'):(4,4),
              Atom('Ag'):(5,5)}

m = MatOptModel(Canv,Atoms)

m.Yi.rules.append(FixedTo(1.0))
GklCoefs = {(Atom('Cu'),Atom('Cu')):3.520,
            (Atom('Cu'),Atom('Ag')):2.112,
            (Atom('Ag'),Atom('Ag')):2.580,
            (Atom('Ag'),Atom('Cu')):3.612}
BEijCoefs = {}
for i in range(len(Canv)):
    CNi = sum(1 for _ in Canv.NeighborhoodIndexes[i] if _ is not None)
    for j in Canv.NeighborhoodIndexes[i]:
        if(j is not None):
            CNj = sum(1 for _ in Canv.NeighborhoodIndexes[j] if _ is not None)
            for k in Atoms:
                for l in Atoms:
                    BEijCoefs[i,j,k,l] = GklCoefs[k,l]*1/sqrt(CNi) + GklCoefs[l,k]*1/sqrt(CNj)
m.addBondsDescriptor('BEij',
                     rules=EqualTo(SumBondTypes(m.Xijkl,coefs=BEijCoefs)),
                     symmetric_bonds=True)

m.addGlobalDescriptor('Ecoh',rules=EqualTo(SumBonds(desc=m.BEij,
                                                    coefs=1.0/(N*sqrt(12)))))
m.addGlobalTypesDescriptor('Composition',bounds=CompBounds,
                           rules=EqualTo(SumSites(desc=m.Yik)))

D = m.maximize(m.Ecoh,tilim=10,trelim=4096,solver='neos-cplex')
if(D is not None):
    D.toPDB('result.pdb')
