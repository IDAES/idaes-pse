
import os
from math import sqrt

from matopt import *

Lat = FCCLattice(IAD=2.7704443686888935)
Canv = Canvas()
Canv.addLocation(np.array([0,0,0],dtype=float))
Canv.addShells(3,Lat.getNeighbors)
Canv.setNeighborsFromFunc(Lat.getNeighbors)
D = Design(Canv,Atom('Pt'))
D.toPDB('canvas.pdb')
# NOTE: If canvas.pdb was already created, we could have simply loaded it:
#Canv = Canvas.fromPDB('canvas.pdb',Lat=Lat)
Atoms = [Atom('Pt')]

N = 13
m = MatOptModel(Canv,Atoms)
m.addSitesDescriptor('CNRi',bounds=(0,sqrt(12)),integer=False,
                     rules=PiecewiseLinear(values=[sqrt(CN) for CN in range(13)],
                                           breakpoints=[CN for CN in range(13)],
                                           input_desc=m.Ci,con_type='UB'))
m.addGlobalDescriptor('Ecoh',rules=EqualTo(SumSites(desc=m.CNRi,
                                                    coefs=(1.0/sqrt(12)*1.0/N))))
m.addGlobalDescriptor('Size',bounds=(N,N),rules=EqualTo(SumSites(desc=m.Yi)))

D = m.maximize(m.Ecoh,tilim=10)
D.toPDB('result.pdb')
