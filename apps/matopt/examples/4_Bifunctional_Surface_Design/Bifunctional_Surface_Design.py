
import numpy as np
from math import sqrt
from copy import deepcopy

from matopt import *


IAD = sqrt(2)*2
Lat = FCCLattice.alignedWith111(IAD)
nUnitCellsOnEdge = 4
nLayers = 4
a = nUnitCellsOnEdge*IAD
b = a
c = nLayers*Lat.FCC111LayerSpacing
alpha = np.pi/2
beta = np.pi/2
gamma = np.pi/3

S = Parallelepiped.fromEdgesAndAngles(a,b,c,alpha,beta,gamma)
S.shift(np.array([-0.01*a,-0.01*b,-0.01*c]))
T = PlanarTiling(S)

Canv = Canvas.fromLatticeAndTilingScan(Lat,T)

D = Design(Canv,Atom('Pt'))
D.toPDB('canvas.pdb')
D.toCFG('canvas.cfg',GS=1.0,BBox=S)

MotifCanvas = Canvas()
MotifCanvas.addLocation(np.array([0,0,0],dtype=float),NNeighbors=12)
MotifCanvas.addShell(Lat.getNeighbors)
Confs = [[None]*len(MotifCanvas.NeighborhoodIndexes[0]) for _ in range(7)]
iToSetNi = [[3,4,5,6,7,8],
            [3,4,5,6],
            [4,5,6,7],
            [5,6,7,8],
            [6,7,8,3],
            [7,8,3,4],
            [8,3,4,5]]
iToSetPt = [[9,10,11],
            [9,10,11],
            [9,10,11],
            [9,10,11],
            [9,10,11],
            [9,10,11],
            [9,10,11]]
for iConf,Conf in enumerate(Confs):
    for i in iToSetNi[iConf]:
        Conf[i] = Atom('Ni')
    for i in iToSetPt[iConf]:
        Conf[i] = Atom('Pt')

TypeAConfs = [0]
TypeBConfs = [1,2,3,4,5,6]
LocsToFixPt = [i for i in range(len(Canv)) if Canv.Points[i][2] < Lat.FCC111LayerSpacing*2.5]
LocsToExcludePt = [i for i in range(len(Canv)) if i not in LocsToFixPt]
CanvTwoBotLayers = [i for i in range(len(Canv)) if Canv.Points[i][2] < Lat.FCC111LayerSpacing*1.5]
CanvMinusTwoBotLayers = [i for i in range(len(Canv)) if i not in CanvTwoBotLayers]
OneLocToFix = [min(LocsToExcludePt)]
TileSizeSquared = nUnitCellsOnEdge**2
CatNorm = TileSizeSquared*6.0
UndefectedSurfE = 0.129758
maxSurfE = 999
CatWeight = 1.0
Atoms = [Atom('Ni'),Atom('Pt')]

m = MatOptModel(Canv,Atoms,Confs)

m.Yik.rules.append(FixedTo(1,sites=LocsToFixPt,site_types=[Atom('Pt')]))
m.Yik.rules.append(FixedTo(0,sites=LocsToExcludePt,site_types=[Atom('Pt')]))

m.Zic.rules.append(FixedTo(1,sites=OneLocToFix,confs=TypeAConfs))
m.Zic.rules.append(Implies(concs=(m.Yik,EqualTo(1,site_types=[Atom('Ni')]))))

SumAConfsExpr = SumConfs(m.Zic,confs_to_sum=TypeAConfs)
SumBConfsExpr = SumConfs(m.Zic,confs_to_sum=TypeBConfs)
m.addBondsDescriptor('SiteCombinations',binary=True,
                     rules=ImpliesSiteCombination(Canv,
                                                  (SumAConfsExpr,GreaterThan(1)),
                                                  (SumBConfsExpr,GreaterThan(1))))
m.addGlobalDescriptor('Activity',
                      rules=EqualTo(SumBonds(m.SiteCombinations,coefs=1/CatNorm)))

EiVals = [0, -0.04293*3+0.41492, -0.04293*10+0.41492, 0.05179*11-0.62148, 0]
EiBPs = [0, 3, 10, 11, 12]
m.addSitesDescriptor('Ei',
                     rules=PiecewiseLinear(values=EiVals,
                                           breakpoints=EiBPs,
                                           input_desc=m.Ci),
                     sites=CanvMinusTwoBotLayers)
m.addGlobalDescriptor('Esurf',
                      rules=EqualTo(SumSites(m.Ei,coefs=1/TileSizeSquared,offset=0.101208)))
m.addGlobalDescriptor('Stability',
                      rules=EqualTo(LinearExpr(m.Esurf,1/UndefectedSurfE)))

m.addGlobalDescriptor('ActAndStab',
                      rules=EqualTo(LinearExpr(descs=[m.Stability,m.Activity],
                                               coefs=[-(1-CatWeight),CatWeight])))


D = m.maximize(m.ActAndStab,tilim=360)
if(D is not None):
    D.toCFG('result.cfg',BBox=S)
    PeriodicD = T.replicateDesign(D,4)
    PeriodicS = deepcopy(S)
    PeriodicS.scale(np.array([4,4,1]))
    PeriodicD.toCFG('periodic_result.cfg',BBox=PeriodicS)
