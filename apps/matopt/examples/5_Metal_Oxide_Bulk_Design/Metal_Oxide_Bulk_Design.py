
from matopt import *

A = 4.0
B = 4.0
C = 4.0
Lat = PerovskiteLattice(A,B,C)

nUnitCellsOnEdge = 2
S = RectPrism(nUnitCellsOnEdge*A,
              nUnitCellsOnEdge*B,
              nUnitCellsOnEdge*C)
S.shift(np.array([-0.01,-0.01,-0.01]))
T = CubicTiling(S)

Canv = Canvas.fromLatticeAndTilingScan(Lat,T)
Atoms = [Atom('Ba'),Atom('Fe'),Atom('In'),Atom('O')]

iDesiredConfs = [394,395,396,397,398,399,400,401,68,69,
                 70,71,162,163,164,165,166,167,168,169]
ConfDesigns = loadFromPDBs([str(i)+'.pdb' for i in iDesiredConfs],folder='./Confs')
Confs = [Conf.Contents for Conf in ConfDesigns]

Sites = [i for i in range(len(Canv))]
ASites = [i for i in Sites if Lat.isASite(Canv.Points[i])]
BSites = [i for i in Sites if Lat.isBSite(Canv.Points[i])]
OSites = [i for i in Sites if Lat.isOSite(Canv.Points[i])]
pctLocalLB,pctLocalUB = 0,1
pctGlobalLB,pctGlobalUB = 0.0,0.3
LocalBounds = {(i,Atom('In')):(round(pctLocalLB*len(Canv.NeighborhoodIndexes[i])),
                               round(pctLocalUB*len(Canv.NeighborhoodIndexes[i]))) for i in OSites}
GlobalLB = round(pctGlobalLB*len(BSites))
GlobalUB = round(pctGlobalUB*len(BSites))

m = MatOptModel(Canv,Atoms,Confs)

m.Yik.rules.append(FixedTo(1,sites=ASites,site_types=[Atom('Ba')]))
m.Yik.rules.append(FixedTo(1,sites=OSites,site_types=[Atom('O')]))
m.Yik.rules.append(FixedTo(0,sites=BSites,site_types=[Atom('Ba'),Atom('O')]))
m.Yi.rules.append(FixedTo(1,sites=BSites))

m.addGlobalDescriptor('Activity',
                      rules=EqualTo(SumSitesAndConfs(m.Zic,coefs=1/len(OSites),sites_to_sum=OSites)))
m.addGlobalTypesDescriptor('GlobalBudget',bounds=(GlobalLB,GlobalUB),
                           rules=EqualTo(SumSites(m.Yik,
                                                  site_types=[Atom('In')],
                                                  sites_to_sum=BSites)))
m.addSitesTypesDescriptor('LocalBudget',bounds=LocalBounds,
                          rules=EqualTo(SumNeighborSites(m.Yik,
                                                         sites=OSites,
                                                         site_types=[Atom('In')])))


D = m.maximize(m.Activity,tilim=360)

if(D is not None):
    for i,c in m.Zic.keys():
        if(m.Zic.values[i,c] > 0.5):
            D.setContent(i,Atom('S'))
    D.toCFG('result.cfg',BBox=S)
