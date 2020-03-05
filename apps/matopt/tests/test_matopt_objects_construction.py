import numpy as np

from apps.matopt.materials import Atom, Canvas, Design, PlanarTiling, CubicTiling
from apps.matopt.materials.geometry import Shape, Cuboctahedron, Parallelepiped, Rhombohedron, RectPrism, Cube
from apps.matopt.materials.lattices import FCCLattice, CubicLattice, PerovskiteLattice
from apps.matopt.materials.transform_func import ShiftFunc, ScaleFunc, RotateFunc, ReflectFunc
from apps.matopt.opt import Coef, LinearExpr, MatOptModel, SumNeighborSites, SumNeighborBonds, SumSites, SumBonds, \
    SumSiteTypes, SumBondTypes, SumSitesAndTypes, SumBondsAndTypes, SumConfs, SumSitesAndConfs, LessThan, EqualTo, \
    GreaterThan, FixedTo, Disallow, PiecewiseLinear, Implies, NegImplies, ImpliesSiteCombination, ImpliesNeighbors


def test_construct_FCCLattice():
    IAD = 2.77
    lattice = FCCLattice(IAD)
    return lattice


def test_construct_CubicLattice():
    IAD = 2.77
    lattice = CubicLattice(IAD)
    return lattice


def test_construct_PerovskiteLattice():
    A, B, C = 4, 4, 4
    lattice = PerovskiteLattice(A, B, C)
    return lattice


def test_construct_Atom():
    atom_zero = Atom(39)
    atom_one = Atom(40)
    atom_two = Atom('Zr')
    assert (atom_one == atom_two)
    return atom_zero, atom_one


def test_construct_Canvas():
    canvas = Canvas()
    return canvas


def test_construct_Design():
    lattice = test_construct_FCCLattice()
    canvas = test_construct_Canvas()
    atom, _ = test_construct_Atom()
    canvas.addLocation(np.array([0, 0, 0], dtype=float))
    canvas.addShell(lattice.getNeighbors)
    canvas.setNeighborsFromFunc(lattice.getNeighbors)
    design = Design(canvas, atom)
    return design


def test_construct_Shape():
    shape = Shape(np.array([0, 0, 0], dtype=float))
    return shape


def test_construct_Cuboctahedron():
    shape = Cuboctahedron(1.0)
    return shape


def test_construct_Parallelepiped():
    shape = Parallelepiped(np.array([1, 0, 0], dtype=float),
                           np.array([0, 1, 0], dtype=float),
                           np.array([0, 0, 1], dtype=float))
    return shape


def test_construct_Rhombohedron():
    shape = Rhombohedron(1.0, 0.5)
    return shape


def test_construct_RectPrism():
    shape = RectPrism(1.0, 1.0, 1.0)
    return shape


def test_construct_Cube():
    shape = Cube(1.0)
    return shape


def test_construct_PlanarTiling():
    shape = test_construct_Parallelepiped()
    tiling = PlanarTiling(shape)
    return tiling


def test_construct_CubicTiling():
    shape = test_construct_Parallelepiped()
    tiling = CubicTiling(shape)
    return tiling


def test_construct_ShiftFunc():
    transformation = ShiftFunc(np.array([1, 0, 0], dtype=float))
    return transformation


def test_construct_ScaleFunc():
    transformation = ScaleFunc(2.0)
    return transformation


def test_construct_RotateFunc():
    transformation = RotateFunc.fromXYZAngles(0.5, 0.5, 0.5)
    return transformation


def test_construct_ReflectFunc():
    transformation = ReflectFunc.acrossX()
    return transformation


def test_construct_Coef():
    coefficient = Coef([0.5, 0.5])
    return coefficient


def test_construct_LinearExpr():
    expression = LinearExpr()
    return expression


def test_construct_SumNeighborSites():
    m = test_construct_MatOptModel()
    expression = SumNeighborSites(m.Yi)
    return expression


def test_construct_SumNeighborBonds():
    m = test_construct_MatOptModel()
    expression = SumNeighborBonds(m.Xij)
    return expression


def test_construct_SumSites():
    m = test_construct_MatOptModel()
    expression = SumSites(m.Yi)
    return expression


def test_construct_SumBonds():
    m = test_construct_MatOptModel()
    expression = SumBonds(m.Xij)
    return expression


def test_construct_SumSiteTypes():
    m = test_construct_MatOptModel()
    expression = SumSiteTypes(m.Yik)
    return expression


def test_construct_SumBondTypes():
    m = test_construct_MatOptModel()
    expression = SumBondTypes(m.Xijkl)
    return expression


def test_construct_SumSitesAndTypes():
    m = test_construct_MatOptModel()
    expression = SumSitesAndTypes(m.Yik)
    return expression


def test_construct_SumBondsAndTypes():
    m = test_construct_MatOptModel()
    expression = SumBondsAndTypes(m.Xijkl)
    return expression


def test_construct_SumConfs():
    lattice = test_construct_FCCLattice()
    canvas = test_construct_Canvas()
    canvas.addLocation(np.array([0, 0, 0], dtype=float))
    canvas.addShell(lattice.getNeighbors)
    canvas.setNeighborsFromFunc(lattice.getNeighbors)
    confs = [[None] * len(canvas.NeighborhoodIndexes[0]) for _ in range(1)]
    a0, a1 = test_construct_Atom()
    atoms = [a0, a1]
    m = MatOptModel(canvas, atoms, confs)
    expression = SumConfs(m.Zic)
    return expression


def test_construct_SumSitesAndConfs():
    lattice = test_construct_FCCLattice()
    canvas = test_construct_Canvas()
    canvas.addLocation(np.array([0, 0, 0], dtype=float))
    canvas.addShell(lattice.getNeighbors)
    canvas.setNeighborsFromFunc(lattice.getNeighbors)
    confs = [[None] * len(canvas.NeighborhoodIndexes[0]) for _ in range(1)]
    a0, a1 = test_construct_Atom()
    atoms = [a0, a1]
    m = MatOptModel(canvas, atoms, confs)
    expression = SumSitesAndConfs(m.Zic)
    return expression


def test_construct_LessThan():
    rule = LessThan(2)
    return rule


def test_construct_EqualTo():
    rule = EqualTo(2)
    return rule


def test_construct_GreaterThan():
    rule = GreaterThan(2)
    return rule


def test_construct_FixedTo():
    rule = FixedTo(2)
    return rule


def test_construct_Disallow():
    design = test_construct_Design()
    rule = Disallow(design)
    return rule


def test_construct_PiecewiseLinear():
    m = test_construct_MatOptModel()
    rule = PiecewiseLinear([x**2 for x in range(3)],
                           [x for x in range(3)],
                           m.Ci)
    return rule


def test_construct_Implies():
    m = test_construct_MatOptModel()
    r = test_construct_LessThan()
    rule = Implies((m.Yi, r))
    return rule


def test_construct_NegImplies():
    m = test_construct_MatOptModel()
    r = test_construct_LessThan()
    rule = NegImplies((m.Yi, r))
    return rule


def test_construct_ImpliesSiteCombination():
    canvas = test_construct_Canvas()
    m = test_construct_MatOptModel()
    r = test_construct_GreaterThan()
    rule = ImpliesSiteCombination(canvas, (m.Yi, r), (m.Yik, r))
    return rule


def test_construct_ImpliesNeighbors():
    m = test_construct_MatOptModel()
    r = test_construct_FixedTo()
    rule = ImpliesNeighbors((m.Yi, r))
    return rule


def test_construct_MatOptModel():
    canvas = test_construct_Canvas()
    a0, a1 = test_construct_Atom()
    atoms = [a0, a1]
    model = MatOptModel(canvas, atoms)
    return model
