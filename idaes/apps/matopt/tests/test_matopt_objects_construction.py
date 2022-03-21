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
import numpy as np
from math import sqrt
from idaes.apps.matopt.materials import (
    Atom,
    Canvas,
    Design,
    LinearTiling,
    PlanarTiling,
    CubicTiling,
)
from idaes.apps.matopt.materials.geometry import (
    Shape,
    Cuboctahedron,
    Parallelepiped,
    Rhombohedron,
    RectPrism,
    Cube,
    Cylinder,
    CylindricalSector,
)
from idaes.apps.matopt.materials.lattices import (
    FCCLattice,
    CubicLattice,
    PerovskiteLattice,
    DiamondLattice,
    WurtziteLattice,
)
from idaes.apps.matopt.materials.transform_func import (
    ShiftFunc,
    ScaleFunc,
    RotateFunc,
    ReflectFunc,
)
from idaes.apps.matopt.opt import (
    Coef,
    LinearExpr,
    MatOptModel,
    SumNeighborSites,
    SumNeighborBonds,
    SumSites,
    SumBonds,
    SumSiteTypes,
    SumBondTypes,
    SumSitesAndTypes,
    SumBondsAndTypes,
    SumConfs,
    SumSitesAndConfs,
    LessThan,
    EqualTo,
    GreaterThan,
    FixedTo,
    Disallow,
    PiecewiseLinear,
    Implies,
    NegImplies,
    ImpliesSiteCombination,
    ImpliesNeighbors,
)
import pytest


@pytest.mark.unit
def test_construct_FCCLattice():
    IAD = sqrt(2) / 2
    lattice = FCCLattice.alignedWith100(IAD)
    lattice = FCCLattice.alignedWith111(IAD)
    lattice = FCCLattice(IAD)
    return lattice


@pytest.mark.unit
def test_construct_CubicLattice():
    IAD = 1.0
    lattice = CubicLattice(IAD)
    return lattice


@pytest.mark.unit
def test_construct_PerovskiteLattice():
    A, B, C = 1.0, 1.0, 1.0
    lattice = PerovskiteLattice(A, B, C)
    return lattice


@pytest.mark.unit
def test_construct_DiamondLattice():
    IAD = sqrt(3) / 4
    lattice = DiamondLattice.alignedWith(IAD, "100")
    lattice = DiamondLattice.alignedWith(IAD, "110")
    lattice = DiamondLattice.alignedWith(IAD, "111")
    lattice = DiamondLattice(IAD)
    return lattice


@pytest.mark.unit
def test_construct_WurtziteLattice():
    IAD = sqrt(3 / 8)
    lattice = WurtziteLattice.alignedWith(IAD, "0001")
    lattice = WurtziteLattice(IAD)
    return lattice


@pytest.mark.unit
def test_construct_Atom():
    atom_zero = Atom(39)
    atom_one = Atom(40)
    atom_two = Atom("Zr")
    assert atom_one == atom_two
    return atom_zero, atom_one


@pytest.mark.unit
def test_construct_Canvas():
    canvas = Canvas()
    return canvas


@pytest.mark.unit
def test_construct_Design():
    lattice = test_construct_FCCLattice()
    canvas = test_construct_Canvas()
    atom, _ = test_construct_Atom()
    canvas.addLocation(np.array([0, 0, 0], dtype=float))
    canvas.addShell(lattice.getNeighbors)
    canvas.setNeighborsFromFunc(lattice.getNeighbors)
    design = Design(canvas, atom)
    return design


@pytest.mark.unit
def test_construct_Shape():
    shape = Shape(np.array([0, 0, 0], dtype=float))
    return shape


@pytest.mark.unit
def test_construct_Cuboctahedron():
    shape = Cuboctahedron(1.0)
    return shape


@pytest.mark.unit
def test_construct_Parallelepiped():
    shape = Parallelepiped(
        np.array([1, 0, 0], dtype=float),
        np.array([0, 1, 0], dtype=float),
        np.array([0, 0, 1], dtype=float),
    )
    return shape


@pytest.mark.unit
def test_construct_Rhombohedron():
    shape = Rhombohedron(1.0, 0.5)
    return shape


@pytest.mark.unit
def test_construct_RectPrism():
    shape = RectPrism(1.0, 1.0, 1.0)
    return shape


@pytest.mark.unit
def test_construct_Cube():
    shape = Cube(1.0)
    return shape


@pytest.mark.unit
def test_construct_Cylinder():
    shape = Cylinder(np.zeros(3, dtype=float), 1.0, 1.0)
    return shape


@pytest.mark.unit
def test_construct_CylindricalSector():
    shape = CylindricalSector(
        np.zeros(3, dtype=float),
        1.0,
        1.0,
        np.array([1, 0, 0], dtype=float),
        np.array([0, 1, 0], dtype=float),
        np.array([0, 0, 1], dtype=float),
    )
    return shape


@pytest.mark.unit
def test_construct_LinearTiling():
    parallelepiped = test_construct_Parallelepiped()
    cylinder = test_construct_Cylinder()
    tiling = LinearTiling.fromParallelepiped(parallelepiped)
    tiling = LinearTiling.fromCylindricalShape(cylinder)
    return tiling


@pytest.mark.unit
def test_construct_PlanarTiling():
    shape = test_construct_Parallelepiped()
    tiling = PlanarTiling(shape)
    return tiling


@pytest.mark.unit
def test_construct_CubicTiling():
    shape = test_construct_Parallelepiped()
    tiling = CubicTiling(shape)
    return tiling


@pytest.mark.unit
def test_construct_ShiftFunc():
    transformation = ShiftFunc(np.array([1, 0, 0], dtype=float))
    return transformation


@pytest.mark.unit
def test_construct_ScaleFunc():
    transformation = ScaleFunc(2.0)
    return transformation


@pytest.mark.unit
def test_construct_RotateFunc():
    transformation = RotateFunc.fromXYZAngles(0.5, 0.5, 0.5)
    return transformation


@pytest.mark.unit
def test_construct_ReflectFunc():
    transformation = ReflectFunc.acrossX()
    return transformation


@pytest.mark.unit
def test_construct_Coef():
    coefficient = Coef([0.5, 0.5])
    return coefficient


@pytest.mark.unit
def test_construct_LinearExpr():
    expression = LinearExpr()
    return expression


@pytest.mark.unit
def test_construct_SumNeighborSites():
    m = test_construct_MatOptModel()
    expression = SumNeighborSites(m.Yi)
    return expression


@pytest.mark.unit
def test_construct_SumNeighborBonds():
    m = test_construct_MatOptModel()
    expression = SumNeighborBonds(m.Xij)
    return expression


@pytest.mark.unit
def test_construct_SumSites():
    m = test_construct_MatOptModel()
    expression = SumSites(m.Yi)
    return expression


@pytest.mark.unit
def test_construct_SumBonds():
    m = test_construct_MatOptModel()
    expression = SumBonds(m.Xij)
    return expression


@pytest.mark.unit
def test_construct_SumSiteTypes():
    m = test_construct_MatOptModel()
    expression = SumSiteTypes(m.Yik)
    return expression


@pytest.mark.unit
def test_construct_SumBondTypes():
    m = test_construct_MatOptModel()
    expression = SumBondTypes(m.Xijkl)
    return expression


@pytest.mark.unit
def test_construct_SumSitesAndTypes():
    m = test_construct_MatOptModel()
    expression = SumSitesAndTypes(m.Yik)
    return expression


@pytest.mark.unit
def test_construct_SumBondsAndTypes():
    m = test_construct_MatOptModel()
    expression = SumBondsAndTypes(m.Xijkl)
    return expression


@pytest.mark.unit
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


@pytest.mark.unit
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


@pytest.mark.unit
def test_construct_LessThan():
    rule = LessThan(2)
    return rule


@pytest.mark.unit
def test_construct_EqualTo():
    rule = EqualTo(2)
    return rule


@pytest.mark.unit
def test_construct_GreaterThan():
    rule = GreaterThan(2)
    return rule


@pytest.mark.unit
def test_construct_FixedTo():
    rule = FixedTo(2)
    return rule


@pytest.mark.unit
def test_construct_Disallow():
    design = test_construct_Design()
    rule = Disallow(design)
    return rule


@pytest.mark.unit
def test_construct_PiecewiseLinear():
    m = test_construct_MatOptModel()
    rule = PiecewiseLinear([x**2 for x in range(3)], [x for x in range(3)], m.Ci)
    return rule


@pytest.mark.unit
def test_construct_Implies():
    m = test_construct_MatOptModel()
    r = test_construct_LessThan()
    rule = Implies((m.Yi, r))
    return rule


@pytest.mark.unit
def test_construct_NegImplies():
    m = test_construct_MatOptModel()
    r = test_construct_LessThan()
    rule = NegImplies((m.Yi, r))
    return rule


@pytest.mark.unit
def test_construct_ImpliesSiteCombination():
    canvas = test_construct_Canvas()
    m = test_construct_MatOptModel()
    r = test_construct_GreaterThan()
    rule = ImpliesSiteCombination(canvas, (m.Yi, r), (m.Yik, r))
    return rule


@pytest.mark.unit
def test_construct_ImpliesNeighbors():
    m = test_construct_MatOptModel()
    r = test_construct_FixedTo()
    rule = ImpliesNeighbors((m.Yi, r))
    return rule


@pytest.mark.unit
def test_construct_MatOptModel():
    canvas = test_construct_Canvas()
    a0, a1 = test_construct_Atom()
    atoms = [a0, a1]
    model = MatOptModel(canvas, atoms)
    return model
