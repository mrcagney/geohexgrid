import pytest

from .context import geohex
from geohex import *

import numpy as np


hgs = HexGridSystem("epsg:2193", 0, 0, 50)


# Test helper functions
def test_axial_to_cartesian():
    get = axial_to_cartesian(0, 0, 5)
    expect = (0, 0)
    assert np.allclose(get, expect)

    R = 5
    get = axial_to_cartesian(1, 0, R)
    expect = (R * 3 / 2, R * sqrt(3) / 2)
    assert np.allclose(get, expect)

    get = axial_to_cartesian(0, 1, R)
    expect = (0, R * sqrt(3))
    assert np.allclose(get, expect)


def test_cartesian_to_axial():
    R = 5
    for x, y in [(1, 1), (-3, -1)]:
        get = axial_to_cartesian(*cartesian_to_axial(x, y, R), R)
        expect = x, y
        assert np.allclose(get, expect)


def test_hexagon_vertex():
    R = 3
    assert np.allclose(hexagon_vertex(0, 0, R, 0), (R, 0))
    assert np.allclose(hexagon_vertex(0, 0, R, 3), (-R, 0))


# Test Cell methods
def test_center():
    assert Cell(hgs, 0, 0).center() == (0, 0)


def test_vertices():
    v = Cell(hgs, 0, 0).vertices()
    assert len(set(v)) == 6


def test_polygon():
    p = Cell(hgs, 1, 2).polygon()
    assert isinstance(p, sg.Polygon)
    # Polygon should have correct area
    assert np.allclose(p.area, 3 * sqrt(3) / 2 * hgs.R**2)


def test_neighbour():
    c = Cell(hgs, 0, 0)
    assert c.neighbour(2) == c.neighbour(-4)
    assert c.neighbour(0) == c.neighbour("ru")
    assert len(set(c.neighbor(i) for i in range(6))) == 6

    with pytest.raises(ValueError):
        c.neighbour("fail")
    with pytest.raises(ValueError):
        c.neighbour(0.8)


# Test HexGridSystem methods
def test_cell_from_axial_point():
    c = hgs.cell_from_axial_point(0.51, 0.1)
    assert isinstance(c, Cell)
    assert (c.a, c.b) == (1, 0)


def test_cell_from_point():
    R = hgs.R
    c = hgs.cell_from_point(1.1 * R, -0.1 * R)
    assert isinstance(c, Cell)
    assert (c.a, c.b) == (1, -1)


def test_grid_from_bbox():
    R = hgs.R

    # One cell
    grid = hgs.grid_from_bbox(0, 0, 0.1 * R, 0.1 * R)
    assert [c.center_axial() for c in grid] == [(0, 0)]

    # Several grid
    grid = hgs.grid_from_bbox(-0.1 * R, -0.1 * R, 0.1 * R, R)
    assert set(c.center_axial() for c in grid) == {(0, 0), (0, 1)}

    grid = hgs.grid_from_bbox(0, 0, 1.1 * R, 0.1 * R)
    assert set(c.center_axial() for c in grid) == {(0, 0), (1, 0)}

    grid = hgs.grid_from_bbox(0, 0, R, R)
    assert set(c.center_axial() for c in grid) == {(0, 0), (1, 0), (0, 1)}

    grid = hgs.grid_from_bbox(0, -0.1 * R, 1.1 * R, 0.1 * R)
    assert set(c.center_axial() for c in grid) == {(0, 0), (1, 0), (1, -1)}

    grid = hgs.grid_from_bbox(0, -0.1 * R, 1.1 * R, 0.1 * R)
    assert set(c.center_axial() for c in grid) == {(0, 0), (1, 0), (1, -1)}

    grid = hgs.grid_from_bbox(0, -0.1 * R, R, R)
    assert set(c.center_axial() for c in grid) == {(1, -1), (0, 0), (1, 0), (0, 1)}

    grid = hgs.grid_from_bbox(-R, 0, R, R)
    assert set(c.center_axial() for c in grid) == {(0, 0), (-1, 1), (1, 0), (0, 1)}

    grid = hgs.grid_from_bbox(-R, -0.1 * R, R, R)
    assert set(c.center_axial() for c in grid) == {
        (0, 0),
        (-1, 1),
        (1, 0),
        (0, 1),
        (-1, 0),
        (1, -1),
    }
