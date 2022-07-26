from .context import geohex
from geohex import *

import numpy as np


hgs = HexGridSystem("epsg:2193", 0, 0, 50)


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


def test_cell_from_axial_point():
    c = hgs.cell_from_axial_point(0.51, 0.1)
    assert isinstance(c, Cell)
    assert (c.a, c.b) == (1, 0)


def test_cell_from_point():
    R = hgs.R
    c = hgs.cell_from_point(1.1 * R, -0.1 * R)
    assert isinstance(c, Cell)
    assert (c.a, c.b) == (1, -1)


def test_cells_from_bbox():
    R = hgs.R

    # One cell
    cells = hgs.cells_from_bbox(0, 0, R / 2, R / 2)
    assert [c.axial_center() for c in cells] == [(0, 0)]

    # Several cells
    cells = hgs.cells_from_bbox(0, -0.1 * R, 1.1 * R, 0.25 * R)
    assert set(c.axial_center() for c in cells) == {(0, 0), (1, 0), (1, -1)}
