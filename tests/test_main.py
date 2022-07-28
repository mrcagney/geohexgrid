import plum
import numpy as np

from .context import geohex, pytest, DATA_DIR
from geohex import *


GHS = GeoHexSystem("epsg:2193", 50, 0, 0)


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
    assert Cell(GHS, 0, 0).center() == (0, 0)


def test_vertices():
    v = Cell(GHS, 0, 0).vertices()
    assert len(set(v)) == 6


def test_polygon():
    p = Cell(GHS, 1, 2).polygon()
    assert isinstance(p, sg.Polygon)
    # Polygon should have correct area
    assert np.allclose(p.area, 3 * sqrt(3) / 2 * GHS.R**2)


def test_neighbour():
    c = Cell(GHS, 0, 0)
    assert c.neighbour(2) == c.neighbour(-4)
    assert c.neighbour(0) == c.neighbour("ru")
    assert len(set(c.neighbor(i) for i in range(6))) == 6

    with pytest.raises(ValueError):
        c.neighbour("fail")
    with pytest.raises(plum.function.NotFoundLookupError):
        c.neighbour(0.8)


# Test GeoHexSystem methods
def test_cell_from_axial_point():
    c = GHS.cell_from_axial_point(0.51, 0.1)
    assert isinstance(c, Cell)
    assert (c.a, c.b) == (1, 0)


def test_cell_from_point():
    R = GHS.R
    c = GHS.cell_from_point(1.1 * R, -0.1 * R)
    assert isinstance(c, Cell)
    assert (c.a, c.b) == (1, -1)


def test_grid_from_bbox():
    R = GHS.R

    # One cell
    grid = GHS.grid_from_bbox(0, 0, 0, 0, as_gdf=False)
    assert [c.center_axial() for c in grid] == [(0, 0)]

    grid = GHS.grid_from_bbox(0, 0, 0.1 * R, 0.1 * R, as_gdf=False)
    assert [c.center_axial() for c in grid] == [(0, 0)]

    # Several grid
    grid = GHS.grid_from_bbox(-0.1 * R, -0.1 * R, 0.1 * R, R, as_gdf=False)
    assert set(c.center_axial() for c in grid) == {(0, 0), (0, 1)}

    grid = GHS.grid_from_bbox(0, 0, 1.1 * R, 0.1 * R, as_gdf=False)
    assert set(c.center_axial() for c in grid) == {(0, 0), (1, 0)}

    grid = GHS.grid_from_bbox(0, 0, R, R, as_gdf=False)
    assert set(c.center_axial() for c in grid) == {(0, 0), (1, 0), (0, 1)}

    grid = GHS.grid_from_bbox(0, -0.1 * R, 1.1 * R, 0.1 * R, as_gdf=False)
    assert set(c.center_axial() for c in grid) == {(0, 0), (1, 0), (1, -1)}

    grid = GHS.grid_from_bbox(0, -0.1 * R, 1.1 * R, 0.1 * R, as_gdf=False)
    assert set(c.center_axial() for c in grid) == {(0, 0), (1, 0), (1, -1)}

    grid = GHS.grid_from_bbox(0, -0.1 * R, R, R, as_gdf=False)
    assert set(c.center_axial() for c in grid) == {(1, -1), (0, 0), (1, 0), (0, 1)}

    grid = GHS.grid_from_bbox(-R, 0, R, R, as_gdf=False)
    assert set(c.center_axial() for c in grid) == {(0, 0), (-1, 1), (1, 0), (0, 1)}

    grid0 = GHS.grid_from_bbox(-R, -0.1 * R, R, R, as_gdf=False)
    assert set(c.center_axial() for c in grid0) == {
        (0, 0),
        (-1, 1),
        (1, 0),
        (0, 1),
        (-1, 0),
        (1, -1),
    }

    grid1 = GHS.grid_from_bbox(-R, -0.1 * R, R, R, as_gdf=True)
    assert grid1.crs == GHS.crs
    assert set(grid1.columns) == {"cell_id", "geometry"}
    assert set(grid1.cell_id) == set(c.id() for c in grid0)


def test_grid_from_gdf():
    ghs = GeoHexSystem("epsg:2193", 250, 0, 0)

    shapes = gpd.read_file(DATA_DIR / "shapes.geojson")

    grid = ghs.grid_from_gdf(shapes, as_gdf=False)
    assert isinstance(grid[0], Cell)

    grid1 = ghs.grid_from_gdf(shapes, as_gdf=True)
    assert grid1.crs == ghs.crs
    assert set(grid1.columns) == {"cell_id", "geometry"}

    grid2 = ghs.grid_from_gdf(shapes, as_gdf=True, intersect=True)
    assert grid2.crs == ghs.crs
    assert grid2.shape[0] <= grid1.shape[0]


def test_geohexgrid():
    shapes = gpd.read_file(DATA_DIR / "shapes.geojson")

    grid1 = geohexgrid(shapes, 0.01)
    assert grid1.crs == shapes.crs

    grid2 = geohexgrid(shapes, 0.01, intersect=True)
    assert grid2.crs == shapes.crs
    assert grid2.shape[0] <= grid1.shape[0]
