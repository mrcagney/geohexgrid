import numpy as np

from .context import geohex, pytest, DATA_DIR
from geohex import *


R = 50
GHS = GeoHexSystem("epsg:2193", R, 0, 0)


# Test helper functions
def test_axial_center_to_id():
    assert axial_center_to_id(26, 530) == "26-530"

def test_cell_id_to_axial_center():
    a, b = 2, 3
    assert id_to_axial_center(axial_center_to_id(a, b)) == (a, b)

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
    assert Cell(0, 0, R).center() == (0, 0)


def test_vertices():
    v = Cell(0, 0, R).vertices()
    assert len(set(v)) == 6


def test_polygon():
    p = Cell(1, 2, R).polygon()
    assert isinstance(p, sg.Polygon)
    # Polygon should have correct area
    assert np.allclose(p.area, 3 * sqrt(3) / 2 * R**2)


def test_neighbour():
    c = Cell(0, 0, R)
    neighbours = ["N", "NW", "SW", "S", "SE", "NE"]
    assert len(set(c.neighbor(i) for i in neighbours))== 6

    with pytest.raises(ValueError):
        c.neighbour("fail")


# Test GeoHexSystem methods
def test_cell_from_axial_point():
    c = GHS.cell_from_axial_point(0.51, 0.1)
    assert isinstance(c, Cell)
    assert (c.a, c.b) == (1, 0)


def test_cell_from_point():
    c = GHS.cell_from_point(1.1 * R, -0.1 * R)
    assert isinstance(c, Cell)
    assert (c.a, c.b) == (1, -1)


def test_grid_from_bbox():
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


# Test main function
def test_make_grid():
    shapes = gpd.read_file(DATA_DIR / "shapes.geojson")

    grid1 = make_grid(shapes, 0.01)
    assert grid1.crs == shapes.crs

    grid2 = make_grid(shapes, 0.01, intersect=True)
    assert grid2.crs == shapes.crs
    assert grid2.shape[0] <= grid1.shape[0]
