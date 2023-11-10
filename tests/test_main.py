import numpy as np

import shapely.geometry as sg
import geopandas as gpd

from .context import geohexgrid, pytest, DATA_DIR
from geohexgrid import main as ghg


def test_axial_round():
    pass


def test_cartesian_to_axial():
    R = 3
    p = (1.7 * R, 0.2 * R)
    assert ghg.cartesian_to_axial(*p, R=R) == (1, 0)


def test_axial_to_double():
    pass


def test_cartesian_to_double():
    R = 3
    p = (1.7 * R, 0.2 * R)
    assert ghg.cartesian_to_double(*p, R) == (1, 1)


def test_make_grid_points():
    nrows = 3
    ncols = 5
    R = 3
    ox, oy = 2, 1
    X, Y = ghg.make_grid_points(nrows, ncols, R=R, ox=ox, oy=oy)
    assert X.shape == (3 * nrows, 2 * ncols)
    assert X[0][1] - X[0][0] == R
    assert Y[1][0] - Y[0][0] == ghg.K * R


def test_make_grid():
    R = 1
    grid = ghg.make_grid(2, 3, R=R, ox=2, oy=1)
    assert set(grid["cell_id"]) == {"0,0", "1,1", "2,0", "0,2", "1,3", "2,2"}
    # Should have no gaps
    assert grid.unary_union.boundary.is_ring


def test_make_grid_from_bounds():
    # A very good example
    rect = gpd.GeoDataFrame(
        geometry=[sg.Polygon([(2.1, -1), (4.9, -1), (4.9, 1.9), (2.1, 1.9)])]
    )
    grid = ghg.make_grid_from_bounds(*rect.total_bounds, R=1, ox=0, oy=0)
    assert set(grid["cell_id"]) == {
        "1,-3",
        "2,-2",
        "3,-3",
        "1,-1",
        "2,0",
        "3,-1",
        "1,1",
        "2,2",
        "3,1",
        "1,3",
        "2,4",
        "3,3",
    }


def test_make_grid_from_gdf():
    shape = gpd.GeoDataFrame(geometry=[sg.Polygon([(1, -1), (3, 1), (0, 3)])])
    grid = ghg.make_grid_from_gdf(shape, R=1, ox=0, oy=0, intersect=True)
    assert set(grid["cell_id"]) == {
        "1,-1",
        "0,0",
        "1,1",
        "2,0",
        "0,2",
        "1,3",
        "2,2",
        "0,4",
    }

    shapes = gpd.read_file(DATA_DIR / "shapes.geojson").to_crs("epsg:2193")
    R = 900
    grid1 = ghg.make_grid_from_gdf(shapes, R=R, intersect=False)
    assert grid1.crs == shapes.crs
    assert set(grid1.columns) == {"cell_id", "geometry"}

    grid2 = ghg.make_grid_from_gdf(shapes, R=R, intersect=True)
    assert grid2.shape[0] <= grid1.shape[0]
    assert grid2.area.sum() >= shapes.area.sum()

    grid3 = ghg.make_grid_from_gdf(shapes, R=R, clip=True)
    assert np.allclose(grid3.area.sum(), shapes.area.sum())
