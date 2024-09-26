import math

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
    x0, y0 = 2, 1
    X, Y = ghg.make_grid_points(nrows, ncols, R=R, x0=x0, y0=y0)
    assert X.shape == (2 * nrows + 2, math.ceil((3 * ncols + 2) / 2))
    assert X[0][1] - X[0][0] == R
    assert Y[1][0] - Y[0][0] == ghg.K * R


def test_make_grid():
    nrows, ncols, R = 2, 4, 1
    a0, b0 = -2, -2
    x0, y0 = ghg.double_to_cartesian(a0, b0, R)
    grid = ghg.make_grid(nrows, ncols, R=R, x0=x0, y0=y0, a0=a0, b0=b0)
    assert set(grid["cell_id"]) == {
        "-2,-2",
        "-1,-1",
        "0,-2",
        "1,-1",
        "-2,0",
        "-1,1",
        "0,0",
        "1,1",
    }
    # Cell areas should be correct
    assert np.allclose(grid.area, 3 * np.sqrt(3) * R**2 / 2)
    # Grid should have no gaps
    assert grid.union_all().boundary.is_ring


def test_make_grid_from_bounds():
    # A good example for checking edges
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
    # Grid should cover rectangle
    grid.union_all().contains(rect.union_all())

    # Two grids with the same origin should have identical cells where they overlap
    rect1 = gpd.GeoDataFrame(geometry=[sg.Polygon([(-2, 1), (3, 1), (3, 5), (-2, 5)])])
    rect2 = rect1.translate(-1, 1)
    R, ox, oy = 1, 0, 0
    grid1 = ghg.make_grid_from_bounds(*rect1.total_bounds, R=R, ox=ox, oy=oy)
    grid2 = ghg.make_grid_from_bounds(*rect2.total_bounds, R=R, ox=ox, oy=oy)
    cell_ids = set(grid1["cell_id"]) & set(grid2["cell_id"])
    g1 = grid1.loc[lambda x: x["cell_id"].isin(cell_ids)].sort_values(
        "cell_id", ignore_index=True
    )
    g2 = grid2.loc[lambda x: x["cell_id"].isin(cell_ids)].sort_values(
        "cell_id", ignore_index=True
    )
    assert g1.geom_equals_exact(g2, tolerance=10e-15).all()

    # Case where one hexagon will just cover
    R = 1
    r = ghg.K * R
    rect = gpd.GeoDataFrame(
        geometry=[sg.Polygon([(-R / 2, -r), (R / 2, -r), (R / 2, r), (-R / 2, r)])]
    )
    grid = ghg.make_grid_from_bounds(*rect.total_bounds, R=R)
    assert grid.union_all().contains(rect.union_all())

    # Same scenario but with origin hexagon shifted to southwest neighbor of covering
    # hexagon
    grid2 = ghg.make_grid_from_bounds(*rect.total_bounds, R=R, ox=-3 * R / 2, oy=-r)
    assert grid2.union_all() == grid.union_all()
    assert grid2.cell_id.tolist() == ["1,1"]

    # Case where all edges aren't initially covered
    R = 1
    r = ghg.K * R
    rect = gpd.GeoDataFrame(
        geometry=[
            sg.Polygon(
                [
                    (-0.6 * R, -0.1),
                    (2.2 * R, -0.1),
                    (2.2 * R, 3.2 * r),
                    (-0.6 * R, 3.2 * r),
                ]
            )
        ]
    )
    grid = ghg.make_grid_from_bounds(*rect.total_bounds, R=R)
    assert grid.union_all().contains(rect.union_all())

    # More edge cases
    rect = gpd.GeoDataFrame(geometry=[sg.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])])
    for grid in [
        # Right edge of rectangle not initially covered
        ghg.make_grid_from_bounds(*rect.total_bounds, R=0.27, ox=0, oy=0),
        # Top edge not
        ghg.make_grid_from_bounds(*rect.total_bounds, R=0.2, ox=0, oy=0),
        # Bottom edge not
        ghg.make_grid_from_bounds(*rect.total_bounds, R=0.5, ox=0, oy=0.1),
        # Left edge not
        ghg.make_grid_from_bounds(*rect.total_bounds, R=1, ox=0.6, oy=0),
    ]:
        assert grid.union_all().contains(rect.union_all())


def test_make_grid_from_gdf():
    shape = gpd.GeoDataFrame(geometry=[sg.Polygon([(1, -1), (3, 1), (0, 3)])])
    grid = ghg.make_grid_from_gdf(shape, R=1, ox=0, oy=0, intersect=True)
    assert set(grid.columns) == {"cell_id", "geometry"}
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

    # Grids with growing circumradius covering horseshoe
    g = gpd.GeoDataFrame(
        {"geometry": [sg.LineString([(0, 0), (1, 0), (1, 1), (0, 1)])]}
    )
    for i in range(0, 20):
        R = (i + 1) / 10
        r = ghg.K * R
        grid = ghg.make_grid_from_gdf(g, R=R, ox=-3 * R / 2 - 0.1, oy=-r - 0.1)
        assert grid.union_all().contains(g.union_all())

    shapes = gpd.read_file(DATA_DIR / "shapes.geojson").to_crs("epsg:2193")
    R = 900
    grid1 = ghg.make_grid_from_gdf(shapes, R=R, intersect=False)
    # CRS should be correct
    assert grid1.crs == shapes.crs
    # Grid should cover shapes
    grid1.dissolve().contains(shapes.dissolve())

    grid2 = ghg.make_grid_from_gdf(shapes, R=R, intersect=True)
    # Intersection should produce fewer cells
    assert grid2.shape[0] <= grid1.shape[0]
    # Intersection should reduce area
    assert shapes.area.sum() <= grid2.area.sum() <= grid1.area.sum()

    # Clipping should produce correct area
    grid3 = ghg.make_grid_from_gdf(shapes, R=R, clip=True)
    assert np.allclose(grid3.area.sum(), shapes.area.sum())
