import marimo

__generated_with = "0.13.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import sys
    import pathlib as pl
    import math
    import marimo as mo

    import numpy as np
    import geopandas as gpd
    import shapely.geometry as sg
    import matplotlib.pyplot as plt

    sys.path.append("../")

    import geohexgrid as ghg

    DATA_DIR = pl.Path("tests/data")
    NZTM = "epsg:2193"  # New Zealand Transverse Mercator CRS
    WGS84 = "epsg:4326"
    return DATA_DIR, NZTM, ghg, gpd, mo, np, plt, sg


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Using `make_grid_points` and `make_grid`

    Hexagons are named (via the `cell_id` column) by their *double coordinates* as shown here 
    and [described in more detail by Red Blob Games](https://www.redblobgames.com/grids/hexagons/#coordinates-doubled) .
    """
    )
    return


@app.cell
def _(ghg, np, plt):
    nrows, ncols, _R = (3, 4, 1)
    a0, b0 = (-2, -2)
    x0, y0 = ghg.double_to_cartesian(a0, b0, _R)
    X, Y = ghg.make_grid_points(nrows, ncols, _R, x0, y0)
    grid = ghg.make_grid(nrows, ncols, _R, x0, y0, a0=a0, b0=b0)
    print(grid.head())
    print(grid["cell_id"].tolist())
    fig, ax = plt.subplots()
    ax.scatter(X, Y, color="black", alpha=0.9)
    grid.plot(ax=ax, color="white", edgecolor="red", alpha=0.9, aspect="equal")
    grid.apply(
        lambda x: ax.annotate(
            text=x["cell_id"],
            xy=x.geometry.centroid.coords[0],
            ha="center",
            color="red",
        ),
        axis=1,
    )
    assert np.allclose(grid.area, 3 * np.sqrt(3) * _R**2 / 2)
    return (grid,)


@app.cell
def _(grid):
    # The grid should be gapless
    p = grid.union_all().boundary
    print(p, p.is_ring)

    # Here's a non-gapless collection of hexagons
    p1 = grid["geometry"].iat[0].buffer(-0.001)
    p2 = grid["geometry"].iat[1]
    q = p1.union(p2).boundary

    print(q, q.is_ring)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Using `make_grid_from_bounds`""")
    return


@app.cell
def _(ghg, gpd, sg):
    rect = gpd.GeoDataFrame(
        geometry=[sg.Polygon([(2.1, -1), (4.9, -1), (4.9, 1.9), (2.1, 1.9)])]
    )
    _R = 1
    grid_1 = ghg.make_grid_from_bounds(*rect.total_bounds, R=_R)
    print(grid_1.head())
    print(grid_1["cell_id"].tolist())
    _base = rect.plot(color="gray", aspect="equal")
    grid_1.plot(ax=_base, color="white", edgecolor="red", alpha=0.5)
    return


@app.cell
def _(ghg, gpd, sg):
    rect1 = gpd.GeoDataFrame(geometry=[sg.Polygon([(-2, 1), (3, 1), (3, 5), (-2, 5)])])
    _base = rect1.plot(alpha=0.5)
    rect2 = rect1.translate(-1, 1)
    rect2.plot(ax=_base, color="orange", alpha=0.5)
    _R = 1
    grid1 = ghg.make_grid_from_bounds(*rect1.total_bounds, R=_R)
    grid2 = ghg.make_grid_from_bounds(*rect2.total_bounds, R=_R)
    _base = grid1.plot(alpha=0.5)
    grid2.plot(ax=_base, alpha=0.5, color="orange")
    cell_ids = set(grid1["cell_id"]) & set(grid2["cell_id"])
    g1 = grid1.loc[lambda x: x["cell_id"].isin(cell_ids)].sort_values(
        "cell_id", ignore_index=True
    )
    g2 = grid2.loc[lambda x: x["cell_id"].isin(cell_ids)].sort_values(
        "cell_id", ignore_index=True
    )
    g1.geom_equals_exact(g2, tolerance=1e-11).all()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Using `make_grid_from_gdf`""")
    return


@app.cell
def _(ghg, gpd, sg):
    shape = gpd.GeoDataFrame(geometry=[sg.Polygon([(1, -1), (3, 1), (0, 3)])])
    grid_2 = ghg.make_grid_from_gdf(shape, R=1)
    print(grid_2)
    _base = shape.plot(color="gray", aspect="equal")
    grid_2.plot(ax=_base, color="white", edgecolor="red", alpha=0.5)
    grid2_1 = ghg.make_grid_from_gdf(shape, R=1, clip=True)
    print(grid2_1)
    _base2 = shape.plot(color="gray", aspect="equal")
    grid2_1.plot(ax=_base2, color="white", edgecolor="red", alpha=0.5)
    return


@app.cell
def _(DATA_DIR, NZTM, ghg, gpd):
    shapes = gpd.read_file(DATA_DIR / "shapes.geojson").to_crs(NZTM)
    _R = 900
    grid_3 = ghg.make_grid_from_gdf(shapes, R=_R)
    print(grid_3.head())
    _base = shapes.plot(color="gray", aspect="equal")
    grid_3.plot(ax=_base, color="white", edgecolor="red", alpha=0.5)
    _base2 = shapes.plot(color="gray", aspect="equal")
    grid2_2 = ghg.make_grid_from_gdf(shapes, R=_R, clip=True)
    assert grid_3.shape[0] == grid2_2.shape[0]
    grid2_2.plot(ax=_base2, color="white", edgecolor="red", alpha=0.5)

    return


@app.cell
def _(DATA_DIR, ghg, gpd):
    nz = gpd.read_file(DATA_DIR / "nz_tas.gpkg")
    print(nz.crs)
    _base = nz.plot(color="gray", aspect="equal", figsize=(10, 10))
    grid_4 = ghg.make_grid_from_gdf(nz, R=10_000)  # 10 km inradius grid
    grid_4.plot(ax=_base, color="white", edgecolor="red", alpha=0.5)

    return


if __name__ == "__main__":
    app.run()
