import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell
def _():
    import sys
    import pathlib as pl
    import math

    import numpy as np
    import geopandas as gpd
    import shapely.geometry as sg
    import matplotlib.pyplot as plt
    from loguru import logger

    sys.path.append("../")

    import geohexgrid as ghg

    # magic command not supported in marimo; please file an issue to add support
    # %load_ext autoreload
    # '%autoreload 2' command supported automatically in marimo

    DATA_DIR = pl.Path("../tests/data")
    NZTM = "epsg:2193"  # New Zealand Transverse Mercator CRS
    WGS84 = "epsg:4326"
    return DATA_DIR, NZTM, ghg, gpd, np, plt, sg


@app.cell
def _(ghg):
    _nrows, _ncols, _R = (1000, 1000, 1)
    _a0, _b0 = (-2, -2)
    _x0, _y0 = ghg.double_to_cartesian(_a0, _b0, _R)
    return


@app.cell
def _(display, ghg, np, plt):
    _nrows, _ncols, _R = (3, 4, 1)
    _a0, _b0 = (-2, -2)
    _x0, _y0 = ghg.double_to_cartesian(_a0, _b0, _R)
    X, Y = ghg.make_grid_points(_nrows, _ncols, _R, _x0, _y0)
    grid = ghg.make_grid(_nrows, _ncols, _R, _x0, _y0, a0=_a0, b0=_b0)
    display(grid.head())
    print(grid["cell_id"].tolist())
    _fig, _ax = plt.subplots()
    _ax.scatter(X, Y, color="black", alpha=0.9)
    grid.plot(ax=_ax, color="white", edgecolor="red", alpha=0.9, aspect="equal")
    grid.apply(
        lambda x: _ax.annotate(
            text=x["cell_id"],
            xy=x.geometry.centroid.coords[0],
            ha="center",
            color="red",
        ),
        axis=1,
    )
    assert np.allclose(grid.area, 3 * np.sqrt(3) * _R**2 / 2)
    return


@app.cell
def _(display, ghg):
    grid_1 = ghg.make_grid(2, 3, R=1, x0=2, y0=1)
    display(grid_1)
    grid_1["cell_id"].tolist()
    p = grid_1.union_all().boundary
    display(p, p.is_ring)
    p1 = grid_1["geometry"].iat[0].buffer(-0.001)
    p2 = grid_1["geometry"].iat[1]
    q = p1.union(p2).boundary
    display(q, q.is_ring)
    return


@app.cell
def _(display, ghg, gpd, sg):
    _R = 1
    _r = ghg.K * _R
    _rect = gpd.GeoDataFrame(
        geometry=[
            sg.Polygon(
                [
                    (-0.6 * _R, -0.1),
                    (2.2 * _R, -0.1),
                    (2.2 * _R, 3.2 * _r),
                    (-0.6 * _R, 3.2 * _r),
                ]
            )
        ]
    )
    grid_2 = ghg.make_grid_from_bounds(*_rect.total_bounds, R=_R)
    display(grid_2.head())
    print(grid_2["cell_id"].tolist())
    _base = _rect.plot(color="gray", aspect="equal")
    grid_2.plot(ax=_base, color="white", edgecolor="red", alpha=0.5)
    return


@app.cell
def _(display, ghg, gpd, sg):
    _R = 1
    _r = ghg.K * _R
    _rect = gpd.GeoDataFrame(
        geometry=[
            sg.Polygon([(-_R / 2, -_r), (_R / 2, -_r), (_R / 2, _r), (-_R / 2, _r)])
        ]
    )
    grid_3 = ghg.make_grid_from_bounds(
        *_rect.total_bounds, R=_R, ox=-3 * _R / 2, oy=-_r
    )
    display(grid_3.head())
    print(grid_3["cell_id"].tolist())
    _base = _rect.plot(color="gray", aspect="equal")
    grid_3.plot(ax=_base, color="white", edgecolor="red", alpha=0.5)
    return


@app.cell
def _(ghg, gpd, plt, sg):
    _rect = gpd.GeoDataFrame(geometry=[sg.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])])
    for grid_4 in [
        ghg.make_grid_from_bounds(*_rect.total_bounds, R=0.27, ox=0, oy=0),
        ghg.make_grid_from_bounds(*_rect.total_bounds, R=0.2, ox=0, oy=0),
        ghg.make_grid_from_bounds(*_rect.total_bounds, R=0.5, ox=0, oy=0.1),
        ghg.make_grid_from_bounds(*_rect.total_bounds, R=1, ox=0.6, oy=0),
    ]:
        _fig, _ax = plt.subplots()
        _rect.plot(ax=_ax, color="red")
        grid_4.plot(ax=_ax, color="blue", alpha=0.5)
        plt.show()
    return


@app.cell
def _(ghg, gpd, plt, sg):
    g = gpd.GeoDataFrame(
        {"geometry": [sg.LineString([(0, 0), (1, 0), (1, 1), (0, 1)])]}
    )
    for i in range(0, 20):
        _R = (i + 1) / 10
        _r = ghg.K * _R
        grid_5 = ghg.make_grid_from_gdf(g, R=_R, ox=-3 * _R / 2 - 0.1, oy=-_r - 0.1)
        assert grid_5.union_all().contains(g.union_all())
        _fig, _ax = plt.subplots()
        g.plot(ax=_ax, color="red")
        grid_5.plot(ax=_ax, color="blue", alpha=0.5)
        plt.show()
    return


@app.cell
def _(display, ghg, gpd, sg):
    shape = gpd.GeoDataFrame(geometry=[sg.Polygon([(1, -1), (3, 1), (0, 3)])])
    grid_6 = ghg.make_grid_from_gdf(shape, R=1)
    display(grid_6.head())
    print(grid_6["cell_id"].tolist())
    _base = shape.plot(color="gray", aspect="equal")
    grid_6.plot(ax=_base, color="white", edgecolor="red", alpha=0.5)
    return


@app.cell
def _(DATA_DIR, NZTM, display, ghg, gpd):
    _shapes = gpd.read_file(DATA_DIR / "shapes.geojson").to_crs(NZTM)
    _R = 900
    grid_7 = ghg.make_grid_from_gdf(_shapes, R=_R)
    display(grid_7.head())
    assert grid_7.union_all().contains(_shapes.union_all())
    _base = _shapes.plot(color="gray", aspect="equal")
    grid_7.plot(ax=_base, color="white", edgecolor="red", alpha=0.5)
    return (grid_7,)


@app.cell
def _(DATA_DIR, display, gpd, grid_7, nz, plt):
    _shapes = gpd.read_file(DATA_DIR / "nz_tas.gpkg")
    display(nz.crs)
    assert grid_7.union_all().contains(_shapes.union_all())
    _base = _shapes.plot(color="black", aspect="equal", figsize=(10, 10))
    grid_7.plot(ax=_base, color="white", edgecolor="red", alpha=0.5)
    plt.axis("off")
    plt.savefig("../nz_10000m.png", bbox_inches="tight")
    return


@app.cell
def _(DATA_DIR, gpd, grid_7):
    nz = gpd.read_file(DATA_DIR / "nz_tas.gpkg")
    akl = nz.loc[lambda x: x["ta2021_name"] == "Auckland"]
    _base = akl.plot(color="black", figsize=(20, 20), aspect="equal")
    grid_7.plot(ax=_base, color="white", edgecolor="red", alpha=0.5)
    assert grid_7.union_all().contains(akl.union_all())
    return (nz,)


if __name__ == "__main__":
    app.run()
