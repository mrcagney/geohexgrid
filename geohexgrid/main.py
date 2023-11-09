from __future__ import annotations
import math

import numpy as np
import geopandas as gpd
import shapely.geometry as sg


K = np.sqrt(3) / 2  # cosine of 30°


def double_to_cartesian(a: float, b: float, R: float) -> tuple[float]:
    """
    Given double coordinates of a hexagon in a flat-top hexagonal
    grid centered at the origin with hexagon circumradius ``R``,
    return the Cartesian coordinates of its center.

    For details on double coordinates see [RBG]_.
    """
    return 3 * a * R / 2, b * K * R


def cartesian_to_double(x: float, y: float, R: float) -> tuple[float]:
    """
    Given Cartesian coordinates of the center of a hexgon in a flat-top hexagon
    grid centered at the origin with hexagon circumradius ``R``,
    return the double coordinates of the hexagon

    For details on double coordinates see [RBG]_.
    """
    return int(2 * x / (3 * R)), int(y / (K * R))


def make_grid_points(nrows, ncols, R: float = 1, ox: float = 0, oy: float = 0):
    """
    Make the vertices and centers of a flat-top hexagon grid with
    ``nrows`` rows, ``ncols`` columns, circumradius ``R``,
    and bottom left hexagon centered at ``(ox, oy)``.
    A row is a left-to-right down-up zig-zag of hexagons and the next row is stacked
    on top of the previous row.
    """
    r = K * R
    x, y = np.meshgrid(
        np.linspace(0, (2 * ncols - 1) * R, 2 * ncols),
        np.linspace(0, (3 * nrows - 1) * r, 3 * nrows),
        sparse=False,
        indexing="xy",
    )
    x[1::2, :] -= R / 2
    return x - R / 2 + ox, y - r + oy


def make_grid(
    nrows,
    ncols,
    R: float = 1,
    ox: float = 0,
    oy: float = 0,
    oa: int = 0,
    ob: int = 0,
) -> gpd.GeoDataFrame:
    """
    Make a grid of flat-top hexagons with ``nrows`` rows, ``ncols`` columns,
    circumradius ``R``, and bottom left hexagon centered at ``(ox, oy)``.
    A row is a left-to-right down-up zig-zag of hexagons and the next row is stacked
    on top of the previous row.

    Use double coordinates (see [RBG]_) for cell IDs starting from the
    bottom-left hexagon/cell with cell ID ``f'oa,ob'``, which defaults to ``'0,0'``.

    NOTES:

    - The cell IDs follow double coordinates, e.g. the cell IDs of the first two rows
      default to
      '0,0', '1,1', '2,0', '3,1', '4,0', '5,1',...
      '0,2', '1,3', '2,2', '3,3', '4,2', '5,3',...
    - Making a 1000 x 1000 grid with this on my computer takes about 12 seconds

    """
    X, Y = make_grid_points(nrows=nrows, ncols=ncols, ox=ox, oy=oy, R=R)
    y = Y[:, 0]
    # Use double coordinates for cell IDs
    cell_id = [f"{oa + j},{ob + i + j%2}" for i in range(nrows) for j in range(ncols)]
    """
    Make hexagons, each of which has vertex order:

       v4  v3

    v5        v6

       v0  v1

    """
    geometry = [
        sg.Polygon(
            [
                [X[j % 2][(3 * j + 1) // 2], y[2 * i + j % 2]],  # v0
                [X[j % 2][math.ceil((3 * j + 2) / 2)], y[2 * i + j % 2]],  # v1
                [X[j % 2 + 1][3 * j // 2 + 2], y[2 * i + j % 2 + 1]],  # v2
                [X[j % 2][math.ceil((3 * j + 2) / 2)], y[2 * i + j % 2 + 2]],  # v3
                [X[j % 2][(3 * j + 1) // 2], y[2 * i + j % 2 + 2]],  # v4
                [X[j % 2 + 1][3 * j // 2], y[2 * i + j % 2 + 1]],  # v5
            ]
        )
        for i in range(nrows)
        for j in range(ncols)
    ]
    return gpd.GeoDataFrame({"cell_id": cell_id, "geometry": geometry})


def make_grid_from_bbox(
    minx: float,
    miny: float,
    maxx: float,
    maxy: float,
    R: float,
    ox: float | None = None,
    oy: float | None = None,
    crs: str | None = None,
) -> gpd.GeoDataFrame:
    """
    Return a grid of flat-top hexagons, each of circumradius ``R``
    (and hence area :math:`\frac{3 \sqrt(3)}{2} R^2`), that minimally covers the
    rectangle  with the given coordinate extrema ``minx``, ``miny``, ``maxx``, ``maxy``.
    Label each hexagon with an ID listed in the column ``'cell_id'``,
    which is the double coordinates of the hexagon relative to a (``'0,0'``) origin
    hexagon centered at point ``(ox, oy)``, which defaults to ``(minx, miny)``.
    For more details on double coordinates, see [RBG]_.

    The grid will lie in the plane of the given CRS (which defaults to ``None``)
    and will use its distance units,
    e.g. metres for the New Zealand Transverse Mercator (NZTM) CRS
    and decimal degrees for the WGS84 CRS.
    """
    # A column of i cells has covering height >= (2 * i - 1) * r, where r = K * R, so
    nrows = math.ceil((maxy - miny) / (2 * K * R) + 1 / 2)
    # A row of j cells has covering width >= (3 * j - 2) * R / 2, so
    ncols = math.ceil(2 * (maxx - minx) / (3 * R) + 2 / 3)

    if ox is None or oy is None:
        # Cover the box with a grid whose origin lies at minx, miny
        grid = make_grid(nrows=nrows, ncols=ncols, ox=minx, oy=miny, R=R)

    else:
        # Cover the box with a grid whose origin lies at p = (ox, oy)
        p = np.array((ox, oy))  # origin of grid
        q = np.array((minx, miny))  # first cell of box covers this point

        # Relative to a hexgrid with origin p, get the cell ID of q
        a, b = cartesian_to_double(*(q - p), R)

        # Cover the box translated to p, then translate it back to q and use cell IDs
        # relative to a p-center
        grid = make_grid(
            nrows=nrows, ncols=ncols, ox=ox, oy=oy, R=R, oa=a, ob=b
        ).assign(geometry=lambda x: x.translate(*(q - p)))
    grid.crs = crs
    return grid


def make_grid_from_gdf(
    g: gpd.GeoDataFrame,
    R: float,
    ox: float | None = None,
    oy: float | None = None,
    *,
    intersect: bool = True,
    clip: bool = False,
) -> gpd.GeoDataFrame:
    """
    Return a grid of flat-top hexagons, each of circumradius ``R``
    (and hence area :math:`\frac{3 \sqrt(3)}{2} R^2`), that minimally covers the
    the bounding box (total bounds) of the given GeoDataFrame ``g``.
    Label each hexagon with an ID listed in the column ``'cell_id'``,
    which is the double coordinates of the hexagon relative to a (``'0,0'``) origin
    hexagon centered at point ``(ox, oy)``, which defaults to the bottom left
    point of the bounding box of ``g``.
    For more details on double coordinates, see [RBG]_.

    The grid will lie in the plane of ``g``'s CRS and will use its distance units,
    e.g. metres for the New Zealand Transverse Mercator (NZTM) CRS
    and decimal degrees for the WGS84 CRS.

    If ``intersect``, then return only the hexagons that intersect ``g``'s features.
    This option computes a spatial join, which can slow things down if the number of
    hexagons is large and the feature set detailed.
    In that case, install pyogrio for a speed up or simplify the features beforehand.

    If ``clip``, then return the grid clipped to ``g``, which may contain fragments
    of hexagons.

    EXAMPLE::

        # Load some New Zealand features and set the CRS to NZTM
        g = gpd.read_file(my_path).to_crs("epsg:2193")

        # Make a hex grid of 250 metre circumradius and keep only the portion
        # that intersects the features of ``g``
        grid = make_grid(g, R=250)

    REFERENCES::

    .. [RBG] Hexagonal Grids, https://www.redblobgames.com/grids/hexagons

    """
    grid = make_grid_from_bbox(*g.total_bounds, R=R, ox=ox, oy=oy, crs=g.crs)
    if intersect:
        grid = (
            grid.sjoin(g)
            .drop_duplicates(subset=["cell_id"])
            .filter(["cell_id", "geometry"])
        )
    if clip:
        grid = grid.clip(g)

    return grid
