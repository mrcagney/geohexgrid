from typing import Optional
from dataclasses import dataclass
from math import sqrt, sin, cos, pi

import pyproj
import geopandas as gpd
import shapely.geometry as sg


@dataclass(frozen=True)
class GeoHexSystem:
    """
    Represents a flat-top hexagon grid system of circumradius ``R`` in the Cartesian
    plane defined by the two axes of the geographic coordinate reference system (CRS)
    ``crs`` and with central hexagon cell centered at the point ``(x, y)``.
    Cell IDs are derived from the axial coordinates of the cell's center.

    For example, ``GeoHexSystem("epsg:2193", 250, 500, 1755147,  5921401)`` creates
    a hexagon grid system in New Zealand Transverse Mercator coordinates with
    hexagons of circumradius 250 metres and central hexagon centered in central Auckland.

    For example, ``GeoHexSystem("epsg:4326", 0.001, 174.74, -36.840556)`` creates
    a hexagon grid system in WGS84 coordinates with hexagons of circumradius 0.001
    decimal degrees and central hexagon centered in central Auckland.

    If you want your hexagon grid system to faithfully represent hexagons on the ground
    in your area, then choose a mostly distance- and area-preserving projection CRS,
    such as NZTM when analysing New Zealand.

    To actually make a grid with this grid system, see e.g. the method
    :meth:`grid_from_gdf`.

    For more about hexagon grids independent of geography, see the excellent website
    [RBG]_.

    REFERENCES::

    .. [RBG] Hexagonal Grids, https://www.redblobgames.com/grids/hexagons

    """

    crs: str
    R: float
    x: float
    y: float

    @staticmethod
    def validate(ghs: "GeoHexSystem") -> "GeoHexSystem":
        """
        Return the given GeoHexSystem if it is valid.
        Otherwise, raise a ValueError.
        """
        try:
            pyproj.CRS.from_user_input(ghs.crs)
        except Exception:
            raise ValueError(f"Invalid CRS {ghs.crs}")

        if ghs.R <= 0:
            raise ValueError(f"Circumradius must be positive but received R={ghs.R}")

        return ghs

    def __post_init__(self):
        GeoHexSystem.validate(self)

    def cell_from_id(self, cell_id: str) -> "Cell":
        """
        Given the ID of a Cell in this GeoHexSystem, as defined by the method
        :meth:`Cell.id`, return the corresponding Cell.
        """
        a, b = axial_center_from_id(cell_id)
        return Cell(int(a), int(b), self.R)

    def cell_from_axial_point(self, a: float, b: float) -> "Cell":
        """
        Given axial coordinates of a point in the plane, return the Cell containing it.
        """
        a_round, b_round = round(a), round(b)
        a, b = a - a_round, b - b_round  # remainders
        if abs(a) >= abs(b):
            center = int(a_round + round(a + 0.5 * b)), int(b_round)
        else:
            center = int(a_round), int(b_round + round(b + 0.5 * a))

        return Cell(center[0], center[1], self.R)

    def cell_from_point(self, x, y) -> "Cell":
        """
        Given Cartesian coordinates of a point in the plane,
        return the Cell containing it.
        """
        return self.cell_from_axial_point(*cartesian_to_axial(x, y, self.R))

    def grid_from_bbox(
        self,
        minx: float,
        miny: float,
        maxx: float,
        maxy: float,
        *,
        as_gdf: bool = True,
    ) -> list["Cell"] | gpd.GeoDataFrame:
        """
        Return a minimal set of Cells that covers the given Cartesian bounding box
        (relative to the CRS of this GeoHexSystem).
        If ``as_gdf``, then return the resulting grid as a GeoDataFrame in the CRS of
        this HexGridStystem, with a ``cell_id`` column.
        """
        # Get the two Cells containing the left-down and right-up corner points of the
        # bounding box
        ld = self.cell_from_point(minx, miny)
        ru = self.cell_from_point(maxx, maxy)

        ld_x, ld_y = ld.center()
        ru_x, ru_y = ru.center()

        # Think in terms of the double coordinate system for hexagons (see [RBG]_)
        # and make the cell cover from left to right (rows) and bottom to top (columns).
        ld_p, ld_q = ld.center_double()
        ru_p, ru_q = ru.center_double()

        # Get horizontal extrema as double coordinate values
        if (ld_x - minx) > self.R / 2 and ld_q != ru_q:
            minp = ld_p - 1
        else:
            minp = ld_p
        if (maxx - ru_x) > self.R / 2 and ld_q != ru_q:
            maxp = ru_p + 1
        else:
            maxp = ru_p

        # Get vertical extrema as double coordinate values
        if miny < ld_y and ld_p != ru_p:
            minq = ld_q - 1
        else:
            minq = ld_q
        if maxy > ru_y and ld_p != ru_p:
            maxq = ru_q + 1
        else:
            maxq = ru_q

        # Make cover
        def make_row(c: "Cell") -> list["Cell"]:
            """
            Make the (double coordinates) row of cells occupied by the given cell.
            Could be an empty list.
            """
            row = []
            while c.center_double()[0] <= maxp:
                row.append(c)
                # Get next cell of row
                c = c.neighbour("NE").neighbour("SE")
            return row

        cover = []
        # Set how to get from the start of one row to the start of the row above:
        # zig-zag or zag-zig
        if minp < ld_p:
            d = {0: "NW", 1: "NE"}
        else:
            d = {0: "NE", 1: "NW"}

        if minq < ld_q:
            # The bottom row lies below cell ld, so build that row first
            row = make_row(ld.neighbour("SE"))
            cover.extend(row)

        # Build rest of rows, working our way upwards
        i = 0
        c = ld
        while c.center_double()[1] <= maxq:
            row = make_row(c)  # Could be empty
            cover.extend(row)
            # Get start cell of next row up
            c = c.neighbour(d[i % 2])
            i += 1

        if as_gdf:
            cover = gpd.GeoDataFrame(
                data={"cell_id": [cell.id() for cell in cover]},
                geometry=[cell.polygon() for cell in cover],
                crs=self.crs,
            )

        return cover

    def grid_from_gdf(
        self, g: gpd.GeoDataFrame, *, as_gdf=True, intersect=False
    ) -> list["Cell"] | gpd.GeoDataFrame:
        """
        Return a minimal set of Cells that covers the bounding box (total bounds) of the
        given GeoDataFrame whose CRS matches that of this GeoHexSystem.
        If ``as_gdf``, then return the resulting grid as a GeoDataFrame in the CRS of
        this GeoHexSystem, with a ``cell_id`` column.
        If ``intersect``, then only keep the cells that cover the features of
        the GeoDataFrame, which can be much fewer than cover its bounding box.
        """
        g = g.to_crs(self.crs)
        if intersect:
            grid = (
                self.grid_from_bbox(*g.total_bounds, as_gdf=True)
                .sjoin(g)
                .drop_duplicates(subset=["cell_id"])
                .filter(["cell_id", "geometry"])
            )
            if not as_gdf:
                grid = [
                    Cell(int(a), int(b), self.R) for a, b in grid.cell_id.str.split(",")
                ]
        else:
            grid = self.grid_from_bbox(*g.total_bounds, as_gdf=as_gdf)

        return grid


@dataclass(frozen=True)
class Cell:
    """
    Represents a flat-top hexagon cell in the Cartesian plane centered at
    axial coordinates ``(a, b)`` and with circumradius ``R``.

    For an explanation of axial coordinates, double coordinates, etc. of hexagon grids,
    see [RBG]_.
    """

    a: int
    b: int
    R: float

    def id(self) -> str:
        """
        Return the (integer) axial coordinates of the center of this Cell in string
        form as a unique identifier for this Cell within the GeoHexSystem
        (but not unique across GeoHexSystems).
        """
        return axial_center_to_id(self.a, self.b)

    def center_axial(self) -> tuple[int]:
        """
        Return the center of this Cell in axial coordinates.
        """
        return self.a, self.b

    def center_double(self) -> tuple[int]:
        """
        Return the center of this Cell in double coordinates.
        """
        return axial_to_double(self.a, self.b)

    def center(self) -> tuple[float]:
        """
        Return the center of this Cell in Cartesian coordinates, that is,
        the coordinates of the CRS of the GeoHexSystem of this Cell.
        """
        return axial_to_cartesian(self.a, self.b, self.R)

    def vertices(self) -> list[tuple[float]]:
        """
        Return the vertices of this Cell in the anticlockwise order pictured here::

               2   1

            3         0

               4   5

        """
        x, y = self.center()
        return [hexagon_vertex(x, y, self.R, i) for i in range(6)]

    def polygon(self) -> sg.Polygon:
        """
        Return this Cell's Polygon with respect to the CRS of the GeoHexSystem of the
        Cell.
        """
        return sg.Polygon(self.vertices() + self.vertices()[:1])

    def neighbour(self, direction: str) -> "Cell":
        """
        Return this Cell's neighbour Cell in the given direction, which must be one of
        the compass direction strings 'N', 'NW', 'SW', 'S', 'SE', 'NE' whose meanings
        are diagramed here with this Cell in the center::

                N
            NW     NE
                .
            SW     SE
                S

        Raise a ValueError when given an invalid direction.
        """
        a, b = self.a, self.b
        d = {
            "NE": (a + 1, b),
            "N": (a, b + 1),
            "NW": (a - 1, b + 1),
            "SW": (a - 1, b),
            "S": (a, b - 1),
            "SE": (a + 1, b - 1),
        }
        try:
            return Cell(*d[direction], self.R)
        except KeyError:
            raise ValueError(f"Invalid direction {direction}")

    def neighbor(self, direction: str) -> "Cell":
        """
        Alias of :method:`neighbour` for the American English spellers.
        """
        return self.neighbour(direction)


# ----------------
# Main functions
# ----------------
def make_grid(
    g: gpd.GeoDataFrame,
    R: float,
    x: Optional[float] = None,
    y: Optional[float] = None,
    *,
    intersect: bool = False,
) -> gpd.GeoDataFrame:
    """
    Return a GeoDataFrame containing a minimal set of flat-top hexagons with
    circumradius ``R`` that covers the bounding box (total bounds) of GeoDataFrame
    ``g``.
    Optionally, center the grid at point ``(x, y)``, which defaults to the bottom left
    corner of ``g``'s bounding box.

    The hexagons will lie in the plane of ``g``'s CRS and will use its distance units,
    e.g. metres for the New Zealand Transverse Mercator (NZTM) CRS
    and decimal degrees for the WGS84 CRS.

    If ``intersect``, then return only the hexagons that intersect ``g``'s features.
    Be warned that this option computes a spatial join, which can slow things down
    if the number of hexagons is large and the feature set detailed.
    In that case, install PyGeos for a speed up or simplify the features beforehand.

    EXAMPLE::

        # Load some New Zealand features and set the CRS to NZTM
        g = gpd.read_file(my_path).to_crs("epsg:2193")

        # Make a hex grid of 250 metre circumradius and keep only the portion
        # that intersects the features of ``g``
        grid = make_grid(g, 250, intersect=True)

    """
    if x is None or y is None:
        x, __, y, __ = g.total_bounds

    return GeoHexSystem(g.crs, R, x, y).grid_from_gdf(g, intersect=intersect)


# -----------------
# Helper functions
# -----------------
def axial_center_to_id(a: int, b: int) -> str:
    return f"{a}-{b}"


def id_to_axial_center(cid: str) -> tuple[int]:
    a, b = cid.split("-")
    return int(a), int(b)


def axial_to_cartesian(a: float, b: float, R: float) -> tuple[float]:
    """
    Given axial coordinates of a point in the plane relative to a flat-top hexagonal
    grid centered at the origin with hexagon circumradius ``R``,
    return its Cartesian coordinates.

    For more details on axial coordinates see [RBG]_.
    """
    return (R * (3 / 2) * a, R * ((sqrt(3) / 2) * a + sqrt(3) * b))


def axial_to_double(a: float, b: float) -> tuple[float]:
    """
    Given axial coordinates of a point in the plane relative to a flat-top hexagonal
    grid centered at the origin, return its double coordinates.
    """
    return a, a + 2 * b


def double_to_axial(r: float, c: float) -> tuple[float]:
    """
    Given double coordinates of a point in the plane relative to a flat-top hexagonal
    grid centered at the origin, return its axial coordinates.
    """
    return r, (c - r) / 2


def cartesian_to_axial(x: float, y: float, R: float) -> tuple[float]:
    """
    Given Cartesian coordinates of a point in the plane, return its axial coordinates
    relative to a flat-top hexagon grid centered at the origin with hexagon
    circumradius ``R``.

    For more details on axial coordinates see [RBG]_.
    """
    return (2 / 3) * x / R, (-x / 3 + (sqrt(3) / 3) * y) / R


def hexagon_vertex(x: float, y: float, R: float, i: int) -> tuple[float]:
    """
    Given a flat-top hexagon centered at Cartesian point ``(x, y)`` with
    circumradius ``R``, return its index ``i`` vertex in counterclockwise order
    according to this diagram::

           2   1

        3         0

           4   5

    Take ``i`` modulo 6 to allow for any integer input.
    """
    θ = i * pi / 3
    return x + R * cos(θ), y + R * sin(θ)
