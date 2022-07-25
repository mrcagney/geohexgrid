"""
REFERENCES::

.. [RBG] Hexagonal Grids, https://www.redblobgames.com/grids/hexagons

"""

import collections
from dataclasses import dataclass
from math import sqrt, sin, cos, pi

import geopandas as gpd
import shapely.geometry as sg


START_ANGLE = 30  # degrees

Point = collections.namedtuple("Point", ["x", "y"])


@dataclass
class HexGridSystem:
    """
    Represents a planar flat-top hexagon grid system in the geographic
    coordinate reference sytem specified by ``crs``.
    """

    crs: str  # CRS string acceptable by GeoPandas, e.g. 'epsg:2193' for NZTM
    x: float  # first Cartesian (CRS) coordinate of system's origin
    y: float  # second Cartesion (CRS) coordinate of system's origin
    R: float  # circumradius of the hexagon grid cells in the CRS's distance units

    def cell_from_axial_point(self, a, b) -> "Cell":
        """
        Given axial coordinates of a point in the plane, return the cell containing it.
        """
        a_round, b_round = round(a), round(b)
        a, b = a - a_round, b - b_round  # remainders
        if abs(a) >= abs(b):
            center = int(a_round + round(a + 0.5 * b)), int(b_round)
        else:
            center = int(a_round), int(b_round + round(b + 0.5 * a))

        return Cell(self, center[0], center[1])

    def cell_from_point(self, x, y) -> "Cell":
        """
        Given Cartesian coordinates of a point in the plane, return the cell containing
        it.
        """
        return self.cell_from_axial_point(*cartesian_to_axial(x, y, self.R))

    def cells_from_bbox(self, minx, miny, maxx, maxy, *, as_gdf=False) -> list["Cell"]:
        ll = self.cell_from_point(minx, miny)
        lr = self.cell_from_point(maxy, miny)
        ul = self.cell_from_point(minx, maxy)
        ur = self.cell_from_point(maxx, maxy)
        # Bottom row
        mina = min(ll.a, ul.a)
        minb = min(ll.b, lr.b)
        maxa = max(lr.a, ur.a)
        maxb = max(ul.b, ur.b)

        result = [
            self.cell_from_axial_point(i, j)
            for i in range(mina, maxa + 1)
            for j in range(minb, maxb + 1)
        ]
        if as_gdf:
            result = gpd.GeoDataFrame(
                data={"cell_id": [cell.id for cell in result]},
                geometry=[cell.polygon() for cell in result],
                crs=self.crs,
            )

        return result

    def cells_from_gdf(self, gdf: gpd.GeoDataFrame) -> list["Cell"]:
        pass


@dataclass
class Cell:
    """
    Represents a hexagon cell within a HexGridSystem.
    """

    hgs: HexGridSystem
    a: int  # first axial coordinate of center of cell
    b: int  # second axial coordinate of center of cell

    def __post_init__(self):
        # Derive some attributes
        self.R = self.hgs.R  # circumradius
        self.id = f"{self.a},{self.b}"

    def center(self) -> Point:
        """
        Return the Cartesian center of this cell (relative to the CRS).
        """
        return Point(*axial_to_cartesian(self.a, self.b, self.R))

    def vertices(self) -> list[tuple[float]]:
        c = self.center()
        return [hexagon_vertex(c.x, c.y, self.R, i) for i in range(6)]

    def polygon(self) -> sg.Polygon:
        """
        Return the Shapely polygon of this cell (relative to the CRS).
        """
        return sg.Polygon(self.vertices() + self.vertices()[:1])


def axial_to_cartesian(a: float, b: float, R: float) -> tuple[float]:
    """
    Given axial coordinates of a point in the plane relative to a flat-top hexagonal
    grid centered at the origin with hexagon circumradius ``R``,
    return its Cartesian coordinates.

    For more details on axial coordinates see [RBG]_.
    """
    return (R * (3 / 2) * a, R * ((sqrt(3) / 2) * a + sqrt(3) * b))


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

        3           0

            4   5

    For ``i`` greater than 5 or negative, then wrap in standard mathematical fashion.
    """
    θ = i * pi / 3
    return x + R * cos(θ), y + R * sin(θ)


def cartesian_point_to_cell_name(x: float, y: float, R: float) -> tuple[int]:
    """
    Given Cartesian coordinates of a point in the plane, return the (integer) axial
    coordinates of the center of the flat-top hexagonal cell containing it relative
    to a hexagon grid of circumradius R.
    """
    return
