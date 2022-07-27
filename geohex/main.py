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


@dataclass(frozen=True)
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

    # TODO: optimize and comment
    def grid_from_bbox(self, minx, miny, maxx, maxy, *, as_gdf=False) -> list["Cell"]:
        """
        Return a minimal cell cover for the bounding box with the given extrema.
        """
        ld = self.cell_from_point(minx, miny)
        ru = self.cell_from_point(maxx, maxy)

        ld_x, ld_y = ld.center()
        ru_x, ru_y = ru.center()

        # Build cover with rows and columns of the double coordinate system.
        ld_p, ld_q = ld.center_double()
        ru_p, ru_q = ru.center_double()

        # Get min and max horizontal coordinates
        if (ld_x - minx) > self.R / 2 and ld_q != ru_q:
            minp = ld_p - 1
        else:
            minp = ld_p
        if (maxx - ru_x) > self.R / 2 and ld_q != ru_q:
            maxp = ru_p + 1
        else:
            maxp = ru_p

        # Get min and max vertical coordinates
        if miny < ld_y and ld_p != ru_p:
            minq = ld_q - 1
        else:
            minq = ld_q
        if maxy > ru_y and ld_p != ru_p:
            maxq = ru_q + 1
        else:
            maxq = ru_q

        # Make cover
        def make_row(start_cell):
            """
            Row could be empty.
            """
            row = []
            c = start_cell
            while c.center_double()[0] <= maxp:
                row.append(c)
                c = c.neighbour("ru").neighbour("rd")
            return row

        cover = []
        if minq < ld_q:
            # Fill row below ld
            row = make_row(ld.neighbour("rd"))
            cover.extend(row)

        # Set how to zig-zag to up neighbour
        if minp < ld_p:
            d = {0: "lu", 1: "ru"}
        else:
            d = {0: "ru", 1: "lu"}

        i = 0
        c = ld
        while c.center_double()[1] <= maxq:
            row = make_row(c)  # Could be empty
            cover.extend(row)
            # Get up neighbour
            c = c.neighbour(d[i % 2])
            i += 1

        if as_gdf:
            cover = gpd.GeoDataFrame(
                data={"cell_id": [cell.id() for cell in cover]},
                geometry=[cell.polygon() for cell in cover],
                crs=self.crs,
            )

        return cover

    def grid_from_gdf(self, gdf: gpd.GeoDataFrame, *, clip=False) -> list["Cell"]:
        pass


@dataclass(frozen=True)
class Cell:
    """
    Represents a hexagon cell within a HexGridSystem.
    """

    hgs: HexGridSystem
    a: int  # first axial coordinate of center of cell
    b: int  # second axial coordinate of center of cell

    def id(self):
        return f"{self.a},{self.b}"

    def center_axial(self) -> tuple[int]:
        return self.a, self.b

    def center_double(self) -> tuple[int]:
        return axial_to_double(self.a, self.b)

    def center(self) -> tuple[float]:
        """
        Return the Cartesian center of this cell (relative to the CRS).
        """
        return axial_to_cartesian(self.a, self.b, self.hgs.R)

    def vertices(self) -> list[tuple[float]]:
        x, y = self.center()
        return [hexagon_vertex(x, y, self.hgs.R, i) for i in range(6)]

    def polygon(self) -> sg.Polygon:
        """
        Return the Shapely polygon of this cell (relative to the CRS).
        """
        return sg.Polygon(self.vertices() + self.vertices()[:1])

    @staticmethod
    def _direction_int_to_str(direction: int) -> "Cell":
        """
        Translate the direction number modulo 6 to a direction string.
        If the direction modulo 6 is not one of the integers
        0, 1, 2, ..., 5, then return ':('.
        Raise a TypeError if taking the direction modulo 6 raises a TypeError,
        e.g. if the direction is a string.
        """
        d = {
            0: "ru",
            1: "u",
            2: "lu",
            3: "ld",
            4: "d",
            5: "rd",
        }
        return d.get(direction % 6, ":(")

    def neighbour(self, direction: str | int) -> "Cell":
        """
        Return the neighbour cell specified by one of the direction strings
        'ru', 'u', 'lu', 'lu', 'd', 'rd' or by an integer that will be taken modulo 6.
        Direction meanings are indicated in the following diagram with this cell
        in the center::

                 u                1
            lu       ru       2       0
                 .                .
            ld       rd       3       5
                 d                4

        """
        try:
            direction = self._direction_int_to_str(direction)
        except TypeError:
            # Already a string
            pass

        a, b = self.a, self.b

        if direction == "ru":
            c = a + 1, b
        elif direction == "u":
            c = a, b + 1
        elif direction == "lu":
            c = a - 1, b + 1
        elif direction == "ld":
            c = a - 1, b
        elif direction == "d":
            c = a, b - 1
        elif direction == "rd":
            c = a + 1, b - 1
        else:
            raise ValueError(f"Invalid direction {direction}")

        return Cell(self.hgs, *c)

    def neighbor(self, direction: str | int) -> "Cell":
        """
        Alias of :method:`neighbour`.
        """
        return self.neighbour(direction)


def axial_to_cartesian(a: float, b: float, R: float) -> tuple[float]:
    """
    Given axial coordinates of a point in the plane relative to a flat-top hexagonal
    grid centered at the origin with hexagon circumradius ``R``,
    return its Cartesian coordinates.

    For more details on axial coordinates see [RBG]_.
    """
    return (R * (3 / 2) * a, R * ((sqrt(3) / 2) * a + sqrt(3) * b))


def axial_to_double(a: float, b: float) -> tuple[float]:
    return a, a + 2 * b


def double_to_axial(r: float, c: float) -> tuple[float]:
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

        3           0

            4   5

    For ``i`` greater than 5 or negative, then wrap in standard mathematical fashion.
    """
    θ = i * pi / 3
    return x + R * cos(θ), y + R * sin(θ)
