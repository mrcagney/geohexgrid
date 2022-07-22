"""
Helpful resource: https://www.redblobgames.com/grids/hexagons
"""

import collections
from dataclasses import dataclass
from math import sqrt


START_ANGLE = 30  # degrees

Point = collections.namedtuple("Point", ["x", "y"])


@dataclass(frozen=True)
class HexSystem:
    """
    Represents a flat top hexagonal grid system of circumradius R and origin...
    """

    pass


@dataclass(frozen=True)
class HexCell:
    """
    Represents a flat top hexagon of circumradius R.
    """

    pass


def cartesian_to_axial(x: float, y: float, R: float) -> tuple(float, float):
    """ """
    return (R * (2 / 3) * x, R * (-x / 3 + (sqrt(3) / 3) * y))


def axial_to_cartesian(q: float, r: float, R: float) -> tuple(float, float):
    """ """
    return (R * (3 / 2) * q, R * ((sqrt(3) / 2) * q + sqrt(3) * r))


def hex_containing_axial_point(q: float, r: float) -> tuple(int, int):
    """
    Given axial coordinates of a point, return the (integer) axial coordinates
    of the center of the hex cell containing it.
    """
    q_grid = round(q)
    r_grid = round(r)
    q -= q_grid
    r -= r_grid  # remainder
    if abs(q) >= abs(r):
        result = int(q_grid + round(q + 0.5 * r)), int(r_grid)
    else:
        result = int(q_grid), int(r_grid + round(r + 0.5 * q))

    return result


def hex_containing_cartesian_point(x: float, y: float, R: float) -> tuple(int, int):
    """
    Return the axial coordinates ``(q, r)`` of the of a hex cell containing the
    Cartesian point ``(x, y)``.
    """
    return hex_containing_axial_point(*cartesian_to_axial(x, y, R))
