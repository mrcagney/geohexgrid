from .context import geohex
from geohex import *

import numpy as np


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
