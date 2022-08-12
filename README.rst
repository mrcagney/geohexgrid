Geohex
******
A Python library to cast evil spells upon the Earth.
Just kidding!
Humanity took care of that already.

Rather, a Python 3.9+ library to make geographically local hexagon grids of any resolution, like those produced by QGIS's `create grid function <https://docs.qgis.org/3.22/en/docs/user_manual/processing_algs/qgis/vectorcreation.html?highlight=create%20grid#create-grid>`_.
Not designed for making `discrete global grid systems <https://en.wikipedia.org/wiki/Discrete_global_grid>`_ like Uber's H3.

.. image:: geohex.png
  :width: 400
  :alt: hexagon grid of 10,000-metre circumradius covering New Zealand


Contributors
============
- Alex Raichev (2014-09), maintainer


Installation
============
Create a Python 3.9+ virtual environment and install from PyPI, e.g. via ``poetry add geohex``.


Examples
=========
See the Jupyter notebook at ``notebooks/examples.ipynb``.


Notes
======
- This project's development status is Alpha.
  Alex uses this project for work and changes it breakingly when it suits his needs.
- This project uses semantic versioning.
- Thanks to `MRCagney <https://mrcagney.com>`_ for periodically funding this project.
- Red Blog Games has a `great write up of non-geographic hexagon grids <https://www.redblobgames.com/grids/hexagons>`_.


Changes
=======

1.0.0, 2022-08-12
-----------------
- First release.