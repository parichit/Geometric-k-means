geokmeans
=========

Fast, eco-friendly **k-means clustering** for Python — a scikit-learn-style
interface to a C++ library implementing seven k-means algorithms, including
**Geometric-k-means**, the bound-free method of Sharma et al. (2026).

All seven algorithms return the *same* Lloyd's k-means solution; they differ
only in how aggressively they prune distance computations, so you get identical
results with far less work.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   quickstart
   api

Citation
--------

If you use this package, please cite:

   Sharma, P., Stanislaw, M., Kurban, H., Kulekci, O., & Dalkilic, M. (2026).
   Geometric-k-means: A Bound Free Approach to Fast and Eco-Friendly k-means.
   *Machine Learning*, 115(2), 30.
   https://doi.org/10.1007/s10994-025-06891-1

Indices
-------

* :ref:`genindex`
* :ref:`search`
