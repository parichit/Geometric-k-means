Quickstart
==========

.. code-block:: python

   import numpy as np
   from geokmeans import GeoKMeans

   X = np.random.default_rng(0).random((1000, 20))

   km = GeoKMeans(n_clusters=8, algorithm="geo", random_state=0).fit(X)

   km.labels_                    # cluster index per point
   km.cluster_centers_           # centroids
   km.n_iter_                    # iterations to convergence
   km.n_distance_calculations_   # distance computations performed

   km.predict(np.random.random((5, 20)))   # assign new points

Choosing an algorithm
---------------------

Pass ``algorithm=`` to select a variant. All produce the same clustering for a
given seed:

``"geo"`` (default), ``"lloyd"``, ``"elkan"``, ``"hamerly"``, ``"annulus"``,
``"exponion"``, ``"ball"``.

Long-running jobs
-----------------

The compute runs in C++ with the GIL released, so a long ``fit`` does not
freeze the interpreter and can run off the main thread:

.. code-block:: python

   from concurrent.futures import ProcessPoolExecutor

   def cluster(X):
       return GeoKMeans(n_clusters=50, random_state=0).fit(X)

   with ProcessPoolExecutor() as pool:
       result = pool.submit(cluster, X).result()

Pass ``verbose=True`` to stream convergence progress.
