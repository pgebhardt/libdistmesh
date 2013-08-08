libDistMesh: A Simple Mesh Generator in C++
==============================================

libDistMesh is a C++ implementation of the original [DistMesh](http://persson.berkeley.edu/distmesh/)
algorithm for generating unstructured triangular and tetrahedral meshes using *signed distance functions*.

Getting Started
---------------

Simply clone the repository, update submodules and make shure all dependencies are installed.
For building the project the [SCons](http://www.scons.org/) build system is used:

    git submodule update --init
    scons
    scons install

Example
-------

* Uniform Mesh on Unit Circle::

```
// bounding box in which the algorithm tries to create points
distmesh::dtype::array<distmesh::dtype::real> bounding_box(2, 2);
bounding_box << -1.0, 1.0, -1.0, 1.0;

// midpoint of unit circle
distmesh::dtype::array<distmesh::dtype::real> midpoint(1, 2);
midpoint.fill(0.0);

// create mesh
auto mesh = distmesh::distmesh(
    distmesh::distance_functions::circular(midpoint, 1.0),
    distmesh::edge_length_functions::uniform(),
    0.02, bounding_box);
```

Dependencies
------------

libDistMesh uses some C++11 features and compiles proper with both clang
and gcc. For linear algebra operations and the delaunay triangulation two
libraries are needed for building and using libDistMesh:

* [Eigen](http://eigen.tuxfamily.org/)
* [QHull](http://www.qhull.org/) (building only)

References
----------

The DistMesh algorithm is described in the following two references.
If you use the algorithm in a program or publication, please
acknowledge its authors by adding a reference to the first paper
below.

* P.-O. Persson, G. Strang, **A Simple Mesh Generator in MATLAB**.
  *SIAM Review*, Volume 46 (2), pp. 329-345, June 2004 ([PDF]
  (http://persson.berkeley.edu/distmesh/persson04mesh.pdf>))

* P.-O. Persson, **Mesh Generation for Implicit Geometries**.
  Ph.D. thesis, *Department of Mathematics, MIT*, Dec 2004 ([PDF]
  (http://persson.berkeley.edu/thesis/persson-thesis-color.pdf))
