libDistMesh: A Simple Mesh Generator in C++
==============================================

libDistMesh is a C++ implementation of the original [DistMesh](http://persson.berkeley.edu/distmesh/)
algorithm for generating unstructured triangular and tetrahedral meshes using *signed distance functions*.

Getting Started
---------------

Simply clone the repository, update submodules and make sure all dependencies are installed.
For building the project the [SCons](http://www.scons.org/) build system is used:

    git clone https://github.com/schansge/libdistmesh.git
    make && make install

Example
-------

* Uniform Mesh on Unit Circle:

```
#include <distmesh/distmesh.h>

int main() {
    // create mesh
    auto mesh = distmesh::distmesh(distmesh::distance_function::circular(1.0), 0.2);

    return 0;
}
```
* Rectangle with circular hole, refined at circle boundary:

```
#include <distmesh/distmesh.h>

int main() {
    // fixed points at the corners of domain to guarantee convergence
    Eigen::ArrayXXd fixed_points(4, 2);
    fixed_points << -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0;

    // create mesh
    auto mesh = distmesh::distmesh(
        distmesh::distance_function::rectangular(rectangle)
            .max(-distmesh::distance_function::circular(0.5)),
        0.05, 0.05 + 0.3 * distmesh::distance_function::circular(0.5),
        distmesh::bounding_box(2), fixed_points);

    return 0;
}
```

Dependencies
------------

libDistMesh uses some C++11 features and compiles properly with both clang
and gcc. For linear algebra operations and the delaunay triangulation two
libraries are needed for building and using libDistMesh:

* [Eigen](http://eigen.tuxfamily.org/) >= 3.2.0
* [QHull](http://www.qhull.org/) >= 2012.1

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
