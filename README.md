libDistMesh: A Simple Mesh Generator in C++
==============================================

libDistMesh is a C++ implementation of the original [DistMesh](http://persson.berkeley.edu/distmesh/)
algorithm for generating unstructured triangular and tetrahedral meshes using *signed distance functions*.

Getting Started
---------------

Simply clone the repository, make sure all dependencies are installed and build it.

    git clone https://github.com/schansge/libdistmesh.git
    cd libdistmesh
    cp Makefile.config.example Makefile.config
    make
    make install

Example
-------

* Uniform Mesh on Unit Circle:

```C++
#include <distmesh/distmesh.h>

int main() {
    // create mesh
    auto mesh = distmesh::distmesh(distmesh::distanceFunction::circular(1.0), 0.2);

    return 0;
}
```
* Rectangle with circular hole, refined at circle boundary:

```C++
#include <distmesh/distmesh.h>

int main() {
    // fixed points at the corners of domain to guarantee convergence
    Eigen::ArrayXXd fixedPoints(4, 2);
    fixedPoints << -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0;

    // create mesh
    auto mesh = distmesh::distmesh(
        distmesh::distanceFunction::rectangular(rectangle)
            .max(-distmesh::distanceFunction::circular(0.5)),
        0.05, 0.05 + 0.3 * distmesh::distanceFunction::circular(0.5),
        distmesh::boundingBox(2), fixedPoints);

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
