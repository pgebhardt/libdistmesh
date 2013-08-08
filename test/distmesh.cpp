#include "gtest/gtest.h"
#include "distmesh/distmesh.h"

TEST(DistmeshTest, Boundedges) {
    // fill points array with data
    auto points =std::make_shared<distmesh::dtype::array<
        distmesh::dtype::real>>(4, 2);
    (*points)(0, 0) = 1.0; (*points)(0, 1) = 0.0;
    (*points)(1, 0) = 0.0; (*points)(1, 1) = 1.0;
    (*points)(2, 0) = -1.0; (*points)(2, 1) = 0.0;
    (*points)(3, 0) = 0.0; (*points)(3, 1) = -1.0;
    auto triangulation = distmesh::triangulation::delaunay(points);

    // calc boundedges
    auto boundary = distmesh::boundedges(triangulation);
    std::cout << "boundary:" << std::endl;
    std::cout << *boundary << std::endl;
};
