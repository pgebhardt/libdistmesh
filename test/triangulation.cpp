#include "gtest/gtest.h"
#include "distmesh/distmesh.h"

TEST(TriangulationTest, Delaunay) {
    // fill points array with data
    auto points =std::make_shared<distmesh::dtype::array<
        distmesh::dtype::real>>(4, 2);
    (*points)(0, 0) = 1.0; (*points)(0, 1) = 0.0;
    (*points)(1, 0) = 0.0; (*points)(1, 1) = 1.0;
    (*points)(2, 0) = -1.0; (*points)(2, 1) = 0.0;
    (*points)(3, 0) = 0.0; (*points)(3, 1) = -1.0;

    // calc delaunay triangulation
    std::shared_ptr<distmesh::dtype::array<distmesh::dtype::index>>
        triangulation = nullptr;
    EXPECT_NO_THROW({
        triangulation = distmesh::triangulation::delaunay(points);
    });

    // check triangulation
    EXPECT_EQ(triangulation->rows(), 2);
    EXPECT_EQ(triangulation->cols(), 3);

    distmesh::dtype::array<distmesh::dtype::index> triangle(1, 3);
    triangle << 3, 1, 0;
    EXPECT_EQ((triangulation->row(0) == triangle).all(), true);
    triangle << 3, 1, 2;
    EXPECT_EQ((triangulation->row(1) == triangle).all(), true);
};
