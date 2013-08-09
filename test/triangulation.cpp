// --------------------------------------------------------------------
// This file is part of libDistMesh.
//
// libDistMesh is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// libDistMesh is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with libDistMesh.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright (C) 2013 Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de
// --------------------------------------------------------------------

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
