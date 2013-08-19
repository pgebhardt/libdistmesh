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

#include <iostream>

TEST(UtilsTest, SelectMaskedArrayElements) {
    distmesh::dtype::array<distmesh::dtype::real> array(4, 2);
    array << 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0;

    distmesh::dtype::array<bool> mask(4, 1);
    mask << true, false, false, true;

    auto masked_array = distmesh::utils::select_masked_array_elements<
        distmesh::dtype::real>(array, mask);

    distmesh::dtype::array<distmesh::dtype::real> test_masked_array(2, 2);
    test_masked_array << 1.0, 2.0, 4.0, 5.0;

    EXPECT_EQ(masked_array.rows(), 2u);
    EXPECT_EQ(masked_array.cols(), 2u);
    EXPECT_TRUE((masked_array == test_masked_array).all());
};

TEST(UtilsTest, SelectMaskedIndexedElements) {
    distmesh::dtype::array<distmesh::dtype::real> array(4, 2);
    array << 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0;

    distmesh::dtype::array<distmesh::dtype::index> indices(2, 1);
    indices << 1, 3;

    auto indexed_array = distmesh::utils::select_indexed_array_elements<
        distmesh::dtype::real>(array, indices);

    distmesh::dtype::array<distmesh::dtype::real> test_indexed_array(2, 2);
    test_indexed_array << 2.0, 3.0, 4.0, 5.0;

    EXPECT_EQ(indexed_array.rows(), 2u);
    EXPECT_EQ(indexed_array.cols(), 2u);
    EXPECT_TRUE((indexed_array == test_indexed_array).all());
};

TEST(UtilsTest, Factorial) {
    EXPECT_EQ(distmesh::utils::factorial(0), 1u);
    EXPECT_EQ(distmesh::utils::factorial(1), 1u);
    EXPECT_EQ(distmesh::utils::factorial(3), 6u);
    EXPECT_EQ(distmesh::utils::factorial(4), 24u);
};

TEST(UtilsTest, NOverK) {
    // test 2 out of 4
    auto combinations = distmesh::utils::n_over_k(4, 2);

    // check combination count
    EXPECT_EQ(combinations.rows(), distmesh::utils::factorial(4) /
        (distmesh::utils::factorial(2) * distmesh::utils::factorial(4 - 2)));

    // check combinations
    distmesh::dtype::array<distmesh::dtype::index> test_combinations(6, 2);
    test_combinations << 0, 1,
                         0, 2,
                         0, 3,
                         1, 2,
                         1, 3,
                         2, 3;
    EXPECT_TRUE((combinations == test_combinations).all());
};

TEST(UtilsTest, PointsInsidePoly) {
    distmesh::dtype::array<distmesh::dtype::real> polygon(4, 2);
    polygon << -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0;

    distmesh::dtype::array<distmesh::dtype::real> points(3, 2);
    points << 0.0, 0.0, -2.0, 0.0, 0.0, 4.0;

    auto inside = distmesh::utils::points_inside_poly(points, polygon);

    distmesh::dtype::array<bool> test_inside(3, 1);
    test_inside << true, false, false;

    EXPECT_TRUE((inside == test_inside).all());
};
