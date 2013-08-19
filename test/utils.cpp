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
