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
// Copyright (C) 2015 Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de
// --------------------------------------------------------------------

#include <set>
#include <array>

#include "distmesh/distmesh.h"
#include "distmesh/constants.h"
#include "distmesh/utils.h"

// calculate factorial recursively
unsigned distmesh::utils::factorial(unsigned const n) {
    if (n <= 1) {
        return 1;
    }
    else {
        return n * factorial(n - 1);
    }
}

// create initial points distribution
Eigen::ArrayXXd distmesh::utils::createInitialPoints(
    Functional const& distanceFunction, double const baseEdgeLength,
    Functional const& edgeLengthFunction, Eigen::Ref<Eigen::ArrayXXd const> const boundingBox,
    Eigen::Ref<Eigen::ArrayXXd const> const fixedPoints) {
    // extract dimension of mesh
    unsigned const dimension = boundingBox.cols();

    // calculate max number of points per dimension and
    // max total point coun and create initial array
    Eigen::ArrayXi maxPointsPerDimension(dimension);
    int maxPointCount = 1;
    for (int dim = 0; dim < dimension; ++dim) {
        maxPointsPerDimension(dim) = ceil((boundingBox(1, dim) - boundingBox(0, dim)) /
            (baseEdgeLength * (dim == 0 ? 1.0 : sqrt(3.0) / 2.0)));
        maxPointCount *= maxPointsPerDimension(dim);
    }

    // fill point list with evenly distributed points
    Eigen::ArrayXXd points(maxPointCount, dimension);
    int sameValueCount = 1;
    for (int dim = 0; dim < dimension; ++dim) {
        for (int point = 0; point < maxPointCount; ++point) {
            points(point, dim) = boundingBox(0, dim) +
                baseEdgeLength * (dim == 0 ? 1.0 : sqrt(3.0) / 2.0) *
                ((point / sameValueCount) % maxPointsPerDimension(dim));
            if (dim > 0) {
                points(point, dim - 1) += point / sameValueCount % 2 != 0 ? baseEdgeLength / 2.0 : 0.0;
            }
        }
        sameValueCount *= maxPointsPerDimension(dim);
    }

    // reject points outside of region defined by distance function
    points = selectMaskedArrayElements<double>(points,
        distanceFunction(points) < constants::geomertyEvaluationTolerance * baseEdgeLength);

    // clear dublicate points
    Eigen::Array<bool, Eigen::Dynamic, 1> isDublicate =
        Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(points.rows(), true);
    for (int i = 0; i < fixedPoints.rows(); ++i)
    for (int j = 0; j < points.rows(); ++j) {
        isDublicate(j) &= !(fixedPoints.row(i) == points.row(j)).all();
    }
    points = selectMaskedArrayElements<double>(points, isDublicate);

    // calculate propability to keep points
    Eigen::ArrayXd propability = 1.0 / edgeLengthFunction(points).pow(dimension);
    propability /= propability.maxCoeff();

    // reject points with wrong propability
    points = selectMaskedArrayElements<double>(points,
        0.5 * (1.0 + Eigen::ArrayXd::Random(points.rows())) < propability);

    // combine fixed and variable points to one array
    Eigen::ArrayXXd finalPoints(points.rows() + fixedPoints.rows(), dimension);
    finalPoints << fixedPoints, points;

    return finalPoints;
}

// create array with all unique combinations n over k
Eigen::ArrayXXi distmesh::utils::nOverK(unsigned const n, unsigned const k) {
    // fill an array with all unique combinations n over k,
    // starting with 0, 1, 2, ...
    Eigen::ArrayXXi combinations = Eigen::ArrayXXi::Zero(factorial(n) /
            (factorial(k) * factorial(n - k)), k);
    combinations.row(0).setLinSpaced(k, 0, k - 1);

    for (int combination = 1; combination < combinations.rows(); ++combination) {
        combinations.row(combination) = combinations.row(combination - 1);
        for (int col = k - 1; col >= 0; --col) {
            combinations.block(combination, col, 1, k - col).row(0)
                .setLinSpaced(k - col, combinations(combination, col) + 1,
                    combinations(combination, col) + k - col);
            if (combinations(combination, k - 1) < (int)n) {
                break;
            }
        }
    }

    return combinations;
}

// find unique bars
Eigen::ArrayXXi distmesh::utils::findUniqueBars(Eigen::Ref<Eigen::ArrayXXi const> const triangulation) {
    // find all unique combinations
    auto combinations = nOverK(triangulation.cols(), 2);

    // find unique bars for all combinations
    std::set<std::array<int, 2>> uniqueBars;
    std::array<int, 2> bar = {{0, 0}};
    for (int combination = 0; combination < combinations.rows(); ++combination)
    for (int triangle = 0; triangle < triangulation.rows(); ++triangle) {
        bar[0] = triangulation(triangle, combinations(combination, 0));
        bar[1] = triangulation(triangle, combinations(combination, 1));

        uniqueBars.insert(bar);
    }

    // copy set to eigen array
    Eigen::ArrayXXi barIndices(uniqueBars.size(), 2);
    int index = 0;
    for (auto const& bar : uniqueBars) {
        barIndices(index, 0) = bar[0];
        barIndices(index, 1) = bar[1];
        index++;
    }

    return barIndices;
}

// project points outside of boundary back to it
void distmesh::utils::projectPointsToFunction(
    Functional const& distanceFunction, double const baseEdgeLength,
    Eigen::Ref<Eigen::ArrayXXd> points) {
    // evaluate distance function at points
    Eigen::ArrayXd distance = distanceFunction(points);

    // check for points outside of boundary
    Eigen::Array<bool, Eigen::Dynamic, 1> outside = distance > 0.0;
    if (outside.any()) {
        // calculate gradient
        Eigen::ArrayXXd gradient(points.rows(), points.cols());
        Eigen::ArrayXXd deltaX(points.rows(), points.cols());
        Eigen::ArrayXXd h;
        deltaX.fill(0.0);
        for (int dim = 0; dim < points.cols(); ++dim) {
            deltaX.col(dim).fill(std::sqrt(std::numeric_limits<double>::epsilon()) * baseEdgeLength);
            h = points + deltaX;
            gradient.col(dim) = (distanceFunction(h) - distance) /
                (std::sqrt(std::numeric_limits<double>::epsilon()) * baseEdgeLength);
            deltaX.col(dim).fill(0.0);
        }

        for (int dim = 0; dim < points.cols(); ++dim) {
            points.col(dim) -= outside.select(
                gradient.col(dim) * distance / gradient.square().rowwise().sum(),
                0.0);
        }
    }
}

// check whether points lies inside or outside of polygon
Eigen::ArrayXXd distmesh::utils::pointsInsidePoly(
    Eigen::Ref<Eigen::ArrayXXd const> const points,
    Eigen::Ref<Eigen::ArrayXXd const> const polygon) {
    Eigen::ArrayXd inside = Eigen::ArrayXd::Zero(points.rows());

    for (int i = 0, j = polygon.rows() - 1; i < polygon.rows(); j = i++) {
        inside = (((points.col(1) < polygon(i, 1)) != (points.col(1) < polygon(j, 1))) &&
            (points.col(0) < (polygon(j, 0) - polygon(i, 0)) * (points.col(1) - polygon(i, 1)) /
            (polygon(j, 1) - polygon(i, 1)) + polygon(i, 0))).select(1.0 - inside, inside);
    }

    return inside;
}
