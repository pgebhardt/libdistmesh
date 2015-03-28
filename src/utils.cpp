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
inline unsigned distmesh::utils::factorial(unsigned const n) {
    if (n <= 1) {
        return 1;
    }
    else {
        return n * factorial(n - 1);
    }
}

// create initial points distribution
Eigen::ArrayXXd distmesh::utils::createInitialPoints(
    Functional const& distanceFunction, double const initialPointDistance,
    Functional const& elementSizeFunction, Eigen::Ref<Eigen::ArrayXXd const> const boundingBox,
    Eigen::Ref<Eigen::ArrayXXd const> const fixedPoints) {
    // extract dimension of mesh
    unsigned const dimension = boundingBox.cols();

    // initially distribute points evenly in complete bounding box
    Eigen::ArrayXi pointsPerDimension(dimension);
    for (int dim = 0; dim < dimension; ++dim) {
        pointsPerDimension(dim) = ceil((boundingBox(1, dim) - boundingBox(0, dim)) /
            (initialPointDistance * (dim == 0 ? 1.0 : sqrt(3.0) / 2.0)));
    }

    Eigen::ArrayXXd points(pointsPerDimension.prod(), dimension);
    for (int point = 0; point < points.rows(); ++point)
    for (int dim = 0; dim < dimension; ++dim) {
        int const pointIndex = (point / std::max(pointsPerDimension.topRows(dim).prod(), 1)) %
            pointsPerDimension(dim);

        points(point, dim) = boundingBox(0, dim) + (double)pointIndex * initialPointDistance *
            (dim == 0 ? 1.0 : sqrt(3.0) / 2.0);

        if (dim > 0) {
            points(point, dim - 1) += pointIndex % 2 != 0 ? initialPointDistance / 2.0 : 0.0;
        }
    }

    // reject points outside of region defined by distance function
    points = selectMaskedArrayElements<double>(points,
        distanceFunction(points) < constants::geometryEvaluationTolerance * initialPointDistance);

    // clear duplicate points
    Eigen::Array<bool, Eigen::Dynamic, 1> isUniquePoint =
        Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(points.rows(), true);
    for (int i = 0; i < fixedPoints.rows(); ++i)
    for (int j = 0; j < points.rows(); ++j) {
        isUniquePoint(j) &= !(fixedPoints.row(i) == points.row(j)).all();
    }
    points = selectMaskedArrayElements<double>(points, isUniquePoint);

    // calculate probability to keep points
    Eigen::ArrayXd probability = 1.0 / elementSizeFunction(points).pow(dimension);
    probability /= probability.maxCoeff();

    // reject points with wrong probability
    points = selectMaskedArrayElements<double>(points,
        0.5 * (1.0 + Eigen::ArrayXd::Random(points.rows())) < probability);

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
    Functional const& distanceFunction, double const initialPointDistance,
    Eigen::Ref<Eigen::ArrayXXd> points) {
    Eigen::ArrayXd distance = distanceFunction(points);

    // check for points outside of boundary
    Eigen::Array<bool, Eigen::Dynamic, 1> outside = distance > 0.0;
    if (outside.any()) {
        // calculate gradient
        Eigen::ArrayXXd gradient(points.rows(), points.cols());
        Eigen::ArrayXXd deltaX = Eigen::ArrayXXd::Zero(points.rows(), points.cols());

        for (int dim = 0; dim < points.cols(); ++dim) {
            deltaX.col(dim).fill(constants::deltaX * initialPointDistance);
            gradient.col(dim) = (distanceFunction(points + deltaX) - distance) /
                (constants::deltaX * initialPointDistance);
            deltaX.col(dim).fill(0.0);
        }

        // project points back to boundary
        points -= outside.replicate(1, points.cols()).select(
            gradient.colwise() * distance / gradient.square().rowwise().sum(), 0.0);
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
