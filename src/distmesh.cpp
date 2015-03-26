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

#include <vector>
#include <set>
#include <algorithm>

#include "distmesh/distmesh.h"
#include "distmesh/constants.h"
#include "distmesh/utils.h"
#include "distmesh/triangulation.h"

// easy creation of n-dimensional bounding box
Eigen::ArrayXXd distmesh::boundingBox(unsigned const dimension) {
    Eigen::ArrayXXd box(2, dimension);
    box.row(0).fill(-1.0);
    box.row(1).fill(1.0);
    return box;
}

// apply the distmesh algorithm
std::tuple<Eigen::ArrayXXd, Eigen::ArrayXXi> distmesh::distmesh(
    Functional const& distanceFunction, double const baseEdgeLength,
    Functional const& edgeLengthFunction, Eigen::Ref<Eigen::ArrayXXd const> const boundingBox,
    Eigen::Ref<Eigen::ArrayXXd const> const fixedPoints) {
    // determine dimension of mesh
    unsigned const dimension = boundingBox.rows();

    // create initial distribution in bounding box
    Eigen::ArrayXXd points = utils::createInitialPoints(distanceFunction,
        baseEdgeLength, edgeLengthFunction, boundingBox, fixedPoints);

    // create initial triangulation
    Eigen::ArrayXXi triangulation = triangulation::delaunay(points);

    // create buffer to store old point locations to calculate
    // retriangulation and stop criterion
    Eigen::ArrayXXd retriangulationCriterionBuffer = Eigen::ArrayXXd::Constant(
        points.rows(), points.cols(), INFINITY);
    Eigen::ArrayXXd stopCriterionBuffer = Eigen::ArrayXXd::Zero(
        points.rows(), points.cols());

    // main distmesh loop
    Eigen::ArrayXXi barIndices;
    for (unsigned step = 0; step < constants::maxSteps; ++step) {
        // retriangulate if point movement is above tolerance
        double const retriangulationCriterion = (points - retriangulationCriterionBuffer)
            .square().rowwise().sum().sqrt().maxCoeff();
        if (retriangulationCriterion > constants::retriangulationTolerance * baseEdgeLength) {
            // update triangulation
            triangulation = triangulation::delaunay(points);

            // calculate circumcenter
            Eigen::ArrayXXd circumcenter = Eigen::ArrayXXd::Zero(
                triangulation.rows(), dimension);
            for (int point = 0; point < triangulation.cols(); ++point) {
                circumcenter += utils::selectIndexedArrayElements<double>(
                    points, triangulation.col(point)) / triangulation.cols();
            }

            // reject triangles with circumcenter outside of the region
            triangulation = utils::selectMaskedArrayElements<int>(triangulation,
                distanceFunction(circumcenter) < -constants::geomertyEvaluationTolerance * baseEdgeLength);

            // find unique bar indices
            barIndices = utils::findUniqueBars(triangulation);

            // store current points positions
            retriangulationCriterionBuffer = points;
        }

        // calculate bar vectors and their length
        Eigen::ArrayXXd barVector = utils::selectIndexedArrayElements<double>(points, barIndices.col(0)) -
            utils::selectIndexedArrayElements<double>(points, barIndices.col(1));
        Eigen::ArrayXd barLength = barVector.square().rowwise().sum().sqrt();

        // evaluate edgeLengthFunction at midpoints of bars
        Eigen::ArrayXd hbars = edgeLengthFunction(0.5 *
            (utils::selectIndexedArrayElements<double>(points, barIndices.col(0)) +
            utils::selectIndexedArrayElements<double>(points, barIndices.col(1))));

        // calculate desired bar length
        Eigen::ArrayXd desiredBarLength = hbars * (1.0 + 0.4 / std::pow(2.0, dimension - 1)) *
            std::pow((barLength.pow(dimension).sum() / hbars.pow(dimension).sum()),
                1.0 / dimension);

        // calculate force vector for each bar
        Eigen::ArrayXXd forceVector = barVector.colwise() *
            ((desiredBarLength - barLength) / barLength).max(0.0);

        // store current points positions
        stopCriterionBuffer = points;

        // move points
        for (int bar = 0; bar < barIndices.rows(); ++bar) {
            if (barIndices(bar, 0) >= fixedPoints.rows()) {
                points.row(barIndices(bar, 0)) += constants::deltaT * forceVector.row(bar);
            }
            if (barIndices(bar, 1) >= fixedPoints.rows()) {
                points.row(barIndices(bar, 1)) -= constants::deltaT * forceVector.row(bar);
            }
        }

        // project points outside of domain to boundary
        utils::projectPointsToFunction(distanceFunction, baseEdgeLength, points);

        // stop criterion
        double const stopCriterion = (points - stopCriterionBuffer).square().rowwise().sum().sqrt().maxCoeff();
        if (stopCriterion < constants::pointsMovementTolerance * baseEdgeLength) {
            break;
        }
    }

    return std::make_tuple(points, triangulation);
}

// determine boundary edges of given triangulation
Eigen::ArrayXXi distmesh::boundEdges(Eigen::Ref<Eigen::ArrayXXi const> const triangulation) {
    std::set<std::vector<int>> uniqueEdges;
    std::vector<std::vector<int>> boundaryEdges;
    std::vector<int> edge(triangulation.cols() - 1);

    // find all possible combinations of edges
    Eigen::ArrayXXi combinations = utils::nOverK(triangulation.cols(), triangulation.cols() - 1);

    // find edges, which only appear once in triangulation
    for (int combination = 0; combination < combinations.rows(); ++combination)
    for (int triangle = 0; triangle < triangulation.rows(); ++triangle) {
        // get current edge
        for (size_t point = 0; point < edge.size(); ++point) {
            edge[point] = triangulation(triangle, combinations(combination, point));
        }

        // insert edge in set to get info about multiple appearance
        if (!std::get<1>(uniqueEdges.insert(edge))) {
            // find edge in vector and delete it
            auto it = std::find(boundaryEdges.begin(), boundaryEdges.end(), edge);
            if (it != boundaryEdges.end()) {
                boundaryEdges.erase(it);
            }
        }
        else {
            boundaryEdges.push_back(edge);
        }
    }

    // convert stl vector to eigen array
    Eigen::ArrayXXi boundary(boundaryEdges.size(), edge.size());
    for (int edge = 0; edge < boundary.rows(); ++edge)
    for (int point = 0; point < boundary.cols(); ++point) {
        boundary(edge, point) = boundaryEdges[edge][point];
    }

    return boundary;
}
