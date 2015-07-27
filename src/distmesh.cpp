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
    Functional const& distanceFunction, double const initialPointDistance,
    Functional const& elementSizeFunction, Eigen::Ref<Eigen::ArrayXXd const> const boundingBox,
    Eigen::Ref<Eigen::ArrayXXd const> const fixedPoints) {
    // determine dimension of mesh
    unsigned const dimension = boundingBox.cols();

    // create initial distribution in bounding box
    Eigen::ArrayXXd points = utils::createInitialPoints(distanceFunction,
        initialPointDistance, elementSizeFunction, boundingBox, fixedPoints);

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
        // retriangulate if point movement is above threshold
        if ((points - retriangulationCriterionBuffer).square().rowwise().sum().sqrt().maxCoeff() >
            constants::retriangulationThreshold * initialPointDistance) {
            // update triangulation
            triangulation = triangulation::delaunay(points);

            // reject triangles with circumcenter outside of the region
            Eigen::ArrayXXd circumcenter = Eigen::ArrayXXd::Zero(triangulation.rows(), dimension);
            for (int point = 0; point < triangulation.cols(); ++point) {
                circumcenter += utils::selectIndexedArrayElements<double>(
                    points, triangulation.col(point)) / triangulation.cols();
            }
            triangulation = utils::selectMaskedArrayElements<int>(triangulation,
                distanceFunction(circumcenter) < -constants::geometryEvaluationThreshold * initialPointDistance);

            // find unique bar indices
            barIndices = utils::findUniqueBars(triangulation);

            // store current points positions
            retriangulationCriterionBuffer = points;
        }

        // calculate bar vectors and their length
        Eigen::ArrayXXd barVector = utils::selectIndexedArrayElements<double>(points, barIndices.col(0)) -
            utils::selectIndexedArrayElements<double>(points, barIndices.col(1));
        Eigen::ArrayXd barLength = barVector.square().rowwise().sum().sqrt();

        // evaluate elementSizeFunction at midpoints of bars
        Eigen::ArrayXd desiredElementSize = elementSizeFunction(0.5 *
            (utils::selectIndexedArrayElements<double>(points, barIndices.col(0)) +
            utils::selectIndexedArrayElements<double>(points, barIndices.col(1))));

        // calculate desired bar length
        Eigen::ArrayXd desiredBarLength = desiredElementSize * (1.0 + 0.4 / std::pow(2.0, dimension - 1)) *
            std::pow((barLength.pow(dimension).sum() / desiredElementSize.pow(dimension).sum()),
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
        utils::projectPointsToBoundary(distanceFunction, initialPointDistance, points);

        // stop, when maximum points movement is below threshold
        if ((points - stopCriterionBuffer).square().rowwise().sum().sqrt().maxCoeff() <
            constants::pointsMovementThreshold * initialPointDistance) {
            break;
        }
    }

    return std::make_tuple(points, triangulation);
}

// determine boundary edges of given triangulation
Eigen::ArrayXi distmesh::boundEdges(Eigen::Ref<Eigen::ArrayXXi const> const triangulation,
    Eigen::Ref<Eigen::ArrayXXi const> const edges) {
    // get edge indices for each triangle in triangulation
    Eigen::ArrayXXi const edgeIndices = [](Eigen::Ref<Eigen::ArrayXXi const> const triangulation,
        Eigen::Ref<Eigen::ArrayXXi const> const edges) {
        if (edges.rows() == 0) {
            return utils::getTriangulationEdgeIndices(triangulation, utils::findUniqueEdges(triangulation));
        }
        else {
            return utils::getTriangulationEdgeIndices(triangulation, edges);            
        }
    }(triangulation, edges);

    // find edges, which only appear once in triangulation
    std::set<int> uniqueEdges;
    std::vector<int> boundaryEdges;
    for (int triangle = 0; triangle < triangulation.rows(); ++triangle)
    for (int edge = 0; edge < triangulation.cols(); ++edge) {
        auto const edgeIndex = edgeIndices(triangle, edge);

        // insert edge in set to get info about multiple appearance
        if (!std::get<1>(uniqueEdges.insert(edgeIndex))) {
            // find edge in vector and delete it
            auto const it = std::find(boundaryEdges.begin(), boundaryEdges.end(), edgeIndex);
            if (it != boundaryEdges.end()) {
                boundaryEdges.erase(it);
            }
        }
        else {
            boundaryEdges.push_back(edgeIndex);
        }
    }

    // convert stl vector to eigen array
    Eigen::ArrayXi boundary(boundaryEdges.size());
    for (int edge = 0; edge < boundary.rows(); ++edge) {
        boundary(edge) = boundaryEdges[edge];
    }

    return boundary;
}
