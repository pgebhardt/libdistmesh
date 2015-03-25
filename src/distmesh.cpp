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

#include "distmesh/distmesh.h"
#include <vector>
#include <set>
#include <algorithm>

// easy creation of n-dimensional bounding_box
Eigen::ArrayXXd distmesh::boundingBox(unsigned const dimension) {
    Eigen::ArrayXXd box(dimension, 2);
    box.col(0).fill(-1.0);
    box.col(1).fill(1.0);
    return box;
}

// apply the distmesh algorithm
std::tuple<Eigen::ArrayXXd, Eigen::ArrayXXi> distmesh::distmesh(
    Functional const& distanceFunction, double const baseEdgeLength,
    Functional const& edgeLengthFunction, Eigen::Ref<Eigen::ArrayXXd const> const boundingBox,
    Eigen::Ref<Eigen::ArrayXXd const> const fixedPoints) {
    // create initial distribution in bounding box
    Eigen::ArrayXXd points = utils::createInitialPoints(distanceFunction,
        baseEdgeLength, edgeLengthFunction, boundingBox, fixedPoints);

    // create initial triangulation
    Eigen::ArrayXXi triangulation = triangulation::delaunay(points);

    // create points buffer for retriangulation and stop criterion
    Eigen::ArrayXXd buffer_retriangulation_criterion(
        points.rows(), points.cols());
    Eigen::ArrayXXd buffer_stop_criterion;
    buffer_retriangulation_criterion.fill(INFINITY);

    // main distmesh loop
    Eigen::ArrayXXi bar_indices;
    for (unsigned step = 0; step < settings::maxSteps; ++step) {
        // retriangulate if point movement is above tolerance
        double retriangulationCriterion = (points - buffer_retriangulation_criterion)
            .square().rowwise().sum().sqrt().maxCoeff();
        if (retriangulationCriterion > settings::retriangulationTolerance * baseEdgeLength) {
            // update triangulation
            triangulation = triangulation::delaunay(points);

            // calculate circumcenter
            Eigen::ArrayXXd circumcenter = Eigen::ArrayXXd::Zero(
                triangulation.rows(), points.cols());
            for (int point = 0; point < triangulation.cols(); ++point) {
                circumcenter += utils::selectIndexedArrayElements<double>(
                    points, triangulation.col(point)) / triangulation.cols();
            }

            // reject triangles with circumcenter outside of the region
            triangulation = utils::selectMaskedArrayElements<int>(triangulation,
                distanceFunction(circumcenter) < -settings::generalPrecision * baseEdgeLength);

            // find unique bar indices
            bar_indices = utils::findUniqueBars(triangulation);

            buffer_retriangulation_criterion = points;
        }

        // calculate bar vectors and their length
        Eigen::ArrayXXd bar_vector = utils::selectIndexedArrayElements<double>(points, bar_indices.col(0)) -
            utils::selectIndexedArrayElements<double>(points, bar_indices.col(1));
        Eigen::ArrayXXd bar_length = bar_vector.square().rowwise().sum().sqrt();

        // evaluate edgeLengthFunction at midpoints of bars
        Eigen::ArrayXXd hbars = edgeLengthFunction(0.5 *
            (utils::selectIndexedArrayElements<double>(points, bar_indices.col(0)) +
            utils::selectIndexedArrayElements<double>(points, bar_indices.col(1))));

        // calculate desired bar length
        Eigen::ArrayXXd desired_bar_length = hbars *
            (1.0 + 0.4 / std::pow(2.0, points.cols() - 1)) *
            std::pow((bar_length.pow(points.cols()).sum() /
            hbars.pow(points.cols()).sum()), 1.0 / points.cols());

        // calculate force vector for each bar
        Eigen::ArrayXXd force_vector = bar_vector.colwise() *
            ((desired_bar_length - bar_length) / bar_length).max(0.0).col(0);

        // move points
        buffer_stop_criterion = points;
        for (int bar = 0; bar < bar_indices.rows(); ++bar) {
            if (bar_indices(bar, 0) >= fixedPoints.rows()) {
                points.row(bar_indices(bar, 0)) += settings::deltaT * force_vector.row(bar);
            }
            if (bar_indices(bar, 1) >= fixedPoints.rows()) {
                points.row(bar_indices(bar, 1)) -= settings::deltaT * force_vector.row(bar);
            }
        }

        // project points outside of domain to boundary
        utils::projectPointsToFunction(distanceFunction,
            baseEdgeLength, points);

        // stop criterion
        double stopCriterion = (points - buffer_stop_criterion).square().rowwise().sum().sqrt().maxCoeff();
        if (stopCriterion < settings::pointMovementTolerance * baseEdgeLength) {
            break;
        }
    }

    return std::make_tuple(points, triangulation);
}

// determine boundary edges of given triangulation
Eigen::ArrayXXi distmesh::boundEdges(Eigen::Ref<Eigen::ArrayXXi const> const triangulation) {
    std::set<std::vector<int>> edge_set;
    std::vector<std::vector<int>> boundary_edges;
    std::vector<int> edge(triangulation.cols() - 1);

    // find all possible combinations of edges
    auto combinations = utils::nOverK(triangulation.cols(), triangulation.cols() - 1);

    // find edges, which only appear once in triangulation
    for (int combination = 0; combination < combinations.rows(); ++combination)
    for (int triangle = 0; triangle < triangulation.rows(); ++triangle) {
        // get current edge
        for (size_t vertex = 0; vertex < edge.size(); ++vertex) {
            edge[vertex] = triangulation(triangle, combinations(combination, vertex));
        }

        // insert edge in set to get info about multiple appearance
        if (!std::get<1>(edge_set.insert(edge))) {
            // find edge in vector and delete it
            auto it = std::find(boundary_edges.begin(), boundary_edges.end(), edge);
            if (it != boundary_edges.end()) {
                boundary_edges.erase(it);
            }
        } else {
            boundary_edges.push_back(edge);
        }
    }

    // convert stl vector to eigen array
    Eigen::ArrayXXi boundary_array(boundary_edges.size(), edge.size());
    for (int edge = 0; edge < boundary_array.rows(); ++edge)
    for (int vertex = 0; vertex < boundary_array.cols(); ++vertex) {
        boundary_array(edge, vertex) = boundary_edges[edge][vertex];
    }

    return boundary_array;
}
