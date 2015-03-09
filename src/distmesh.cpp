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

#include "distmesh/distmesh.h"
#include <vector>
#include <set>
#include <algorithm>

// easy creation of n-dimensional bounding_box
Eigen::ArrayXXd distmesh::bounding_box(unsigned dimension) {
    Eigen::ArrayXXd box(dimension, 2);
    box.col(0).fill(-1.0);
    box.col(1).fill(1.0);
    return box;
}

// apply the distmesh algorithm
std::tuple<Eigen::ArrayXXd, Eigen::ArrayXXi> distmesh::distmesh(
    Functional distance_function, double edge_length_base,
    Functional edge_length_function, Eigen::Ref<const Eigen::ArrayXXd> bounding_box,
    Eigen::Ref<const Eigen::ArrayXXd> fixed_points) {
    // create initial distribution in bounding_box
    Eigen::ArrayXXd points = utils::create_point_list(distance_function,
        edge_length_base, edge_length_function, bounding_box, fixed_points);

    // create initial triangulation
    Eigen::ArrayXXi triangulation = triangulation::delaunay(points);

    // create points buffer for retriangulation and stop criterion
    Eigen::ArrayXXd buffer_retriangulation_criterion(
        points.rows(), points.cols());
    Eigen::ArrayXXd buffer_stop_criterion;
    buffer_retriangulation_criterion.fill(INFINITY);

    // main distmesh loop
    Eigen::ArrayXXi bar_indices;
    for (unsigned step = 0; step < settings::max_steps; ++step) {
        // retriangulate if point movement is above tolerance
        double retriangulation_criterion =
            ((points - buffer_retriangulation_criterion).square().rowwise().sum().sqrt() /
            edge_length_base).maxCoeff();
        if (retriangulation_criterion > settings::retriangulation_tolerance) {
            // update triangulation
            triangulation = triangulation::delaunay(points);

            // calculate circumcenter
            Eigen::ArrayXXd circumcenter = Eigen::ArrayXXd::Zero(
                triangulation.rows(), points.cols());
            for (int point = 0; point < triangulation.cols(); ++point) {
                circumcenter += utils::select_indexed_array_elements<double>(
                    points, triangulation.col(point)) / triangulation.cols();
            }

            // reject triangles with circumcenter outside of the region
            triangulation = utils::select_masked_array_elements<int>(triangulation,
                distance_function(circumcenter) < -settings::general_precision * edge_length_base);

            // find unique bar indices
            bar_indices = utils::find_unique_bars(triangulation);

            buffer_retriangulation_criterion = points;
        }

        // calculate bar vectors and their length
        Eigen::ArrayXXd bar_vector = utils::select_indexed_array_elements<double>(points, bar_indices.col(0)) -
            utils::select_indexed_array_elements<double>(points, bar_indices.col(1));
        Eigen::ArrayXXd bar_length = bar_vector.square().rowwise().sum().sqrt();

        // evaluate edge_length_function at midpoints of bars
        Eigen::ArrayXXd hbars = edge_length_function(0.5 *
            (utils::select_indexed_array_elements<double>(points, bar_indices.col(0)) +
            utils::select_indexed_array_elements<double>(points, bar_indices.col(1))));

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
            if (bar_indices(bar, 0) >= fixed_points.rows()) {
                points.row(bar_indices(bar, 0)) += settings::deltaT * force_vector.row(bar);
            }
            if (bar_indices(bar, 1) >= fixed_points.rows()) {
                points.row(bar_indices(bar, 1)) -= settings::deltaT * force_vector.row(bar);
            }
        }

        // project points outside of domain to boundary
        utils::project_points_to_function(distance_function,
            edge_length_base, points);

        // stop criterion
        double stop_criterion = ((points - buffer_stop_criterion)
            .square().rowwise().sum().sqrt() / edge_length_base)
            .maxCoeff();
        if (stop_criterion < settings::point_movement_tolerance) {
            break;
        }
    }

    return std::make_tuple(points, triangulation);
}

// determine boundary edges of given triangulation
Eigen::ArrayXXi distmesh::boundedges(Eigen::Ref<const Eigen::ArrayXXi> triangulation) {
    std::set<std::vector<int>> edge_set;
    std::vector<std::vector<int>> boundary_edges;
    std::vector<int> edge(triangulation.cols() - 1);

    // find all possible combinations of edges
    auto combinations = utils::n_over_k(triangulation.cols(), triangulation.cols() - 1);

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
