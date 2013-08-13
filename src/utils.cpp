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
#include <random>
#include <chrono>
#include <limits>
#include <set>
#include <array>

// create point list
distmesh::dtype::array<distmesh::dtype::real> distmesh::utils::create_point_list(
    functional::function_t distance_function,
    functional::function_t edge_length_function,
    dtype::real edge_length_base,
    dtype::array<dtype::real> bounding_box,
    dtype::array<dtype::real> fixed_points) {
    // calculate max number of points per dimension and
    // max total point coun and create initial array
    dtype::array<dtype::index> max_points_per_dimension(
        bounding_box.rows(), 1);
    dtype::index max_point_count = 1;
    for (dtype::index dim = 0; dim < bounding_box.rows(); ++dim) {
        max_points_per_dimension(dim, 0) = 1 +
            (bounding_box(dim, 1) - bounding_box(dim, 0)) /
            edge_length_base;
        max_point_count *= max_points_per_dimension(dim, 0);
    }
    dtype::array<dtype::real> initial_points(
        max_point_count, bounding_box.rows());

    // fill point list with evenly distributed points
    dtype::index same_value_count = 1;
    for (dtype::index dim = 0; dim < bounding_box.rows(); ++dim) {
        for (dtype::index point = 0; point < max_point_count; ++point) {
            initial_points(point, dim) = bounding_box(dim, 0) +
                edge_length_base * ((point / same_value_count) %
                max_points_per_dimension(dim, 0));
        }
        same_value_count *= max_points_per_dimension(dim, 0);
    }

    // reject points outside of region defined by distance_function
    dtype::array<bool> inside = distance_function(initial_points) < settings::general_precision * edge_length_base;
    dtype::array<dtype::real> inside_points = select_masked_array_elements<dtype::real>(
        initial_points, inside);

    // initialize random number generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine random_generator(seed);
    std::uniform_real_distribution<dtype::real> random_distribution(0.0, 1.0);

    // calculate propability to keep point in point list based on
    // edge_length_function
    dtype::array<dtype::real> propability = edge_length_function(inside_points);

    // add fixed points to final list first
    dtype::array<dtype::real> final_points(inside_points.rows() + fixed_points.rows(),
        bounding_box.rows());
    for (dtype::index fixed_point = 0; fixed_point < fixed_points.rows(); ++fixed_point) {
        final_points.row(fixed_point) = fixed_points.row(fixed_point);
    }

    // reject points with wrong propability
    dtype::real propability_norm = propability.minCoeff();
    dtype::index final_point_count = fixed_points.rows();
    for (dtype::index point = 0; point < propability.rows(); ++point) {
        if (random_distribution(random_generator) <
            std::pow(propability_norm / propability(point, 0),
                bounding_box.rows())) {
            final_points.row(final_point_count + fixed_points.rows()) = inside_points.row(point);
            final_point_count++;
        }
    }
    final_points.conservativeResize(final_point_count + fixed_points.rows(),
        bounding_box.rows());

    return final_points;
}

// find unique bars
distmesh::dtype::array<distmesh::dtype::index> distmesh::utils::find_unique_bars(
    const Eigen::Ref<dtype::array<dtype::real>>& points,
    const Eigen::Ref<dtype::array<dtype::index>>& triangulation) {
    // fill set of sorted bar indices
    std::set<std::array<dtype::index, 2>> bar_indices_set;
    std::array<dtype::index, 2> bar = {{0, 0}};
    for (dtype::index triangle = 0; triangle < triangulation.rows(); ++triangle)
    for (dtype::index i = 0; i < points.cols() + 1; ++i) {
        if (triangulation(triangle, i) > triangulation(triangle, (i + 1) % (points.cols() + 1))) {
            bar[0] = triangulation(triangle, i);
            bar[1] = triangulation(triangle, (i + 1) % (points.cols() + 1));
        } else {
            bar[0] = triangulation(triangle, (i + 1) % (points.cols() + 1));
            bar[1] = triangulation(triangle, i);
        }
        bar_indices_set.insert(bar);
    }

    // copy set to eigen array
    dtype::array<dtype::index> bar_indices(bar_indices_set.size(), 2);
    dtype::index bar_index = 0;
    for (auto & bar : bar_indices_set) {
        bar_indices(bar_index, 0) = bar[0];
        bar_indices(bar_index, 1) = bar[1];
        bar_index++;
    }

    return bar_indices;
}

// project points outside of boundary back to it
void distmesh::utils::project_points_to_function(
    functional::function_t distance_function,
    dtype::real edge_length_base,
    Eigen::Ref<dtype::array<dtype::real>> points) {
    // evaluate distance function at points
    dtype::array<dtype::real> distance = distance_function(points);

    // check for points outside of boundary
    dtype::array<bool> outside = distance > 0.0;
    if (outside.any()) {
        // calculate gradient
        dtype::array<dtype::real> gradient(points.rows(), points.cols());
        dtype::array<dtype::real> deltaX(points.rows(), points.cols());
        dtype::array<dtype::real> h;
        deltaX.fill(0.0);
        for (dtype::index dim = 0; dim < points.cols(); ++dim) {
            deltaX.col(dim).fill(std::sqrt(std::numeric_limits<dtype::real>::epsilon()) * edge_length_base);
            h = points + deltaX;
            gradient.col(dim) = (distance_function(h) - distance) /
                (std::sqrt(std::numeric_limits<dtype::real>::epsilon()) * edge_length_base);
            deltaX.col(dim).fill(0.0);
        }

        for (dtype::index dim = 0; dim < points.cols(); ++dim) {
            points.col(dim) -= outside.select(
                gradient.col(dim) * distance.col(0) / gradient.square().rowwise().sum(),
                0.0);
        }
    }
}
