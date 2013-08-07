// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#include "distmesh/distmesh.h"
#include <random>

// create point list
std::shared_ptr<distmesh::dtype::array<distmesh::dtype::real>>
    distmesh::meshgen::create_point_list(
    std::function<dtype::real(dtype::array<dtype::real>)> distance_function,
    std::function<dtype::real(dtype::array<dtype::real>)> edge_length_function,
    dtype::real initial_edge_length, dtype::array<dtype::real> bounding_box) {
    // calculate max number of points per dimension and max total point count
    // and create initial array
    dtype::array<dtype::index> max_points_per_dimension(bounding_box.rows(), 1);
    dtype::index max_point_count = 1;
    for (dtype::index dim = 0; dim < bounding_box.rows(); ++dim) {
        max_points_per_dimension(dim, 0) = 1 +
            (bounding_box(dim, 1) - bounding_box(dim, 0)) / initial_edge_length;
        max_point_count *= max_points_per_dimension(dim, 0);
    }
    auto initial_points = std::make_shared<dtype::array<dtype::real>>(
        max_point_count, bounding_box.rows());

    // fill point list with evenly distributed points
    dtype::index same_value_count = 1;
    for (dtype::index dim = 0; dim < bounding_box.rows(); ++dim) {
        for (dtype::index point = 0; point < max_point_count; ++point) {
            (*initial_points)(point, dim) = bounding_box(dim, 0) +
                initial_edge_length * ((point / same_value_count) %
                max_points_per_dimension(dim, 0));
        }
        same_value_count *= max_points_per_dimension(dim, 0);
    }

    // reject points outside of region defined by distance_function
    auto inside_points = std::make_shared<dtype::array<dtype::real>>(
        max_point_count, bounding_box.rows());
    dtype::index inside_point_count = 0;
    for (dtype::index point = 0; point < initial_points->rows(); ++point) {
        if (distance_function(initial_points->row(point)) <
            settings::general_precision * initial_edge_length) {
            inside_points->row(inside_point_count) = initial_points->row(point);
            inside_point_count++;
        }
    }
    inside_points->conservativeResize(inside_point_count, bounding_box.rows());

    // initialize random number generator
    std::default_random_engine random_generator;
    std::uniform_real_distribution<dtype::real> random_distribution(0.0, 1.0);

    // calculate propability to keep point in point list based on
    // edge_length_function
    auto propability = std::make_shared<dtype::array<dtype::real>>(
        inside_points->rows(), 1);
    for (dtype::index point = 0; point < propability->rows(); ++point) {
        (*propability)(point, 0) = edge_length_function(inside_points->row(point));
    }

    // reject points with wrong propability
    auto final_points = std::make_shared<dtype::array<dtype::real>>(
        inside_point_count, bounding_box.rows());
    auto propability_norm = propability->minCoeff();
    dtype::index final_point_count = 0;
    for (dtype::index point = 0; point < propability->rows(); ++point) {
        if (random_distribution(random_generator) <
            std::pow(propability_norm / (*propability)(point, 0),
                bounding_box.rows())) {
            final_points->row(final_point_count) = inside_points->row(point);
            final_point_count++;
        }
    }
    final_points->conservativeResize(final_point_count, bounding_box.rows());

    return final_points;
}

// find unique bars
std::shared_ptr<distmesh::dtype::array<distmesh::dtype::index>>
    distmesh::meshgen::find_unique_bars(
    std::shared_ptr<dtype::array<dtype::real>> points,
    std::shared_ptr<dtype::array<dtype::index>> triangulation) {
    // create initial list of bar indices
    dtype::array<dtype::index> initial_bar_indices(
        (points->cols() + 1) * triangulation->rows(), 2);
    for (dtype::index triangle = 0; triangle < triangulation->rows(); ++triangle)
    for (dtype::index bar = 0; bar < points->cols() + 1; ++bar) {
        initial_bar_indices(bar + triangle * (points->cols() + 1), 0) =
            (*triangulation)(triangle, bar);
        initial_bar_indices(bar + triangle * (points->cols() + 1), 1) =
            (*triangulation)(triangle, (bar + 1) % (points->cols() + 1));
    }

    // sort bar indices rowwise
    for (dtype::index bar = 0; bar < initial_bar_indices.rows(); ++bar) {
        if (initial_bar_indices(bar, 0) > initial_bar_indices(bar, 1)) {
            dtype::real temp = initial_bar_indices(bar, 0);
            initial_bar_indices(bar, 0) = initial_bar_indices(bar, 1);
            initial_bar_indices(bar, 1) = temp;
        }
    }

    // reject duplicated bars
    auto bar_indices = std::make_shared<dtype::array<dtype::index>>(
        initial_bar_indices.rows(), initial_bar_indices.cols());
    bar_indices->fill(0);
    dtype::index bar_count = 0;
    for (dtype::index bar = 0; bar < initial_bar_indices.rows(); ++bar) {
        bool unique = true;
        for (dtype::index unique_bar = 0; unique_bar < bar_count; ++unique_bar) {
            if ((initial_bar_indices.row(bar) == bar_indices->row(unique_bar)).all()) {
                unique = false;
                break;
            }
        }
        if (unique == true) {
            bar_indices->row(bar_count) = initial_bar_indices.row(bar);
            bar_count++;
        }
    }
    bar_indices->conservativeResize(bar_count, bar_indices->cols());

    return bar_indices;
}
