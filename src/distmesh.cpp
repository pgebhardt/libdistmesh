// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#include "distmesh/distmesh.h"

// apply the distmesh algorithm
std::tuple<std::shared_ptr<distmesh::dtype::array<distmesh::dtype::real>>,
    std::shared_ptr<distmesh::dtype::array<distmesh::dtype::index>>>
    distmesh::distmesh(
    std::function<dtype::real(dtype::array<dtype::real>)> distance_function,
    std::function<dtype::real(dtype::array<dtype::real>)> edge_length_function,
    dtype::real initial_edge_length, dtype::array<dtype::real> bounding_box) {
    // create initial distribution in bounding_box
    auto points = meshgen::create_point_list(distance_function,
        edge_length_function, initial_edge_length, bounding_box);

    // create initial triangulation
    auto triangulation = triangulation::delaunay(points);

    // create array of points of last iteration, but initialize with inf for
    // first iteration
    auto old_points = std::make_shared<dtype::array<dtype::real>>(
        points->rows(), points->cols());
    // old_points->fill(INFINITY);

    // main distmesh loop
    while (1) {
        // retriangulate if point movement is above tolerance
        auto retriangulation_criterium =
            (((*points) - (*old_points)).square().rowwise().sum() /
            initial_edge_length).maxCoeff();
        if (retriangulation_criterium > settings::retriangulation_tolerance) {

        }
        break;
    }

    return std::make_tuple(points, triangulation);
}
