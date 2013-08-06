// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#include "distmesh/distmesh.h"

std::shared_ptr<distmesh::dtype::matrix<distmesh::dtype::real>>
    distmesh::meshgen::even_node_list(
    std::shared_ptr<dtype::matrix<dtype::real>> boundaries, dtype::real distance) {
    // calculate number of nodes per dimension
    dtype::matrix<dtype::index> nodes_per_dimension(boundaries->rows(), 1);
    dtype::index total_node_count = 1;
    for (dtype::index dim = 0; dim < boundaries->rows(); ++dim) {
        nodes_per_dimension(dim, 1) =
            ((*boundaries)(dim, 1) - (*boundaries)(dim, 0)) / distance;
        total_node_count *= nodes_per_dimension(dim, 0);
    }

    // initialize empty node list
    auto node_list = std::make_shared<dtype::matrix<dtype::real>>(
        total_node_count, boundaries->rows());

    // fill node_list
    dtype::index same_value_count = 1;
    for (dtype::index dim = 0; dim < boundaries->rows(); ++dim) {
        for (dtype::index node = 0; node < total_node_count; ++node) {
            (*node_list)(node, dim) = (*boundaries)(dim, 0) +
                distance * ((node / same_value_count) %
                nodes_per_dimension(dim, 0));
        }
        same_value_count *= nodes_per_dimension(dim, 0);
    }

    return node_list;
}
