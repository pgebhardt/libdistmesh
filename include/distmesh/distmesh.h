// libdistmesh
//
// Copyright (C) 2013 Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#ifndef LIBDISTMESH_INCLUDE_DISTMESH_H
#define LIBDISTMESH_INCLUDE_DISTMESH_H

#include "common.h"
#include "dtype.h"
#include "settings.h"
#include "triangulation.h"
#include "utils.h"
#include "distance_functions.h"

// namespace distmesh
namespace distmesh {
    // apply the distmesh algorithm
    std::tuple<std::shared_ptr<dtype::array<dtype::real>>,
        std::shared_ptr<dtype::array<dtype::index>>> distmesh(
        std::function<dtype::real(dtype::array<dtype::real>)> distance_function,
        std::function<dtype::real(dtype::array<dtype::real>)> edge_length_function,
        dtype::real initial_edge_length, dtype::array<dtype::real> bounding_box);

    // determine boundary edges of given triangulation
    std::shared_ptr<dtype::array<dtype::index>> boundedges(
        std::shared_ptr<dtype::array<dtype::index>> triangulation);
}

#endif
