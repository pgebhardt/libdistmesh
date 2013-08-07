// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#ifndef LIBDISTMESH_INCLUDE_DISTANCE_FUNCTIONS_H
#define LIBDISTMESH_INCLUDE_DISTANCE_FUNCTIONS_H

// namespace distmesh::distance_functions
namespace distmesh {
namespace distance_functions {
    // creates distance function for circular domains
    std::function<dtype::real(dtype::array<dtype::real>)> circular(
        dtype::array<dtype::real> midpoint, dtype::real radius);
}
}

#endif
