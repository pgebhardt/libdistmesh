// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#ifndef LIBDISTMESH_INCLUDE_DISTANCE_FUNCTIONS_H
#define LIBDISTMESH_INCLUDE_DISTANCE_FUNCTIONS_H

// namespace distmesh::distance_functions
namespace distmesh {
namespace distance_functions {
    // generate new distance function with difference of two ones
    std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> diff(
        std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> function1,
        std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> function2);

    // creates distance function of rectangular domain
    std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> rectangular(
        dtype::array<dtype::real> rectangle);

    // creates distance function for circular domains
    std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> circular(
        dtype::array<dtype::real> midpoint, dtype::real radius);
}
}

#endif
