// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#ifndef LIBDISTMESH_INCLUDE_EDGE_LENGTH_FUNCTIONS_H
#define LIBDISTMESH_INCLUDE_EDGE_LENGTH_FUNCTIONS_H

// namespace distmesh::edge_length_functions
namespace distmesh {
namespace edge_length_functions {
    // uniform edge length
    std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> uniform();
}
}

#endif
