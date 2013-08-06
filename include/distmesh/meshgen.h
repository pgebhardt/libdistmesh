// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#ifndef LIBDISTMESH_INCLUDE_MESHGEN_H
#define LIBDISTMESH_INCLUDE_MESHGEN_H

// namespace distmesh::meshgen
namespace distmesh {
namespace meshgen {
    // create evenly distributed node list
    std::shared_ptr<dtype::matrix<dtype::real>> even_node_list(
        std::shared_ptr<dtype::matrix<dtype::real>> boundaries,
        dtype::real distance);
}
}

#endif
