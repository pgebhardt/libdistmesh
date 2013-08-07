// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#ifndef LIBDISTMESH_INCLUDE_TRIANGULATION_H
#define LIBDISTMESH_INCLUDE_TRIANGULATION_H

// namespace distmesh::triangulation
namespace distmesh {
namespace triangulation {
    // create delaunay triangulation from points array
    std::shared_ptr<dtype::array<dtype::index>> delaunay(
        std::shared_ptr<dtype::array<dtype::real>> points);
}
}

#endif
