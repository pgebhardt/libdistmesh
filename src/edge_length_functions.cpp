// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#include "distmesh/distmesh.h"

// uniform edge length
std::function<distmesh::dtype::array<distmesh::dtype::real>(
    distmesh::dtype::array<distmesh::dtype::real>&)>
    distmesh::edge_length_functions::uniform() {
    return [=](dtype::array<dtype::real>& points) {
        return dtype::array<dtype::real>::Ones(points.rows(), 1).eval();
    };
}
