// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#include "distmesh/distmesh.h"

// creates distance function for circular domains
std::function<distmesh::dtype::array<distmesh::dtype::real>(
    distmesh::dtype::array<distmesh::dtype::real>&)>
    distmesh::edge_length_functions::uniform() {
    return [=](dtype::array<dtype::real>& points) {
        return Eigen::MatrixXd::Ones(points.rows(), 1).array();
    };
}
