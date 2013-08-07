// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#include "distmesh/distmesh.h"

// creates distance function for circular domains
std::function<distmesh::dtype::real(
    distmesh::dtype::array<distmesh::dtype::real>)>
    distmesh::distance_functions::circular(
    dtype::array<dtype::real> midpoint, dtype::real radius) {
    return [=](dtype::array<dtype::real> point) {
        return std::sqrt((point - midpoint).square().sum()) - radius;
    };
}
