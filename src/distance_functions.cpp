// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#include "distmesh/distmesh.h"

// generate new distance function with difference of two ones
std::function<distmesh::dtype::real(
    distmesh::dtype::array<distmesh::dtype::real>)>
    distmesh::distance_functions::diff(
    std::function<dtype::real(dtype::array<dtype::real>)> function1,
    std::function<dtype::real(dtype::array<dtype::real>)> function2) {
    return [=](dtype::array<dtype::real> point) {
        return std::max(function1(point), -function2(point));
    };
}

// creates distance function of rectangular domain
std::function<distmesh::dtype::real(
    distmesh::dtype::array<distmesh::dtype::real>)>
    distmesh::distance_functions::rectangular(
    dtype::array<dtype::real> corners) {
    return [=](dtype::array<dtype::real> point) {
        return -std::min(std::min(
            std::min(-corners(0, 1) + point(1), corners(1, 1) - point(1)),
            -corners(0, 0) + point(0)), corners(1, 0) - point(0));
    };
}

// creates distance function for circular domains
std::function<distmesh::dtype::real(
    distmesh::dtype::array<distmesh::dtype::real>)>
    distmesh::distance_functions::circular(
    dtype::array<dtype::real> midpoint, dtype::real radius) {
    return [=](dtype::array<dtype::real> point) {
        return std::sqrt((point - midpoint).square().sum()) - radius;
    };
}
