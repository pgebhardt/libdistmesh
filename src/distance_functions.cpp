// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#include "distmesh/distmesh.h"

// generate new distance function with difference of two ones
std::function<distmesh::dtype::array<distmesh::dtype::real>(
    distmesh::dtype::array<distmesh::dtype::real>&)>
    distmesh::distance_functions::diff(
    std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> function1,
    std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> function2) {
    return [=](dtype::array<dtype::real>& points) {
        dtype::array<dtype::real> result(points.rows(), 1);
        result = function1(points).max(-function2(points));
        return result;
    };
}

// creates distance function of rectangular domain
std::function<distmesh::dtype::array<distmesh::dtype::real>(
    distmesh::dtype::array<distmesh::dtype::real>&)>
    distmesh::distance_functions::rectangular(
    dtype::array<dtype::real> corners) {
    return [=](dtype::array<dtype::real>& points) {
        dtype::array<dtype::real> result(points.rows(), 1);
        result = (-corners(0, 1) + points.col(1));
        result = result.max(corners(1, 1) - points.col(1));
        result = result.max(-corners(0, 0) + points.col(0));
        result = result.max(corners(1, 0) - points.col(0));
        return result;
    };
}

// creates distance function for circular domains
std::function<distmesh::dtype::array<distmesh::dtype::real>(
    distmesh::dtype::array<distmesh::dtype::real>&)>
    distmesh::distance_functions::circular(
    dtype::array<dtype::real> midpoint, dtype::real radius) {
    return [=](dtype::array<dtype::real>& points) {
        dtype::array<dtype::real> norm_points(points.rows(), points.cols());
        dtype::array<dtype::real> result(points.rows(), 1);
        for (dtype::index dim = 0; dim < points.cols(); ++dim) {
            norm_points.col(dim) = points.col(dim) - midpoint(0, dim);
        }
        result = norm_points.square().rowwise().sum().sqrt() - radius;
        return result;
    };
}
