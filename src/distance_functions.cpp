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
    dtype::array<dtype::real> rectangle) {
    return [=](dtype::array<dtype::real>& points) {
        dtype::array<dtype::real> result(points.rows(), 1);
        result = (points.col(0) - rectangle(0, 0))
            .min(rectangle(0, 1) - points.col(0));
        for (dtype::index dim = 1; dim < points.cols(); ++dim) {
            result = result
                .min((points.col(dim) - rectangle(dim, 0)))
                .min(rectangle(dim, 1) - points.col(dim));
        }
        result = -result;
        return result;
    };
}

// creates distance function for circular domains
std::function<distmesh::dtype::array<distmesh::dtype::real>(
    distmesh::dtype::array<distmesh::dtype::real>&)>
    distmesh::distance_functions::circular(
    dtype::real radius, dtype::array<dtype::real> midpoint) {
    return [=](dtype::array<dtype::real>& points) {
        // move points towards midpoint
        dtype::array<dtype::real> norm_points(points.rows(), points.cols());
        if (midpoint.cols() == points.cols()) {
            for (dtype::index dim = 0; dim < points.cols(); ++dim) {
                norm_points.col(dim) = points.col(dim) - midpoint(0, dim);
            }
        } else {
            norm_points = points;
        }

        // apply circle equation
        dtype::array<dtype::real> result(points.rows(), 1);
        result = norm_points.square().rowwise().sum().sqrt() - radius;
        return result;
    };
}
