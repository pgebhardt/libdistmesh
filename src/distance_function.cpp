// --------------------------------------------------------------------
// This file is part of libDistMesh.
//
// libDistMesh is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// libDistMesh is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with libDistMesh.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright (C) 2013 Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de
// --------------------------------------------------------------------

#include "distmesh/distmesh.h"

// generate new distance function with difference of two ones
distmesh::distance_function::function_t
    distmesh::distance_function::diff(function_t function1,
    function_t function2) {
    return DISTMESH_DISTANCE_FUNCTION({
        return function1(points).max(-function2(points)).eval();
    });
}

// creates distance function of rectangular domain
distmesh::distance_function::function_t
    distmesh::distance_function::rectangular(
    dtype::array<dtype::real> rectangle) {
    return DISTMESH_DISTANCE_FUNCTION({
        dtype::array<dtype::real> result(points.rows(), 1);
        result = (points.col(0) - rectangle(0, 0))
            .min(rectangle(0, 1) - points.col(0));
        for (dtype::index dim = 1; dim < points.cols(); ++dim) {
            result = result
                .min((points.col(dim) - rectangle(dim, 0)))
                .min(rectangle(dim, 1) - points.col(dim));
        }
        return (-result).eval();
    });
}

// creates distance function for circular domains
distmesh::distance_function::function_t
    distmesh::distance_function::circular(
    dtype::real radius, dtype::array<dtype::real> midpoint) {
    return DISTMESH_DISTANCE_FUNCTION({
        // move points towards midpoint
        dtype::array<dtype::real> norm_points;
        if (midpoint.cols() == points.cols()) {
            norm_points = points.rowwise() - midpoint.row(0);
        } else {
            norm_points = points;
        }

        // apply circle equation
        return (norm_points.square().rowwise().sum().sqrt() - radius)
            .eval();
    });
}
