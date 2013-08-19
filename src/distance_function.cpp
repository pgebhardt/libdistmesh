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

// creates distance function of rectangular domain
distmesh::Functional distmesh::distance_function::rectangular(
    dtype::array<dtype::real> rectangle) {
    return DISTMESH_FUNCTIONAL({
        dtype::array<dtype::real> result =
            (points.col(0) - rectangle(0, 0))
            .min(rectangle(0, 1) - points.col(0));
        for (dtype::index dim = 1; dim < points.cols(); ++dim) {
            result = result
                .min((points.col(dim) - rectangle(dim, 0)))
                .min(rectangle(dim, 1) - points.col(dim));
        }
        return -result;
    });
}

// creates distance function for elliptical domains
// Note: not a real distance function but a level function,
// which is sufficient
distmesh::Functional distmesh::distance_function::elliptical(
    dtype::array<dtype::real> radii, dtype::array<dtype::real> midpoint) {
    return DISTMESH_FUNCTIONAL({
        if (midpoint.cols() == points.cols()) {
            if (radii.cols() == points.cols()) {
                return ((points.rowwise() - midpoint.row(0)).rowwise() / radii.row(0))
                    .square().rowwise().sum().sqrt() - 1.0;
            } else {
                return (points.rowwise() - midpoint.row(0))
                    .square().rowwise().sum().sqrt() - 1.0;
            }
        } else {
            if (radii.cols() == points.cols()) {
                return (points.rowwise() / radii.row(0))
                    .square().rowwise().sum().sqrt() - 1.0;
            } else {
                return points.square().rowwise().sum().sqrt() - 1.0;
            }
        }
    });
}

// creates distance function for circular domains
distmesh::Functional
    distmesh::distance_function::circular(
    dtype::real radius, dtype::array<dtype::real> midpoint) {
    return DISTMESH_FUNCTIONAL({
        if (midpoint.cols() == points.cols()) {
            return (points.rowwise() - midpoint.row(0))
                .square().rowwise().sum().sqrt() - radius;
        } else {
            return points.square().rowwise().sum().sqrt() - radius;
        }
    });
}

// creates distance function for domain described by polygon
distmesh::Functional distmesh::distance_function::polygon(
    const Eigen::Ref<dtype::array<dtype::real>>& polygon) {
    return DISTMESH_FUNCTIONAL({
        dtype::array<dtype::real> v(points.rows(), 2);
        dtype::array<dtype::real> w(points.rows(), 2);
        dtype::array<dtype::real> c1(points.rows(), 1);
        dtype::array<dtype::real> c2(points.rows(), 1);
        dtype::array<dtype::real> distance(points.rows(), polygon.rows());

        for (dtype::index i = 0, j = polygon.rows() - 1;
            i < polygon.rows(); j = i++) {
            v.rowwise() = polygon.row(i) - polygon.row(j);
            w = points.rowwise() - polygon.row(j);

            c1 = v.col(0) * w.col(0) + v.col(1) * w.col(1);
            c2 = v.col(0) * v.col(0) + v.col(1) * v.col(1);

            distance.col(i) =
                (c1 <= 0.0).select(
                    (points.rowwise() - polygon.row(j)).square().rowwise().sum().sqrt(),
                (c1 >= c2).select(
                    (points.rowwise() - polygon.row(i)).square().rowwise().sum().sqrt(),
                    (points - ((v.colwise() * (c1 / c2).col(0)).rowwise() + polygon.row(j)))
                        .square().rowwise().sum().sqrt()
                ));
        }

        return (1.0 - 2.0 * utils::points_inside_poly(points, polygon)) *
            distance.rowwise().minCoeff();
    });
}
