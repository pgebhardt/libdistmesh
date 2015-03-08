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

#ifndef LIBDISTMESH_INCLUDE_UTILS_H
#define LIBDISTMESH_INCLUDE_UTILS_H

// namespace distmesh::utils
namespace distmesh {
namespace utils {
    // select array elements based on mask
    template <
        class type
    >
    Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> select_masked_array_elements(
        Eigen::Ref<const Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic>> array,
        Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>> mask) {
        Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> result(mask.count(), array.cols());

        int result_count = 0;
        for (int row = 0; row < array.rows(); ++row) {
            if (mask(row, 0)) {
                result.row(result_count) = array.row(row);
                result_count++;
            }
        }

        return result;
    }

    // select array elements based on indices
    template <
        class type
    >
    Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> select_indexed_array_elements(
        Eigen::Ref<const Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic>> array,
        Eigen::Ref<const Eigen::ArrayXXi> indices) {
        Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> result(indices.rows(), array.cols());

        for (int row = 0; row < indices.rows(); ++row) {
            result.row(row) = array.row(indices(row, 0));
        }

        return result;
    }

    // calculate factorial recursively
    unsigned factorial(unsigned n);

    // create point list
    Eigen::ArrayXXd create_point_list(Functional distance_function,
        double edge_length_base, Functional edge_length_function,
        Eigen::Ref<const Eigen::ArrayXXd> bounding_box,
        Eigen::Ref<const Eigen::ArrayXXd> fixed_points);

    // create array with all unique combinations n over k
    Eigen::ArrayXXi n_over_k(unsigned n, unsigned k);

    // find unique bars
    Eigen::ArrayXXi find_unique_bars(Eigen::Ref<const Eigen::ArrayXXi> triangulation);

    // project points outside of boundary back to it
    void project_points_to_function(
        Functional distance_function, double edge_length_base,
        Eigen::Ref<Eigen::ArrayXXd> points);

    // check whether points lies inside or outside of polygon
    Eigen::ArrayXXd points_inside_poly(
        Eigen::Ref<const Eigen::ArrayXXd> points,
        Eigen::Ref<const Eigen::ArrayXXd> polygon);
}
}

#endif
