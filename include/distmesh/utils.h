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
    dtype::array<type> select_masked_array_elements(
        const Eigen::Ref<dtype::array<type>> array,
        const Eigen::Ref<dtype::array<bool>> mask) {
        dtype::array<type> result(array.rows(), array.cols());
        dtype::index result_count = 0;
        for (dtype::index row = 0; row < array.rows(); ++row) {
            if (mask(row, 0)) {
                result.row(result_count) = array.row(row);
                result_count++;
            }
        }
        result.conservativeResize(result_count, array.cols());

        return result;
    }

    // select array elements based on indices
    template <
        class type
    >
    dtype::array<type> select_indiced_array_elements(
        const Eigen::Ref<dtype::array<type>> array,
        const Eigen::Ref<dtype::array<dtype::index>> indices) {
        dtype::array<type> result(indices.rows(), array.cols());
        for (dtype::index row = 0; row < indices.rows(); ++row) {
            result.row(row) = array.row(indices(row, 0));
        }

        return result;
    }


    // create point list
    std::shared_ptr<dtype::array<dtype::real>> create_point_list(
        distance_function::function_t distance_function,
        edge_length_function::function_t edge_length_function,
        dtype::real initial_edge_length,
        dtype::array<dtype::real> bounding_box,
        dtype::array<dtype::real> fixed_points);

    // find unique bars
    std::shared_ptr<dtype::array<dtype::index>> find_unique_bars(
        std::shared_ptr<dtype::array<dtype::real>> points,
        std::shared_ptr<dtype::array<dtype::index>> triangulation);

    // project points outside of boundary back to it
    void project_points_to_function(
        distance_function::function_t distance_function,
        dtype::real initial_edge_length,
        std::shared_ptr<dtype::array<dtype::real>> points);
}
}

#endif
