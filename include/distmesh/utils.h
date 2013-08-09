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
    // create point list
    std::shared_ptr<dtype::array<dtype::real>> create_point_list(
        std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> distance_function,
        std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> edge_length_function,
        dtype::real initial_edge_length, dtype::array<dtype::real> bounding_box,
        dtype::array<dtype::real> fixed_points);

    // find unique bars
    std::shared_ptr<dtype::array<dtype::index>> find_unique_bars(
        std::shared_ptr<dtype::array<dtype::real>> points,
        std::shared_ptr<dtype::array<dtype::index>> triangulation);

    // project points outside of boundary back to it
    void project_points_to_function(
        std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> distance_function,
        dtype::real initial_edge_length, std::shared_ptr<dtype::array<dtype::real>> points);
}
}

#endif
