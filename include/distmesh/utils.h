// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#ifndef LIBDISTMESH_INCLUDE_UTILS_H
#define LIBDISTMESH_INCLUDE_UTILS_H

// namespace distmesh::utils
namespace distmesh {
namespace utils {
    // create point list
    std::shared_ptr<dtype::array<dtype::real>> create_point_list(
        std::function<dtype::real(dtype::array<dtype::real>)> distance_function,
        std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> edge_length_function,
        dtype::real initial_edge_length, dtype::array<dtype::real> bounding_box);

    // find unique bars
    std::shared_ptr<dtype::array<dtype::index>> find_unique_bars(
        std::shared_ptr<dtype::array<dtype::real>> points,
        std::shared_ptr<dtype::array<dtype::index>> triangulation);

    // project points outside of boundary back to it
    void project_points_to_function(
        std::function<dtype::real(dtype::array<dtype::real>)> distance_function,
        dtype::real initial_edge_length, std::shared_ptr<dtype::array<dtype::real>> points);
}
}

#endif
