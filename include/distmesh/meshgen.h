// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#ifndef LIBDISTMESH_INCLUDE_MESHGEN_H
#define LIBDISTMESH_INCLUDE_MESHGEN_H

// namespace distmesh::meshgen
namespace distmesh {
namespace meshgen {
    // create point list
    std::shared_ptr<dtype::array<dtype::real>> create_point_list(
        std::function<dtype::real(dtype::array<dtype::real>)> distance_function,
        std::function<
            dtype::real(dtype::array<dtype::real>)> edge_length_function,
        dtype::real initial_edge_length, dtype::array<dtype::real> bounding_box);

    // find unique bars
    std::shared_ptr<dtype::array<dtype::index>> find_unique_bars(
        std::shared_ptr<dtype::array<dtype::real>> points,
        std::shared_ptr<dtype::array<dtype::index>> triangulation);

    // calculate forces on each point based on lengths of each bar
    std::shared_ptr<dtype::array<dtype::real>> calculate_forced(
        std::shared_ptr<dtype::array<dtype::real>> points,
        std::shared_ptr<dtype::array<dtype::index>> bars_vector,
        std::function<
            dtype::real(dtype::array<dtype::real>)
            > edge_length_function);
}
}

#endif
