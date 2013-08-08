// libdistmesh
//
// Copyright (C) 2013  Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de

#ifndef LIBDISTMESH_INCLUDE_SETTINGS_H
#define LIBDISTMESH_INCLUDE_SETTINGS_H

// namespace libdistmesh::settings
namespace distmesh {
namespace settings {
    static const dtype::real retriangulation_tolerance = 1e-1;
    static const dtype::real point_movement_tolerance = 1e-3;
    static const dtype::real general_precision = 1e-1;
    static const dtype::real deltaT = 1e-1;
    static const dtype::index max_steps = 10000;
}
}

#endif
