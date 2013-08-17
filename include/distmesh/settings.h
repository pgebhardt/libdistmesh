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

#ifndef LIBDISTMESH_INCLUDE_SETTINGS_H
#define LIBDISTMESH_INCLUDE_SETTINGS_H

// namespace distmesh::settings
namespace distmesh {
namespace settings {
    static const dtype::real retriangulation_tolerance = 1e-1;
    static const dtype::real point_movement_tolerance = 1e-3;
    static const dtype::real general_precision = 1e-1;
    static const dtype::real deltaT = 1e-1;
    static const dtype::index max_steps = 10000;
    static const unsigned random_seed = std::chrono::system_clock::now().time_since_epoch().count();
}
}

#endif
