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

#ifndef LIBDISTMESH_INCLUDE_DISTANCE_FUNCTIONS_H
#define LIBDISTMESH_INCLUDE_DISTANCE_FUNCTIONS_H

// namespace distmesh::distance_functions
namespace distmesh {
namespace distance_functions {
    // generate new distance function with difference of two ones
    std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> diff(
        std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> function1,
        std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> function2);

    // creates distance function of rectangular domain
    std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> rectangular(
        dtype::array<dtype::real> rectangle);

    // creates distance function for circular domains
    std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> circular(
        dtype::real radius=1.0,
        dtype::array<dtype::real> midpoint=dtype::array<dtype::real>());
}
}

#endif
