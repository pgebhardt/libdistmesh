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
// Copyright (C) 2015 Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de
// --------------------------------------------------------------------

#ifndef LIBDISTMESH_INCLUDE_TRIANGULATION_H
#define LIBDISTMESH_INCLUDE_TRIANGULATION_H

namespace distmesh {
namespace triangulation {
    // create delaunay triangulation from points array
    Eigen::ArrayXXi delaunay(Eigen::Ref<Eigen::ArrayXXd const> const points);
}
}

#endif
