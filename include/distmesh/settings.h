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

#ifndef LIBDISTMESH_INCLUDE_SETTINGS_H
#define LIBDISTMESH_INCLUDE_SETTINGS_H

namespace distmesh {
namespace settings {
    static const double retriangulationTolerance = 1e-1;
    static const double pointMovementTolerance = 1e-3;
    static const double generalPrecision = 1e-3;
    static const double deltaT = 1e-1;
    static const unsigned maxSteps = 100000;
    static const unsigned randomSeed = std::chrono::system_clock::now().time_since_epoch().count();
}
}

#endif
