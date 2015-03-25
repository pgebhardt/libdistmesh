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

#ifndef _57ca4052_ba33_4235_9a5f_84154336d924
#define _57ca4052_ba33_4235_9a5f_84154336d924

namespace distmesh {
namespace settings {
    static double const retriangulationTolerance = 1e-1;
    static double const pointMovementTolerance = 1e-3;
    static double const generalPrecision = 1e-3;
    static double const deltaT = 1e-1;
    static unsigned const maxSteps = 10000;
}
}

#endif
