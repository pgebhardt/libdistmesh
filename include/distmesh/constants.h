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
namespace constants {
    // algorithm stops, when maximum relative points movement is below tolerance
    static double const pointsMovementTolerance = 1e-3;

    // triangulation is updated, when maximum relative points movement is above tolerance
    static double const retriangulationTolerance = 1e-1;

    // time step for updating points positions with Euler's method
    static double const deltaT = 1e-1;

    // relative tolerance in the geometry evaluations
    static double const geomertyEvaluationTolerance = 1e-3;

    // algorithm will be terminated after the maximum number of iterations,
    // when no convergence can be achieved
    static unsigned const maxSteps = 10000;
}
}

#endif
