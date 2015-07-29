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

#ifndef _c7492357_ec3f_4dbf_b941_9175e9f79ab0
#define _c7492357_ec3f_4dbf_b941_9175e9f79ab0

// standard c++ lib
#include <functional>
#include <tuple>

// Eigen lib for array handling
#include <Eigen/Core>

// libdistmesh includes
#include "functional.h"
#include "distance_function.h"
#include "utils.h"

namespace distmesh {
    // apply the distmesh algorithm
    std::tuple<Eigen::ArrayXXd, Eigen::ArrayXXi> distmesh(
        Functional const& distanceFunction, double const initialPointDistance,
        Functional const& elementSizeFunction=1.0,
        Eigen::Ref<Eigen::ArrayXXd const> const boundingBox=utils::boundingBox(2),
        Eigen::Ref<Eigen::ArrayXXd const> const fixedPoints=Eigen::ArrayXXd());
}

#endif
