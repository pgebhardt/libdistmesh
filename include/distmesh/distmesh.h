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

#ifndef LIBDISTMESH_INCLUDE_DISTMESH_H
#define LIBDISTMESH_INCLUDE_DISTMESH_H

// std c++ lib for handling functors
#include <functional>

// Eigen lib for array handling
#include <Eigen/Dense>

// libdistmesh includes
#include "functional.h"
#include "distance_function.h"

namespace distmesh {
    // easy creation of n-dimensional bounding_box
    Eigen::ArrayXXd boundingBox(unsigned const dimension);

    // apply the distmesh algorithm
    std::tuple<Eigen::ArrayXXd, Eigen::ArrayXXi> distmesh(
        Functional const& distanceFunction, double const baseEdgeLength,
        Functional const& edge_length_function=1.0,
        Eigen::Ref<Eigen::ArrayXXd const> const boundingBox=distmesh::boundingBox(2),
        Eigen::Ref<Eigen::ArrayXXd const> const fixedPoints=Eigen::ArrayXXd());

    // determine boundary edges of given triangulation
    Eigen::ArrayXXi boundEdges(Eigen::Ref<Eigen::ArrayXXi const> const triangulation);
}

#endif
