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

#ifndef LIBDISTMESH_INCLUDE_DISTANCE_FUNCTION_H
#define LIBDISTMESH_INCLUDE_DISTANCE_FUNCTION_H

// namespace distmesh::distance_function
namespace distmesh {
namespace distance_function {
    // creates distance function of rectangular domain
    Functional rectangular(Eigen::Ref<Eigen::ArrayXXd const> const rectangle);

    // creates distance function for elliptical domains
    // Note: not a real distance function but a level function,
    // which is sufficient
    Functional elliptical(Eigen::Ref<Eigen::ArrayXd const> const radii=Eigen::ArrayXd(),
        Eigen::Ref<Eigen::ArrayXd const> const midpoint=Eigen::ArrayXd());

    // creates distance function for circular domains
    Functional circular(double const radius=1.0,
        Eigen::Ref<Eigen::ArrayXd const> const midpoint=Eigen::ArrayXd());

    // creates distance function for domain described by polygon
    // Only 2D domains supported
    Functional polygon(Eigen::Ref<Eigen::ArrayXXd const> const polygon);
}
}

#endif
