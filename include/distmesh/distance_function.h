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
// along with libDistMesh. If not, see <http://www.gnu.org/licenses/>.
//
// Copyright (C) 2015 Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de
// --------------------------------------------------------------------

#ifndef _7f312bb5_67a5_4eb7_bde7_a4baaed8a573
#define _7f312bb5_67a5_4eb7_bde7_a4baaed8a573

namespace distmesh {
namespace distanceFunction {
    // creates distance function for a nd rectangular domain
    // Attention: Not a real distance function at the corners of domainm
    // you have to give the corners as fixed points to distmesh algorithm
    Functional rectangular(Eigen::Ref<Eigen::ArrayXXd const> const rectangle);

    // creates the true distance function for a 2d rectangular domain
    Functional rectangle(Eigen::Ref<Eigen::ArrayXXd const> const rectangle);

    // creates distance function for elliptical domains
    // Note: not a real distance function but a level function,
    // which is sufficient
    Functional elliptical(Eigen::Ref<Eigen::ArrayXd const> const radii=Eigen::ArrayXd(),
        Eigen::Ref<Eigen::ArrayXd const> const midpoint=Eigen::ArrayXd());

    // creates the true distance function for circular domains
    Functional circular(double const radius=1.0,
        Eigen::Ref<Eigen::ArrayXd const> const midpoint=Eigen::ArrayXd());

    // creates distance function for a 2d domain described by polygon
    // Attention: Not a real distance function at the corners of domainm
    // you have to give the corners as fixed points to distmesh algorithm
    Functional polygon(Eigen::Ref<Eigen::ArrayXXd const> const polygon);
}
}

#endif
