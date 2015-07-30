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

#include <iostream>
#include <distmesh/distmesh.h>
#include "helper.h"

int main() {
    distmesh::helper::HighPrecisionTime time;

    // create mesh
    Eigen::ArrayXXd points;
    Eigen::ArrayXXi elements;

    std::tie(points, elements) = distmesh::distmesh(
        distmesh::distanceFunction::rectangle(distmesh::utils::boundingBox(2))
            .max(-distmesh::distanceFunction::circular(0.5)),
        0.05, 0.05 + 0.3 * distmesh::distanceFunction::circular(0.5),
        distmesh::utils::boundingBox(2));

    // print mesh properties and elapsed time
    std::cout << "Created mesh with " << points.rows() << " points and " << elements.rows() <<
        " elements in " << time.elapsed() * 1e3 << " ms." << std::endl;

    // save mesh to file
    distmesh::helper::savetxt<double>(points, "points.txt");
    distmesh::helper::savetxt<int>(elements, "triangulation.txt");

    // plot mesh using python
    return system("python plot_mesh.py");
}
