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

#include <distmesh/distmesh.h>
#include <fstream>

int main() {
    // fixed points at the corners of domain to guarantee convergence
    Eigen::ArrayXXd fixed_points(4, 2);
    fixed_points << -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0;

    // midpoints of eyes, mouth and nose
    Eigen::ArrayXXd midpoints(4, 2);
    midpoints << -0.5, 0.5, 0.5, 0.5, 0.0, -0.5, 0.0, 0.1;

    // radii of eliptical mouth and nose
    Eigen::ArrayXXd radii(2, 2);
    radii << 0.75, 0.75 / std::sqrt(10.0), 0.15 / std::sqrt(3.0), 0.15;

    auto distance_function = distmesh::distance_function::rectangular(
        distmesh::bounding_box(2))
        .max(-distmesh::distance_function::circular(0.25, midpoints.row(0)))
        .max(-distmesh::distance_function::circular(0.25, midpoints.row(1)))
        .max(-0.75 * distmesh::distance_function::elliptical(radii.row(0), midpoints.row(2)))
        .max(-0.15 * distmesh::distance_function::elliptical(radii.row(1), midpoints.row(3)));

    // create mesh
    auto mesh = distmesh::distmesh(distance_function, 0.02,
        0.2 - distance_function, distmesh::bounding_box(2),
        fixed_points);

    // plot mesh
    std::ofstream points_file;
    std::ofstream triangulation_file;
    points_file.open("points.txt");
    triangulation_file.open("triangulation.txt");

    points_file << std::get<0>(mesh);
    triangulation_file << std::get<1>(mesh);

    points_file.close();
    triangulation_file.close();

    system("python plot_mesh.py");

    return 0;
}
