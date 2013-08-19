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
    distmesh::dtype::array<distmesh::dtype::real> bounding_box(2, 2);
    bounding_box << 0.0, 1.0, 0.0, 1.0;

    // fixed points at corners of domain to guarantee convergence
    distmesh::dtype::array<distmesh::dtype::real> fixed_points(4, 2);
    fixed_points << 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0;

    // distance function for rectangular domain
    auto distance_function =
        distmesh::distance_function::rectangular(bounding_box);

    // corner points of polygon
    distmesh::dtype::array<distmesh::dtype::real> poly(2, 2);
    poly << 0.3, 0.7, 0.7, 0.5;

    // edge size function
    auto size_function =
        (0.01 + 0.3 * distmesh::distance_function::circular(0.0).abs())
        .min(0.025 + 0.3 * distmesh::distance_function::polygon(poly).abs())
        .min(0.15);

    // create mesh
    auto mesh = distmesh::distmesh(distance_function, 0.01,
        size_function, bounding_box, fixed_points);

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
