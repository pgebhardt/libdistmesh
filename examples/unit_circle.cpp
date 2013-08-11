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
    // bounding box in which the algorithm tries to create points
    distmesh::dtype::array<distmesh::dtype::real> bounding_box(2, 2);
    bounding_box << -1.0, 1.0, -1.0, 1.0;

    auto mesh = distmesh::distmesh(
        distmesh::distance_function::circular(1.0),
        distmesh::edge_length_function::uniform(),
        0.2, bounding_box);

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
